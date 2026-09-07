"""Per-ply stress from ``*EL PRINT, S``, corrected into the material frame.

ccx will not report a reaction moment for a composite shell, and it will not
report co-rotated stress either. What it does report is measured in
``cases/cantilever_ansys/README.md``; the parts that shape this module:

- **The printed tensor does not co-rotate.** It is the Cauchy stress expressed
  in a basis frozen at t=0 -- the ply's ``*ORIENTATION`` frame for a bare ``S``
  card, or the global axes under ``GLOBAL=YES``. Under NLGEOM the material
  rotates away from that basis, so at 89 degrees the raw output claims 392 MPa
  of through-thickness normal stress on a strip that physically carries none.
  Reading ``sxx`` as "fibre-direction stress" at large deflection produces a
  plausible number that is badly wrong, which is the whole reason this module
  exists.
- **Undoing it is possible because the curvature is developable.** ``geometry``
  only supports single curvature, so the material rotates about the width axis
  alone. That axis' normal stress passes through untouched, the
  (axial, through-thickness) block rotates within itself, and the
  (axial-width, normal-width) shear pair rotates as a vector.
- **The correction is checkable.** A thin shell carries no through-thickness
  normal stress, so ``sigma_33`` after un-rotation is a residual that is near
  zero when the angle and its sign are right and large when they are not.
  ``material_frame`` returns it rather than asserting it away.

The rotation angle comes from the deformed geometry (``frd.deformed_midsurface``),
not from the stress. The stress can supply its own estimate -- the principal
direction of the rotating block -- but using that and then checking the
through-thickness residual would be circular, so it is offered separately as a
cross-check.

``GLOBAL=YES`` is the frame to ask for. A bare ``S`` card puts every ply in its
own t=0 orientation frame, so plies at different angles arrive in different
bases and cannot be compared until each is rotated back; the global card puts
them all in one basis for the cost of a ply-angle rotation we have to do anyway.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from compfea.run import SolveError

# " stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set BLADE and time  0.1E+01"
_STRESS_HEADER = re.compile(
    r"^\s*stresses \(elem, integ\.pnt\.,.*\) for set (?P<elset>\S+) "
    r"and time\s+(?P<time>\S+)\s*$"
)
_TENSOR = ("sxx", "syy", "szz", "sxy", "sxz", "syz")
STRESS_COLUMNS = ("increment", "time", "elset", "elem", "ip", *_TENSOR, "ori", "frame")

#: ccx expands a composite shell into one solid layer per ply, 8 points each.
IPS_PER_PLY = 8


def parse_dat_stress(path: str | Path) -> pd.DataFrame:
    """Parse ``*EL PRINT ... S`` blocks from a ``.dat``.

    One row per (element, integration point, time). ``elem`` is the **shell**
    element id from the deck, not an expanded solid id -- verified on a 16
    element deck that printed elements 1-16. This is the opposite convention to
    the ``.frd``, whose node ids are the expanded mesh; the two outputs do not
    share a numbering scheme and must never be joined on id.

    The frame is detected from the row width rather than trusted from the card:
    a local row carries a trailing ``*ORIENTATION`` name and has 9 fields, a
    ``GLOBAL=YES`` row has 8 and no name.
    """
    path = Path(path)
    if not path.is_file():
        raise SolveError(f"ccx wrote no results file at {path}")
    lines = path.read_text().splitlines()
    rows: list[dict] = []
    for i, line in enumerate(lines):
        header = _STRESS_HEADER.match(line)
        if header is None:
            continue
        time = float(header["time"])
        elset = header["elset"].lower()
        # A blank line separates the header from the first record.
        for row in _rows_after(lines, i + 1):
            rows.append({"time": time, "elset": elset, **row})
    if not rows:
        raise SolveError(
            f"no stress blocks in {path}; the deck needs "
            f"*EL PRINT, ELSET=..., GLOBAL=YES with the S label"
        )
    frame = pd.DataFrame(rows)
    frame.insert(0, "increment", frame["time"].rank(method="dense").astype(int))
    return frame[list(STRESS_COLUMNS)]


def _rows_after(lines: list[str], start: int) -> list[dict]:
    """Records following a header, skipping the blank line ccx puts first."""
    out: list[dict] = []
    for line in lines[start:]:
        fields = line.split()
        if not fields:
            if out:  # blank line after data ends the block
                break
            continue  # the blank line ccx writes right after the header
        if len(fields) not in (8, 9) or not fields[0].isdigit():
            break
        try:
            values = [float(v) for v in fields[2:8]]
        except ValueError:
            break
        out.append(
            {
                "elem": int(fields[0]),
                "ip": int(fields[1]),
                **dict(zip(_TENSOR, values, strict=True)),
                "ori": fields[8] if len(fields) == 9 else None,
                "frame": "local" if len(fields) == 9 else "global",
            }
        )
    return out


def ips_per_ply_from(max_ip: int, n_plies: int) -> int:
    """Integration points per ply, measured from an element's own point count.

    ccx uses 8 per expanded layer for an S8R and **6 for an S6** -- measured, by
    solving a 4-ply composite of each. Assuming 8 on a triangle is worse than
    wrong: a 4-ply S6 tops out at ip 24, so ``ply_of_ip``'s ``ply >= n_plies``
    guard never fires, plies are silently mis-assigned, and the top ply is never
    reached at all. ``geometry`` puts S6 elements in the same ELSET as S8R, so a
    mesh with a few stray triangles would report the peak in the wrong ply with
    no diagnostic.
    """
    if n_plies < 1:
        raise ValueError(f"a stack needs at least one ply, got {n_plies}")
    ips, remainder = divmod(max_ip, n_plies)
    if remainder or ips < 1:
        raise ValueError(
            f"{max_ip} integration points do not divide into {n_plies} plies; "
            "the element type or the stack is not what this code assumes"
        )
    return ips


def ply_of_ip(ip: int, n_plies: int, ips_per_ply: int = IPS_PER_PLY) -> int:
    """Zero-based ply index for a 1-based integration point.

    Thickness is the slowest index and **ply 0 is the -z ply** -- the first line
    of the ``*SHELL SECTION, COMPOSITE`` card. Getting this backwards silently
    mirrors every unsymmetric stack, which is the failure the smoke case's
    reversed-stack check exists to catch.
    """
    if ip < 1:
        raise ValueError(f"integration points are 1-based, got {ip}")
    ply = (ip - 1) // ips_per_ply
    if ply >= n_plies:
        raise ValueError(
            f"ip {ip} implies ply {ply}, but the stack has {n_plies} plies; "
            f"ips_per_ply={ips_per_ply} is probably wrong for this element type"
        )
    return ply


def _axes(long_axis: str) -> tuple[int, int, int, float]:
    """(width, axial, normal) indices and the sign of the width basis vector.

    The fibre transform needs a **right-handed** (axial, width, normal) triple,
    and one of the two span directions does not give one from the raw axes.
    With the span along x the triple is (x, y, z) and x_hat x y_hat = +z_hat, so
    the width axis is +y. With the span along y it is (y, x, z) and
    y_hat x x_hat = -z_hat: left-handed, so the width basis vector is **-x**.

    That sign matters only for the shear term, which is why a cross-ply layup
    cannot detect it -- tau_12 is ~0 there, and sigma_11/sigma_22 are unchanged
    at 0 and 90 degrees. On a +/-45 stack it swaps sigma_11 with sigma_22, and
    it flips the sign of tau_12 for every ply at every angle. CLAUDE.md flags
    exactly this class of mirrored convention for *ORIENTATION; it applies here
    for the same reason.
    """
    if long_axis == "y":
        return 0, 1, 2, -1.0
    if long_axis == "x":
        return 1, 0, 2, +1.0
    raise ValueError(f"long_axis must be 'x' or 'y', not {long_axis!r}")


def _tensor_matrix(row) -> np.ndarray:
    return np.array(
        [
            [row.sxx, row.sxy, row.sxz],
            [row.sxy, row.syy, row.syz],
            [row.sxz, row.syz, row.szz],
        ]
    )


def unrotate(sigma: np.ndarray, phi_rad: float, width_axis: int) -> np.ndarray:
    """Undo a rigid rotation of ``phi`` about ``width_axis``.

    Returns ``R.T @ sigma @ R``: the same Cauchy stress expressed in the frame
    that has rotated with the material, rather than the frozen t=0 frame ccx
    prints in.
    """
    c, s = math.cos(phi_rad), math.sin(phi_rad)
    r = np.eye(3)
    a, b = [i for i in range(3) if i != width_axis]
    r[a, a] = c
    r[a, b] = -s
    r[b, a] = s
    r[b, b] = c
    return r.T @ sigma @ r


def phi_from_stress(row, long_axis: str = "y") -> float:
    """Rotation angle implied by the stress itself, in radians.

    The principal direction of the rotating (axial, normal) block. Independent
    of the geometry, and therefore useful as a cross-check -- but not as the
    input to the through-thickness residual test, which it would make circular.
    """
    w, a, n, _ = _axes(long_axis)
    sigma = _tensor_matrix(row)
    return 0.5 * math.atan2(2.0 * sigma[a, n], sigma[a, a] - sigma[n, n])


def in_plane_to_fibre(
    s_axial: float, s_width: float, t_aw: float, angle_deg: float
) -> tuple[float, float, float]:
    """Rotate in-plane stress into fibre coordinates for a ply at ``angle_deg``.

    A 0-degree ply runs along the part's long axis (see CLAUDE.md), so the
    laminate axes are (axial, width) and the rotation is the standard plane
    transformation about the shell normal.
    """
    t = math.radians(angle_deg)
    c, s = math.cos(t), math.sin(t)
    s11 = s_axial * c * c + s_width * s * s + 2.0 * t_aw * c * s
    s22 = s_axial * s * s + s_width * c * c - 2.0 * t_aw * c * s
    t12 = (s_width - s_axial) * c * s + t_aw * (c * c - s * s)
    return s11, s22, t12


def material_frame(
    df: pd.DataFrame,
    phi_of_elem: dict[int, float] | float,
    ply_angles: Sequence[float],
    *,
    long_axis: str = "y",
) -> pd.DataFrame:
    """Un-rotate each row and express it in its ply's fibre coordinates.

    ``phi_of_elem`` is the local material rotation in radians, per shell element
    (or one value for all). Adds ``ply``, ``ply_angle_deg``, ``s11``, ``s22``,
    ``t12``, and ``s33_residual`` -- the through-thickness normal stress a shell
    cannot carry, as a fraction of the **peak** in-plane stress in this block.

    That normalisation is deliberate. Scaling each row by its own in-plane
    magnitude looks more natural and is useless: near the neutral axis both the
    residual and the divisor go to zero, and the ratio blows up to tens of
    percent on rows carrying 10 MPa out of a 950 MPa field. Measured against the
    peak, the residual answers the question actually worth asking -- how big is
    the impossible stress next to the stresses that matter.
    """
    w, a, n, w_sign = _axes(long_axis)
    n_plies = len(ply_angles)
    # Measured per element: S8R gives 8 points per ply, S6 gives 6, and a mesh
    # may carry both in one ELSET.
    ips = {
        int(e): ips_per_ply_from(int(g["ip"].max()), n_plies)
        for e, g in df.groupby("elem")
    }
    out = []
    for row in df.itertuples():
        phi = (
            float(phi_of_elem)
            if isinstance(phi_of_elem, (int, float))
            else phi_of_elem[row.elem]
        )
        mat = unrotate(_tensor_matrix(row), phi, w)
        ply = ply_of_ip(row.ip, n_plies, ips[row.elem])
        angle = float(ply_angles[ply])
        # w_sign makes (axial, width, normal) right-handed; sigma_ww is
        # unaffected by flipping the width basis, the shear is not.
        s11, s22, t12 = in_plane_to_fibre(
            mat[a, a], mat[w, w], w_sign * mat[w, a], angle
        )
        out.append(
            {
                "ply": ply,
                "ply_angle_deg": angle,
                "s11": s11,
                "s22": s22,
                "t12": t12,
                "_s33": abs(mat[n, n]),
            }
        )
    frame = pd.concat([df.reset_index(drop=True), pd.DataFrame(out)], axis=1)
    peak = max(frame[["s11", "s22"]].abs().to_numpy().max(), 1e-30)
    frame["s33_residual"] = frame.pop("_s33") / peak
    return frame


def stress_summary(df: pd.DataFrame) -> dict:
    """Peak tensile and compressive stress per component, and where it sits.

    Tension and compression are reported separately because CFRP's tensile and
    compressive strengths differ by roughly a factor of two, so a single
    magnitude would hide which limit a design is actually near.
    """
    out: dict = {}
    if df.empty:
        return out
    for comp in ("s11", "s22"):
        for tag, idx in (("t", df[comp].idxmax()), ("c", df[comp].idxmin())):
            row = df.loc[idx]
            out[f"{comp}_max_{tag}"] = float(row[comp])
            out[f"{comp}_max_{tag}_ply"] = int(row["ply"])
            out[f"{comp}_max_{tag}_angle"] = float(row["ply_angle_deg"])
            out[f"{comp}_max_{tag}_elem"] = int(row["elem"])
    idx = df["t12"].abs().idxmax()
    out["t12_max_abs"] = float(df.loc[idx, "t12"])
    out["t12_max_abs_ply"] = int(df.loc[idx, "ply"])
    out["t12_max_abs_angle"] = float(df.loc[idx, "ply_angle_deg"])
    out["t12_max_abs_elem"] = int(df.loc[idx, "elem"])
    out["s33_residual_max"] = float(df["s33_residual"].max())
    return out


_ELEMENT_HEADER = re.compile(r"^\s*\*ELEMENT\s*(,|$)", re.IGNORECASE)
_NODE_HEADER = re.compile(r"^\s*\*NODE\s*(,|$)", re.IGNORECASE)


def element_spans(inp_path: str | Path, *, long_axis: str = "y") -> dict[int, float]:
    """Undeformed span coordinate of each element's centroid, keyed by element id.

    Parsed from the deck rather than taken from a ``Mesh`` because the ``.dat``
    reports **shell element ids**, and only the deck ties an id to a position.
    Corner nodes only: midside nodes would bias the centroid of a curved edge
    and the span station does not need that precision.
    """
    axis = 0 if long_axis == "x" else 1
    nodes: dict[int, tuple[float, float, float]] = {}
    elements: dict[int, list[int]] = {}
    mode = None
    for line in Path(inp_path).read_text().splitlines():
        if line.lstrip().startswith("*"):
            mode = (
                "node" if _NODE_HEADER.match(line)
                else "elem" if _ELEMENT_HEADER.match(line)
                else None
            )
            continue
        if mode is None:
            continue
        fields = [f.strip() for f in line.split(",") if f.strip()]
        if not fields or not fields[0].lstrip("-").isdigit():
            continue
        if mode == "node" and len(fields) >= 4:
            nodes[int(fields[0])] = tuple(float(v) for v in fields[1:4])
        elif mode == "elem" and len(fields) >= 5:
            elements[int(fields[0])] = [int(v) for v in fields[1:5]]
    if not elements:
        raise ValueError(f"no *ELEMENT block in {inp_path}")
    return {
        eid: sum(nodes[n][axis] for n in conn) / len(conn)
        for eid, conn in elements.items()
        if all(n in nodes for n in conn)
    }


def bend_angle_by_span(mid: pd.DataFrame):
    """Interpolator for the local material rotation from a deformed mid-surface.

    ``phi`` is the tangent angle of the deformed centreline, so it comes from the
    geometry and owes nothing to the stress. That independence is the point: an
    angle fitted to the stress tensor would make the through-thickness residual
    check circular, since the residual would be the minor principal value by
    construction and near zero whatever the true rotation was.

    Takes the frame from ``frd.deformed_midsurface`` and returns
    ``phi(span0_mm) -> radians``, keyed on the **undeformed** span so an element
    can be located by its position in the deck.
    """
    ordered = mid.sort_values("span0_mm")
    s0 = ordered["span0_mm"].to_numpy(dtype=float)
    ds = np.gradient(ordered["span_mm"].to_numpy(dtype=float))
    dz = np.gradient(ordered["z_mm"].to_numpy(dtype=float))
    phi = np.unwrap(np.arctan2(dz, ds))

    def at(span: float) -> float:
        return float(np.interp(float(span), s0, phi))

    return at
