"""Flex and kick-point metrics for a solved tip-weight buckling run.

Post-processing only: no solver, no gmsh. The tip-weight bench
(``cases/fin_test_3/README.md``) clamps the blade at the heel, holds its axis
along a hung weight's line of action, and lets the weight compress it past its
Euler load until the tip tangent lies over. ``run_tipweight.py`` drives a tip
clip patch axially under displacement control and reads the equivalent weight
``W`` as the axial reaction at the clamp; the tip angle ``theta`` is an output.

This module reduces such a run to the two numbers a layup is designed to:

- **flex** -- ``W`` at tip tangent ``theta = 90 deg``, in newtons (~5 N
  supersoft, ~15-20 N hard).
- **kick point** -- where along the span the blade bends most: the arc-length
  location of the peak of a smoothed ``|kappa(s)|``, where ``kappa`` is the rate
  of change with undeformed span of the local material tangent rotation (each
  band segment's deformed chord measured off its own undeformed chord, so the
  blade's built-in camber is not projected away). Reported in mm from the clamp,
  as a fraction of free length, and bucketed ``heel`` (< 1/3) / ``mid`` /
  ``tip`` (> 2/3).

The centreline comes from a row of mid-chord ``*NODE PRINT`` U stations along
the span, read from the ``.dat`` (deck nodes) -- **not** the ``.frd``: ccx's
expanded solid mesh fans one planform station into many wherever the shell
normal tilts, so a midsurface slope fit there is a chord across curvature, not a
tangent (it inflated an earlier fin ``W(90 deg)`` from ~9 to ~13 N).

The flat-strip closed form this is checked against is
``cases/fin_test_3/compressed_elastica.py`` (``curvature_profile`` /
``rotation_median_s``): a uniform column's curvature is maximal at the clamp, so
its kick point is the heel; a stiffness taper that thins the tip moves the peak
outboard.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from compfea.geometry import Mesh
from compfea.ubend import tip_length_mm

BAND_PREFIX = "band_"


# --------------------------------------------------------------------------
# Span / clamp geometry -- shared with cases/fin_test_3/run_tipweight.py, which
# imports these rather than keeping its own copies, so the clamp face and axis
# sense are defined in exactly one place.
# --------------------------------------------------------------------------


def axis_index(long_axis: str) -> int:
    """0 for ``x``, 1 for ``y``; the blade's long axis in a node tuple."""
    if long_axis not in ("x", "y"):
        raise ValueError(f"long_axis must be 'x' or 'y', not {long_axis!r}")
    return 0 if long_axis == "x" else 1


def _station(mesh: Mesh, nset: str, axis: int) -> float:
    ids = mesh.nsets[nset]
    return sum(mesh.nodes[n][axis] for n in ids) / len(ids)


def clamp_face(
    mesh: Mesh, long_axis: str = "y", *, root_nset: str = "fixed_end"
) -> tuple[float, float]:
    """``(s0, sign)``: the clamp edge nearest the tip and the tip-ward direction.

    The CAD origin need not sit at the heel, so the tip can be at the high or the
    low end of the span. ``s0`` is the clamp station a free length is measured
    from; ``sign`` is +1 if the tip is at higher ``long_axis`` than the clamp,
    -1 otherwise. Matches ``ubend.tip_length_mm``'s sense.
    """
    axis = axis_index(long_axis)
    root = [mesh.nodes[n][axis] for n in mesh.nsets[root_nset]]
    tip_mid = _station(mesh, "far_face", axis)
    if tip_mid >= 0.5 * (min(root) + max(root)):
        return max(root), 1.0
    return min(root), -1.0


# --------------------------------------------------------------------------
# .dat U reader -- one implementation, imported by run_tipweight.py.
# --------------------------------------------------------------------------

_DISP_HEADER = re.compile(
    r"^\s*displacements \(vx,vy,vz\) for set (?P<nset>\S+) and time\s+(?P<time>\S+)\s*$"
)


def read_sets_disp(
    dat_path: str | Path, nsets: Sequence[str]
) -> dict[str, dict[float, tuple[float, float, float]]]:
    """Mean ``(ux, uy, uz)`` per printed time for each of ``nsets``, one file pass.

    A ``.dat`` from a many-increment solve reaches tens of MB; a row of band
    stations must not read it once per station. Keyed by the name as passed
    (ccx prints set names upper-case regardless of deck spelling).
    """
    wanted = {n.upper(): n for n in nsets}
    acc: dict[str, dict[float, list[tuple[float, float, float]]]] = {
        n: {} for n in nsets
    }
    key: str | None = None
    t: float | None = None
    for line in Path(dat_path).read_text().splitlines():
        header = _DISP_HEADER.match(line)
        if header:
            key = wanted.get(header["nset"])
            t = float(header["time"]) if key is not None else None
            if key is not None:
                acc[key].setdefault(t, [])
            continue
        if key is None or not line.strip():
            continue
        fields = line.split()
        if len(fields) == 4 and fields[0].isdigit():
            acc[key][t].append(tuple(float(v) for v in fields[1:]))
        else:
            key = None
    return {
        name: {
            time: tuple(sum(v[i] for v in rows) / len(rows) for i in range(3))
            for time, rows in per_time.items()
            if rows
        }
        for name, per_time in acc.items()
    }


def read_set_disp(
    dat_path: str | Path, nset: str
) -> dict[float, tuple[float, float, float]]:
    """Mean ``(ux, uy, uz)`` over ``nset`` per printed time, from ``*NODE PRINT`` U."""
    return read_sets_disp(dat_path, [nset])[nset]


# --------------------------------------------------------------------------
# Band node sets along the span
# --------------------------------------------------------------------------


def span_band_nsets(
    mesh: Mesh,
    *,
    long_axis: str = "y",
    n_bands: int = 21,
    prefix: str = BAND_PREFIX,
    chord_frac: float = 0.5,
    k: int = 5,
    clamp_margin_mm: float = 8.0,
    tip_margin_mm: float = 8.0,
    exclude: frozenset[int] | set[int] | None = None,
) -> dict[str, tuple[int, ...]]:
    """``n_bands`` mid-chord node clusters evenly spaced clamp -> tip.

    Each cluster is up to ``k`` free interior nodes near the chord midpoint
    within a spanwise window half a band-spacing wide, ranked by distance to the
    band station first and chord offset second. Only nodes strictly between
    ``clamp_margin_mm`` outboard of the clamp face and ``tip_margin_mm`` short of
    the tip edge are eligible -- a wide window on a long blade would otherwise
    reach clamped nodes (U == 0, which pulls the centreline's zero) or the driven
    tip clip. Pass the clamp / tip / clip node ids in ``exclude`` to drop them
    outright. Returned names are ``{prefix}00 .. {prefix}{n_bands-1}`` in
    clamp -> tip order. Generalises ``run_tipweight.py:tip_band`` (a single
    window) to a row.
    """
    if n_bands < 5:
        raise ValueError(
            f"n_bands must be >= 5 (tangent_and_curvature needs >= 4 stations "
            f"after the segment midpoints), got {n_bands}"
        )
    axis = axis_index(long_axis)
    across = 1 - axis
    s0, sign = clamp_face(mesh, long_axis)
    free_len = tip_length_mm(mesh, long_axis=long_axis)
    if free_len <= clamp_margin_mm + tip_margin_mm:
        raise ValueError(
            f"free length {free_len:g} mm too short for the "
            f"{clamp_margin_mm:g}+{tip_margin_mm:g} mm band margins"
        )
    drop = frozenset(exclude or ())

    chord = [xyz[across] for xyz in mesh.nodes.values()]
    chord_mid = 0.5 * (min(chord) + max(chord))
    chord_half = chord_frac * 0.5 * (max(chord) - min(chord))
    s_lo, s_hi = clamp_margin_mm, free_len - tip_margin_mm

    targets = np.linspace(s_lo, s_hi, n_bands)
    spacing = float(targets[1] - targets[0])

    # eligible node id -> (spanwise distance from clamp face, chord offset)
    span_of = {
        nid: (sign * (xyz[axis] - s0), abs(xyz[across] - chord_mid))
        for nid, xyz in mesh.nodes.items()
        if nid not in drop
        and s_lo <= sign * (xyz[axis] - s0) <= s_hi
        and abs(xyz[across] - chord_mid) <= chord_half
    }
    # Assign each eligible node to its single nearest band station -- disjoint by
    # construction, and a sparse region (a tapered tip) still fills from whatever
    # is closest rather than widening into a neighbour's territory.
    owned: dict[int, list[int]] = {j: [] for j in range(n_bands)}
    for nid, (s, _dc) in span_of.items():
        j = int(np.argmin(np.abs(targets - s)))
        if abs(s - targets[j]) <= spacing:  # else it is nobody's -- a real gap
            owned[j].append(nid)

    bands: dict[str, tuple[int, ...]] = {}
    for j, s_t in enumerate(targets):
        cand = owned[j]
        if len(cand) < 2:
            raise ValueError(
                f"band {j} at s={s_t:.1f} mm owns {len(cand)} eligible nodes; "
                f"refine the mesh or lower n_bands"
            )
        # Chord offset first: a smooth mid-chord centreline matters more than
        # hitting the exact station (the slab is already ~one spacing wide, and
        # clamp / tip / clip nodes are out of `span_of` entirely, so this cannot
        # pull in a constrained node the way it could before they were excluded).
        cand.sort(
            key=lambda nid: (span_of[nid][1], abs(span_of[nid][0] - s_t))
        )
        picked = cand[: max(2, k)]
        bands[f"{prefix}{j:02d}"] = tuple(sorted(picked))
    return bands


# --------------------------------------------------------------------------
# Deformed centreline -> tangent angle -> curvature
# --------------------------------------------------------------------------


def centerline_from_dat(
    dat_path: str | Path,
    mesh: Mesh,
    *,
    long_axis: str = "y",
    prefix: str = BAND_PREFIX,
) -> pd.DataFrame:
    """Undeformed and deformed mid-chord centreline per band per printed time.

    Columns: ``time``, ``band``, ``s0_mm`` (**undeformed arc length** from the
    clamp face along the band centroids, tip-ward positive -- the axial gap to
    the first band plus the running chord length between centroids, so blade
    camber is not projected away), ``x0_mm``/``y0_mm``/``z0_mm`` (undeformed band
    centroid) and ``x_mm``/``y_mm``/``z_mm`` (that centroid plus the band's mean
    U at that time). Sorted by ``time`` then ``s0_mm``.
    """
    axis = axis_index(long_axis)
    s0_clamp, sign = clamp_face(mesh, long_axis)
    names = sorted(n for n in mesh.nsets if n.startswith(prefix))
    if not names:
        raise ValueError(f"no band nsets with prefix {prefix!r} on the mesh")

    # Band centroids in clamp -> tip order, then developed span length along
    # them: the axial gap from the clamp face to the first centroid, then
    # cumulative length in the (long_axis, z) bend plane -- chord (x) wander from
    # which nodes each band happened to pick is not part of the span, and the
    # repo models only single developable curvature. Reduces to the axial
    # projection on a flat blade.
    c0 = {name: np.mean([mesh.nodes[n] for n in mesh.nsets[name]], axis=0)
          for name in names}
    order = sorted(names, key=lambda n: sign * (c0[n][axis] - s0_clamp))
    s_arc = {order[0]: sign * (c0[order[0]][axis] - s0_clamp)}
    for a, b in zip(order[:-1], order[1:], strict=True):
        d = c0[b] - c0[a]
        s_arc[b] = s_arc[a] + float(math.hypot(d[axis], d[2]))

    disp = read_sets_disp(dat_path, names)
    empty = [n for n in names if not disp[n]]
    if empty:
        raise ValueError(
            f"{len(empty)} band nset(s) have no U block in {dat_path} "
            f"(first: {empty[0]}); the deck was built with a different "
            f"--kick-bands than this mesh"
        )
    rows = []
    for name in names:
        cen = c0[name]
        for t, u in disp[name].items():
            defp = cen + np.asarray(u, float)
            rows.append(
                {
                    "time": t,
                    "band": name,
                    "s0_mm": float(s_arc[name]),
                    "x0_mm": float(cen[0]),
                    "y0_mm": float(cen[1]),
                    "z0_mm": float(cen[2]),
                    "x_mm": float(defp[0]),
                    "y_mm": float(defp[1]),
                    "z_mm": float(defp[2]),
                }
            )
    if not rows:
        raise ValueError(f"no U prints for any band in {dat_path}")
    return (
        pd.DataFrame(rows)
        .sort_values(["time", "s0_mm"])
        .reset_index(drop=True)
    )


def _odd_window(n_points: int, requested: int | None) -> int:
    """A small symmetric (odd) smoothing window.

    Fixed at 3 by default, not a fraction of ``n_points``: the smoothing is only
    there to kill single-station noise in ``|kappa|``, and a window that grows
    with the station count would reshape the profile and walk the peak toward
    whichever side carries more of its mass -- more bands should mean *finer*
    resolution, not heavier smoothing.
    """
    w = int(requested) if requested is not None else 3
    if w % 2 == 0:
        w += 1
    return max(1, min(w, n_points if n_points % 2 else n_points - 1))


def _smooth(y: np.ndarray, window: int) -> np.ndarray:
    """Centred moving average; edges shrink the window so no phase shift."""
    if window <= 1:
        return y.astype(float)
    half = window // 2
    out = np.empty_like(y, dtype=float)
    for i in range(y.size):
        lo = max(0, i - half)
        hi = min(y.size, i + half + 1)
        out[i] = y[lo:hi].mean()
    return out


def tangent_and_curvature(
    centerline_at_t: pd.DataFrame, *, long_axis: str = "y"
) -> pd.DataFrame:
    """Tangent angle and curvature of one time slice of the centreline.

    For each station-to-station segment, ``phi`` is the rotation of **that
    segment alone**: the angle from its undeformed chord to its deformed chord,
    projected onto the ``(long_axis, z)`` bend plane (the repo models single
    developable curvature). Measuring off the segment's own undeformed chord,
    not the global axis, cancels the blade's built-in camber. Same 2-D
    ``atan2`` form as ``run_tipweight.py:station_tangent_deg``, one segment at a
    time -- so ``phi`` is the local material tangent angle at the segment
    midpoint, **not** a running sum (summing absolute slopes would run away).
    The whole array's sign is normalised so the largest rotation is positive,
    matching the bench's 0..180 deg convention whichever way the blade buckled;
    a genuine S-curve still shows as a ``kappa`` sign change along the span.

    ``kappa = dphi/ds`` from the raw (already smooth, integral-quantity) ``phi``
    by ``np.gradient`` -- central differences inside, a plain one-sided slope at
    each end. Pre-smoothing ``phi`` would flatten a real curvature peak sitting
    against the clamp; peak finding does its own smoothing on ``|kappa|`` in
    ``kick_point``. Rows are the ``len - 1`` segment midpoints. Columns:
    ``s0_mm`` (midpoint arc position), ``phi_rad``, ``kappa_1pmm``.
    """
    df = centerline_at_t.sort_values("s0_mm").reset_index(drop=True)
    if len(df) < 4:
        raise ValueError(f"need >= 4 bands for a curvature, got {len(df)}")
    axis = axis_index(long_axis)
    p0 = df[["x0_mm", "y0_mm", "z0_mm"]].to_numpy()
    pdef = df[["x_mm", "y_mm", "z_mm"]].to_numpy()
    s0 = df["s0_mm"].to_numpy()

    s_mid = 0.5 * (s0[:-1] + s0[1:])
    phi = np.zeros(len(df) - 1)
    for i in range(len(df) - 1):
        b = p0[i + 1] - p0[i]
        a = pdef[i + 1] - pdef[i]
        b2 = np.array([b[axis], b[2]])
        a2 = np.array([a[axis], a[2]])
        phi[i] = math.atan2(
            b2[0] * a2[1] - b2[1] * a2[0], float(np.dot(b2, a2))
        )

    if phi[int(np.argmax(np.abs(phi)))] < 0.0:
        phi = -phi

    kappa = np.gradient(phi, s_mid, edge_order=1)
    return pd.DataFrame({"s0_mm": s_mid, "phi_rad": phi, "kappa_1pmm": kappa})


# --------------------------------------------------------------------------
# Kick point
# --------------------------------------------------------------------------


def _bucket(frac: float) -> str:
    if frac < 1.0 / 3.0:
        return "heel"
    if frac < 2.0 / 3.0:
        return "mid"
    return "tip"


@dataclass(frozen=True)
class KickPoint:
    """Where the blade bends most, at one deformation state."""

    s_kick_mm: float
    s_kick_frac: float
    bucket: str
    s_peak_raw_mm: float
    s_rotmed_mm: float
    s_rotmed_frac: float
    plateau_frac: float
    theta_deg: float
    warning: str | None


def kick_point(curv_df: pd.DataFrame, *, free_len_mm: float) -> KickPoint:
    """Reduce a ``tangent_and_curvature`` frame to a kick point.

    ``s_kick`` is the peak of a lightly smoothed ``|kappa(s)|`` -- "where you
    see the most bend" -- refined by a 3-point parabola vertex for sub-station
    resolution. ``s_rotmed`` (station splitting the tip rotation in half) is
    reported as context, not as an agreement check: for a monotone curvature
    (uniform column) the two legitimately differ. ``warning`` fires only when
    the peak is genuinely not localisable: a wide interior plateau, or two
    comparable maxima far apart.
    """
    interior = curv_df.dropna(subset=["kappa_1pmm"]).reset_index(drop=True)
    if len(interior) < 4:
        raise ValueError("need >= 4 curvature points")
    s = interior["s0_mm"].to_numpy()
    absk = np.abs(interior["kappa_1pmm"].to_numpy())

    j_raw = int(np.argmax(absk))
    s_peak_raw = float(s[j_raw])
    absk_s = _smooth(absk, _odd_window(len(absk), None))
    j = int(np.argmax(absk_s))
    s_kick = float(s[j])
    if j == 0 or j == len(s) - 1:
        # Peak sits against a boundary (a monotone-ish taper): the metric's
        # resolution is one band spacing; report the edge band, no parabola.
        pass
    else:
        y0, y1, y2 = absk_s[j - 1], absk_s[j], absk_s[j + 1]
        denom = y0 - 2.0 * y1 + y2
        if denom != 0.0:
            # Parabola vertex; a real maximum keeps the vertex within +/- half a
            # sample of j. Clip there so a near-flat triple cannot jump a station.
            offset = float(np.clip(0.5 * (y0 - y2) / denom, -0.5, 0.5))
            step = s[j + 1] - s[j] if offset >= 0 else s[j] - s[j - 1]
            s_kick = float(np.clip(s[j] + offset * step, s[0], s[-1]))
    frac = s_kick / free_len_mm

    # phi is already the absolute material tangent angle at each station, so the
    # station splitting the tip rotation in half is where |phi| reaches half its
    # tip value. maximum.accumulate guards the interp against a non-monotone
    # wiggle near the clamp where phi ~ 0.
    phi_abs = np.abs(curv_df["phi_rad"].to_numpy())
    s_all = curv_df["s0_mm"].to_numpy()
    tip_rot = phi_abs[-1]
    s_rotmed = (
        float(np.interp(0.5 * tip_rot, np.maximum.accumulate(phi_abs), s_all))
        if tip_rot > 0
        else float(s_all[0])
    )

    peak = absk_s[j]
    wide = absk_s >= 0.8 * peak
    lo = j
    while lo - 1 >= 0 and wide[lo - 1]:
        lo -= 1
    hi = j
    while hi + 1 < len(wide) and wide[hi + 1]:
        hi += 1
    plateau_frac = float(s[hi] - s[lo]) / free_len_mm
    at_boundary = frac < 0.1 or frac > 0.9

    warning: str | None = None
    if not at_boundary and plateau_frac > 0.4:
        warning = (
            f"curvature is within 80% of its peak over {plateau_frac:.0%} of the "
            "span -- the kick point is a broad hinge, not a point"
        )
    else:
        # A second local max >= 0.7 * peak, separated by >= 0.25 L with a
        # valley <= 0.5 * peak between it and the main peak.
        for m in range(1, len(absk_s) - 1):
            if m == j or absk_s[m] < 0.7 * peak:
                continue
            if absk_s[m] < absk_s[m - 1] or absk_s[m] < absk_s[m + 1]:
                continue
            lo, hi = sorted((j, m))
            far = s[hi] - s[lo] >= 0.25 * free_len_mm
            valley = absk_s[lo : hi + 1].min() <= 0.5 * peak
            if far and valley:
                warning = (
                    "two comparable curvature maxima far apart -- 'kick point' "
                    "is ill-posed for this blade"
                )
                break

    theta_deg = float(math.degrees(curv_df["phi_rad"].to_numpy()[-1]))
    return KickPoint(
        s_kick_mm=s_kick,
        s_kick_frac=frac,
        bucket=_bucket(frac),
        s_peak_raw_mm=s_peak_raw,
        s_rotmed_mm=s_rotmed,
        s_rotmed_frac=s_rotmed / free_len_mm,
        plateau_frac=plateau_frac,
        theta_deg=theta_deg,
        warning=warning,
    )


# --------------------------------------------------------------------------
# Flex + kick, tied to a solved run
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class FlexKick:
    """The design numbers for one solved tip-weight run.

    ``flex_n`` is ``W`` at ``theta = 90 deg`` as the caller measured it from the
    canonical clip -> tip_band angle (``run_tipweight.collect``); ``None`` if the
    solve never reached 90 deg. The kick points are evaluated at the ``.dat``
    increments nearest the caller's headline / migration times, so they sit at
    the same load states the flex number is read at. ``KickPoint.theta_deg`` is
    the band-centreline tip tangent at that increment -- reported for
    transparency, a few degrees stiffer than the clip angle because the
    outermost band segment curls past the clip -> tip_band chord.
    """

    flex_n: float | None
    free_len_mm: float
    headline: KickPoint
    migration: KickPoint
    migration_delta_mm: float

    def as_dict(self) -> dict:
        d = asdict(self)
        d["kick_bucket"] = self.headline.bucket
        d["kick_s_mm"] = self.headline.s_kick_mm
        d["kick_s_frac"] = self.headline.s_kick_frac
        return d


def tip_tangent_series(
    dat_path: str | Path,
    mesh: Mesh,
    *,
    long_axis: str = "y",
    prefix: str = BAND_PREFIX,
) -> pd.DataFrame:
    """Per printed time: the band-centreline tip tangent, ``[time, theta_deg]``.

    ``theta`` is the last segment's ``phi`` (already sign-normalised positive).
    A few degrees stiffer than ``run_tipweight.collect``'s clip -> tip_band
    angle, because the outermost band segment curls past that chord -- use
    ``collect`` for the number of record, this only to locate load states.
    """
    cl = centerline_from_dat(dat_path, mesh, long_axis=long_axis, prefix=prefix)
    rows = [
        {
            "time": float(t),
            "theta_deg": math.degrees(
                tangent_and_curvature(sub, long_axis=long_axis)["phi_rad"].iloc[-1]
            ),
        }
        for t, sub in cl.groupby("time")
    ]
    return pd.DataFrame(rows).sort_values("time").reset_index(drop=True)


def flex_kick(
    dat_path: str | Path,
    mesh: Mesh,
    *,
    time_headline: float,
    time_migration: float,
    flex_n: float | None,
    long_axis: str = "y",
    prefix: str = BAND_PREFIX,
    time_min: float = 0.0,
    out_dir: str | Path | None = None,
) -> FlexKick:
    """Kick point at two load states, bundled with a caller-supplied flex number.

    ``dat_path`` is the run's ``job.dat`` from a **converged** solve (this reads
    the ``.dat`` directly and does not re-check the ``.sta`` -- the caller runs
    ``run.solve``, which raises on a partial run). The deck must carry a row of
    ``span_band_nsets`` U prints (``run_tipweight.py --kick-bands`` adds them).
    ``time_headline`` / ``time_migration`` are ccx TOT TIMEs -- the caller maps
    them from its canonical ``theta -> time`` curve (e.g. 90 deg and 15 deg), so
    the kick point is evaluated at the same increment the flex number is. The
    nearest available band-print time to each is used. ``time_min`` drops a
    gravity-settle pre-step (pass 1.0 when ``run_tipweight`` was run with
    ``--gravity``). Writes ``kick_kappa.csv``, ``kick_kappa.svg`` (|kappa| vs
    span) and ``kick_shape.svg`` (the deformed centreline, coloured by |kappa|,
    kick point ringed) into ``out_dir`` (default: the run directory,
    ``dat_path``'s grandparent).
    """
    dat_path = Path(dat_path)
    out = Path(out_dir) if out_dir is not None else dat_path.parents[1]
    free_len = tip_length_mm(mesh, long_axis=long_axis)

    cl = centerline_from_dat(dat_path, mesh, long_axis=long_axis, prefix=prefix)
    pos_by_time = {
        float(t): sub.sort_values("s0_mm").reset_index(drop=True)
        for t, sub in cl.groupby("time")
        if float(t) > time_min + 1e-9
    }
    per_time = {
        t: tangent_and_curvature(sub, long_axis=long_axis)
        for t, sub in pos_by_time.items()
    }
    if len(per_time) < 2:
        raise ValueError(
            f"only {len(per_time)} band-print increments after time_min="
            f"{time_min:g}; nothing to interpolate"
        )
    times = np.array(sorted(per_time))

    def _at(target_time: float) -> tuple[pd.DataFrame, pd.DataFrame, KickPoint]:
        t = float(times[int(np.abs(times - target_time).argmin())])
        curv = per_time[t]
        return pos_by_time[t], curv, kick_point(curv, free_len_mm=free_len)

    head_pos, head_curv, headline = _at(time_headline)
    mig_pos, mig_curv, migration = _at(time_migration)

    _write_outputs(
        out, long_axis,
        head_pos, head_curv, headline,
        mig_pos, mig_curv, migration,
    )

    return FlexKick(
        flex_n=flex_n,
        free_len_mm=free_len,
        headline=headline,
        migration=migration,
        migration_delta_mm=headline.s_kick_mm - migration.s_kick_mm,
    )


def _write_outputs(
    out: Path,
    long_axis: str,
    head_pos: pd.DataFrame,
    head_curv: pd.DataFrame,
    headline: KickPoint,
    mig_pos: pd.DataFrame,
    mig_curv: pd.DataFrame,
    migration: KickPoint,
) -> None:
    out.mkdir(parents=True, exist_ok=True)

    hi = head_curv.rename(
        columns={"phi_rad": "phi_headline", "kappa_1pmm": "kappa_headline"}
    )
    mig = mig_curv.rename(
        columns={"phi_rad": "phi_migration", "kappa_1pmm": "kappa_migration"}
    )
    merged = hi.merge(mig, on="s0_mm", how="outer").sort_values("s0_mm")
    merged.to_csv(out / "kick_kappa.csv", index=False)

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
    except Exception:  # pragma: no cover - plotting is optional
        return

    # --- kick_kappa.svg: |kappa| against span --------------------------------
    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    ax.plot(
        hi["s0_mm"], hi["kappa_headline"].abs(), "o-", lw=1.4, ms=3,
        label=f"|kappa| at theta~{headline.theta_deg:.0f} deg",
    )
    ax.plot(
        mig["s0_mm"], mig["kappa_migration"].abs(), "s--", lw=1.1, ms=3,
        label=f"|kappa| at theta~{migration.theta_deg:.0f} deg",
    )
    ax.axvline(headline.s_kick_mm, color="C3", lw=1.4, label="kick point")
    ax.axvline(
        headline.s_rotmed_mm, color="0.6", lw=1.0, ls=":", label="rotation median"
    )
    ax.set_xlabel("arc length from clamp s (mm)")
    ax.set_ylabel("curvature |kappa| (1/mm)")
    ax.set_title(
        f"kick point: {headline.bucket.upper()}  (f = {headline.s_kick_frac:.2f})"
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "kick_kappa.svg", format="svg")
    plt.close(fig)

    # --- kick_shape.svg: the deformed centreline, coloured by |kappa| -------
    axis = axis_index(long_axis)
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(
        head_pos[f"{'xyz'[axis]}0_mm"], head_pos["z0_mm"],
        "-", color="0.75", lw=2.0, label="undeformed", zorder=1,
    )
    kmax = float(np.abs(head_curv["kappa_1pmm"]).max()) or 1.0
    for pos, curv, kp, mk in (
        (mig_pos, mig_curv, migration, "s"),
        (head_pos, head_curv, headline, "o"),
    ):
        p = np.column_stack([pos[f"{'xyz'[axis]}_mm"], pos["z_mm"]])
        segs = np.stack([p[:-1], p[1:]], axis=1)
        lc = LineCollection(
            segs, cmap="viridis", norm=plt.Normalize(0.0, kmax), zorder=2
        )
        lc.set_array(np.abs(curv["kappa_1pmm"].to_numpy()))
        lc.set_linewidth(4.0)
        ax.add_collection(lc)
        j = int(np.abs(pos["s0_mm"].to_numpy() - kp.s_kick_mm).argmin())
        ax.plot(
            p[j, 0], p[j, 1], mk, color="C3", ms=13, mfc="none", mew=2.5,
            zorder=3,
            label=f"theta~{kp.theta_deg:.0f} deg  kick f={kp.s_kick_frac:.2f} "
            f"({kp.bucket})",
        )
    ax.plot(
        head_pos[f"{'xyz'[axis]}0_mm"].iloc[0], head_pos["z0_mm"].iloc[0],
        "s", color="k", ms=9, zorder=4,
    )
    fig.colorbar(lc, ax=ax, label="|kappa| (1/mm)")
    ax.set_aspect("equal")
    ax.set_xlabel(f"span {long_axis} (mm)")
    ax.set_ylabel("rise z (mm)")
    ax.set_title("deformed mid-chord centreline, kick point ringed")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "kick_shape.svg", format="svg")
    plt.close(fig)


def write_flex_kick_json(fk: FlexKick, path: str | Path) -> None:
    """Serialise a :class:`FlexKick` (dataclasses -> plain dict) to JSON."""
    Path(path).write_text(json.dumps(fk.as_dict(), indent=2, sort_keys=True))
