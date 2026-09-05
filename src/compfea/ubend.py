"""Strip U-bend: tip-edge prescribed U on a circular-arc path.

Tip tangent θ from 0 to 180° maps to an undeformed-length circular arc.
Bench force uses energy, not tip |RF|::

    M = 2U / θ
    F = M / arm

See ``compfea.metrics.f_spring`` and ``cases/u_bend_path/README.md``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from compfea.deck import StaticStep, assemble, tip_u_clamp_body
from compfea.geometry import Mesh
from compfea.layup import Layup
from compfea.metrics import f_spring

# Default path that reached 180° on the placeholder-CFRP strip.
DEFAULT_STEP_DEG = 5.0
DEFAULT_START_DEG = 5.0
DEFAULT_END_DEG = 180.0
DEFAULT_STATIC_LINE = "0.005, 1.0, 1.E-8, 0.02"
DEFAULT_INC = 8000


def tip_length_mm(
    mesh: Mesh,
    *,
    tip_nset: str = "far_face",
    root_nset: str = "fixed_end",
    long_axis: str = "y",
) -> float:
    """Free length: tip station minus outboard clamp station along ``long_axis``.

    For a bonded clamp patch (e.g. HEAL mask), the arm is tip minus *max*
    root-set coordinate along the span, not the inboard edge.
    """
    if long_axis not in ("x", "y"):
        raise ValueError(f"long_axis must be 'x' or 'y', not {long_axis!r}")
    axis = 0 if long_axis == "x" else 1
    tip = [mesh.nodes[n][axis] for n in mesh.nsets[tip_nset]]
    root = [mesh.nodes[n][axis] for n in mesh.nsets[root_nset]]
    return max(tip) - max(root)


def tip_displacements(
    mesh: Mesh,
    theta_rad: float,
    *,
    tip_nset: str = "far_face",
    length_mm: float | None = None,
    long_axis: str = "y",
) -> dict[int, tuple[float, float, float]]:
    """Prescribed U for each tip node at tip-tangent ``theta_rad``.

    Tip edge translates as a rigid line on a circular arc of radius ``L / θ``
    in the (long_axis, z) plane. The transverse in-plane coordinate is held.
    """
    if theta_rad <= 0:
        raise ValueError("theta must be > 0")
    if long_axis not in ("x", "y"):
        raise ValueError(f"long_axis must be 'x' or 'y', not {long_axis!r}")
    length = (
        tip_length_mm(mesh, tip_nset=tip_nset, long_axis=long_axis)
        if length_mm is None
        else length_mm
    )
    if length <= 0:
        raise ValueError(f"length must be > 0, got {length}")
    r = length / theta_rad
    s_t = r * math.sin(theta_rad)  # along undeformed long axis from clamp edge
    z_lift = r * (1.0 - math.cos(theta_rad))
    axis = 0 if long_axis == "x" else 1
    # Clamp outboard station; tip target along axis = clamp + s_t
    root_axis = [mesh.nodes[n][axis] for n in mesh.nsets["fixed_end"]]
    s0 = max(root_axis)
    target_s = s0 + s_t
    out: dict[int, tuple[float, float, float]] = {}
    for nid in mesh.nsets[tip_nset]:
        x0, y0, z0 = mesh.nodes[nid]
        if long_axis == "y":
            # keep x; drive y,z (strip convention)
            out[nid] = (0.0, target_s - y0, z0 + z_lift - z0)
        else:
            # keep y; drive x,z (fin span along +x)
            out[nid] = (target_s - x0, 0.0, z0 + z_lift - z0)
    return out


def theta_grid_deg(
    *,
    step_deg: float = DEFAULT_STEP_DEG,
    start_deg: float = DEFAULT_START_DEG,
    end_deg: float = DEFAULT_END_DEG,
) -> list[float]:
    """Inclusive multi-step angles from ``start_deg`` to ``end_deg``."""
    if step_deg <= 0:
        raise ValueError("step_deg must be > 0")
    if end_deg < start_deg:
        raise ValueError("end_deg must be >= start_deg")
    n = int(round((end_deg - start_deg) / step_deg))
    angles = [start_deg + i * step_deg for i in range(n + 1)]
    angles[-1] = end_deg
    return angles


def build_deck(
    mesh: Mesh,
    layup: Layup,
    angles_deg: Sequence[float],
    *,
    long_axis: str = "y",
    heading: str = "",
    static_line: str = DEFAULT_STATIC_LINE,
    inc: int = DEFAULT_INC,
    file_deg: Sequence[float] | None = (90.0, 180.0),
) -> str:
    """Multi-step NLGEOM deck: root clamp + tip U at each θ.

    ``file_deg`` angles get ``*NODE FILE`` / U (DISP in the ``.frd``). Default
    is 90° and 180° only so FRDs stay small; pass ``()`` to disable.
    """
    if not angles_deg:
        raise ValueError("angles_deg must not be empty")
    file_set = {float(d) for d in (file_deg or ())}
    steps = [
        StaticStep(
            tip_u_clamp_body(
                tip_displacements(
                    mesh, math.radians(float(deg)), long_axis=long_axis
                ),
                node_file=any(abs(float(deg) - d) <= 1e-9 for d in file_set),
            ),
            inc=inc,
            static_line=static_line,
        )
        for deg in angles_deg
    ]
    default_heading = (
        f"U-bend tip-U path: {len(angles_deg)} steps, "
        f"{float(angles_deg[0]):g}→{float(angles_deg[-1]):g} deg; "
        "metric F=M/L from ELSE"
    )
    return assemble(
        mesh_inp=mesh.to_inp(),
        layup=layup,
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=steps,
        heading=heading or default_heading,
    )


def final_time_for(angles_deg: Sequence[float]) -> float:
    """ccx TOT TIME after N unit-period static steps."""
    return float(len(angles_deg))


def force_at_theta(
    energy_rows,
    *,
    theta_deg: float,
    step_index: int,
    arm_mm: float,
    elset: str = "blade",
) -> tuple[float, float, float]:
    """Return (U, M, F) at the end of step ``step_index`` (1-based).

    ``step_index`` matches ccx TOT TIME when each *STATIC period is 1.0.
    """
    target_time = float(step_index)
    rows = energy_rows[
        (energy_rows["elset"] == elset.lower())
        & ((energy_rows["time"] - target_time).abs() <= 1e-6)
    ]
    if rows.empty:
        raise LookupError(
            f"no ELSE energy for elset={elset!r} at time={target_time:g} "
            f"(θ={theta_deg:g}°)"
        )
    u = float(rows["energy"].iloc[-1])
    m, f = f_spring(u, math.radians(theta_deg), arm_mm)
    return u, m, f


def step_index_for(angles_deg: Sequence[float], theta_deg: float) -> int:
    """1-based step index whose angle equals ``theta_deg``."""
    for i, deg in enumerate(angles_deg, start=1):
        if abs(float(deg) - theta_deg) <= 1e-9:
            return i
    raise LookupError(f"θ={theta_deg:g}° is not in the angle grid")
