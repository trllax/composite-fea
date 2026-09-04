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
) -> float:
    """Undeformed tip station minus root station along y (long axis)."""
    tip_ys = [mesh.nodes[n][1] for n in mesh.nsets[tip_nset]]
    root_ys = [mesh.nodes[n][1] for n in mesh.nsets[root_nset]]
    return max(tip_ys) - min(root_ys)


def tip_displacements(
    mesh: Mesh,
    theta_rad: float,
    *,
    tip_nset: str = "far_face",
    length_mm: float | None = None,
) -> dict[int, tuple[float, float, float]]:
    """Prescribed U (ux, uy, uz) for each tip node at tip-tangent ``theta_rad``.

    The tip edge translates as a rigid width line on a circular arc of radius
    ``L / θ``. Each node keeps its undeformed x.
    """
    if theta_rad <= 0:
        raise ValueError("theta must be > 0")
    length = tip_length_mm(mesh) if length_mm is None else length_mm
    if length <= 0:
        raise ValueError(f"length must be > 0, got {length}")
    r = length / theta_rad
    y_t = r * math.sin(theta_rad)
    z_t = r * (1.0 - math.cos(theta_rad))
    out: dict[int, tuple[float, float, float]] = {}
    for nid in mesh.nsets[tip_nset]:
        x0, y0, z0 = mesh.nodes[nid]
        out[nid] = (0.0, y_t - y0, z_t - z0)
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
    heading: str = "",
    static_line: str = DEFAULT_STATIC_LINE,
    inc: int = DEFAULT_INC,
) -> str:
    """Multi-step NLGEOM deck: root clamp + tip U at each θ."""
    if not angles_deg:
        raise ValueError("angles_deg must not be empty")
    steps = [
        StaticStep(
            tip_u_clamp_body(tip_displacements(mesh, math.radians(float(deg)))),
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
