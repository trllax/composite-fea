#!/usr/bin/env python3
"""Build U-clamp decks at tip tangent θ and solve with compfea.run."""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

from compfea.run import solve

CASE = Path(__file__).resolve().parent
ROOT = CASE.parents[1]
BASE = ROOT / "cases" / "cantilever_ansys" / "cantilever_89deg.inp"
RUNS = CASE / "runs"

# Mesh geometry from cantilever_89deg (y = long axis, L = 100 mm).
L_MM = 100.0
FAR_FACE = (49, 50, 51)  # x = 0, 10, 20 at y = L


def tip_targets(theta_rad: float) -> dict[int, tuple[float, float, float]]:
    """Map far_face node -> target (x,y,z) for circular-arc tip tangent θ."""
    if theta_rad <= 0:
        raise ValueError("theta must be > 0")
    r = L_MM / theta_rad
    y_t = r * math.sin(theta_rad)
    z_t = r * (1.0 - math.cos(theta_rad))
    # Undeformed tip nodes keep their x; edge translates as a rigid width line.
    undeformed_x = {49: 0.0, 50: 10.0, 51: 20.0}
    return {n: (undeformed_x[n], y_t, z_t) for n in FAR_FACE}


def undeformed_xyz() -> dict[int, tuple[float, float, float]]:
    text = BASE.read_text()
    nodes = {}
    for line in text.splitlines():
        m = re.match(
            r"^(\d+)\s*,\s*([-+eE0-9.]+)\s*,\s*([-+eE0-9.]+)\s*,\s*([-+eE0-9.]+)\s*$",
            line.strip(),
        )
        if m:
            nodes[int(m.group(1))] = tuple(float(m.group(i)) for i in range(2, 5))
    return {n: nodes[n] for n in FAR_FACE}


def build_deck(theta_deg: float) -> Path:
    theta = math.radians(theta_deg)
    targets = tip_targets(theta)
    origin = undeformed_xyz()
    base = BASE.read_text()
    # Drop the old step (rotation drive + prints) and append U-clamp step.
    head, _, _ = base.partition("*STEP")
    lines = [
        head.rstrip(),
        "** tip_clamp_u_drive: prescribed U on far_face for tip tangent "
        f"theta={theta_deg:g} deg (circular arc pose)",
        f"** R = L/theta = {L_MM / theta:.6f} mm",
        "*STEP, NLGEOM, INC=5000",
        "*STATIC",
        "0.01, 1.0, 1.E-6, 0.02",
        "*BOUNDARY",
    ]
    for n in FAR_FACE:
        x0, y0, z0 = origin[n]
        xt, yt, zt = targets[n]
        ux, uy, uz = xt - x0, yt - y0, zt - z0
        lines.append(f"{n}, 1, 1, {ux:.10f}")
        lines.append(f"{n}, 2, 2, {uy:.10f}")
        lines.append(f"{n}, 3, 3, {uz:.10f}")
    lines += [
        "*NODE PRINT, NSET=far_face, TOTALS=YES",
        "RF",
        "*NODE PRINT, NSET=far_face",
        "U",
        "*NODE PRINT, NSET=fixed_end, TOTALS=YES",
        "RF",
        "*END STEP",
        "",
    ]
    out = CASE / f"clamp_{int(theta_deg)}deg.inp"
    out.write_text("\n".join(lines))
    return out


def main(argv: list[str]) -> int:
    angles = [float(a) for a in (argv[1:] or ["90", "180"])]
    RUNS.mkdir(parents=True, exist_ok=True)
    print(f"base deck: {BASE}")
    for deg in angles:
        deck = build_deck(deg)
        run_dir = RUNS / f"theta_{int(deg)}"
        print(f"\n=== θ = {deg:g}°  deck={deck.name} ===")
        try:
            result = solve(deck, run_dir, job_name="job", timeout_s=600.0)
        except Exception as exc:
            print(f"FAIL: {type(exc).__name__}: {exc}")
            continue
        tip = result.final[result.final["nset"] == "far_face"].iloc[0]
        root = result.final[result.final["nset"] == "fixed_end"].iloc[0]
        print(
            f"OK increments={result.increments} wall={result.wall_time_s:.2f}s\n"
            f"  tip  RF fx={tip.fx:.4f} fy={tip.fy:.4f} fz={tip.fz:.4f} N\n"
            f"  root RF fx={root.fx:.4f} fy={root.fy:.4f} fz={root.fz:.4f} N\n"
            f"  |F_tip|={math.sqrt(tip.fx**2 + tip.fy**2 + tip.fz**2):.4f} N"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
