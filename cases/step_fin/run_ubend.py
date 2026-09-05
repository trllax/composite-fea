#!/usr/bin/env python3
"""Tip U-clamp path on test_fin_2: HEAL mask fixed, tip edge driven."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from compfea.deck import assemble
from compfea.geometry import Mesh
from compfea.layup import (
    Ply,
    coverages_from_mesh,
    layup_from_coverage,
    mesh_elsets_for_stacks,
)
from compfea.run import parse_dat_energy, solve
from compfea.step_mesh import mesh_step
from compfea.ubend import (
    build_deck,
    final_time_for,
    force_at_theta,
    step_index_for,
    tip_length_mm,
    theta_grid_deg,
)

ROOT = Path(__file__).resolve().parents[2]
STEP = ROOT / "test_fin_2.step"
LONG_AXIS = "x"

# The fin needs its own increment settings; the strip's calibrated default
# (max 0.25) diverges here. Measured on this STEP at 40 mm, 1-degree steps,
# varying only the max increment with the minimum held at 1.E-10:
#
#   max 0.25 -> diverges at ~88% of the first step
#   max 0.1  -> converges, 21.5 s, 44 increments   <- adopted
#   max 0.01 -> converges, 51.8 s, 211 increments  (the old hardcoded value)
#
# So there is a 2.4x saving here, not the strip's 5.4x. The wall is in step 1,
# where a flat blade takes its first bend, and it is not caused by the S6
# elements: none of the 46 diverging nodes belongs to a triangle, and the
# all-quad 16 mm mesh diverges at the same point. The minimum increment matters
# as much as the maximum -- 1.E-10 is what lets ccx grind through that step.
FIN_STATIC_LINE = "0.001, 1.0, 1.E-10, 0.1"


def default_plies():
    return [
        Ply(0.15, 0.0, coverage="FULL"),
        Ply(0.15, 90.0, coverage="FULL"),
        Ply(0.10, 0.0, coverage="z_3_4ths"),
        Ply(0.10, 90.0, coverage="HALF"),
        Ply(0.20, 0.0, coverage="TIP"),
        Ply(0.25, 45.0, coverage="HEAL"),
    ]


def build_mesh_and_layup(size_mm: float):
    raw = mesh_step(STEP, size_mm=size_mm)
    cov = coverages_from_mesh(raw.elsets)
    layup, stacks = layup_from_coverage(default_plies(), cov, long_axis=LONG_AXIS)
    mesh = Mesh(
        nodes=raw.nodes,
        elements=raw.elements,
        nsets=raw.nsets,
        elsets=mesh_elsets_for_stacks(stacks, all_elements=raw.elsets["blade"]),
        heading=raw.heading,
    )
    return mesh, layup


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--end-deg", type=float, default=180.0)
    p.add_argument("--start-deg", type=float, default=1.0)
    p.add_argument("--step-deg", type=float, default=1.0)
    # 40 mm leaves triangles on one tile and mesh_step refuses that now.
    # 17.0 is the coarsest size that recombines to all quad8 on this STEP,
    # but 17.25 already fails; 16.0 keeps margin to that cliff.
    p.add_argument("--size-mm", type=float, default=16.0)
    p.add_argument("--timeout-s", type=float, default=3600.0)
    p.add_argument("--static-line", default=FIN_STATIC_LINE)
    p.add_argument("--threads", type=int, default=1,
                   help="OMP_NUM_THREADS for this solve (default 1)")
    p.add_argument(
        "--run-dir",
        type=Path,
        default=ROOT / "results" / "fin_ubend",
    )
    args = p.parse_args(argv)

    if not STEP.is_file():
        print(f"missing {STEP}", file=sys.stderr)
        return 2

    mesh, layup = build_mesh_and_layup(args.size_mm)
    arm = tip_length_mm(mesh, long_axis=LONG_AXIS)
    angles = theta_grid_deg(step_deg=args.step_deg, start_deg=args.start_deg, end_deg=args.end_deg)
    deck = build_deck(
        mesh,
        layup,
        angles,
        long_axis=LONG_AXIS,
        static_line=args.static_line,
        heading=(
            f"fin U-bend: HEAL clamp, tip edge U, L={arm:.3f} mm, "
            f"{angles[0]:g}→{angles[-1]:g} deg"
        ),
    )
    args.run_dir.mkdir(parents=True, exist_ok=True)
    deck_path = args.run_dir / "deck.inp"
    deck_path.write_text(deck)
    print(f"deck {deck_path}  elems={len(mesh.elements)}  arm={arm:.2f} mm  steps={len(angles)}")

    result = solve(
        deck_path,
        args.run_dir / "ccx",
        job_name="job",
        timeout_s=args.timeout_s,
        final_time=final_time_for(angles),
        threads=args.threads,
    )
    energy = parse_dat_energy(args.run_dir / "ccx" / "job.dat")
    print(f"OK wall={result.wall_time_s:.1f}s increments={result.increments}")
    for deg in (90.0, 180.0):
        if deg > args.end_deg + 1e-9:
            continue
        idx = step_index_for(angles, deg)
        u, m, f = force_at_theta(
            energy, theta_deg=deg, step_index=idx, arm_mm=arm
        )
        print(f"  θ={deg:g}°  U={u:.4g} N.mm  M={m:.4g} N.mm  F={f:.4g} N")
    if args.end_deg >= 180.0:
        # ratio if both present
        try:
            _, _, f90 = force_at_theta(
                energy, theta_deg=90.0, step_index=step_index_for(angles, 90.0), arm_mm=arm
            )
            _, _, f180 = force_at_theta(
                energy, theta_deg=180.0, step_index=step_index_for(angles, 180.0), arm_mm=arm
            )
            print(f"  F_180/F_90 = {f180/f90:.3f}")
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
