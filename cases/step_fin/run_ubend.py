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
    p.add_argument("--size-mm", type=float, default=40.0)
    p.add_argument("--timeout-s", type=float, default=3600.0)
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
        static_line="0.001, 1.0, 1.E-10, 0.01",
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
