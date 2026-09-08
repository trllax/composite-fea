#!/usr/bin/env python3
"""Tip U-clamp on a plybook-built test_fin_2 (build → ubend, does not call solve in build)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from compfea.build import build
from compfea.geometry import Mesh
from compfea.layup import coverages_from_mesh, mesh_elsets_for_stacks
from compfea.materials import load_materials
from compfea.plybook import load_plybook, resolve
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
FIN_STATIC_LINE = "0.001, 1.0, 1.E-10, 0.1"


def mesh_and_layup(
    *,
    plybook: Path,
    materials: Path,
    size_mm: float,
    out: Path,
):
    """Same resolve path as compfea-build; also writes reports under ``out``."""
    build(
        step=STEP,
        plybook=plybook,
        materials=materials,
        long_axis=LONG_AXIS,
        out=out,
        size_mm=size_mm,
        clamp_coverage="HEAL",
    )
    # rebuild live objects for the solve (build only returns a manifest)
    library = load_materials(materials)
    rows = load_plybook(plybook)
    raw = mesh_step(STEP, size_mm=size_mm, clamp_coverage="HEAL", long_axis=LONG_AXIS)
    coverages = coverages_from_mesh(raw.elsets)
    layup, stacks = resolve(
        rows, library, coverages, long_axis=LONG_AXIS, elsets=raw.elsets
    )
    mesh = Mesh(
        nodes=raw.nodes,
        elements=raw.elements,
        nsets=raw.nsets,
        elsets=mesh_elsets_for_stacks(stacks, all_elements=tuple(raw.elements)),
        heading=raw.heading,
    )
    return mesh, layup


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--plybook",
        type=Path,
        default=ROOT / "cases/fin_zoned/plybook_shop.csv",
    )
    p.add_argument(
        "--materials",
        type=Path,
        default=ROOT / "materials/from_shop.csv",
    )
    p.add_argument("--end-deg", type=float, default=90.0)
    p.add_argument("--start-deg", type=float, default=1.0)
    p.add_argument("--step-deg", type=float, default=1.0)
    p.add_argument("--size-mm", type=float, default=40.0,
                   help="40 mm is faster for layup hunting; confirm winners at 16")
    p.add_argument("--timeout-s", type=float, default=18000.0)
    p.add_argument("--threads", type=int, default=5)
    p.add_argument("--static-line", default=FIN_STATIC_LINE)
    p.add_argument(
        "--run-dir",
        type=Path,
        default=ROOT / "results" / "fin_shop_demo",
    )
    args = p.parse_args(argv)

    if not STEP.is_file():
        print(f"missing {STEP}", file=sys.stderr)
        return 2

    args.run_dir.mkdir(parents=True, exist_ok=True)
    mesh, layup = mesh_and_layup(
        plybook=args.plybook,
        materials=args.materials,
        size_mm=args.size_mm,
        out=args.run_dir / "build",
    )
    arm = tip_length_mm(mesh, long_axis=LONG_AXIS)
    angles = theta_grid_deg(
        step_deg=args.step_deg, start_deg=args.start_deg, end_deg=args.end_deg
    )
    deck = build_deck(
        mesh,
        layup,
        angles,
        long_axis=LONG_AXIS,
        static_line=args.static_line,
        file_deg=(90.0,) if args.end_deg >= 90 else (),
        heading=(
            f"fin_zoned U-bend plybook={args.plybook.name} "
            f"L={arm:.3f} mm {angles[0]:g}→{angles[-1]:g} deg"
        ),
    )
    deck_path = args.run_dir / "deck.inp"
    deck_path.write_text(deck)
    print(
        f"deck {deck_path} elems={len(mesh.elements)} arm={arm:.2f} mm "
        f"steps={len(angles)} mats={[m.name for m in layup.materials]}"
    )

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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
