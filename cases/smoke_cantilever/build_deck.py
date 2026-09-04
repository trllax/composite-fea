#!/usr/bin/env python3
"""Build the smoke_cantilever decks through layup.py and deck.py.

The case exists to gate those two modules, so the deck is assembled by them
rather than checked in: a hand-written .inp would still pass after a refactor
scrambled the ply order, which is the failure this case is here to catch.

Units: mm, N, tonne, MPa, s. See cases/smoke_cantilever/README.md for the
pinned numbers and where they come from.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from compfea.deck import StaticStep, assemble
from compfea.layup import PLACEHOLDER_CFRP, Layup, Ply

CASE = Path(__file__).resolve().parent
MESH = CASE / "mesh_strip_32el.inp"

LENGTH_MM = 100.0
WIDTH_MM = 20.0
PLY_MM = 0.25

# x is the width and y is the long axis in this mesh, so a 0-degree ply runs
# along y. Stated here and in the README because layup.py has no default.
LONG_AXIS = "y"

# Bottom (-z) ply first, always. Four stacks, four 0.25 mm plies each, so every
# one has the same mass and the same thickness and only the z stations differ.
#
#   0_90_s / 90_0_s   symmetric, B = 0. Same plies at swapped z stations, 4.8x
#                     apart in bending stiffness. These carry the CLPT checks.
#   0_0_90_90 and     unsymmetric and each other's reverse, so they have exactly
#   90_90_0_0         the same D and opposite B: identical stiffness, opposite
#                     bend-extension coupling. Reversal is invisible in the
#                     reaction force and visible in the tip draw-in, which is the
#                     only way to pin "first ply line = -z" through the solver.
STACKS: dict[str, tuple[float, ...]] = {
    "0_90_s": (0.0, 90.0, 90.0, 0.0),
    "90_0_s": (90.0, 0.0, 0.0, 90.0),
    "0_0_90_90": (0.0, 0.0, 90.0, 90.0),
    "90_90_0_0": (90.0, 90.0, 0.0, 0.0),
}


def layup_for(stack: str) -> Layup:
    if stack not in STACKS:
        raise KeyError(f"stack must be one of {sorted(STACKS)}, not {stack!r}")
    return Layup.uniform(
        [Ply(PLY_MM, angle) for angle in STACKS[stack]],
        long_axis=LONG_AXIS,
        elset="blade",
        material=PLACEHOLDER_CFRP,
    )


def build(stack: str, delta_mm: float) -> str:
    """Deck text: tip edge driven to ``delta_mm`` in +z, free to draw in.

    Displacement control per CLAUDE.md. DOF 1 and 2 are left free at the tip so
    the strip shortens its projection as it bends, which makes the reaction a
    pure transverse force and the closed-form elastica the right reference. The
    root is clamped in all six DOF and is where the reaction is read.
    """
    body = "\n".join(
        [
            "*BOUNDARY",
            f"far_face, 3, 3, {delta_mm:.10f}",
            "*NODE PRINT, NSET=fixed_end, TOTALS=YES",
            "RF",
            # Per-node, not TOTALS: the tip draw-in is what shows bend-extension
            # coupling, and summing it over the edge would only scale it.
            "*NODE PRINT, NSET=far_face",
            "U",
        ]
    )
    # Increment size scaled to the drive: 1 mm is a linear-range nudge, 30 mm is
    # 26 degrees of tip rotation and needs the path resolved.
    coarse = delta_mm <= 0.05 * LENGTH_MM
    step = StaticStep(
        body,
        inc=2000,
        static_line="0.25, 1.0, 1.E-6, 1.0" if coarse else "0.05, 1.0, 1.E-6, 0.1",
    )
    return assemble(
        mesh_inp=MESH.read_text(),
        layup=layup_for(stack),
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=[step],
        heading=(
            f"smoke_cantilever: [{stack}] strip, tip edge driven to "
            f"uz = {delta_mm:g} mm (delta/L = {delta_mm / LENGTH_MM:g})\n"
            f"Built by cases/smoke_cantilever/build_deck.py -- do not edit by hand."
        ),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack", default="0_90_s", choices=sorted(STACKS))
    parser.add_argument("--delta", type=float, default=1.0, help="tip uz, mm")
    parser.add_argument("--print", action="store_true", dest="show")
    parser.add_argument("--out", type=Path, help="write the deck here")
    args = parser.parse_args(argv)

    deck = build(args.stack, args.delta)
    if args.out:
        args.out.write_text(deck)
        print(f"wrote {args.out}")
    if args.show or not args.out:
        sys.stdout.write(deck)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
