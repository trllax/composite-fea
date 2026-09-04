"""The fast gate: cases/smoke_cantilever, solved with ccx.

CLAUDE.md requires this case to pass after any change to deck.py, layup.py or
geometry.py. tests/test_layup_deck.py proves the deck says what it was asked to
say; this proves ccx reads it as the physics requires, against references that
share no code with compfea: classical lamination theory (clpt.py) in the linear
range, and the closed-form elastica (elastica.py) at 26 degrees of tip rotation.

Like tests/test_run.py there is deliberately no skipif on a missing ccx. A gate
that quietly skips is not a gate.

Tolerances come from a measured mesh study, not from what would pass. See
cases/smoke_cantilever/README.md.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import pytest

from compfea.run import solve

CASE = Path(__file__).resolve().parents[1] / "cases" / "smoke_cantilever"


def _case_module(name: str):
    """Import a module out of the case directory without installing it."""
    spec = importlib.util.spec_from_file_location(f"smoke_{name}", CASE / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


if str(CASE) not in sys.path:
    sys.path.insert(0, str(CASE))
build_deck = _case_module("build_deck")
clpt = _case_module("clpt")
elastica = _case_module("elastica")

MATERIAL = dict(
    e1=build_deck.PLACEHOLDER_CFRP.e1,
    e2=build_deck.PLACEHOLDER_CFRP.e2,
    nu12=build_deck.PLACEHOLDER_CFRP.nu12,
    g12=build_deck.PLACEHOLDER_CFRP.g12,
)
LENGTH = build_deck.LENGTH_MM
WIDTH = build_deck.WIDTH_MM
PLY = build_deck.PLY_MM

LINEAR_DELTA_MM = 1.0  # delta/L = 1%, small-deflection range
LARGE_DELTA_MM = 30.0  # delta/L = 30%, 26.3 deg of tip rotation

# Hand-computed and written down so a change in clpt.py cannot move the target
# quietly. EI is the free-edge narrow-strip reduction b(D11 - D12^2/D22), which
# is what this strip is; D11*b would be 0.247% stiffer. Derivation in README.md.
PINNED = {
    "0_90_s": {"d11": 9997.4849, "ei": 199455.826, "p_1mm": 0.59836748},
    "90_0_s": {"d11": 2074.9497, "ei": 41396.492, "p_1mm": 0.12418948},
    "0_0_90_90": {"d11": 6036.2173, "ei": 120554.577, "b11_over_d11": +0.21875},
    "90_90_0_0": {"d11": 6036.2173, "ei": 120554.577, "b11_over_d11": -0.21875},
}

# Measured on the shipped 1x32 mesh: +0.00% and +0.24% against CLPT, -0.03%
# against the elastica. The residual is discretization, not a solver offset: on
# a 1x8 mesh [90/0]s is +0.97% and it falls monotonically with refinement. See
# the mesh study in README.md.
CLPT_TOL = 0.005
ELASTICA_TOL = 0.003
# The two references for the ply-order check are exact statements about D and B,
# so they get tighter bands: D is invariant under reversal to machine precision
# in theory and to the solver's repeatability in practice, while the draw-in
# prediction is small-deflection beam theory and carries its idealisation error.
REVERSAL_FORCE_TOL = 0.001
REVERSAL_DRAW_IN_TOL = 0.01

_DISPLACEMENTS = re.compile(
    r"^\s*displacements \(vx,vy,vz\) for set (?P<nset>\S+) and time\s+(?P<time>\S+)\s*$"
)


def tip_displacement(result, nset: str = "FAR_FACE") -> tuple[float, float, float]:
    """Mean (ux, uy, uz) over ``nset`` at the last time in the .dat.

    run.py parses reactions and deliberately stops there. The draw-in is read
    here, in the case that needs it, from the same converged .dat -- solve() has
    already refused anything that did not reach the final time, so there is no
    path from a partial solve to a number.
    """
    lines = (result.job_dir / f"{result.job_name}.dat").read_text().splitlines()
    blocks: dict[float, list[tuple[float, float, float]]] = {}
    current: float | None = None
    for line in lines:
        header = _DISPLACEMENTS.match(line)
        if header:
            current = float(header["time"]) if header["nset"] == nset else None
            if current is not None:
                blocks[current] = []
            continue
        if current is None or not line.strip():
            continue
        fields = line.split()
        if len(fields) == 4 and fields[0].isdigit():
            blocks[current].append(tuple(float(v) for v in fields[1:]))
        else:
            current = None
    assert blocks, f"no displacement block for {nset} in {result.job_name}.dat"
    rows = blocks[max(blocks)]
    return tuple(sum(v[i] for v in rows) / len(rows) for i in range(3))


@pytest.fixture(scope="module")
def solved(tmp_path_factory):
    """Every deck this module needs, solved once. Under two seconds in total."""
    root = tmp_path_factory.mktemp("smoke_cantilever")
    results = {}
    for stack, delta in (
        ("0_90_s", LINEAR_DELTA_MM),
        ("90_0_s", LINEAR_DELTA_MM),
        ("0_0_90_90", LINEAR_DELTA_MM),
        ("90_90_0_0", LINEAR_DELTA_MM),
        ("0_90_s", LARGE_DELTA_MM),
    ):
        run_dir = root / f"{stack}_{delta:g}mm"
        run_dir.mkdir()
        deck = run_dir / "deck.inp"
        deck.write_text(build_deck.build(stack, delta))
        results[(stack, delta)] = solve(deck, run_dir)
    return results


def _root_reaction(solved, stack: str, delta: float) -> float:
    """Transverse reaction at the clamped end, N, sign dropped."""
    return abs(solved[(stack, delta)].total_force("fixed_end", "fz"))


def _ei(stack: str) -> float:
    """Narrow-strip EI, cross-checked against the pinned value."""
    stack_angles = build_deck.STACKS[stack]
    assert clpt.d11(stack_angles, PLY, **MATERIAL) == pytest.approx(
        PINNED[stack]["d11"], rel=1e-8
    )
    ei = clpt.ei_narrow_strip(stack_angles, PLY, width=WIDTH, **MATERIAL)
    assert ei == pytest.approx(PINNED[stack]["ei"], rel=1e-8)
    return ei


@pytest.mark.parametrize("stack", ["0_90_s", "90_0_s"])
def test_linear_stiffness_matches_hand_clpt(solved, stack):
    """delta/L = 1%: the root reaction is 3 EI delta / L^3 from lamination theory."""
    expected = clpt.tip_load_small_deflection(
        LINEAR_DELTA_MM, ei=_ei(stack), length=LENGTH
    )
    assert expected == pytest.approx(PINNED[stack]["p_1mm"], rel=1e-7)
    measured = _root_reaction(solved, stack, LINEAR_DELTA_MM)
    assert measured == pytest.approx(expected, rel=CLPT_TOL)


def test_which_ply_sits_outside_reaches_the_solver(solved):
    """[0/90]s and [90/0]s are the same plies at swapped z stations.

    Same ply count, thickness, mass and material: everything a deck diff would
    compare is equal, and the bending stiffness differs 4.8x because the stiff
    plies move from the outside to the middle. Any bug that loses which ply sits
    where -- sorting the stack, deduplicating it, hanging the wrong section on
    the wrong ELSET -- collapses this ratio toward 1.

    It also pins the long axis: built with long_axis="x" the 0-degree plies point
    across the width and each stack returns the other one's stiffness.

    What it cannot see is a reversal; that is test_reversing_the_stack below.
    """
    stiff = _root_reaction(solved, "0_90_s", LINEAR_DELTA_MM)
    soft = _root_reaction(solved, "90_0_s", LINEAR_DELTA_MM)
    expected = PINNED["0_90_s"]["ei"] / PINNED["90_0_s"]["ei"]
    assert expected == pytest.approx(4.8182, rel=1e-4)
    assert stiff / soft == pytest.approx(expected, rel=CLPT_TOL)
    # Blunt version of the same statement, immune to any tolerance argument.
    assert stiff / soft > 4.0


def test_reversing_the_stack_is_invisible_in_the_force(solved):
    """D is exactly invariant under reversal, so the two must agree.

    Ply k spans [z_{k-1}, z_k] and its mirror spans [-z_k, -z_{k-1}]; both
    contribute the same (z^3 - z^3)/3 to D. This is why no stiffness check can
    ever catch a reversed stack, and why the next test exists. Asserted rather
    than assumed: if these two ever disagree, the deck is not what it claims.
    """
    a = _root_reaction(solved, "0_0_90_90", LINEAR_DELTA_MM)
    b = _root_reaction(solved, "90_90_0_0", LINEAR_DELTA_MM)
    assert a == pytest.approx(b, rel=REVERSAL_FORCE_TOL)


def test_reversing_the_stack_flips_the_draw_in(solved):
    """The first ply line is the -z ply, checked through the solver.

    B changes sign when a stack is reversed and D does not, so [0/0/90/90] and
    [90/90/0/0] bend identically and draw in differently. Under pure bending of a
    narrow strip the mid-plane axial strain is (b11/d11) * kappa, so integrating
    along the beam the tip's axial position moves by (b11/d11) * theta, and the
    two stacks differ by twice that.

    The sign is not fitted to ccx. Pushing the tip to +z puts the -z face in
    tension; in [0/0/90/90] the stiff 0-degree plies are the ones at -z, so
    holding the axial resultant at zero drives the mid-plane into extra
    compression and that stack must draw in *more* than its reverse. In the plate
    curvature convention (kappa = -w'') that is du_y = -(b11/d11) * theta.

    Nothing else in the repo checks this. Verified by sabotage: emitting the
    stack top-to-bottom leaves every other test in this module green.
    """
    ratio = clpt.axial_coupling_ratio(build_deck.STACKS["0_0_90_90"], PLY, **MATERIAL)
    assert ratio == pytest.approx(PINNED["0_0_90_90"]["b11_over_d11"], rel=1e-6)
    assert clpt.axial_coupling_ratio(
        build_deck.STACKS["90_90_0_0"], PLY, **MATERIAL
    ) == pytest.approx(-ratio, rel=1e-9)
    # A symmetric stack has no coupling at all -- the control on the whole idea.
    assert clpt.axial_coupling_ratio(
        build_deck.STACKS["0_90_s"], PLY, **MATERIAL
    ) == pytest.approx(0.0, abs=1e-15)

    theta = clpt.tip_angle_small_deflection(LINEAR_DELTA_MM, length=LENGTH)
    expected = -2.0 * ratio * theta
    stiff_at_bottom = tip_displacement(solved[("0_0_90_90", LINEAR_DELTA_MM)])[1]
    stiff_at_top = tip_displacement(solved[("90_90_0_0", LINEAR_DELTA_MM)])[1]
    assert stiff_at_bottom - stiff_at_top == pytest.approx(
        expected, rel=REVERSAL_DRAW_IN_TOL
    )
    # Both still draw in: the coupling modulates the geometric shortening.
    assert stiff_at_bottom < stiff_at_top < 0.0


def test_large_deflection_matches_the_closed_form_elastica(solved):
    """delta/L = 30%, 26.3 deg of tip rotation, against the elliptic-integral form."""
    expected = elastica.tip_load(
        LARGE_DELTA_MM / LENGTH, ei=_ei("0_90_s"), length=LENGTH
    )
    measured = _root_reaction(solved, "0_90_s", LARGE_DELTA_MM)
    assert measured == pytest.approx(expected, rel=ELASTICA_TOL)


def test_nlgeom_is_actually_on(solved):
    """The large-deflection point must not agree with small-deflection theory.

    ccx accepts a deck whose *STEP lost NLGEOM without complaint, and the answer
    stays plausible -- it just stops being large-deflection. At delta/L = 30% the
    linear result is 10.3% low, far outside anything discretization can explain.
    """
    linear = clpt.tip_load_small_deflection(
        LARGE_DELTA_MM, ei=_ei("0_90_s"), length=LENGTH
    )
    measured = _root_reaction(solved, "0_90_s", LARGE_DELTA_MM)
    assert measured / linear > 1.05
    assert "NLGEOM" in build_deck.build("0_90_s", LARGE_DELTA_MM)


def test_the_case_is_fast(solved):
    """CLAUDE.md calls this the fast gate. Keep it one."""
    total = sum(r.wall_time_s for r in solved.values())
    assert total < 5.0, f"smoke_cantilever took {total:.2f} s of solver time"
