"""Unit tests for layup and deck assembly (no ccx).

These are string-level tests: they prove the deck says what it was asked to say.
Whether ccx then reads it the way the physics needs is the job of
`tests/test_smoke_cantilever.py`, which solves and checks against hand CLPT.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

from compfea.deck import StaticStep, assemble, tip_u_clamp_body
from compfea.layup import (
    EngineeringConstants,
    Layup,
    Ply,
    ZoneLayup,
    UD_CFRP_GENERIC,
    cross_ply_symmetric,
    orientation_card,
    orientation_name,
    woven_from_ud,
)
from compfea.metrics import f_spring

MAT = EngineeringConstants(
    135000, 9000, 9000, 0.3, 0.3, 0.45, 4500, 4500, 3000, name="t700"
)


def _cosines(card: str) -> str:
    """The direction-cosine line of an *ORIENTATION card, without its name."""
    return card.splitlines()[1]


def _ply_lines(text: str, elset: str) -> list[str]:
    """The ply lines of one *SHELL SECTION block, in the order they appear."""
    lines = text.splitlines()
    start = lines.index(f"*SHELL SECTION, COMPOSITE, ELSET={elset}") + 1
    out = []
    for line in lines[start:]:
        if line.startswith("*"):
            break
        out.append(line)
    return out


def test_orientation_names():
    assert orientation_name(0) == "ori_p0"
    assert orientation_name(-0.0) == "ori_p0"
    assert orientation_name(90) == "ori_p90"
    assert orientation_name(-45) == "ori_m45"
    assert orientation_name(22.5) == "ori_p22p5"


def test_angles_a_fraction_of_a_degree_apart_get_distinct_names():
    """Two angles must never collapse onto one *ORIENTATION name.

    A collision emits two cards with the same NAME; ccx keeps one and both plies
    silently take a fibre direction nobody asked for.
    """
    assert orientation_name(30.02) != orientation_name(30.04)


def test_one_orientation_card_per_distinct_angle():
    angles = (0.0, 45.0, -45.0, 22.5, 30.02, 30.04)
    layup = Layup.uniform(
        [Ply(0.1, a) for a in angles], long_axis="y", material=MAT
    )
    text = layup.to_inp()
    assert text.count("*ORIENTATION") == len(angles)
    names = [
        line.split("NAME=")[1].split(",")[0]
        for line in text.splitlines()
        if line.startswith("*ORIENTATION")
    ]
    assert len(set(names)) == len(angles)
    # Every ply line refers to a card that the deck actually defines.
    for line in _ply_lines(text, "blade"):
        assert line.rsplit(", ", 1)[1] in names


def test_swapping_which_ply_is_outside_swaps_the_deck():
    """[0/90]s and [90/0]s must not emit the same section.

    Both stacks are their own reverse, so this says nothing about reversal --
    test_unsymmetric_stack_keeps_its_order does that. What it pins is that the
    two are distinguishable at all, which is what cases/smoke_cantilever then
    turns into a 4.8x stiffness difference through the solver.
    """
    t = 0.25
    down = Layup.uniform(
        [Ply(t, 0), Ply(t, 90), Ply(t, 90), Ply(t, 0)], long_axis="y", material=MAT
    )
    up = Layup.uniform(
        [Ply(t, 90), Ply(t, 0), Ply(t, 0), Ply(t, 90)], long_axis="y", material=MAT
    )
    orientations = [
        line.rsplit(", ", 1)[1] for line in _ply_lines(down.to_inp(), "blade")
    ]
    assert orientations == ["ori_p0", "ori_p90", "ori_p90", "ori_p0"]
    flipped = [line.rsplit(", ", 1)[1] for line in _ply_lines(up.to_inp(), "blade")]
    assert flipped == ["ori_p90", "ori_p0", "ori_p0", "ori_p90"]


def test_unsymmetric_stack_keeps_its_order():
    """A stack with no symmetry to hide behind, emitted in the order given."""
    plies = [Ply(0.1, a) for a in (0, 45, 90, -45)]
    layup = Layup.uniform(plies, long_axis="x", material=MAT)
    emitted = [
        line.rsplit(", ", 1)[1] for line in _ply_lines(layup.to_inp(), "blade")
    ]
    assert emitted == ["ori_p0", "ori_p45", "ori_p90", "ori_m45"]


def test_long_axis_y_is_a_rotation_not_a_mirror():
    """The y branch must rotate the same way about +z as the x branch.

    The mirrored form is exactly angle negation, so a 0/90 stack cannot detect it
    (+90 and -90 are the same fibre direction) and every other stack shows it only
    in the sign of bend-twist coupling, never in a force magnitude. No solver
    check in this repo can catch it -- see cases/smoke_cantilever/README.md -- so
    it is pinned here, on the cosines themselves.

    The invariant: both branches measure the same physical angle from their own
    0-degree axis, and the y axis is the x axis rotated +90 deg about +z, so
    y at theta must equal x at theta + 90 for every theta. The mirrored form
    fails this at any angle that is not a multiple of 90.
    """
    for angle in (0, 10, 30, 45, 60, 90, -25):
        # Cosines only: the two cards carry different names, by construction.
        assert _cosines(orientation_card(angle, long_axis="y")) == _cosines(
            orientation_card(angle + 90, long_axis="x")
        )
    assert orientation_card(45, long_axis="y") != orientation_card(45, long_axis="x")
    assert orientation_card(30, long_axis="y") != orientation_card(30, long_axis="x")


@pytest.mark.parametrize("long_axis,angle,expected", [
    ("x", 0, (1.0, 0.0)),
    ("x", 90, (0.0, 1.0)),
    ("y", 0, (0.0, 1.0)),
    ("y", 90, (-1.0, 0.0)),
])
def test_zero_degree_ply_runs_along_the_long_axis(long_axis, angle, expected):
    """0 deg is the long axis; positive angle turns counter-clockwise about +z."""
    cosines = orientation_card(angle, long_axis=long_axis).splitlines()[1]
    a1, a2 = (float(v) for v in cosines.split(",")[:2])
    assert (a1, a2) == pytest.approx(expected, abs=1e-9)


def test_orientation_frame_is_right_handed():
    for long_axis in ("x", "y"):
        for angle in (0, 15, 45, 90, 137.5, -60):
            v = [
                float(x)
                for x in orientation_card(angle, long_axis=long_axis)
                .splitlines()[1]
                .split(",")
            ]
            a, b = v[:3], v[3:]
            cross_z = a[0] * b[1] - a[1] * b[0]
            assert cross_z == pytest.approx(1.0, abs=1e-9)  # card carries 10 places


def test_long_axis_is_required():
    """No default: CLAUDE.md says +x, every case deck so far uses +y."""
    with pytest.raises(TypeError):
        cross_ply_symmetric(0.1)  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        Layup(materials=(MAT,), zones=(ZoneLayup("z", (Ply(0.1, 0, "t700"),)),))  # type: ignore[call-arg]
    with pytest.raises(ValueError):
        Layup(
            materials=(MAT,),
            zones=(ZoneLayup("z", (Ply(0.1, 0, "t700"),)),),
            long_axis="z",
        )


def test_cross_ply_emits_orientation_names_not_angles():
    layup = cross_ply_symmetric(0.1, long_axis="y")
    text = layup.to_inp()
    assert "*SHELL SECTION, COMPOSITE, ELSET=blade" in text
    assert "ori_p0" in text and "ori_p90" in text
    section = text.split("COMPOSITE")[1]
    assert ", 0," not in section
    assert text.count("0.1, , cfrp,") == 4
    assert "ENGINEERING CONSTANTS" in text


def test_complex_zone_layup():
    layup = Layup(
        materials=(MAT,),
        zones=(
            ZoneLayup(
                "zone_root",
                (
                    Ply(0.1, 0, "t700"),
                    Ply(0.1, 45, "t700"),
                    Ply(0.1, -45, "t700"),
                    Ply(0.1, 90, "t700"),
                    Ply(0.1, 90, "t700"),
                    Ply(0.1, -45, "t700"),
                    Ply(0.1, 45, "t700"),
                    Ply(0.1, 0, "t700"),
                ),
            ),
            ZoneLayup(
                "zone_tip",
                (Ply(0.1, 0, "t700"), Ply(0.1, 90, "t700"), Ply(0.1, 0, "t700")),
            ),
        ),
        long_axis="y",
    )
    text = layup.to_inp()
    assert "zone_root" in text and "zone_tip" in text
    assert "ori_p45" in text and "ori_m45" in text
    assert text.count("*SHELL SECTION, COMPOSITE") == 2
    # Ply drops are longer stacks inboard: the zones differ in ply count.
    assert len(_ply_lines(text, "zone_root")) == 8
    assert len(_ply_lines(text, "zone_tip")) == 3


def test_a_non_finite_ply_angle_is_refused():
    """A sweep hands this module computed numbers; nan must not reach a deck."""
    for bad in (float("nan"), float("inf")):
        with pytest.raises(ValueError, match="must be finite"):
            Layup.uniform([Ply(0.25, bad)], long_axis="y", material=MAT).to_inp()


def test_a_ply_material_must_be_declared():
    with pytest.raises(ValueError, match="not in materials"):
        Layup(
            materials=(MAT,),
            zones=(ZoneLayup("z", (Ply(0.1, 0, "unobtainium"),)),),
            long_axis="y",
        )


def test_deck_assemble_orders_blocks():
    layup = cross_ply_symmetric(0.1, long_axis="y")
    mesh = (
        "*NODE\n1, 0, 0, 0\n"
        "*ELEMENT, TYPE=S8R, ELSET=blade\n"
        "1, 1, 1, 1, 1, 1, 1, 1, 1"
    )
    step = StaticStep(tip_u_clamp_body({10: (0.0, -1.0, 2.0)}), inc=100)
    deck = assemble(
        mesh_inp=mesh,
        layup=layup,
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=[step],
        heading="unit test deck",
    )
    assert deck.index("*NODE") < deck.index("*MATERIAL")
    assert deck.index("*MATERIAL") < deck.index("*ORIENTATION")
    assert deck.index("*ORIENTATION") < deck.index("*SHELL SECTION")
    assert deck.index("*SHELL SECTION") < deck.index("*STEP")
    assert "*EL PRINT, ELSET=blade, TOTALS=ONLY" in deck
    assert "ELSE" in deck
    assert "10, 3, 3, 2.0000000000" in deck


def test_deck_has_no_trailing_comments():
    """ccx fails to parse a card with `**` after data on the same line."""
    layup = cross_ply_symmetric(0.1, long_axis="y")
    deck = assemble(
        mesh_inp="*NODE\n1, 0, 0, 0",
        layup=layup,
        steps=[StaticStep("*BOUNDARY\n1, 3, 3, 1.0")],
        heading="a heading\nover two lines",
    )
    for line in deck.splitlines():
        assert "**" not in line or line.startswith("**")


def test_f_spring_doubles_with_theta():
    c = 100.0
    arm = 100.0
    th = math.pi / 2
    u90 = 0.5 * c * th**2
    u180 = 0.5 * c * (2 * th) ** 2
    _m90, f90 = f_spring(u90, th, arm)
    _m180, f180 = f_spring(u180, 2 * th, arm)
    assert f180 / f90 == pytest.approx(2.0)


# --- woven derivation -------------------------------------------------------
# Checked against cases/smoke_cantilever/clpt.py, which shares no code with
# compfea. woven_from_ud claims a membrane identity; these pin exactly how far
# that claim reaches, so nobody can quietly widen it into a bending one.

_CASE = Path(__file__).resolve().parents[1] / "cases" / "smoke_cantilever"


def _clpt():
    spec = importlib.util.spec_from_file_location("layup_clpt", _CASE / "clpt.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _kw(ec):
    return dict(e1=ec.e1, e2=ec.e2, nu12=ec.nu12, g12=ec.g12)


WOVEN = woven_from_ud(UD_CFRP_GENERIC)
TOTAL_MM = 1.0


def test_woven_is_balanced_and_lands_on_real_weave_values():
    assert WOVEN.e1 == pytest.approx(WOVEN.e2)
    assert WOVEN.g13 == pytest.approx(WOVEN.g23)
    assert WOVEN.nu13 == pytest.approx(WOVEN.nu23)
    # out-of-plane pairs are averaged across the 0 and 90 plies, like g13/g23
    assert WOVEN.nu13 == pytest.approx(0.375)
    # Derived, not typed in -- but it has to come out somewhere a real 2x2
    # twill actually lives, or the derivation is wrong.
    assert WOVEN.e1 == pytest.approx(72332.7, rel=1e-4)
    assert WOVEN.nu12 == pytest.approx(0.0375, rel=1e-4)
    assert WOVEN.g12 == pytest.approx(UD_CFRP_GENERIC.g12)


def test_woven_membrane_matches_a_ud_cross_ply_pair():
    """The claim: A is identical to a [0/90] UD pair of the same thickness."""
    clpt = _clpt()
    ud = clpt.abd([0.0, 90.0], TOTAL_MM / 2, **_kw(UD_CFRP_GENERIC))
    wov = clpt.abd([0.0], TOTAL_MM, **_kw(WOVEN))
    for i in range(3):
        for j in range(3):
            assert ud[i, j] == pytest.approx(wov[i, j], rel=1e-9, abs=1e-9)


def test_woven_has_no_extension_bending_coupling_but_the_ud_pair_does():
    clpt = _clpt()
    ud = clpt.abd([0.0, 90.0], TOTAL_MM / 2, **_kw(UD_CFRP_GENERIC))
    wov = clpt.abd([0.0], TOTAL_MM, **_kw(WOVEN))
    assert abs(wov[0, 3]) < 1e-9
    assert abs(ud[0, 3]) > 1e3


def test_woven_bending_matches_only_the_interleaved_stack():
    """D is equal for [0/90]_n and far apart for the symmetric [0/90]_ns.

    sweep.py builds the symmetric one, so this gap is what a UD-vs-woven sweep
    point actually measures. If this test ever goes quiet, the two fibre types
    have collapsed onto the same bending stiffness and the axis is dead.
    """
    clpt = _clpt()
    wov2 = clpt.abd([0.0, 0.0], TOTAL_MM / 2, **_kw(WOVEN))

    interleaved = clpt.abd([0.0, 90.0, 0.0, 90.0], TOTAL_MM / 4, **_kw(UD_CFRP_GENERIC))
    assert interleaved[3, 3] == pytest.approx(wov2[3, 3], rel=1e-9)

    symmetric = clpt.abd([0.0, 90.0, 90.0, 0.0], TOTAL_MM / 4, **_kw(UD_CFRP_GENERIC))
    assert symmetric[3, 3] / wov2[3, 3] == pytest.approx(1.65625, rel=1e-6)


def test_crimp_factor_scales_in_plane_only():
    soft = woven_from_ud(UD_CFRP_GENERIC, crimp_factor=0.9)
    assert soft.e1 == pytest.approx(0.9 * WOVEN.e1)
    assert soft.g12 == pytest.approx(0.9 * WOVEN.g12)
    assert soft.e3 == pytest.approx(WOVEN.e3)
    assert soft.g13 == pytest.approx(WOVEN.g13)
    assert soft.nu12 == pytest.approx(WOVEN.nu12)
    with pytest.raises(ValueError):
        woven_from_ud(UD_CFRP_GENERIC, crimp_factor=0.0)
    with pytest.raises(ValueError):
        woven_from_ud(UD_CFRP_GENERIC, crimp_factor=1.5)
