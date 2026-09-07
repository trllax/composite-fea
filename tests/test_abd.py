"""Classical lamination theory: A, B, D, and the bend-twist sign.

Validated against values this repo already pins independently -- the four
reference stacks in cases/smoke_cantilever, whose D11, narrow-strip EI and
inverted-ABD coupling ratio are hand-computed and written down there.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from compfea.abd import (
    abd,
    beam_ei_free_edge,
    bend_twist_sign,
    is_balanced,
    is_symmetric,
    laminate_summary,
    ply_z_stations,
    qbar,
)
from compfea.layup import (
    UD_CFRP_GENERIC,
    Layup,
    Ply,
    ZoneLayup,
    woven_from_ud,
)

ROOT = Path(__file__).resolve().parents[1]
MATS = {"cfrp": UD_CFRP_GENERIC}
PLY_MM = 0.25
WIDTH_MM = 20.0

#: Hand-computed in cases/smoke_cantilever/README.md and pinned in
#: tests/test_smoke_cantilever.py. Reproduced here from an independent
#: implementation: (angles, D11, narrow-strip EI, inverted-ABD b11/d11).
REFERENCE = {
    "[0/90]s": ([0, 90, 90, 0], 9997.4849, 199455.826, 0.0),
    "[90/0]s": ([90, 0, 0, 90], 2074.9497, 41396.492, 0.0),
    "[0/0/90/90]": ([0, 0, 90, 90], 6036.2173, 120554.577, +0.21875),
    "[90/90/0/0]": ([90, 90, 0, 0], 6036.2173, 120554.577, -0.21875),
}


def zone(angles, *, thickness=PLY_MM, material="cfrp", name="blade") -> ZoneLayup:
    return ZoneLayup(name, tuple(Ply(thickness, a, material) for a in angles))


def coupling_ratio(a, b, d) -> float:
    """b11/d11 from the inverted 6x6 ABD, mm -- clpt.axial_coupling_ratio."""
    compliance = np.linalg.inv(np.block([[a, b], [b, d]]))
    return float(compliance[0, 3] / compliance[3, 3])


# --------------------------------------------------------------------------
# against the pinned reference


@pytest.mark.parametrize("label", sorted(REFERENCE))
def test_d11_and_ei_match_the_pinned_hand_computation(label):
    angles, d11, ei, _ = REFERENCE[label]
    _a, _b, d = abd(zone(angles), MATS)
    assert d[0, 0] == pytest.approx(d11, rel=1e-7)
    assert beam_ei_free_edge(d, WIDTH_MM) == pytest.approx(ei, rel=1e-7)


@pytest.mark.parametrize("label", sorted(REFERENCE))
def test_the_coupling_ratio_matches_and_pins_B(label):
    """The only reference here that exercises B rather than D."""
    angles, _d11, _ei, ratio = REFERENCE[label]
    a, b, d = abd(zone(angles), MATS)
    assert coupling_ratio(a, b, d) == pytest.approx(ratio, abs=1e-9)


def test_this_agrees_with_the_independent_clpt_module():
    """Two implementations, written apart; they must not have drifted."""
    sys.path.insert(0, str(ROOT / "cases" / "smoke_cantilever"))
    try:
        import clpt
    finally:
        sys.path.pop(0)
    material = dict(e1=135000.0, e2=9000.0, nu12=0.30, g12=4500.0)
    for angles, *_ in REFERENCE.values():
        _a, _b, d = abd(zone(angles), MATS)
        assert d[0, 0] == pytest.approx(clpt.d11(angles, PLY_MM, **material), rel=1e-12)
        assert beam_ei_free_edge(d, WIDTH_MM) == pytest.approx(
            clpt.ei_narrow_strip(angles, PLY_MM, width=WIDTH_MM, **material),
            rel=1e-12,
        )


def test_reversing_a_stack_leaves_D_and_flips_B():
    """The identity the -z convention rests on, stated on the matrices."""
    a1, b1, d1 = abd(zone([0, 0, 90, 90]), MATS)
    a2, b2, d2 = abd(zone([90, 90, 0, 0]), MATS)
    assert np.allclose(d1, d2, rtol=1e-12)
    assert np.allclose(a1, a2, rtol=1e-12)
    assert np.allclose(b1, -b2, rtol=1e-12)
    assert not np.allclose(b1, np.zeros((3, 3)))


def test_free_edge_and_wide_plate_differ_by_the_documented_0_247_percent():
    _a, _b, d = abd(zone([0, 90, 90, 0]), MATS)
    wide = d[0, 0] * WIDTH_MM
    narrow = beam_ei_free_edge(d, WIDTH_MM)
    assert (wide - narrow) / narrow == pytest.approx(0.00247, abs=5e-5)


# --------------------------------------------------------------------------
# the bend-twist detector


def test_only_a_cross_ply_stack_is_truly_blind_to_the_mirror():
    """D16 vanishes exactly when every angle is its own negation: 0 and 90."""
    for angles in ([0, 90, 90, 0], [0, 90], [90, 0], [0, 0, 0], [90, 90]):
        _a, _b, d = abd(zone(angles), MATS)
        assert bend_twist_sign(d) == 0, f"{angles} should be blind to the mirror"


def test_a_balanced_symmetric_angle_ply_is_caught_by_D16_alone():
    """The stack CLAUDE.md says no force check in this repo can distinguish.

    [45/-45]s is balanced (A16 = 0) and symmetric (B = 0), so its reaction is
    identical to its mirror's -- which is why layup.py says the angle
    convention cannot be caught by a force. It is *not* free of bend-twist
    coupling: the +45 and -45 plies sit at different z, so they cancel in
    membrane and not in bending. D16 is 1980.6 N.mm here, and it changes sign
    under the mirror. This is the gap the module was built to close.
    """
    angles = [45, -45, -45, 45]
    a, b, d = abd(zone(angles), MATS)
    assert a[0, 2] == pytest.approx(0.0, abs=1e-9), "balanced: A16 must vanish"
    assert is_symmetric(b, a), "symmetric: B must vanish"
    assert d[0, 2] == pytest.approx(1980.63380282, rel=1e-8)

    _a2, _b2, dm = abd(zone([-x for x in angles]), MATS)
    assert bend_twist_sign(d) == +1
    assert bend_twist_sign(dm) == -1


def test_bend_twist_sign_flips_when_every_angle_is_mirrored():
    """The mirror is angle negation, and D16 is what sees it.

    layup.py states that writing the mirrored *ORIENTATION form "flips the sign
    of bend-twist coupling" and that no force check in this repo can catch it.
    This is that check.
    """
    angles = [45, 45, 0, 90]           # unbalanced on purpose
    _a, _b, d = abd(zone(angles), MATS)
    _a2, _b2, d2 = abd(zone([-x for x in angles]), MATS)
    assert bend_twist_sign(d) != 0, "the probe stack must be sensitive"
    assert bend_twist_sign(d) == -bend_twist_sign(d2)
    assert d[0, 2] == pytest.approx(-d2[0, 2], rel=1e-12)


def test_a_plus_45_stack_has_a_specific_D16_sign():
    """Pinned, not printed. A sign that can flip silently is not a check."""
    _a, _b, d = abd(zone([45, 45, 45, 45]), MATS)
    assert d[0, 2] > 0.0
    assert bend_twist_sign(d) == +1
    _a, _b, dm = abd(zone([-45, -45, -45, -45]), MATS)
    assert bend_twist_sign(dm) == -1


# --------------------------------------------------------------------------
# mechanics of the module


def test_ply_z_stations_put_ply_one_at_minus_half_thickness():
    z = ply_z_stations(zone([0, 90, 90, 0]))
    assert z[0][0] == pytest.approx(-0.5)
    assert z[-1][1] == pytest.approx(+0.5)
    assert all(b == pytest.approx(a) for (_, b), (a, _) in zip(z, z[1:], strict=False))


def test_unequal_plies_put_ply_one_thickness_at_the_bottom():
    """Uniform stacks cannot see a reversed z assignment; unequal ones can.

    Every other stack in this file has one ply thickness, which makes reversing
    the z stations a no-op -- a mutation doing exactly that passed all of them.
    A real ply book has mixed thicknesses (the fin runs 0.15/0.10/0.20/0.25),
    so this is the case that pins ply 1 to -z.
    """
    stack = ZoneLayup("b", (Ply(0.1, 0.0), Ply(0.3, 90.0)))
    assert ply_z_stations(stack) == [
        pytest.approx((-0.2, -0.1)),
        pytest.approx((-0.1, 0.2)),
    ]


def test_the_stiff_ply_at_the_bottom_gives_B_the_opposite_sign():
    """Unequal thicknesses, so the check cannot be satisfied by symmetry alone."""
    low = ZoneLayup("b", (Ply(0.3, 0.0), Ply(0.1, 90.0)))     # stiff ply at -z
    high = ZoneLayup("b", (Ply(0.1, 90.0), Ply(0.3, 0.0)))    # stiff ply at +z
    _a, b_low, _d = abd(low, MATS)
    _a, b_high, _d = abd(high, MATS)
    assert b_low[0, 0] < 0.0 < b_high[0, 0]
    assert b_low[0, 0] == pytest.approx(-b_high[0, 0], rel=1e-12)


def test_qbar_at_zero_is_the_on_axis_reduced_stiffness():
    q = qbar(UD_CFRP_GENERIC, 0.0)
    nu21 = 0.30 * 9000.0 / 135000.0
    denom = 1.0 - 0.30 * nu21
    assert q[0, 0] == pytest.approx(135000.0 / denom)
    assert q[1, 1] == pytest.approx(9000.0 / denom)
    assert q[2, 2] == pytest.approx(4500.0)
    assert q[0, 2] == pytest.approx(0.0, abs=1e-9)


def test_qbar_at_ninety_swaps_the_axial_and_transverse_terms():
    q0, q90 = qbar(UD_CFRP_GENERIC, 0.0), qbar(UD_CFRP_GENERIC, 90.0)
    assert q90[0, 0] == pytest.approx(q0[1, 1])
    assert q90[1, 1] == pytest.approx(q0[0, 0])
    assert q90[0, 2] == pytest.approx(0.0, abs=1e-9)


def test_qbar_is_symmetric_at_every_angle():
    for angle in (-67.5, -45.0, 0.0, 22.5, 45.0, 90.0, 135.0):
        q = qbar(UD_CFRP_GENERIC, angle)
        assert np.allclose(q, q.T, rtol=1e-12)


def test_thickness_scaling_is_cubic_in_D_and_linear_in_A():
    a1, _b1, d1 = abd(zone([0, 90, 90, 0], thickness=0.25), MATS)
    a2, _b2, d2 = abd(zone([0, 90, 90, 0], thickness=0.50), MATS)
    assert a2[0, 0] == pytest.approx(2.0 * a1[0, 0], rel=1e-12)
    assert d2[0, 0] == pytest.approx(8.0 * d1[0, 0], rel=1e-12)


def test_symmetric_and_balanced_flags():
    a, b, _d = abd(zone([0, 90, 90, 0]), MATS)
    assert is_symmetric(b, a) and is_balanced(zone([0, 90, 90, 0]))
    a, b, _d = abd(zone([0, 0, 90, 90]), MATS)
    assert not is_symmetric(b, a)
    assert is_balanced(zone([45, -45, -45, 45]))
    assert not is_balanced(zone([45, 45, -45, 0]))
    # 0 and 90 are their own opposites and never unbalance a stack.
    assert is_balanced(zone([0, 90, 0, 90]))


def test_an_unknown_ply_material_names_the_zone():
    with pytest.raises(KeyError, match="zone 'tip'"):
        abd(zone([0, 90], material="basalt", name="tip"), MATS)


# --------------------------------------------------------------------------
# the report


def test_laminate_summary_reports_every_zone_with_flat_keys():
    layup = Layup(
        materials=(UD_CFRP_GENERIC, woven_from_ud(UD_CFRP_GENERIC)),
        zones=(
            zone([0, 90, 90, 0], name="root"),
            ZoneLayup("tip", (Ply(0.2, 45.0, "cfrp_woven"), Ply(0.2, -45.0, "cfrp"))),
        ),
        long_axis="x",
    )
    rows = {r["zone"]: r for r in laminate_summary(layup)}
    assert set(rows) == {"root", "tip"}
    assert rows["root"]["stack"] == "[0/90/90/0]"
    assert rows["root"]["n_plies"] == 4
    assert rows["root"]["thickness_mm"] == pytest.approx(1.0)
    assert rows["root"]["symmetric"] is True
    assert rows["tip"]["symmetric"] is False
    for key in ("a11", "a12", "a16", "b11", "d11", "d16", "d66"):
        assert isinstance(rows["root"][key], float)
    assert rows["root"]["d11"] == pytest.approx(REFERENCE["[0/90]s"][1], rel=1e-7)
    assert rows["root"]["ef_x_mpa"] > rows["root"]["ef_y_mpa"]


def test_flexural_modulus_of_a_single_ply_is_its_own_modulus():
    """A homogeneous 0-degree laminate must bend at E1, near enough."""
    layup = Layup(
        materials=(UD_CFRP_GENERIC,), zones=(zone([0, 0, 0, 0]),), long_axis="x"
    )
    row = laminate_summary(layup)[0]
    assert row["ef_x_mpa"] == pytest.approx(135000.0, rel=1e-6)
