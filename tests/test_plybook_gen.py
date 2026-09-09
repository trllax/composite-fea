"""Design vector -> ply book: stacking order, symmetry, inventory thickness.

Pure: no solver, no gmsh. The generator turns a ``TaperDesign`` into the repo's
existing ``plybook.PlyRow`` tuple such that the laminate is mirror-symmetric
about the mid-plane by construction (angles unchanged across the mirror, *not*
sign-flipped -- see CLAUDE.md on ``D16``), its thicknesses come from the shop
inventory and not a guess, and it round-trips through ``plybook.load_plybook``
unchanged.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from compfea.abd import COUPLING_REL_TOL, laminate_summary
from compfea.plybook import check_against_library, load_plybook, resolve
from compfea.plybook_gen import (
    PlyBookGenError,
    PlySpec,
    TaperDesign,
    design_fingerprint,
    design_from_dict,
    design_rows,
    design_to_dict,
    expand_grid,
    write_plybook,
)
from compfea.shop_inventory import load_shop_inventory, materials_library

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "materials" / "shop_inventory.csv"

TWILL = "hexcel_im7_twill_205"
UD = "hexcel_im2_uni_193"
BIAX = "hexcel_himax_biax_100"


@pytest.fixture(scope="module")
def skus():
    return load_shop_inventory(INVENTORY)


def _shop_like_design() -> TaperDesign:
    """A FIN_TEST_3-shaped design: twill skin, UD+biax core, inboard pad pairs."""
    return TaperDesign(
        through_zone="FULL",
        skins=(PlySpec(TWILL, 0.0),),
        core=(PlySpec(UD, 0.0), PlySpec(BIAX, 0.0)),
        pads=(
            ("HEAL", (PlySpec(TWILL, 0.0), PlySpec(UD, 0.0))),
            ("QUARTER", (PlySpec(TWILL, 0.0),)),
        ),
        pad_order=("TIP", "MID", "QUARTER", "HEAL"),
    )


def _unbalanced_design() -> TaperDesign:
    """One off-axis unbalanced UD ply in the core -- exercises the mirror."""
    return TaperDesign(
        through_zone="FULL",
        skins=(PlySpec(TWILL, 0.0),),
        core=(PlySpec(UD, 30.0), PlySpec(BIAX, 0.0)),
        core_center=PlySpec(UD, 90.0),
        pads=(("HEAL", (PlySpec(UD, 45.0),)),),
        pad_order=("HEAL",),
    )


def _coverages():
    """Synthetic element coverage: one FULL-only, one +HEAL, one +QUARTER."""
    return {
        1: frozenset({"FULL"}),
        2: frozenset({"FULL", "HEAL"}),
        3: frozenset({"FULL", "QUARTER"}),
    }


def _zone_thickness_for_element(layup, elsets, eid):
    """Resolved-zone thickness for the element with the given id."""
    name = next(n for n, eids in elsets.items() if eid in eids)
    return next(z.thickness for z in layup.zones if z.elset == name)


# --------------------------------------------------------------------------
# symmetry -- B must cancel to float noise for every zone


@pytest.mark.parametrize("design_fn", [_shop_like_design, _unbalanced_design])
def test_design_has_no_bending_extension_coupling(design_fn, skus):
    rows = design_rows(design_fn(), skus=skus)
    library = materials_library(skus)
    layup, _ = resolve(rows, library, _coverages(), long_axis="y")

    summary = laminate_summary(layup)
    assert summary, "no zones resolved"
    for zone in summary:
        assert zone["symmetric"] is True, f"{zone['zone']} is not symmetric"
        scale = abs(zone["a11"]) + abs(zone["a22"])
        for term in ("b11", "b22", "b12", "b16", "b26", "b66"):
            assert abs(zone[term]) <= COUPLING_REL_TOL * scale, term


def test_mirror_keeps_angles_not_sign_flips_them(skus):
    """The one bug D16 catches: a negated mirror. Palindrome, no sign flip."""
    rows = design_rows(_unbalanced_design(), skus=skus)
    full = [r for r in rows if r.zone == "FULL"]
    angles = [r.angle_deg for r in full]
    # skin, core(30, biax0), centre(90), core reversed, skin
    assert angles == [0.0, 30.0, 0.0, 90.0, 0.0, 30.0, 0.0]
    assert angles == list(reversed(angles))  # palindrome
    assert angles != [-a for a in reversed(angles)]  # NOT an antisymmetric mirror
    materials = [r.material for r in full]
    assert materials == list(reversed(materials))


def test_core_center_sits_once_at_the_mid_plane(skus):
    rows = design_rows(_unbalanced_design(), skus=skus)
    full = [(r.material, r.angle_deg) for r in rows if r.zone == "FULL"]
    assert len(full) % 2 == 1  # odd: a centre ply exists
    assert full[len(full) // 2] == (UD, 90.0)
    assert [(r.material, r.angle_deg) for r in rows if r.zone == "FULL"].count(
        (UD, 90.0)
    ) == 1


def test_stacking_order_is_mirrored_about_the_mid_plane(skus):
    rows = design_rows(_shop_like_design(), skus=skus)
    full = [r for r in rows if r.zone == "FULL"]
    seq = [(r.material, r.angle_deg) for r in full]
    assert seq == list(reversed(seq))
    assert seq[0] == (TWILL, 0.0)
    assert [m for m, _ in seq] == [TWILL, UD, BIAX, BIAX, UD, TWILL]


# --------------------------------------------------------------------------
# round-trip through the existing loader


def test_write_plybook_round_trips_through_load_plybook(skus, tmp_path):
    design = _unbalanced_design()  # carries a non-integer-free angle path
    path = write_plybook(design, tmp_path / "gen.csv", skus=skus, note="unit test")
    loaded = load_plybook(path)
    generated = design_rows(design, skus=skus)

    assert len(loaded) == len(generated)
    for got, want in zip(loaded, generated, strict=True):
        assert got.ply == want.ply
        assert got.zone == want.zone
        assert got.material == want.material
        assert got.angle_deg == want.angle_deg  # exact: written at 6 dp
        assert got.thickness_mm == want.thickness_mm  # exact: written at 9 dp
        assert got.kind == want.kind

    check_against_library(loaded, materials_library(skus))


def test_write_plybook_is_exact_for_a_fractional_angle(skus, tmp_path):
    design = TaperDesign(
        through_zone="FULL", core=(PlySpec(UD, 22.512345),)
    )
    path = write_plybook(design, tmp_path / "frac.csv", skus=skus)
    loaded = load_plybook(path)
    assert loaded[0].angle_deg == pytest.approx(22.512345, abs=1e-9)


# --------------------------------------------------------------------------
# thickness provenance and the taper


def test_every_ply_thickness_comes_from_the_inventory(skus):
    rows = design_rows(_shop_like_design(), skus=skus)
    for r in rows:
        assert r.thickness_mm == pytest.approx(skus[r.material].thickness_mm)


def test_pads_thicken_the_zone_they_name(skus):
    rows = design_rows(_shop_like_design(), skus=skus)
    library = materials_library(skus)
    layup, elsets = resolve(rows, library, _coverages(), long_axis="y")

    t_full = _zone_thickness_for_element(layup, elsets, 1)
    t_heal = _zone_thickness_for_element(layup, elsets, 2)
    t_quarter = _zone_thickness_for_element(layup, elsets, 3)
    # HEAL adds a pad pair (2x twill + 2x UD), QUARTER one pair (2x twill).
    assert t_heal == pytest.approx(
        t_full + 2 * skus[TWILL].thickness_mm + 2 * skus[UD].thickness_mm
    )
    assert t_quarter == pytest.approx(t_full + 2 * skus[TWILL].thickness_mm)
    assert t_full < t_quarter < t_heal


# --------------------------------------------------------------------------
# the generator refuses a design nobody can build


def test_unstocked_material_raises(skus):
    design = TaperDesign(through_zone="FULL", core=(PlySpec("unobtanium_ud", 0.0),))
    with pytest.raises(PlyBookGenError, match="not stocked"):
        design_rows(design, skus=skus)


def test_non_finite_angle_raises(skus):
    design = TaperDesign(through_zone="FULL", core=(PlySpec(UD, math.nan),))
    with pytest.raises(PlyBookGenError, match="finite"):
        design_rows(design, skus=skus)


def test_nothing_covering_the_through_zone_raises(skus):
    design = TaperDesign(
        through_zone="FULL",
        pads=(("HEAL", (PlySpec(UD, 0.0),)),),
        pad_order=("HEAL",),
    )
    with pytest.raises(PlyBookGenError, match="bare"):
        design_rows(design, skus=skus)


def test_pad_zone_missing_from_pad_order_raises_on_construction():
    # Direct construction, not just design_from_dict: a pad silently dropped
    # here would fingerprint-collide with a no-pad design in the sweep cache.
    with pytest.raises(PlyBookGenError, match="pad_order"):
        TaperDesign(
            through_zone="FULL",
            core=(PlySpec(UD, 0.0),),
            pads=(("HEAL", (PlySpec(TWILL, 0.0),)),),
            pad_order=(),
        )


def test_duplicate_zone_in_pad_order_raises():
    with pytest.raises(PlyBookGenError, match="duplicate"):
        TaperDesign(
            through_zone="FULL",
            core=(PlySpec(UD, 0.0),),
            pads=(("HEAL", (PlySpec(TWILL, 0.0),)),),
            pad_order=("HEAL", "HEAL"),
        )


def test_duplicate_pad_zone_in_sequence_form_raises_cleanly():
    # sequence-form pads with a repeated zone must be a PlyBookGenError, not a
    # raw TypeError from sorting two PlySpec tuples.
    with pytest.raises(PlyBookGenError, match="more than once"):
        TaperDesign(
            through_zone="FULL",
            core=(PlySpec(UD, 0.0),),
            pads=(
                ("HEAL", (PlySpec(TWILL, 0.0),)),
                ("HEAL", (PlySpec(UD, 0.0),)),
            ),
            pad_order=("HEAL",),
        )


@pytest.mark.parametrize("name", ["all", "*", "", "  "])
def test_whole_part_zone_name_is_rejected(name):
    with pytest.raises(PlyBookGenError, match="named zone"):
        TaperDesign(through_zone=name, core=(PlySpec(UD, 0.0),))


# --------------------------------------------------------------------------
# dict round-trip and the sweep cache key


def test_dict_round_trip_is_identity():
    for design in (_shop_like_design(), _unbalanced_design()):
        assert design_from_dict(design_to_dict(design)) == design


def test_dict_round_trip_normalises_pad_order_independent_of_input_order():
    a = TaperDesign(
        through_zone="FULL",
        core=(PlySpec(UD, 0.0),),
        pads=(("QUARTER", (PlySpec(UD, 0.0),)), ("HEAL", (PlySpec(UD, 0.0),))),
        pad_order=("QUARTER", "HEAL"),
    )
    b = design_from_dict(design_to_dict(a))
    assert a == b  # pads normalised to sorted items on both paths


def test_fingerprint_is_stable_and_moves_with_a_pad(skus):
    base = _shop_like_design()
    fp = design_fingerprint(base, skus=skus)
    assert fp == design_fingerprint(design_from_dict(design_to_dict(base)), skus=skus)
    assert len(fp) == 16

    heavier = design_from_dict(
        {
            **design_to_dict(base),
            "pads": {
                **design_to_dict(base)["pads"],
                "HEAL": [[TWILL, 0], [UD, 0], [BIAX, 0]],
            },
        }
    )
    assert design_fingerprint(heavier, skus=skus) != fp


def test_fingerprint_ignores_pad_mapping_order(skus):
    base = design_to_dict(_shop_like_design())
    reordered = {**base, "pads": dict(reversed(list(base["pads"].items())))}
    assert design_fingerprint(
        design_from_dict(base), skus=skus
    ) == design_fingerprint(design_from_dict(reordered), skus=skus)


def test_fingerprint_moves_when_a_lamina_modulus_changes(skus):
    """A modulus edit under a stable SKU name must miss the cache."""
    import dataclasses

    base = _shop_like_design()
    fp = design_fingerprint(base, skus=skus)

    stiffer = dict(skus)
    sku = stiffer[UD]
    bumped = dataclasses.replace(
        sku.record.constants, e1=sku.record.constants.e1 * 1.1
    )
    stiffer[UD] = dataclasses.replace(
        sku, record=dataclasses.replace(sku.record, constants=bumped)
    )
    assert design_fingerprint(base, skus=stiffer) != fp


# --------------------------------------------------------------------------
# grid expansion for the sweep


def test_expand_grid_cardinality_and_validity(skus):
    base = design_to_dict(_shop_like_design())
    grid = expand_grid(
        base,
        {
            "core": [
                [[UD, 0]],
                [[UD, 0], [BIAX, 0]],
                [[UD, 0], [UD, 0], [BIAX, 0]],
            ],
            "pads.HEAL": [[], [[TWILL, 0]]],
        },
    )
    assert len(grid) == 6
    fps = {design_fingerprint(d, skus=skus) for d in grid}
    assert len(fps) == 6  # every grid point is a distinct laminate


def test_expand_grid_tolerates_a_null_pads_field(skus):
    base = {
        **design_to_dict(_shop_like_design()),
        "pads": None,
        "pad_order": ["HEAL"],
    }
    grid = expand_grid(base, {"pads.HEAL": [[[TWILL, 0]]]})
    assert len(grid) == 1
    assert grid[0].pads_map["HEAL"] == (PlySpec(TWILL, 0.0),)
