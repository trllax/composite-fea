"""Ply book CSV -> Layup: ordering, zone binding, material resolution."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from compfea.layup import Layup, Ply, ZoneLayup, layup_from_coverage
from compfea.materials import load_materials
from compfea.plybook import (
    PlyBookError,
    check_against_library,
    check_against_mesh,
    load_plybook,
    plies_from_rows,
    resolve,
    stack_label,
    stack_table,
    unused_zones,
)

ROOT = Path(__file__).resolve().parents[1]
GENERIC = ROOT / "materials" / "generic.csv"
FIN2 = ROOT / "test_fin_2.step"
EQUIV = ROOT / "cases" / "step_fin" / "plybook_equivalent.csv"

HEADER = "ply,zone,material,angle_deg,thickness_mm,kind"
ROWS = [
    "1,FULL,cfrp,0,0.125,ud",
    "2,FULL,cfrp,45,0.125,ud",
    "3,SPAR,cfrp,0,0.250,ud",
    "4,FULL,cfrp,-45,0.125,ud",
]


def write(tmp_path: Path, *lines: str, name: str = "plies.csv") -> Path:
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------------
# ordering -- the -z convention is the bug class this repo exists to catch


def test_rows_come_back_in_stacking_order_ply_one_first(tmp_path):
    rows = load_plybook(write(tmp_path, HEADER, *reversed(ROWS)))
    assert [r.ply for r in rows] == [1, 2, 3, 4]
    assert [r.angle_deg for r in rows] == [0.0, 45.0, 0.0, -45.0]


def test_ply_one_is_the_minus_z_ply(tmp_path):
    """First row -> first *SHELL SECTION line -> the -z ply, per CLAUDE.md."""
    rows = load_plybook(write(tmp_path, HEADER, *ROWS))
    layup = Layup(
        materials=(load_materials(GENERIC)["cfrp"].constants,),
        zones=(ZoneLayup("blade", plies_from_rows(rows)),),
        long_axis="x",
    )
    lines = [
        ln for ln in layup.to_inp().splitlines() if ln.startswith(("0.125", "0.25"))
    ]
    assert lines[0].endswith("ori_p0"), "ply 1 must be the first section line"
    assert lines[-1].endswith("ori_m45")
    table = stack_table(layup)
    assert table[0]["z_bot_mm"] == pytest.approx(-0.5 * 0.625)
    assert table[0]["ply"] == 1 and table[0]["angle_deg"] == 0.0
    assert table[-1]["z_top_mm"] == pytest.approx(+0.5 * 0.625)


def test_a_gap_in_the_ply_numbers_is_refused(tmp_path):
    """A gap means a row was deleted; the stack is not the one that was drawn."""
    bad = [r for r in ROWS if not r.startswith("3,")]
    with pytest.raises(PlyBookError, match="no gaps or duplicates"):
        load_plybook(write(tmp_path, HEADER, *bad))


def test_a_duplicate_ply_number_is_refused(tmp_path):
    dup = ROWS[:3] + ["3,FULL,cfrp,-45,0.125,ud"]
    with pytest.raises(PlyBookError, match="duplicated: \\[3\\]"):
        load_plybook(write(tmp_path, HEADER, *dup))


def test_thickness_and_angle_are_validated(tmp_path):
    with pytest.raises(PlyBookError, match="must be finite and > 0"):
        load_plybook(write(tmp_path, HEADER, "1,FULL,cfrp,0,0,ud"))
    with pytest.raises(PlyBookError, match="must be finite and > 0"):
        load_plybook(write(tmp_path, HEADER, "1,FULL,cfrp,0,-0.1,ud"))
    with pytest.raises(PlyBookError, match="angle_deg"):
        load_plybook(write(tmp_path, HEADER, "1,FULL,cfrp,nan,0.1,ud"))


def test_missing_columns_are_refused(tmp_path):
    with pytest.raises(PlyBookError, match="missing required columns"):
        load_plybook(write(tmp_path, "ply,zone,material,angle_deg", "1,FULL,cfrp,0"))


# --------------------------------------------------------------------------
# zone binding


def test_zone_names_are_sanitized_the_way_the_mesh_sanitizes_them(tmp_path):
    """You write 3_4ths; ccx and the mesh both call it z_3_4ths."""
    rows = load_plybook(write(tmp_path, HEADER, "1,3_4ths,cfrp,0,0.1,ud"))
    assert rows[0].zone == "z_3_4ths"


def test_all_and_blank_mean_every_element(tmp_path):
    rows = load_plybook(
        write(tmp_path, HEADER, "1,ALL,cfrp,0,0.1,ud", "2,,cfrp,90,0.1,ud")
    )
    assert [r.zone for r in rows] == [None, None]
    assert [p.coverage for p in plies_from_rows(rows)] == [None, None]


def test_a_zone_the_mesh_does_not_have_is_refused(tmp_path):
    rows = load_plybook(write(tmp_path, HEADER, *ROWS))
    with pytest.raises(PlyBookError, match="does not have"):
        check_against_mesh(rows, {"blade": (1,), "FULL": (1,)})


def test_a_nested_zone_with_no_ply_of_its_own_is_legal(tmp_path):
    """QUARTER on the real fin: covered by the plies above it, needs none itself.

    Refusing this would refuse a perfectly good laminate. The hazard it looks
    like -- an element with an empty stack -- is caught in layup_from_coverage.
    """
    rows = load_plybook(write(tmp_path, HEADER, *ROWS))
    elsets = {"blade": (1, 2), "FULL": (1, 2), "SPAR": (1,), "TIP": (2,)}
    check_against_mesh(rows, elsets)
    assert unused_zones(rows, elsets) == ["TIP"]


def test_an_element_with_no_plies_at_all_is_still_refused(tmp_path):
    rows = load_plybook(write(tmp_path, HEADER, *ROWS))
    with pytest.raises(ValueError, match="matches no plies"):
        resolve(
            rows, load_materials(GENERIC),
            {1: frozenset({"FULL", "SPAR"}), 2: frozenset({"TIP"})},
            long_axis="x",
        )


def test_a_whole_part_ply_covers_every_zone(tmp_path):
    """coverage=None reaches everything, so no zone can be left bare."""
    rows = load_plybook(write(tmp_path, HEADER, "1,ALL,cfrp,0,0.1,ud"))
    check_against_mesh(rows, {"blade": (1, 2), "FULL": (1,), "TIP": (2,)})


# --------------------------------------------------------------------------
# materials


def test_an_unknown_material_is_refused(tmp_path):
    rows = load_plybook(write(tmp_path, HEADER, "1,FULL,basalt_ud,0,0.1,ud"))
    with pytest.raises(PlyBookError, match="not in the materials library"):
        check_against_library(rows, load_materials(GENERIC))


def test_a_kind_that_contradicts_the_library_is_refused(tmp_path):
    """Naming the UD card where the woven one was meant is silent otherwise."""
    rows = load_plybook(write(tmp_path, HEADER, "1,FULL,cfrp,0,0.1,woven"))
    with pytest.raises(PlyBookError, match="is 'ud' in the library"):
        check_against_library(rows, load_materials(GENERIC))


def test_kind_is_optional(tmp_path):
    rows = load_plybook(write(tmp_path, "ply,zone,material,angle_deg,thickness_mm",
                              "1,FULL,cfrp,0,0.1"))
    assert rows[0].kind is None
    check_against_library(rows, load_materials(GENERIC))


def test_two_materials_reach_one_deck(tmp_path):
    """Structurally supported forever; nothing had ever produced one."""
    rows = load_plybook(
        write(tmp_path, HEADER,
              "1,ALL,cfrp,0,0.1,ud",
              "2,ALL,cfrp_woven,45,0.2,woven")
    )
    layup, _ = resolve(
        rows, load_materials(GENERIC), {1: frozenset(), 2: frozenset()}, long_axis="x"
    )
    text = layup.to_inp()
    assert text.count("*MATERIAL, NAME=") == 2
    assert "*MATERIAL, NAME=cfrp\n" in text
    assert "*MATERIAL, NAME=cfrp_woven\n" in text
    assert ", cfrp, ori_p0" in text
    assert ", cfrp_woven, ori_p45" in text


def test_only_referenced_materials_are_emitted(tmp_path):
    """An unused library row must not drift into the deck, or a cache key."""
    rows = load_plybook(write(tmp_path, HEADER, "1,ALL,cfrp,0,0.1,ud"))
    layup, _ = resolve(rows, load_materials(GENERIC), {1: frozenset()}, long_axis="x")
    assert [m.name for m in layup.materials] == ["cfrp"]


def test_material_names_are_never_rewritten_on_the_library_path():
    """layup_from_coverage(material=...) rewrites 'cfrp'; materials=... must not."""
    lib = load_materials(GENERIC)
    plies = (Ply(0.1, 0.0, "cfrp_woven"),)
    layup, _ = layup_from_coverage(
        plies, {1: frozenset()}, long_axis="x",
        materials=(lib["cfrp"].constants, lib["cfrp_woven"].constants),
    )
    assert layup.zones[0].plies[0].material == "cfrp_woven"
    assert [m.name for m in layup.materials] == ["cfrp_woven"]


def test_material_and_materials_together_are_refused():
    lib = load_materials(GENERIC)
    with pytest.raises(ValueError, match="not both"):
        layup_from_coverage(
            (Ply(0.1, 0.0),), {1: frozenset()}, long_axis="x",
            material=lib["cfrp"].constants,
            materials=(lib["cfrp"].constants,),
        )


# --------------------------------------------------------------------------
# the equivalence gate


@pytest.mark.skipif(not FIN2.is_file(), reason="test_fin_2.step not in repo root")
def test_csv_reproduces_the_hardcoded_fin_layup():
    """The strongest check in this phase, and it costs no solve.

    cases/step_fin/plybook_equivalent.csv is the CSV form of run_ubend's
    default_plies(). If the two build the same deck text, then the loader, the
    zone binding, the material resolution and the ply ordering are all right --
    measured against a layup this repo already solves, not against my own
    restatement of it.
    """
    sys.path.insert(0, str(ROOT / "cases" / "step_fin"))
    try:
        from run_ubend import default_plies
    finally:
        sys.path.pop(0)

    from compfea.layup import coverages_from_mesh
    from compfea.step_mesh import mesh_step

    mesh = mesh_step(FIN2, size_mm=40.0)
    coverages = coverages_from_mesh(mesh.elsets)

    from_py, elsets_py = layup_from_coverage(
        default_plies(), coverages, long_axis="x"
    )
    rows = load_plybook(EQUIV)
    from_csv, elsets_csv = resolve(
        rows, load_materials(GENERIC), coverages, long_axis="x", elsets=mesh.elsets
    )

    assert plies_from_rows(rows) == tuple(default_plies())
    assert elsets_csv == elsets_py
    assert from_csv.to_inp() == from_py.to_inp()


@pytest.mark.skipif(not FIN2.is_file(), reason="test_fin_2.step not in repo root")
def test_the_equivalence_book_names_every_zone_the_fin_mesh_has():
    from compfea.step_mesh import mesh_step

    mesh = mesh_step(FIN2, size_mm=40.0)
    check_against_mesh(load_plybook(EQUIV), mesh.elsets)


def test_stack_label_reads_bottom_to_top(tmp_path):
    rows = load_plybook(write(tmp_path, HEADER, *ROWS))
    layup, _ = resolve(
        rows, load_materials(GENERIC),
        {1: frozenset({"FULL"}), 2: frozenset({"FULL", "SPAR"})},
        long_axis="x",
    )
    labels = sorted(stack_label(z) for z in layup.zones)
    assert labels == ["[0/45/-45]", "[0/45/0/-45]"]
