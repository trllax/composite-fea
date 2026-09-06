"""ACP-style coverage plies -> exclusive COMPOSITE ELSETs."""

from __future__ import annotations

import pytest

from compfea.deck import assemble
from compfea.layup import (
    Ply,
    coverages_from_mesh,
    layup_from_coverage,
    mesh_elsets_for_stacks,
)


def test_ply_coverage_defaults_to_everywhere():
    assert Ply(0.1, 0.0).coverage is None


def test_coverages_from_mesh_inverts_overlapping_masks():
    elsets = {
        "blade": (1, 2, 3, 4),
        "FULL": (1, 2, 3, 4),
        "HALF": (1, 2),
        "TIP": (3, 4),
    }
    cov = coverages_from_mesh(elsets)
    assert cov[1] == frozenset({"FULL", "HALF"})
    assert cov[3] == frozenset({"FULL", "TIP"})


def test_layup_from_coverage_partitions_by_stack():
    # elems 1,2: FULL+HALF; 3,4: FULL+TIP
    element_coverages = {
        1: frozenset({"FULL", "HALF"}),
        2: frozenset({"FULL", "HALF"}),
        3: frozenset({"FULL", "TIP"}),
        4: frozenset({"FULL", "TIP"}),
    }
    plies = [
        Ply(0.2, 0.0, coverage="FULL"),
        Ply(0.2, 90.0, coverage="FULL"),
        Ply(0.1, 0.0, coverage="HALF"),
        Ply(0.15, 45.0, coverage="TIP"),
    ]
    layup, stacks = layup_from_coverage(plies, element_coverages, long_axis="x")
    assert set(stacks) == {"cov_1", "cov_2"}
    # first group is elems 1,2 (sorted eid order -> first stack key seen)
    assert set(stacks["cov_1"]) == {1, 2}
    assert set(stacks["cov_2"]) == {3, 4}
    # HALF stack: FULL plies + HALF ply (3 plies); TIP stack: FULL + TIP (3 plies)
    z1 = next(z for z in layup.zones if z.elset == "cov_1")
    z2 = next(z for z in layup.zones if z.elset == "cov_2")
    assert [p.angle_deg for p in z1.plies] == [0.0, 90.0, 0.0]
    assert [p.angle_deg for p in z2.plies] == [0.0, 90.0, 45.0]
    text = layup.to_inp()
    assert text.count("*SHELL SECTION, COMPOSITE") == 2
    assert "ELSET=cov_1" in text and "ELSET=cov_2" in text


def test_unknown_coverage_name_raises():
    with pytest.raises(ValueError, match="not in mesh masks"):
        layup_from_coverage(
            [Ply(0.1, 0.0, coverage="NOPE")],
            {1: frozenset({"FULL"})},
            long_axis="y",
        )


def test_mesh_elsets_for_stacks_partitions_blade():
    stacks = {"cov_1": (1, 2), "cov_2": (3,)}
    elsets = mesh_elsets_for_stacks(stacks, all_elements=(1, 2, 3))
    assert elsets["blade"] == (1, 2, 3)
    assert set(elsets["cov_1"]) | set(elsets["cov_2"]) == {1, 2, 3}


def test_coverage_deck_assembles():
    element_coverages = {
        1: frozenset({"FULL"}),
        2: frozenset({"FULL", "TIP"}),
    }
    layup, stacks = layup_from_coverage(
        [Ply(0.1, 0.0, coverage="FULL"), Ply(0.1, 90.0, coverage="TIP")],
        element_coverages,
        long_axis="x",
    )
    mesh_inp = (
        "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n"
        "*ELEMENT, TYPE=S8R, ELSET=blade\n"
        "1, 1, 1, 1, 1, 1, 1, 1, 1\n"
        "2, 2, 2, 2, 2, 2, 2, 2, 2\n"
    )
    for name, eids in stacks.items():
        mesh_inp += f"*ELSET, ELSET={name}\n" + ", ".join(str(e) for e in eids) + "\n"
    deck = assemble(mesh_inp=mesh_inp, layup=layup, steps=[])
    assert deck.count("*SHELL SECTION, COMPOSITE") == 2
