"""STEP import with nested ACP coverage masks."""

from __future__ import annotations

from pathlib import Path

import pytest

from compfea.step_mesh import mesh_step, shell_x_ranges_mm

ROOT = Path(__file__).resolve().parents[1]
FIN2 = ROOT / "test_fin_2.step"

pytestmark = pytest.mark.skipif(not FIN2.is_file(), reason="test_fin_2.step not in repo root")


def test_shell_x_ranges_named_products():
    ranges = shell_x_ranges_mm(FIN2)
    assert set(ranges) == {"FULL", "3_4ths", "HEAL", "QUARTER", "HALF", "TIP"}
    assert ranges["FULL"][0] == pytest.approx(0.0, abs=1e-3)
    assert ranges["FULL"][1] == pytest.approx(1174.9, rel=1e-3)
    assert ranges["HEAL"][1] < ranges["HALF"][1] <= ranges["3_4ths"][1] < ranges["FULL"][1]
    assert ranges["TIP"][0] == pytest.approx(ranges["HALF"][1], rel=1e-3)


def test_mesh_step_imprints_and_tags_nested_coverage():
    mesh = mesh_step(FIN2, size_mm=40.0)
    assert len(mesh.elements) > 100
    assert "blade" in mesh.elsets
    assert set(mesh.nsets) >= {"fixed_end", "far_face"}
    assert len(mesh.nsets["fixed_end"]) >= 2
    assert len(mesh.nsets["far_face"]) >= 2

    full = set(mesh.elsets["FULL"])
    assert full == set(mesh.elsets["blade"])
    assert set(mesh.elsets["HEAL"]).issubset(full)
    assert set(mesh.elsets["TIP"]).issubset(full)
    assert set(mesh.elsets["HALF"]).issubset(full)
    # 3_4ths starts with a digit -> sanitized
    assert "z_3_4ths" in mesh.elsets
    assert len(mesh.elsets["z_3_4ths"]) >= len(mesh.elsets["HALF"])

    # Fragment imprint: more than one unique element x-station near known drops
    xs = sorted(
        {
            round(sum(mesh.nodes[n][0] for n in conn[:4]) / 4.0, 0)
            for conn in mesh.elements.values()
        }
    )
    assert len(xs) >= 5

    inp = mesh.to_inp()
    assert "*ELEMENT, TYPE=S8R, ELSET=blade" in inp
    assert "*ELSET, ELSET=FULL" in inp
    assert "*NSET, NSET=far_face" in inp
