"""STEP import with nested ACP coverage masks."""

from __future__ import annotations

from pathlib import Path

import pytest

from compfea.geometry import GeometryError
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


# 17.0 is the coarsest that recombines to all quad8; 17.25 already leaves
# 6-node triangles. 16.0 is the default because it keeps margin to that edge.
CLEAN_SIZE_MM = 16.0


def _boundary_loops(mesh) -> list[int]:
    """Sizes of the connected components of the free-edge graph.

    One component means the quad mesh is watertight -- its only free edges are
    the outline. Two means there is a hole somewhere inside the part.
    """
    from collections import defaultdict, deque

    used: dict[frozenset, int] = defaultdict(int)
    for conn in mesh.elements.values():
        corners = conn[:4]
        for i in range(4):
            used[frozenset((corners[i], corners[(i + 1) % 4]))] += 1
    adjacent: dict[int, set[int]] = defaultdict(set)
    for edge, count in used.items():
        if count == 1:
            a, b = tuple(edge)
            adjacent[a].add(b)
            adjacent[b].add(a)
    seen: set[int] = set()
    sizes: list[int] = []
    for start in adjacent:
        if start in seen:
            continue
        queue = deque([start])
        seen.add(start)
        size = 0
        while queue:
            node = queue.popleft()
            size += 1
            for other in adjacent[node]:
                if other not in seen:
                    seen.add(other)
                    queue.append(other)
        sizes.append(size)
    return sorted(sizes)


def test_a_mesh_with_triangles_is_refused_not_silently_holed():
    """Only quad8 is read back, so a stray triangle would vanish.

    At 40 mm -- the size the fin U-bend runs used -- one tile recombines to
    quad8 plus 6-node triangles. Dropping those leaves a hole that ccx solves
    without complaint, and no node is orphaned by it, so nothing downstream
    notices. It has to fail here.
    """
    with pytest.raises(GeometryError, match="non-quad8"):
        mesh_step(FIN2, size_mm=40.0)


def test_the_clean_mesh_is_watertight():
    mesh = mesh_step(FIN2, size_mm=CLEAN_SIZE_MM)
    loops = _boundary_loops(mesh)
    assert len(loops) == 1, f"interior hole: free-edge components {loops}"


def test_mesh_step_imprints_and_tags_nested_coverage():
    mesh = mesh_step(FIN2, size_mm=CLEAN_SIZE_MM)
    assert len(mesh.elements) > 100
    assert "blade" in mesh.elsets
    assert set(mesh.nsets) >= {"fixed_end", "far_face"}
    assert len(mesh.nsets["fixed_end"]) >= 2
    assert len(mesh.nsets["far_face"]) >= 2
    # Default clamp is every node under HEAL, not just the xmin edge.
    heal_nodes = {n for eid in mesh.elsets["HEAL"] for n in mesh.elements[eid]}
    assert set(mesh.nsets["fixed_end"]) == heal_nodes
    assert len(mesh.nsets["fixed_end"]) > len(mesh.nsets["far_face"])

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
