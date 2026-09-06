"""STEP import with nested ACP coverage masks."""

from __future__ import annotations

from pathlib import Path

import pytest

from compfea.geometry import GeometryError, check_watertight
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


def test_stray_triangles_are_kept_as_s6_not_dropped():
    """At 40 mm one tile recombines to quad8 plus two tri6.

    Dropping them leaves a hole that ccx solves without complaint, orphans no
    node, and is invisible in gmsh (which still holds the elements the deck
    lost). On the strip, deleting one interior element of 128 moved the
    reported force 2.5% against 0.8% of area. So they are read, as S6.
    """
    mesh = mesh_step(FIN2, size_mm=40.0)
    tris = [c for c in mesh.elements.values() if len(c) == 6]
    quads = [c for c in mesh.elements.values() if len(c) == 8]
    assert len(tris) == 2 and len(quads) == 335
    inp = mesh.to_inp()
    # both types, one ELSET, so a single COMPOSITE section covers them
    assert "*ELEMENT, TYPE=S8R, ELSET=blade" in inp
    assert "*ELEMENT, TYPE=S6, ELSET=blade" in inp


def test_a_triangle_dominated_mesh_is_refused():
    """S6 is a stiffer bending element; a few is a tolerance, mostly is not."""
    with pytest.raises(GeometryError, match="below the"):
        mesh_step(FIN2, size_mm=40.0, quad_floor=0.999)


def test_the_coarse_mesh_is_watertight_too():
    """The property, not the element type, is what keeps holes out."""
    check_watertight(mesh_step(FIN2, size_mm=40.0))


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
