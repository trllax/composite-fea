"""STEP import with nested ACP coverage masks."""

from __future__ import annotations

from pathlib import Path

import pytest

from compfea.geometry import GeometryError, check_watertight
from compfea.step_mesh import (
    mesh_step,
    shell_faces,
    shell_hulls_mm,
    shell_x_ranges_mm,
    zone_report,
)

ROOT = Path(__file__).resolve().parents[1]
FIN2 = ROOT / "test_fin_2.step"

pytestmark = pytest.mark.skipif(
    not FIN2.is_file(), reason="test_fin_2.step not in repo root"
)


NAMES = {"FULL", "3_4ths", "HEAL", "QUARTER", "HALF", "TIP"}

#: OCC face areas of each named shell's tiles, mm^2. Measured once from
#: occ.getMass on test_fin_2.step; these are CAD numbers, not mesh numbers.
CAD_AREA_MM2 = {
    "FULL": 336223.808,
    "z_3_4ths": 193746.626,
    "HALF": 118746.626,
    "QUARTER": 73746.626,
    "TIP": 60000.000,
    "HEAL": 37502.639,
}


def test_shell_x_ranges_are_a_hull_not_a_mask():
    """The control-point hull overhangs the real face, and by a lot.

    This function is kept because it is how the names are read, but membership
    must never come from it: HALF and QUARTER have the *same* hull x-range here
    while their true extents are 450 and 300 mm, and TIP's hull starts 300 mm
    inboard of the real TIP face. Asserting the overhang keeps the distinction
    from quietly being forgotten again.
    """
    ranges = shell_x_ranges_mm(FIN2)
    assert set(ranges) == NAMES
    assert ranges["FULL"][1] == pytest.approx(1174.9, rel=1e-3)
    # The bug the fragment map replaced: these two hulls are indistinguishable.
    assert ranges["HALF"] == pytest.approx(ranges["QUARTER"])
    assert ranges["TIP"][0] == pytest.approx(ranges["HALF"][1], rel=1e-3)


def test_shell_faces_are_read_in_step_order():
    """21 faces over 6 shells, contiguous ids -- what the mapping rests on."""
    faces = shell_faces(FIN2)
    assert set(faces) == NAMES
    assert {n: len(f) for n, f in faces.items()} == {
        "FULL": 7, "3_4ths": 5, "HEAL": 1, "QUARTER": 3, "HALF": 4, "TIP": 1
    }
    flat = sorted(f for ids in faces.values() for f in ids)
    assert len(flat) == 21
    assert flat == list(range(flat[0], flat[0] + 21)), "face ids must be contiguous"
    # Each shell owns a contiguous block, root-of-the-file first.
    for ids in faces.values():
        assert list(ids) == list(range(ids[0], ids[0] + len(ids)))


def test_shell_hulls_contain_every_true_zone():
    """The containment the alignment check relies on, stated directly."""
    hulls = shell_hulls_mm(FIN2)
    mesh = mesh_step(FIN2, size_mm=CLEAN_SIZE_MM)
    sane = {"3_4ths": "z_3_4ths"}
    for name, hull in hulls.items():
        eids = mesh.elsets[sane.get(name, name)]
        pts = [mesh.nodes[n] for e in eids for n in mesh.elements[e]]
        for axis in (0, 1, 2):
            assert min(p[axis] for p in pts) >= hull[axis] - 1e-6
            assert max(p[axis] for p in pts) <= hull[axis + 3] + 1e-6


def test_zones_match_their_cad_faces_not_a_bounding_box():
    """The regression this phase exists for.

    Zone membership used to be "element centroid x inside the shell's hull".
    That gave HALF and QUARTER *identical* ELSETs and made TIP 2.5x too large.
    Membership now follows the OCC fragment map, so each zone's meshed area
    must reproduce the CAD area of the faces it came from.
    """
    mesh = mesh_step(FIN2, size_mm=CLEAN_SIZE_MM)
    rows = {r["zone"]: r for r in zone_report(mesh)}
    assert set(rows) == set(CAD_AREA_MM2) | {"blade"}

    for zone, cad in CAD_AREA_MM2.items():
        got = rows[zone]["area_mm2"]
        # Chordal error under-measures a curved band; worst here is 0.007%.
        assert got == pytest.approx(cad, rel=1e-3), f"{zone}: {got} vs CAD {cad}"

    # The two that used to collapse onto one mask are now distinct...
    assert set(mesh.elsets["HALF"]) != set(mesh.elsets["QUARTER"])
    assert set(mesh.elsets["QUARTER"]) < set(mesh.elsets["HALF"])
    # ...and their true extents differ by 150 mm, which the hull could not see.
    assert rows["HALF"]["x_max_mm"] == pytest.approx(450.0, abs=1e-6)
    assert rows["QUARTER"]["x_max_mm"] == pytest.approx(300.0, abs=1e-6)
    # TIP starts at its real face, 300 mm outboard of where the hull put it.
    assert rows["TIP"]["x_min_mm"] == pytest.approx(974.92, rel=1e-4)


def test_zones_nest_the_way_the_cad_does():
    mesh = mesh_step(FIN2, size_mm=CLEAN_SIZE_MM)
    z = {n: set(mesh.elsets[n]) for n in mesh.elsets}
    assert z["HEAL"] < z["QUARTER"] < z["HALF"] < z["z_3_4ths"] < z["FULL"]
    assert z["FULL"] == z["blade"]
    # TIP is outboard of everything else and shares no element with 3_4ths.
    assert not (z["TIP"] & z["z_3_4ths"])


def test_coverage_tol_mm_is_dead_and_says_so():
    with pytest.deprecated_call():
        mesh_step(FIN2, size_mm=40.0, coverage_tol_mm=2.0)


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
    assert len(mesh.elsets["z_3_4ths"]) > len(mesh.elsets["HALF"])

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


# --------------------------------------------------------------------------
# the alignment guard itself
#
# Everything above pins the *outcome* on this STEP, which a hardcoded
# CAD_AREA_MM2 would catch even with the guard removed. These exercise the
# guard, i.e. the behaviour that protects a STEP nobody has tried yet, or this
# one after a gmsh upgrade. Deleting _check_alignment used to leave the whole
# suite green.


def alignment_inputs():
    """Everything _check_alignment needs, taken from the real STEP."""
    import gmsh

    from compfea.step_mesh import _sanitize_elset, _tiles_by_shell, face_hulls_mm

    faces = shell_faces(FIN2)
    by_id = sorted((f, n) for n, ids in faces.items() for f in ids)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.model.add("align")
        gmsh.option.setNumber("Geometry.OCCImportLabels", 1)
        gmsh.model.occ.importShapes(str(FIN2.resolve()))
        gmsh.model.occ.synchronize()
        surfs = gmsh.model.getEntities(2)
        tags = [t for _d, t in surfs]
        labels = {
            t: gmsh.model.getEntityName(2, t).rsplit("/", 1)[-1] for t in tags
        }
        _out, out_map = gmsh.model.occ.fragment([surfs[0]], list(surfs[1:]))
        gmsh.model.occ.synchronize()
        boxes = {
            t: gmsh.model.getBoundingBox(2, t) for _d, t in gmsh.model.getEntities(2)
        }
    finally:
        gmsh.finalize()
    names = [n for _f, n in by_id]
    return {
        "face_ids": [f for f, _n in by_id],
        "face_hulls": face_hulls_mm(FIN2),
        "tile_bbox": boxes,
        "labels": labels,
        "face_order": names,
        "surf_tags": tags,
        "out_map": out_map,
        "tiles_by_shell": _tiles_by_shell(
            [_sanitize_elset(n) for n in names], out_map
        ),
    }


def test_check_alignment_accepts_the_step_face_order():
    from compfea.step_mesh import _check_alignment

    _check_alignment(**alignment_inputs())   # must not raise


@pytest.mark.parametrize(
    "mutate,label",
    [
        (lambda f: f[::-1], "reversed"),
        (lambda f: f[1:] + f[:1], "rotated by one"),
        # HEAL's only face against TIP's only face: equal face counts, so no
        # count check could see it.
        (lambda f: f[:12] + [f[20]] + f[13:20] + [f[12]], "HEAL/TIP faces"),
        # A HALF face against a QUARTER face. Those two shells have
        # byte-identical control-point hulls, so a per-shell test is blind
        # here -- and this is the pair the old x-range mask collapsed.
        (lambda f: f[:13] + [f[16]] + f[14:16] + [f[13]] + f[17:], "HALF/QUARTER"),
    ],
)
def test_check_alignment_refuses_a_misordered_import(mutate, label):
    from compfea.step_mesh import _check_alignment

    inputs = alignment_inputs()
    inputs["face_ids"] = mutate(inputs["face_ids"])
    with pytest.raises(GeometryError, match="control-point hull|not within"):
        _check_alignment(**inputs)


def test_half_and_quarter_share_a_hull_so_per_shell_containment_is_blind():
    """Why the guard is per face. If this ever stops being true, say so."""
    hulls = shell_hulls_mm(FIN2)
    assert hulls["HALF"] == hulls["QUARTER"]
    faces = shell_faces(FIN2)
    from compfea.step_mesh import face_hulls_mm

    fh = face_hulls_mm(FIN2)
    assert {fh[f] for f in faces["HALF"]} != {fh[f] for f in faces["QUARTER"]}
