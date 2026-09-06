"""Unit tests for geometry.py, plus two ccx round trips for the curvature map.

The physical validation of this module is `cases/smoke_cantilever`, which now
builds its mesh here: hand CLPT on two stacks, the stiffness ratio, the
closed-form elastica and the reversed-stack draw-in all run on a geometry.py
mesh. What is left for this file is the things a reaction force cannot see --
element type, node ordering, normal direction, set membership, reproducibility --
and the arc-length property of the camber map, which a reaction force *can* see
and which gets a solve of its own.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from compfea.deck import StaticStep, assemble
from compfea.geometry import (
    check_watertight,
    CircularCamber,
    CurvatureProfile,
    GeometryError,
    Mesh,
    Outline,
    mesh_outline,
)
from compfea.layup import Layup, Ply
from compfea.run import solve

CHORD = 20.0
SPAN = 100.0
N_CHORD = 1
N_SPAN = 32


@pytest.fixture(scope="module")
def strip():
    return mesh_outline(
        Outline.rectangle(CHORD, SPAN), n_chord=N_CHORD, n_span=N_SPAN
    )


# --------------------------------------------------------------------------
# Element type, node ordering, normals
# --------------------------------------------------------------------------


def test_every_element_is_an_eight_node_quad(strip):
    assert len(strip.elements) == N_CHORD * N_SPAN
    assert all(len(conn) == 8 for conn in strip.elements.values())
    assert all(len(set(conn)) == 8 for conn in strip.elements.values())
    assert "*ELEMENT, TYPE=S8R, ELSET=blade" in strip.to_inp()


def test_gmsh_node_order_is_calculix_s8r_order(strip):
    """Mid-side node k sits between corners k and k+1.

    gmsh and CalculiX happen to agree, so connectivity passes straight through
    with no permutation. That is an assumption worth failing loudly on rather
    than discovering through a deck that solves and is wrong.
    """
    for conn in strip.elements.values():
        corners = [np.array(strip.nodes[n]) for n in conn[:4]]
        midsides = [np.array(strip.nodes[n]) for n in conn[4:]]
        for k, mid in enumerate(midsides):
            expected = (corners[k] + corners[(k + 1) % 4]) / 2.0
            assert mid == pytest.approx(expected, abs=1e-9)


def test_every_element_normal_points_at_plus_z(strip):
    """The normal decides which face is +z, and so which ply the first section
    line lands on. A -z normal inverts every unsymmetric stack in the model."""
    for normal in strip.element_normals().values():
        assert normal == pytest.approx((0.0, 0.0, 1.0), abs=1e-12)


def test_a_clockwise_outline_is_reoriented_rather_than_inverted():
    """The same strip, listed the other way round, must still mesh +z.

    root and tip keep their identity through the flip -- they are what fixed_end
    and far_face are built from -- so the reoriented outline still clamps the
    y = 0 edge.
    """
    clockwise = Outline(
        root=((CHORD, 0.0), (0.0, 0.0)),
        leading=((0.0, 0.0), (0.0, SPAN)),
        tip=((0.0, SPAN), (CHORD, SPAN)),
        trailing=((CHORD, SPAN), (CHORD, 0.0)),
    )
    assert clockwise.signed_area() > 0.0  # normalised on construction
    assert all(y == 0.0 for _, y in clockwise.root)
    assert all(y == SPAN for _, y in clockwise.tip)
    mesh = mesh_outline(clockwise, n_chord=N_CHORD, n_span=N_SPAN)
    for normal in mesh.element_normals().values():
        assert normal == pytest.approx((0.0, 0.0, 1.0), abs=1e-12)


def test_a_mixed_element_mesh_is_refused(monkeypatch):
    """A stray triangle would be dropped, not written.

    Only quad8 elements are read back, so a mixed mesh loses its triangles and
    leaves a hole in the part with nodes attached to nothing -- which ccx solves.
    Recombination has produced pure quads on every outline tried here, so the
    guard is exercised by faking gmsh's answer rather than by hunting for an
    outline that defeats the recombiner.
    """
    import compfea.geometry as geometry

    real = geometry.gmsh.model.mesh.getElements

    def fake(*args, **kwargs):
        types, tags, nodes = real(*args, **kwargs)
        return ([*list(types), 9], list(tags) + [[]], list(nodes) + [[]])

    monkeypatch.setattr(geometry.gmsh.model.mesh, "getElements", fake)
    with pytest.raises(GeometryError, match="element types"):
        mesh_outline(Outline.rectangle(CHORD, SPAN), n_chord=1, n_span=4)


def test_the_normal_check_fires_when_an_outline_slips_through_clockwise():
    """The guard itself, exercised directly.

    Outline normalisation means no ordinary path can reach _check_normals with a
    clockwise loop, so the guard was previously unpinned -- making it `return`
    left the whole suite green. Here the outline is flipped back to clockwise
    after construction, which is what a bug in the normalisation would do.
    """
    outline = Outline.rectangle(CHORD, SPAN)
    flipped = {
        "root": tuple(reversed(outline.root)),
        "leading": tuple(reversed(outline.trailing)),
        "tip": tuple(reversed(outline.tip)),
        "trailing": tuple(reversed(outline.leading)),
    }
    for role, edge in flipped.items():
        object.__setattr__(outline, role, edge)
    assert outline.signed_area() < 0.0
    with pytest.raises(GeometryError, match="faces away from the laminate"):
        mesh_outline(outline, n_chord=1, n_span=4)


def test_a_clockwise_outline_splines_the_edge_it_was_told_to():
    """`splines` names roles, and reorientation moves the edges between roles.

    The same physical blade listed both ways must mesh to the same part. Leaving
    `splines` unmapped chords the curved edge and bulges the straight one, which
    changes the area by over 1% with no error anywhere.
    """
    curve = ((24.0, 0.0), (26.0, 30.0), (22.0, 70.0), (10.0, 100.0))
    counter_clockwise = Outline(
        root=((0.0, 0.0), (24.0, 0.0)),
        leading=curve,
        tip=((10.0, 100.0), (0.0, 100.0)),
        trailing=((0.0, 100.0), (0.0, 0.0)),
        splines=frozenset({"leading"}),
    )
    clockwise = Outline(
        root=((24.0, 0.0), (0.0, 0.0)),
        leading=((0.0, 0.0), (0.0, 100.0)),
        tip=((0.0, 100.0), (10.0, 100.0)),
        trailing=tuple(reversed(curve)),
        splines=frozenset({"trailing"}),
    )
    areas = [
        sum(mesh_outline(o, size=6.0).element_areas().values())
        for o in (counter_clockwise, clockwise)
    ]
    assert areas[1] == pytest.approx(areas[0], rel=1e-9)


def test_a_splined_tip_keeps_every_node_of_its_edge():
    """far_face is the driven set in every case deck.

    Finding edge nodes by distance to the control points drops exactly the nodes
    that bulge away from the chord -- and never reports an empty set, because the
    control points themselves always match. That turns an edge drive into a
    four-point drive.
    """
    outline = Outline(
        root=((0.0, 0.0), (30.0, 0.0)),
        leading=((30.0, 0.0), (30.0, 100.0)),
        tip=((30.0, 100.0), (20.0, 106.0), (10.0, 108.0), (0.0, 100.0)),
        trailing=((0.0, 100.0), (0.0, 0.0)),
        splines=frozenset({"tip"}),
    )
    mesh = mesh_outline(outline, size=6.0)
    far = set(mesh.nsets["far_face"])
    assert len(far) > len(outline.tip)
    # Every node of every element edge that lies wholly in far_face is in it:
    # the set is a complete boundary, not a sample of one.
    boundary_ys = [mesh.nodes[n][1] for n in far]
    assert min(boundary_ys) == pytest.approx(100.0, abs=1e-9)
    # The spline overshoots its control points -- 108.29 against a highest
    # control point of 108.0 -- which is exactly why measuring distance to the
    # control polyline dropped these nodes.
    assert max(boundary_ys) > 108.0


def test_an_outline_that_does_not_start_at_the_origin_still_zones_correctly():
    """Spanwise position is measured from the root station, not from y = 0."""
    offset = Outline(
        root=((0.0, 20.0), (CHORD, 20.0)),
        leading=((CHORD, 20.0), (CHORD, 120.0)),
        tip=((CHORD, 120.0), (0.0, 120.0)),
        trailing=((0.0, 120.0), (0.0, 20.0)),
    )
    mesh = mesh_outline(offset, n_chord=1, n_span=10, zones=(0.5,))
    assert len(mesh.elsets["zone_1"]) == len(mesh.elsets["zone_2"]) == 5
    # And the camber map must not run off the end of its integration grid, which
    # would clamp every outboard node onto the tip and collapse whole elements.
    cambered = mesh_outline(offset, n_chord=1, n_span=10, camber=CircularCamber(200.0))
    assert sum(cambered.element_areas().values()) == pytest.approx(
        CHORD * SPAN, rel=1e-3
    )


def test_the_root_must_be_the_inboard_end():
    """A blade listed with its root outboard would zone and camber backwards."""
    with pytest.raises(GeometryError, match="root edge must reach the inboard"):
        Outline(
            root=((CHORD, SPAN), (0.0, SPAN)),
            leading=((0.0, SPAN), (0.0, 0.0)),
            tip=((0.0, 0.0), (CHORD, 0.0)),
            trailing=((CHORD, 0.0), (CHORD, SPAN)),
        )


def test_a_zone_thinner_than_one_element_row_is_refused():
    """An empty ELSET writes a bare *ELSET card and ccx solves the deck anyway,
    with the ply drop simply missing and the reaction looking fine."""
    with pytest.raises(GeometryError, match="caught no elements"):
        mesh_outline(
            Outline.rectangle(CHORD, SPAN), n_chord=1, n_span=8, zones=(0.0625,)
        )


def test_a_size_coarser_than_the_outline_is_refused_with_advice():
    """gmsh raises a bare Exception here; a narrow tip chord is the usual cause."""
    with pytest.raises(GeometryError, match="could not mesh this outline"):
        mesh_outline(Outline.tapered(40.0, 0.6, SPAN), size=5.0)


# --------------------------------------------------------------------------
# The grid, areas, node sets
# --------------------------------------------------------------------------


def test_structured_strip_is_the_documented_grid(strip):
    """Corner rows carry 2*n_chord+1 nodes, mid-side rows n_chord+1."""
    corner_rows = N_SPAN + 1
    mid_rows = N_SPAN
    expected = corner_rows * (2 * N_CHORD + 1) + mid_rows * (N_CHORD + 1)
    assert len(strip.nodes) == expected == 163
    ys = sorted({round(y, 9) for _, y, _ in strip.nodes.values()})
    assert len(ys) == corner_rows + mid_rows
    assert ys[0] == 0.0 and ys[-1] == SPAN


@pytest.mark.parametrize(
    "outline,expected",
    [
        (Outline.rectangle(CHORD, SPAN), CHORD * SPAN),
        (Outline.tapered(30.0, 10.0, SPAN), 0.5 * (30.0 + 10.0) * SPAN),
        (Outline.tapered(30.0, 10.0, SPAN, sweep=8.0), 0.5 * (30.0 + 10.0) * SPAN),
    ],
)
def test_meshed_area_matches_the_analytic_planform(outline, expected):
    """Element areas must sum to the shape asked for -- a mesher that quietly
    trims or overshoots the outline changes the part."""
    mesh = mesh_outline(outline, n_chord=4, n_span=16)
    assert sum(mesh.element_areas().values()) == pytest.approx(expected, rel=1e-9)


def test_a_spline_edge_meshes_and_keeps_its_area():
    """Curved leading edge, area checked against a densely sampled polygon."""
    control = ((24.0, 0.0), (26.0, 30.0), (24.0, 65.0), (18.0, 95.0), (10.0, 110.0))
    outline = Outline(
        root=((0.0, 0.0), (24.0, 0.0)),
        leading=control,
        tip=((10.0, 110.0), (0.0, 110.0)),
        trailing=((0.0, 110.0), (0.0, 0.0)),
        splines=frozenset({"leading"}),
    )
    mesh = mesh_outline(outline, size=5.0)
    assert len(mesh.elements) > 50
    for normal in mesh.element_normals().values():
        assert normal[2] > 0.0
    # The spline bulges outside the control polygon, so the polygon area is a
    # lower bound; the two agree to a few percent on this shape.
    polygon_area = outline.signed_area()
    assert sum(mesh.element_areas().values()) == pytest.approx(polygon_area, rel=0.05)


def test_node_sets_hold_exactly_the_root_and_tip_edges(strip):
    fixed = strip.nsets["fixed_end"]
    far = strip.nsets["far_face"]
    assert len(fixed) == len(far) == 2 * N_CHORD + 1
    assert all(strip.nodes[n][1] == pytest.approx(0.0, abs=1e-9) for n in fixed)
    assert all(strip.nodes[n][1] == pytest.approx(SPAN, abs=1e-9) for n in far)
    assert not set(fixed) & set(far)


def test_the_inp_declares_the_sets_the_decks_reference(strip):
    text = strip.to_inp()
    for card in (
        "*NODE, NSET=all_nodes",
        "*ELEMENT, TYPE=S8R, ELSET=blade",
        "*NSET, NSET=fixed_end",
        "*NSET, NSET=far_face",
    ):
        assert card in text
    # ccx cannot parse a comment that follows data on the same line.
    for line in text.splitlines():
        assert "**" not in line or line.startswith("**")


# --------------------------------------------------------------------------
# Zones
# --------------------------------------------------------------------------


def test_zones_partition_every_element_exactly_once():
    mesh = mesh_outline(
        Outline.tapered(30.0, 10.0, SPAN), n_chord=2, n_span=20, zones=(0.4, 0.7)
    )
    zone_names = [n for n in mesh.elsets if n.startswith("zone_")]
    assert zone_names == ["zone_1", "zone_2", "zone_3"]
    members = [set(mesh.elsets[n]) for n in zone_names]
    assert set().union(*members) == set(mesh.elements)
    assert sum(len(m) for m in members) == len(mesh.elements)


def test_zone_boundaries_land_where_asked():
    mesh = mesh_outline(
        Outline.rectangle(CHORD, SPAN), n_chord=1, n_span=32, zones=(0.4, 0.7)
    )

    def band(name):
        ys = [
            sum(mesh.nodes[n][1] for n in mesh.elements[e][:4]) / 4.0
            for e in mesh.elsets[name]
        ]
        return min(ys) / SPAN, max(ys) / SPAN

    assert band("zone_1")[0] < 0.02 and band("zone_1")[1] < 0.4
    assert band("zone_2")[0] >= 0.4 and band("zone_2")[1] < 0.7
    assert band("zone_3")[0] >= 0.7 and band("zone_3")[1] > 0.98


def test_zone_fractions_are_validated():
    for bad in ((0.7, 0.4), (0.0, 0.5), (0.5, 1.0), (-0.1,), (0.4, 0.4)):
        with pytest.raises(GeometryError, match="zone fractions"):
            mesh_outline(
                Outline.rectangle(CHORD, SPAN), n_chord=1, n_span=8, zones=bad
            )


# --------------------------------------------------------------------------
# Reproducibility and argument validation
# --------------------------------------------------------------------------


def test_two_identical_calls_produce_identical_text():
    """A sweep caches on the design vector, so the same vector must give the
    same deck.

    On this gmsh build the raw tags happen to be stable too, so the deterministic
    renumbering in `_generate` is defence in depth rather than the thing making
    this pass -- verified by swapping the sort key for the raw tags, which leaves
    the suite green.
    """
    kwargs = dict(n_chord=2, n_span=9, zones=(0.5,))
    first = mesh_outline(Outline.tapered(30.0, 10.0, SPAN), **kwargs).to_inp()
    second = mesh_outline(Outline.tapered(30.0, 10.0, SPAN), **kwargs).to_inp()
    assert first == second


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"n_chord": 4},
        {"n_span": 4},
        {"n_chord": 0, "n_span": 4},
        {"size": 0.0},
        {"size": float("inf")},
    ],
)
def test_bad_mesh_requests_are_refused(kwargs):
    with pytest.raises(GeometryError):
        mesh_outline(Outline.rectangle(CHORD, SPAN), **kwargs)


def test_bad_outlines_are_refused():
    with pytest.raises(GeometryError, match="join end to end"):
        Outline(
            root=((0.0, 0.0), (CHORD, 0.0)),
            leading=((CHORD, 1.0), (CHORD, SPAN)),  # gap
            tip=((CHORD, SPAN), (0.0, SPAN)),
            trailing=((0.0, SPAN), (0.0, 0.0)),
        )
    with pytest.raises(GeometryError, match="non-finite"):
        Outline.rectangle(float("nan"), SPAN)
    with pytest.raises(GeometryError, match="at least two points"):
        Outline(
            root=((0.0, 0.0),),
            leading=((0.0, 0.0), (0.0, SPAN)),
            tip=((0.0, SPAN), (CHORD, SPAN)),
            trailing=((CHORD, SPAN), (0.0, 0.0)),
        )
    with pytest.raises(GeometryError, match="unknown spline"):
        Outline(
            root=((0.0, 0.0), (CHORD, 0.0)),
            leading=((CHORD, 0.0), (CHORD, SPAN)),
            tip=((CHORD, SPAN), (0.0, SPAN)),
            trailing=((0.0, SPAN), (0.0, 0.0)),
            splines=frozenset({"nose"}),
        )


# --------------------------------------------------------------------------
# Curvature
# --------------------------------------------------------------------------


def test_circular_camber_puts_every_node_on_the_arc():
    radius = 200.0
    mesh = mesh_outline(
        Outline.rectangle(CHORD, SPAN),
        n_chord=1,
        n_span=32,
        camber=CircularCamber(radius),
    )
    pts = np.array(list(mesh.nodes.values()))
    # Centre of curvature is (y, z) = (0, R).
    assert np.hypot(pts[:, 1], pts[:, 2] - radius) == pytest.approx(radius, abs=1e-5)
    assert pts[:, 0].min() == 0.0 and pts[:, 0].max() == CHORD  # chord untouched


def test_camber_preserves_developed_span():
    """The map is an isometry: a node at flat y sits at arc length y on the arc.

    This is the property that keeps the blade the length it was laid up at. Get
    it wrong -- lift with z = f(y) instead -- and the span grows with the bend,
    changing the stiffness as the cube of the length.
    """
    radius = 200.0
    mesh = mesh_outline(
        Outline.rectangle(CHORD, SPAN),
        n_chord=1,
        n_span=32,
        camber=CircularCamber(radius),
    )
    flat = mesh_outline(Outline.rectangle(CHORD, SPAN), n_chord=1, n_span=32)
    for nid, (_, y, z) in mesh.nodes.items():
        arc_length = radius * math.atan2(y, radius - z)
        assert arc_length == pytest.approx(flat.nodes[nid][1], abs=1e-5)


def test_circular_camber_tip_matches_the_closed_form():
    radius = 200.0
    mesh = mesh_outline(
        Outline.rectangle(CHORD, SPAN),
        n_chord=1,
        n_span=32,
        camber=CircularCamber(radius),
    )
    tip = [mesh.nodes[n] for n in mesh.nsets["far_face"]]
    angle = SPAN / radius
    for _, y, z in tip:
        assert y == pytest.approx(radius * math.sin(angle), abs=1e-6)
        assert z == pytest.approx(radius * (1.0 - math.cos(angle)), abs=1e-6)


def test_a_curvature_profile_integrates_to_the_same_shape_as_a_constant_one():
    """A flat kappa profile is a circular arc, by a different code path."""
    radius = 250.0
    a = mesh_outline(
        Outline.rectangle(CHORD, SPAN),
        n_chord=1,
        n_span=16,
        camber=CircularCamber(radius),
    )
    b = mesh_outline(
        Outline.rectangle(CHORD, SPAN),
        n_chord=1,
        n_span=16,
        camber=CurvatureProfile(((0.0, 1.0 / radius), (1.0, 1.0 / radius))),
    )
    for nid in a.nodes:
        assert b.nodes[nid] == pytest.approx(a.nodes[nid], abs=1e-6)


def test_camber_arguments_are_validated():
    with pytest.raises(GeometryError, match="finite and nonzero"):
        CircularCamber(0.0)
    with pytest.raises(GeometryError, match="at least two points"):
        CurvatureProfile(((0.0, 0.01),))
    with pytest.raises(GeometryError, match="increasing"):
        CurvatureProfile(((0.6, 0.01), (0.2, 0.02)))
    with pytest.raises(GeometryError, match="lie in"):
        CurvatureProfile(((0.0, 0.01), (1.4, 0.02)))
    with pytest.raises(GeometryError, match="finite"):
        CurvatureProfile(((0.0, 0.01), (1.0, float("nan"))))


# --------------------------------------------------------------------------
# Curvature, through ccx
# --------------------------------------------------------------------------


def _strip_deck(camber, delta_mm: float = 1.0) -> str:
    mesh = mesh_outline(
        Outline.rectangle(CHORD, SPAN), n_chord=1, n_span=32, camber=camber
    )
    layup = Layup.uniform(
        [Ply(0.25, angle) for angle in (0.0, 90.0, 90.0, 0.0)], long_axis="y"
    )
    body = "\n".join(
        [
            "*BOUNDARY",
            f"far_face, 3, 3, {delta_mm:.10f}",
            "*NODE PRINT, NSET=fixed_end, TOTALS=YES",
            "RF",
        ]
    )
    return assemble(
        mesh_inp=mesh.to_inp(),
        layup=layup,
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=[StaticStep(body, inc=2000, static_line="0.25, 1.0, 1.E-6, 1.0")],
    )


@pytest.fixture(scope="module")
def curved_reactions(tmp_path_factory):
    root = tmp_path_factory.mktemp("geometry_camber")
    out = {}
    for name, camber in (
        ("flat", None),
        ("near_flat", CircularCamber(1.0e7)),
        ("curved", CircularCamber(200.0)),
    ):
        run_dir = root / name
        run_dir.mkdir()
        deck = run_dir / "deck.inp"
        deck.write_text(_strip_deck(camber))
        out[name] = abs(solve(deck, run_dir).total_force("fixed_end", "fz"))
    return out


def test_a_nearly_flat_blade_returns_the_flat_answer(curved_reactions):
    """R = 1e7 mm over a 100 mm span is 1e-5 rad of bend: physically flat.

    Catches a camber map that misbehaves in the small-curvature limit -- a scale
    factor, a sign, a unit -- where the answer must converge back on the flat
    strip. It does *not* catch a map that stretches the span: verified by
    sabotage, replacing the isometry with z = f(y) leaves this test green, since
    at 1e-5 rad the two maps agree. Arc-length preservation is pinned
    geometrically, by test_camber_preserves_developed_span.
    """
    assert curved_reactions["near_flat"] == pytest.approx(
        curved_reactions["flat"], rel=1e-3
    )


def test_a_curved_blade_solves_and_is_not_the_flat_one(curved_reactions):
    """R = 200 mm is a 28.6 degree arc. Measured 12.8% stiffer than flat.

    The point is not the number -- there is no closed form for it -- but that
    curvature reaches the solver at all. A camber dropped between the mesh and
    the deck lands back on the flat answer; a camber that stretches the span
    instead of bending it also moves this ratio, and does fail here.
    """
    ratio = curved_reactions["curved"] / curved_reactions["flat"]
    assert ratio > 1.05


# --- holes ------------------------------------------------------------------
# A hole reaches the deck whenever elements go missing between gmsh and the
# .inp: a tile that recombined to quad8 plus triangles (only quad8 is read
# back), a tile that meshed to nothing, or a seam whose nodes are not shared.
# Checking those causes one at a time is how the first one survived. These pin
# the property instead.


def _interior_elements(mesh):
    from collections import Counter

    shared = Counter()
    for conn in mesh.elements.values():
        q = conn[:4]
        for i in range(4):
            shared[frozenset((q[i], q[(i + 1) % 4]))] += 1
    return [
        eid
        for eid, conn in mesh.elements.items()
        if all(
            shared[frozenset((conn[:4][i], conn[:4][(i + 1) % 4]))] == 2
            for i in range(4)
        )
    ]


def _strip():
    return mesh_outline(Outline.rectangle(chord=20.0, span=100.0), n_chord=4, n_span=8)


def test_a_sound_mesh_is_watertight():
    check_watertight(_strip())


def test_a_missing_interior_element_is_caught():
    mesh = _strip()
    victim = _interior_elements(mesh)[0]
    holed = Mesh(
        nodes=mesh.nodes,
        elements={k: v for k, v in mesh.elements.items() if k != victim},
        nsets=mesh.nsets,
        elsets=mesh.elsets,
        heading=mesh.heading,
    )
    with pytest.raises(GeometryError, match="hole"):
        check_watertight(holed)


def test_counting_orphan_nodes_would_not_have_caught_it():
    """Why the check is on free edges and not on nodes.

    An interior element shares every one of its nodes with a neighbour, so
    deleting it orphans nothing. A node-based check reports a clean mesh.
    """
    mesh = _strip()
    victim = _interior_elements(mesh)[0]
    remaining = {k: v for k, v in mesh.elements.items() if k != victim}
    used = {n for conn in remaining.values() for n in conn}
    assert set(mesh.nodes) - used == set()


def test_a_non_manifold_edge_is_caught():
    mesh = _strip()
    doubled = dict(mesh.elements)
    doubled[max(doubled) + 1] = mesh.elements[_interior_elements(mesh)[0]]
    with pytest.raises(GeometryError, match="more than two elements"):
        check_watertight(
            Mesh(
                nodes=mesh.nodes,
                elements=doubled,
                nsets=mesh.nsets,
                elsets=mesh.elsets,
                heading=mesh.heading,
            )
        )
