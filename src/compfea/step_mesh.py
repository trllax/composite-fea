"""STEP mid-surface import with nested ACP-style coverage masks.

Onshape/ANSYS-style overlapping named shells (FULL, HALF, TIP, ...) are
coverage masks, not exclusive tiles. This module:

1. Imports all shells via OpenCASCADE
2. ``fragment``s them so ply-drop boundaries become mesh edges (vertices)
3. Meshes the unique tiles as second-order incomplete quads (S8R)
4. Tags each element with every shell whose faces became its tile (exact)
5. Builds root/tip NSETs from the min/max-x boundary nodes

Units: STEP SI metres are scaled to mm to match the rest of the repo.
"""

from __future__ import annotations

import math
import re
import warnings
from collections import defaultdict
from pathlib import Path

import gmsh
import numpy as np

from compfea.geometry import (
    _GMSH_QUAD8,
    GeometryError,
    Mesh,
    _snap,
    check_quad_fraction,
    check_watertight,
)

_GMSH_TRI6 = 9

# Onshape STEP AP242 length unit is metre; CalculiX decks here are mm.
_STEP_M_TO_MM = 1000.0
_SHELL_MODEL = "SHELL_BASED_SURFACE_MODEL"
_CARTESIAN_POINT = "CARTESIAN_POINT"

#: gmsh pads every OCC bounding box by a **fixed, absolute** 1e-7 model units --
#: measured identical on rectangles of 1, 10, 100 and 1000 mm. So a purely
#: relative budget is wrong: 1e-9 * span alone passes only where a face spans
#: more than 100 mm, and a correctly ordered STEP of a part ten times smaller
#: than test_fin_2 is refused. The absolute floor is that pad with a decade of
#: margin; the relative term carries large parts where float noise scales.
_HULL_PAD_MM = 1e-6
_HULL_TOL_REL = 1e-9

#: Entity types whose extent is a radius rather than a CARTESIAN_POINT, with the
#: number of leading numeric arguments that are lengths. A face bounded by one
#: of these reaches up to that radius away from any point recorded for it, so
#: the point hull alone is NOT a superset -- see _collect_extent.
_RADIUS_ARGS = {
    "CIRCLE": 1,
    "ELLIPSE": 2,
    "CYLINDRICAL_SURFACE": 1,
    "CONICAL_SURFACE": 1,
    "SPHERICAL_SURFACE": 1,
    "TOROIDAL_SURFACE": 2,
    "DEGENERATE_TOROIDAL_SURFACE": 2,
}


def _parse_step_entities(text: str) -> dict[int, tuple[str, str]]:
    ents: dict[int, tuple[str, str]] = {}
    for m in re.finditer(r"#(\d+)=([A-Z0-9_]+)\((.*?)\);", text, re.S):
        ents[int(m.group(1))] = (m.group(2), m.group(3))
    return ents


def _collect_points(
    ents: dict[int, tuple[str, str]], eid: int, seen: set[int] | None = None
) -> list[list[float]]:
    if seen is None:
        seen = set()
    if eid in seen or eid not in ents:
        return []
    seen.add(eid)
    typ, arg = ents[eid]
    pts: list[list[float]] = []
    if typ == _CARTESIAN_POINT:
        m = re.search(r"\(([^#)]+)\)\s*$", arg.replace("\n", ""))
        if m:
            nums = [float(x) for x in m.group(1).split(",") if x.strip()]
            if len(nums) >= 3:
                pts.append(nums[:3])
    for ref in re.findall(r"#(\d+)", arg):
        pts.extend(_collect_points(ents, int(ref), seen))
    return pts


def _collect_radii(
    ents: dict[int, tuple[str, str]], eid: int, seen: set[int] | None = None
) -> float:
    """Largest radius-like length anywhere under ``eid``, in STEP units.

    ``_collect_points`` sees only ``CARTESIAN_POINT``s, and for a circle,
    ellipse or cylinder the extent lives in a **radius**, not a point. A face
    bounded by a full circle records its centre and nothing else, so its point
    hull is a degenerate segment -- tighter than the face, not a superset.
    Inflating every hull by the largest radius under it restores the superset
    property conservatively.
    """
    if seen is None:
        seen = set()
    if eid in seen or eid not in ents:
        return 0.0
    seen.add(eid)
    typ, arg = ents[eid]
    radius = 0.0
    count = _RADIUS_ARGS.get(typ)
    if count:
        numbers = re.findall(r"(?<![#\w.])(-?\d+\.\d*(?:[eE][-+]?\d+)?)", arg)
        for value in numbers[:count]:
            radius = max(radius, abs(float(value)))
    for ref in re.findall(r"#(\d+)", arg):
        radius = max(radius, _collect_radii(ents, int(ref), seen))
    return radius


def shell_faces(step_path: str | Path) -> dict[str, tuple[int, ...]]:
    """Named ``SHELL_BASED_SURFACE_MODEL`` -> its ``ADVANCED_FACE`` ids, in order.

    A named shell is a *collection* of faces, not one face: in
    ``test_fin_2.step`` the six names hold 7/5/1/3/4/1 faces for 21 in total,
    and OCC imports exactly those 21 surfaces in declaration order. That
    correspondence is what lets a name reach an element exactly, with no
    geometric tolerance -- see ``cases/step_fin/README.md``.
    """
    text = Path(step_path).read_text()
    ents = _parse_step_entities(text)
    out: dict[str, tuple[int, ...]] = {}
    for _eid, (typ, arg) in ents.items():
        if typ != _SHELL_MODEL:
            continue
        name_m = re.search(r"'([^']+)'", arg)
        if name_m is None:
            continue
        faces: list[int] = []
        for shell_ref in (int(r) for r in re.findall(r"#(\d+)", arg)):
            if shell_ref not in ents:
                raise GeometryError(
                    f"STEP shell {name_m.group(1)!r} references missing #{shell_ref}"
                )
            faces.extend(int(r) for r in re.findall(r"#(\d+)", ents[shell_ref][1]))
        if not faces:
            raise GeometryError(f"STEP shell {name_m.group(1)!r} holds no faces")
        out[name_m.group(1)] = tuple(faces)
    if not out:
        raise GeometryError(
            f"no {_SHELL_MODEL} names in {step_path}; export overlapping "
            "named surfaces/bodies from Onshape"
        )
    return out


def shell_hulls_mm(step_path: str | Path) -> dict[str, tuple[float, ...]]:
    """Named shell -> ``(xmin, ymin, zmin, xmax, ymax, zmax)`` of its points, mm.

    These are ``CARTESIAN_POINT``s, which for a B-spline face are **control
    points**. By the convex-hull property the trimmed face lies inside them, so
    this box is a superset of the real surface and never a description of it.

    That distinction is the whole reason this function is not a mask. Zone
    membership used to be "is the element centroid inside this shell's x-range",
    and on ``test_fin_2.step`` the hull overhangs the true face by up to 225 mm:
    ``HALF`` and ``QUARTER`` came out with identical ELSETs and ``TIP`` covered
    2.5x its true area. Use it only as a containment check.
    """
    text = Path(step_path).read_text()
    ents = _parse_step_entities(text)
    out: dict[str, tuple[float, ...]] = {}
    for eid, (typ, arg) in ents.items():
        if typ != _SHELL_MODEL:
            continue
        name_m = re.search(r"'([^']+)'", arg)
        if name_m is None:
            continue
        pts = _collect_points(ents, eid)
        if not pts:
            raise GeometryError(
                f"STEP shell {name_m.group(1)!r} has no CARTESIAN_POINT data"
            )
        cols = [[p[k] * _STEP_M_TO_MM for p in pts] for k in (0, 1, 2)]
        out[name_m.group(1)] = (
            min(cols[0]), min(cols[1]), min(cols[2]),
            max(cols[0]), max(cols[1]), max(cols[2]),
        )
    if not out:
        raise GeometryError(
            f"no {_SHELL_MODEL} names in {step_path}; export overlapping "
            "named surfaces/bodies from Onshape"
        )
    return out


def face_hulls_mm(step_path: str | Path) -> dict[int, tuple[float, ...]]:
    """``ADVANCED_FACE`` id -> its control-point bounding box in mm.

    Per **face**, not per shell, and that distinction carries the alignment
    check. Two shells can share a hull exactly -- ``HALF`` and ``QUARTER`` in
    ``test_fin_2.step`` do, to the last digit -- so a per-shell containment test
    has no power between them, which is the pair that caused the original zone
    bug.

    Their individual faces do not all differ either: three of ``QUARTER``'s
    hulls match three of ``HALF``'s exactly. What makes the per-face check
    sufficient is not that every face is distinguishable, but that **every
    swap it cannot see moves no element**: of the 210 pairwise face swaps on
    this STEP, 27 pass, and in all 27 the two faces fragment to identical tile
    sets, so no ELSET changes. Zero harmful misses, and 0/2000 random
    permutations missed. ``tests/test_step_mesh.py`` asserts that property
    rather than the distinguishability it does not have.

    The box is inflated by any radius under the face -- see ``_collect_radii``,
    without which this is not a superset and refuses valid geometry.
    """
    ents = _parse_step_entities(Path(step_path).read_text())
    wanted: set[int] = set()
    for _eid, (typ, arg) in ents.items():
        if typ != _SHELL_MODEL:
            continue
        for shell_ref in (int(r) for r in re.findall(r"#(\d+)", arg)):
            if shell_ref in ents:
                wanted.update(int(r) for r in re.findall(r"#(\d+)", ents[shell_ref][1]))
    out: dict[int, tuple[float, ...]] = {}
    for fid in wanted:
        pts = _collect_points(ents, fid)
        if not pts:
            raise GeometryError(f"STEP face #{fid} has no CARTESIAN_POINT data")
        cols = [[p[k] * _STEP_M_TO_MM for p in pts] for k in (0, 1, 2)]
        grow = _collect_radii(ents, fid) * _STEP_M_TO_MM
        out[fid] = (
            min(cols[0]) - grow, min(cols[1]) - grow, min(cols[2]) - grow,
            max(cols[0]) + grow, max(cols[1]) + grow, max(cols[2]) + grow,
        )
    return out


def shell_x_ranges_mm(step_path: str | Path) -> dict[str, tuple[float, float]]:
    """Named shell -> ``(xmin, xmax)`` of its control-point hull, mm.

    Kept because it is how the names are read and reported. It is a **hull**,
    not a mask: see ``shell_hulls_mm``. Nothing decides zone membership from it.
    """
    return {n: (h[0], h[3]) for n, h in shell_hulls_mm(step_path).items()}


def _tiles_by_shell(
    face_order: list[str],
    out_map: list[list[tuple[int, int]]],
) -> dict[str, set[int]]:
    """Shell name -> the fragment tiles its faces became.

    ``out_map`` is what ``occ.fragment`` returns alongside the entity list: it
    is parallel to ``objects + tools``, and entry *i* lists the tiles that input
    surface *i* was cut into. ``face_order`` is the shell name of each input
    surface, in the same order.
    """
    tiles: dict[str, set[int]] = {}
    for index, name in enumerate(face_order):
        for _dim, tag in out_map[index]:
            tiles.setdefault(name, set()).add(tag)
    return tiles


def _check_alignment(
    face_ids: list[int],
    face_hulls: dict[int, tuple[float, ...]],
    tile_bbox: dict[int, tuple[float, ...]],
    labels: dict[int, str],
    face_order: list[str],
    surf_tags: list[int],
    out_map: list[list[tuple[int, int]]],
    tiles_by_shell: dict[str, set[int]],
) -> None:
    """Prove the STEP-face -> OCC-surface correspondence, or refuse.

    Exactly one thing here is a guess. The name of a shell and the faces it owns
    are read straight out of ``SHELL_BASED_SURFACE_MODEL`` -> ``OPEN_SHELL``,
    with nothing to get wrong. The guess is that **OCC imports surfaces in STEP
    face order**, so that imported surface *i* is ``face_ids[i]``. That is
    deterministic, but assuming it silently is how this module got its zones
    wrong the first time.

    The check is **per face**, not per shell. A shell's tiles must lie inside
    that shell's own control-point box is far too weak: ``HALF`` and ``QUARTER``
    in ``test_fin_2.step`` have byte-identical hulls, so a per-shell test has no
    power between exactly the pair whose ELSETs the old x-range mask collapsed.
    Per face it separates them by 122 mm. Measured on this STEP, per face:

        correct order                 0.0000 mm overhang
        HALF/QUARTER faces swapped  122.4899
        rotated by one              674.9239
        reversed                    997.4138
        HEAL/TIP faces swapped     1002.2619
        200 random permutations     200/200 caught

    Containment is a superset property, but only once the hull is inflated by
    any radius under the face: ``_collect_points`` sees ``CARTESIAN_POINT``s
    only, and a circle or cylinder carries its extent in a radius, so the raw
    point hull can be *tighter* than the face. ``face_hulls_mm`` inflates for
    that. Note this STEP has **no B-splines at all** -- 17 planes and 4
    cylinders, bounded by 69 lines and 8 ellipses -- so the convex-hull property
    of a B-spline, which an earlier version of this docstring cited, is not what
    makes the check sound here.

    The tolerance is ``max(rel * span, _HULL_PAD_MM)`` because gmsh pads every
    bounding box by a fixed absolute 1e-7; a purely relative budget refuses a
    correctly ordered STEP of a part under 100 mm.

    Second, weaker check: OCC labels, where they exist. ``OCCImportLabels``
    carries only 2 of this file's 6 names and pins them to coincident faces from
    other shells too, so it is read one way only -- a labelled tag's tiles must
    sit inside the tiles of the shell it names.

    Deliberately absent: **area conservation across ``fragment``**. It reads
    like the natural invariant and it is a tautology -- ``fragment`` conserves
    area, so every block partition passes it, including a HEAL/TIP swap, at
    ``rel = 0.000e+00``. Do not add it back believing it tests something.
    """
    for index, fid in enumerate(face_ids):
        hull = face_hulls[fid]
        span = max(
            hull[3] - hull[0], hull[4] - hull[1], hull[5] - hull[2], 1.0
        )
        overhang = 0.0
        for _dim, tag in out_map[index]:
            box = tile_bbox[tag]
            for axis in (0, 1, 2):
                overhang = max(
                    overhang, hull[axis] - box[axis], box[axis + 3] - hull[axis + 3]
                )
        if overhang > max(_HULL_TOL_REL * span, _HULL_PAD_MM):
            raise GeometryError(
                f"imported surface {surf_tags[index]} was matched to STEP face "
                f"#{fid} (shell {face_order[index]!r}), but its tiles stick "
                f"{overhang:.4g} mm outside that face's own control-point hull. "
                "OCC did not import the faces in STEP order, so every zone in "
                "this mesh is suspect. Re-export the STEP, or map the names by "
                "hand."
            )

    for index, tag in enumerate(surf_tags):
        label = labels.get(tag, "")
        if not label or label not in tiles_by_shell:
            continue
        got = {t for _dim, t in out_map[index]}
        if not got <= tiles_by_shell[label]:
            raise GeometryError(
                f"OCC labelled surface {tag} as {label!r}, but STEP face order "
                f"assigned it to {face_order[index]!r}; its tiles "
                f"{sorted(got)} are not within {label!r}'s "
                f"{sorted(tiles_by_shell[label])}"
            )


#: A meshed zone under-measures a curved CAD face by chordal error -- on
#: test_fin_2.step the worst is 0.007%, all of it in one 4.8 mm curved band. 1%
#: leaves that alone while catching a tile the mesher failed to fill.
_ZONE_AREA_TOL_REL = 0.01


def _check_zone_areas(
    elsets: dict[str, tuple[int, ...]],
    cad_area: dict[str, float],
    element_area: dict[int, float],
) -> None:
    """Elements must cover the tiles their zone was built from.

    **This does not check the name-to-tile assignment**, and an earlier version
    of this docstring claimed it did. Both sides are keyed off the same tile
    set: ``cad_area[z]`` sums the tiles assigned to ``z`` and ``elsets[z]`` is
    the elements meshed on those same tiles, so permuting the assignment moves
    both together and the test passes. That is the same tautology as area
    conservation across ``fragment``, reproduced one level up. Alignment is
    checked per face in ``_check_alignment``; that is the only thing that
    checks it.

    What this does catch is a mesh that failed to fill its geometry -- a tile
    that produced no elements, or far too few -- which is a real failure mode
    here and one ccx would solve without complaint. Keep it for that, and do
    not read it as a statement about zones.
    """
    for name, area in sorted(cad_area.items()):
        meshed = sum(element_area[e] for e in elsets.get(name, ()))
        if area <= 0.0:
            raise GeometryError(f"zone {name!r} has non-positive CAD area {area:g}")
        rel = abs(meshed - area) / area
        if rel > _ZONE_AREA_TOL_REL:
            raise GeometryError(
                f"zone {name!r} meshes to {meshed:.4g} mm^2 but its CAD faces "
                f"cover {area:.4g} mm^2 ({rel:.2%} off). The elements assigned "
                "to this zone are not the ones the named shell covers"
            )


def zone_report(mesh: Mesh) -> list[dict[str, float | int | str]]:
    """Per-zone element count, area and extent -- look before you solve.

    The point is that a zone you drew in CAD and a zone the deck actually
    carries are different objects, and the only cheap way to know they agree is
    to read the areas off. ``blade`` is included: every other zone is a subset
    of it.
    """
    areas = mesh.element_areas()
    rows: list[dict[str, float | int | str]] = []
    for name in sorted(mesh.elsets, key=lambda n: (-len(mesh.elsets[n]), n)):
        eids = mesh.elsets[name]
        pts = [mesh.nodes[n] for e in eids for n in mesh.elements[e]]
        rows.append(
            {
                "zone": name,
                "elements": len(eids),
                "area_mm2": sum(areas[e] for e in eids),
                "x_min_mm": min(p[0] for p in pts),
                "x_max_mm": max(p[0] for p in pts),
                "y_min_mm": min(p[1] for p in pts),
                "y_max_mm": max(p[1] for p in pts),
            }
        )
    return rows


def _sanitize_elset(name: str) -> str:
    """CalculiX set names: letters, digits, underscore."""
    cleaned = re.sub(r"[^A-Za-z0-9_]", "_", name.strip())
    if not cleaned:
        raise GeometryError(f"empty coverage name from {name!r}")
    if cleaned[0].isdigit():
        cleaned = f"z_{cleaned}"
    return cleaned


def mesh_step(
    step_path: str | Path,
    *,
    size_mm: float = 40.0,
    coverage_tol_mm: float | None = None,
    clamp_coverage: str | None = "HEAL",
    long_axis: str = "x",
    quad_floor: float = 0.98,
    heading: str = "",
) -> Mesh:
    """Import a STEP of nested coverage shells -> a shell ``Mesh``.

    Mostly S8R, with any stray tri6 kept as S6 rather than dropped.
    ``quad_floor`` is how triangle-heavy the result may get before it is
    refused outright -- a few strays cost a little local bending accuracy, a
    mesh that is mostly triangles is a different model.

    ``fragment`` imprints overlapping shell boundaries onto a single set of
    tiles so ply drops get mesh edges. Coverage ELSETs are ACP-style masks:
    an element may belong to several (e.g. FULL and HALF and HEAL).

    Membership is **exact**: each named shell's faces are followed through the
    ``fragment`` map to the tiles they became, and an element belongs to a zone
    if its tile does. There is no geometric tolerance in that decision. It used
    to be a centroid-in-x-range test against the shell's control-point hull,
    which on ``test_fin_2.step`` gave ``HALF`` and ``QUARTER`` identical ELSETs
    and made ``TIP`` 2.5x too large -- see ``cases/step_fin/README.md``.
    ``coverage_tol_mm`` is therefore dead and warns if passed.

    Tip drive set (``far_face``) is the free tip edge: the mesh end along
    ``long_axis`` opposite the clamp. Clamp set (``fixed_end``) defaults to
    every node under the ``HEAL`` coverage mask; pass ``clamp_coverage=None``
    to fall back to the low end of ``long_axis`` only. ``long_axis`` must match
    the ply book / build flag (``x`` for ``test_fin_2``, ``y`` for ``FIN_TEST_3``).
    """
    step_path = Path(step_path)
    if not step_path.is_file():
        raise FileNotFoundError(step_path)
    if not (math.isfinite(size_mm) and size_mm > 0):
        raise GeometryError(f"size_mm must be positive, got {size_mm}")
    if long_axis not in ("x", "y"):
        raise GeometryError(f"long_axis must be 'x' or 'y', got {long_axis!r}")

    if coverage_tol_mm is not None:
        warnings.warn(
            "coverage_tol_mm no longer does anything: zone membership comes "
            "from the OCC fragment map, not from a centroid-in-range test",
            DeprecationWarning,
            stacklevel=2,
        )

    faces_by_shell = shell_faces(step_path)
    face_hulls = face_hulls_mm(step_path)
    coverages = {name: _sanitize_elset(name) for name in faces_by_shell}
    if len(set(coverages.values())) != len(coverages):
        raise GeometryError(
            f"shell names collide after sanitizing for ccx: {coverages}"
        )
    # Faces in ascending STEP entity id, which is the order OCC imports them in
    # -- the one assumption in this module, checked in _check_alignment.
    _by_id = sorted(
        (fid, name) for name, fids in faces_by_shell.items() for fid in fids
    )
    face_ids = [fid for fid, _name in _by_id]
    face_order = [name for _fid, name in _by_id]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.model.add("step_mesh")
        # Labels are partial -- OCC carries only some STEP names and pins them
        # to coincident faces from other shells too -- so they are a check in
        # _check_alignment, never the mapping.
        gmsh.option.setNumber("Geometry.OCCImportLabels", 1)
        gmsh.model.occ.importShapes(str(step_path.resolve()))
        gmsh.model.occ.synchronize()
        surfs = gmsh.model.getEntities(2)
        if len(surfs) < 1:
            raise GeometryError(f"no surfaces in {step_path}")
        if len(surfs) != len(face_order):
            raise GeometryError(
                f"{step_path.name} declares {len(face_order)} faces across "
                f"{len(faces_by_shell)} named shells but OCC imported "
                f"{len(surfs)} surfaces; the face-to-surface correspondence "
                "zone membership depends on cannot be established"
            )
        surf_tags = [tag for _dim, tag in surfs]
        labels = {
            tag: gmsh.model.getEntityName(2, tag).rsplit("/", 1)[-1]
            for tag in surf_tags
        }

        # Imprint overlapping shells: creates vertices along ply-drop curves.
        # out_map is parallel to objects+tools and carries the face -> tile
        # correspondence that every zone ELSET is built from.
        if len(surfs) > 1:
            _out, out_map = gmsh.model.occ.fragment([surfs[0]], list(surfs[1:]))
            gmsh.model.occ.synchronize()
        else:
            out_map = [[surfs[0]]]

        tiles = gmsh.model.getEntities(2)
        if not tiles:
            raise GeometryError("fragment left no surfaces to mesh")

        tiles_by_shell = _tiles_by_shell(face_order, out_map)
        _check_alignment(
            face_ids,
            face_hulls,
            {tag: gmsh.model.getBoundingBox(2, tag) for _dim, tag in tiles},
            labels,
            face_order,
            surf_tags,
            out_map,
            tiles_by_shell,
        )
        shells_by_tile: dict[int, set[str]] = {}
        for name, tile_tags in tiles_by_shell.items():
            for tag in tile_tags:
                shells_by_tile.setdefault(tag, set()).add(coverages[name])

        for _dim, _tag in tiles:
            gmsh.model.mesh.setRecombine(2, _tag)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
        gmsh.option.setNumber("Mesh.ElementOrder", 2)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", size_mm * 0.4)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", size_mm)
        gmsh.model.mesh.generate(2)

        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        raw_xyz = {
            int(t): _snap(coords[3 * i : 3 * i + 3]) for i, t in enumerate(node_tags)
        }

        # Elements per tile entity (for coverage), then global renumber.
        tile_elems: dict[int, list[tuple[int, ...]]] = {}
        for dim, tag in tiles:
            # Read quad8 AND tri6. Recombination does not always give all
            # quads on a fragmented tile, and dropping the strays is not a
            # small error: deleting one interior element of 128 moved the
            # reported force 2.5%, against 0.8% of area, because it severs
            # load path. ccx takes S6 in a *SHELL SECTION, COMPOSITE and in
            # the same ELSET as S8R. Anything that is neither is refused --
            # it would still vanish.
            types, type_tags, _ = gmsh.model.mesh.getElements(dim, tag)
            bad = [
                int(t)
                for t, tags in zip(types, type_tags, strict=True)
                if int(t) not in (_GMSH_QUAD8, _GMSH_TRI6) and len(tags)
            ]
            if bad:
                raise GeometryError(
                    f"surface {tag} meshed as unsupported types {bad}; only "
                    "quad8 (S8R) and tri6 (S6) are read back, so those "
                    "elements would vanish and leave a hole"
                )
            rows: list[tuple[int, ...]] = []
            for gmsh_type, n_nodes in ((_GMSH_QUAD8, 8), (_GMSH_TRI6, 6)):
                etags, flat = gmsh.model.mesh.getElementsByType(gmsh_type, tag)
                if len(etags) == 0:
                    continue
                conn = np.array(flat, dtype=np.int64).reshape(len(etags), n_nodes)
                rows.extend(tuple(int(n) for n in row) for row in conn)
            if rows:
                tile_elems[tag] = rows

        if not tile_elems:
            raise GeometryError(
                f"no shell elements in {step_path}; try a smaller size_mm"
            )

        order = sorted(raw_xyz, key=lambda t: (raw_xyz[t][0], raw_xyz[t][1], t))
        renumber = {old: new for new, old in enumerate(order, start=1)}
        nodes = {renumber[t]: raw_xyz[t] for t in order}

        elements: dict[int, tuple[int, ...]] = {}
        elem_coverage: dict[int, set[str]] = {}
        next_eid = 1
        for tag, rows in sorted(tile_elems.items()):
            for row in rows:
                mapped = tuple(renumber[n] for n in row)
                # Exact: the element inherits the zones of the tile it was
                # meshed on. No centroid, no tolerance, no bounding box.
                elements[next_eid] = mapped
                elem_coverage[next_eid] = set(shells_by_tile.get(tag, ()))
                next_eid += 1

        elsets: dict[str, tuple[int, ...]] = {
            "blade": tuple(sorted(elements)),
        }
        by_cov: dict[str, list[int]] = defaultdict(list)
        for eid, names in elem_coverage.items():
            for name in names:
                by_cov[name].append(eid)
        for name, members in sorted(by_cov.items()):
            elsets[name] = tuple(sorted(members))

        # Sanity: every named shell should catch at least one element.
        missing = [n for n in coverages.values() if n not in elsets]
        if missing:
            raise GeometryError(
                f"coverage ELSETs empty for {missing}; the named shells "
                "reached no mesh tile, so those zones would silently carry no "
                "plies"
            )
        _check_zone_areas(
            elsets,
            {
                coverages[name]: sum(
                    gmsh.model.occ.getMass(2, t) for t in tile_tags
                )
                for name, tile_tags in tiles_by_shell.items()
            },
            {
                eid: area
                for eid, area in Mesh(
                    nodes=nodes, elements=elements, nsets={}, elsets={}
                )
                .element_areas()
                .items()
            },
        )
        orphans = sorted(e for e, names in elem_coverage.items() if not names)
        if orphans:
            raise GeometryError(
                f"{len(orphans)} elements (e.g. {orphans[:5]}) lie on a tile "
                "that no named shell covers, so no ply reaches them. Every "
                "element must be inside at least the outermost shell"
            )

        # Clamp = HEAL mask (default) or low end of long_axis.
        # Tip = mesh end along long_axis opposite the clamp (not always xmax).
        axis = 0 if long_axis == "x" else 1
        span = [p[axis] for p in nodes.values()]
        smin, smax = min(span), max(span)
        tip_tol = max(1e-6, 1e-4 * (smax - smin))

        def on_station(nid: int, target: float) -> bool:
            return abs(nodes[nid][axis] - target) <= tip_tol

        clamp_name = None
        if clamp_coverage is not None:
            clamp_name = _sanitize_elset(clamp_coverage)
            if clamp_name not in elsets:
                raise GeometryError(
                    f"clamp_coverage {clamp_coverage!r} -> {clamp_name!r} "
                    f"not in mesh elsets {sorted(elsets)}"
                )
            clamp_nodes: set[int] = set()
            for eid in elsets[clamp_name]:
                clamp_nodes.update(elements[eid])
            root_nodes = tuple(sorted(clamp_nodes))
        else:
            root_nodes = tuple(
                sorted(nid for nid in nodes if on_station(nid, smin))
            )
        if not root_nodes:
            raise GeometryError("clamp NSET (fixed_end) is empty")

        clamp_span = [nodes[nid][axis] for nid in root_nodes]
        clamp_mid = 0.5 * (min(clamp_span) + max(clamp_span))
        # Tip station = global end farther from the clamp patch.
        if abs(smax - clamp_mid) >= abs(smin - clamp_mid):
            tip_station = smax
            tip_side = f"{long_axis}max"
        else:
            tip_station = smin
            tip_side = f"{long_axis}min"
        tip_nodes = tuple(
            sorted(nid for nid in nodes if on_station(nid, tip_station))
        )
        if not tip_nodes:
            raise GeometryError(
                f"no tip-edge nodes at {long_axis}={tip_station:g} "
                f"(long_axis={long_axis!r})"
            )
        # Tip must sit outside the clamp along the span -- origin anywhere is fine.
        tip_span = [nodes[nid][axis] for nid in tip_nodes]
        if max(tip_span) < min(clamp_span) - tip_tol:
            pass  # tip entirely below clamp
        elif min(tip_span) > max(clamp_span) + tip_tol:
            pass  # tip entirely above clamp
        else:
            raise GeometryError(
                f"tip edge at {long_axis}={tip_station:g} overlaps clamp "
                f"[{min(clamp_span):g}, {max(clamp_span):g}]; check long_axis "
                f"({long_axis!r}) and clamp_coverage"
            )

        default_heading = (
            f"STEP import {step_path.name}: {len(elements)} S8R, "
            f"coverages {sorted(coverages.values())}; "
            f"clamp={'mask '+clamp_name if clamp_name else long_axis + 'min edge'}, "
            f"tip={tip_side} edge; zones from the OCC fragment map (exact)"
        )
        mesh = Mesh(
            nodes=nodes,
            elements=elements,
            nsets={"fixed_end": root_nodes, "far_face": tip_nodes},
            elsets=elsets,
            heading=heading or default_heading,
        )
    finally:
        gmsh.finalize()

    return _orient_normals_outward(mesh, quad_floor)


def _flip_quad8(conn: tuple[int, ...]) -> tuple[int, ...]:
    """Reverse a serendipity quad so the normal flips."""
    a, b, c, d, e, f, g, h = conn
    return (a, d, c, b, h, g, f, e)


def _flip_tri6(conn: tuple[int, ...]) -> tuple[int, ...]:
    """Reverse a 6-node triangle. Midsides follow their edges: 1-2, 2-3, 3-1."""
    a, b, c, d, e, f = conn
    return (a, c, b, f, e, d)


def _flip(conn: tuple[int, ...]) -> tuple[int, ...]:
    return _flip_quad8(conn) if len(conn) == 8 else _flip_tri6(conn)


def _orient_normals_outward(mesh: Mesh, quad_floor: float = 0.98) -> Mesh:
    """Force shell normals toward +z on the flat skin (first ply = -z face).

    Onshape surfaces in this export face -z, so any element with a clear -z
    normal is flipped. The flip alone is a heuristic and is not the guarantee:
    the check afterwards is. Without it an element whose normal sits in the
    ``[-0.1, 0)`` band is neither flipped nor caught, and every unsymmetric
    stack on it is silently upside down -- the deck still solves and the
    reaction still looks reasonable. geometry.py raises on exactly this; so
    does this now.
    """
    flipped = dict(mesh.elements)
    n_flip = 0
    for eid, normal in mesh.element_normals().items():
        if normal[2] < -0.1:
            flipped[eid] = _flip(mesh.elements[eid])
            n_flip += 1
    out = (
        mesh
        if n_flip == 0
        else Mesh(
            nodes=mesh.nodes,
            elements=flipped,
            nsets=mesh.nsets,
            elsets=mesh.elsets,
            heading=mesh.heading,
        )
    )
    _check_normals_up(out)
    check_watertight(out, what="STEP mesh")
    check_quad_fraction(out, quad_floor, what="STEP mesh")
    return out


def _check_normals_up(mesh: Mesh) -> None:
    """Every element normal must have a positive +z component after flipping.

    Written as ``not (nz > 0)`` rather than ``nz <= 0`` because a degenerate
    element gives a nan normal and every comparison against nan is False; the
    other spelling waves those through. Same reasoning as
    ``geometry._check_normals``.
    """
    for eid, normal in mesh.element_normals().items():
        nz = float(normal[2])
        if not nz > 0.0:
            raise GeometryError(
                f"element {eid} has normal {normal} after reorientation, which "
                "does not face the laminate's +z side; the first ply line would "
                "land on the wrong face and every unsymmetric stack on that "
                "element would be inverted"
            )
