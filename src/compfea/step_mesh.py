"""STEP mid-surface import with nested ACP-style coverage masks.

Onshape/ANSYS-style overlapping named shells (FULL, HALF, TIP, ...) are
coverage masks, not exclusive tiles. This module:

1. Imports all shells via OpenCASCADE
2. ``fragment``s them so ply-drop boundaries become mesh edges (vertices)
3. Meshes the unique tiles as second-order incomplete quads (S8R)
4. Tags each element with every shell whose x-span contains its centroid
5. Builds root/tip NSETs from the min/max-x boundary nodes

Units: STEP SI metres are scaled to mm to match the rest of the repo.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from pathlib import Path

import gmsh
import numpy as np

from compfea.geometry import GeometryError, Mesh, _GMSH_QUAD8, _snap

# Onshape STEP AP242 length unit is metre; CalculiX decks here are mm.
_STEP_M_TO_MM = 1000.0
_SHELL_MODEL = "SHELL_BASED_SURFACE_MODEL"
_CARTESIAN_POINT = "CARTESIAN_POINT"


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


def shell_x_ranges_mm(step_path: str | Path) -> dict[str, tuple[float, float]]:
    """Named SHELL_BASED_SURFACE_MODEL -> (xmin, xmax) in mm."""
    text = Path(step_path).read_text()
    ents = _parse_step_entities(text)
    out: dict[str, tuple[float, float]] = {}
    for eid, (typ, arg) in ents.items():
        if typ != _SHELL_MODEL:
            continue
        name_m = re.search(r"'([^']+)'", arg)
        if name_m is None:
            continue
        name = name_m.group(1)
        pts = _collect_points(ents, eid)
        if not pts:
            raise GeometryError(f"STEP shell {name!r} has no CARTESIAN_POINT data")
        xs = [p[0] * _STEP_M_TO_MM for p in pts]
        out[name] = (min(xs), max(xs))
    if not out:
        raise GeometryError(
            f"no {_SHELL_MODEL} names in {step_path}; export overlapping "
            "named surfaces/bodies from Onshape"
        )
    return out


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
    coverage_tol_mm: float = 2.0,
    clamp_coverage: str | None = "HEAL",
    heading: str = "",
) -> Mesh:
    """Import a STEP of nested coverage shells -> S8R ``Mesh``.

    ``fragment`` imprints overlapping shell boundaries onto a single set of
    tiles so ply drops get mesh edges. Coverage ELSETs are ACP-style masks:
    an element may belong to several (e.g. FULL and HALF and HEAL).

    Tip drive set (``far_face``) is the maximum-x boundary edge (fin span
    along +x). Clamp set (``fixed_end``) defaults to every node under the
    ``HEAL`` coverage mask; pass ``clamp_coverage=None`` to fall back to the
    minimum-x edge only.
    """
    step_path = Path(step_path)
    if not step_path.is_file():
        raise FileNotFoundError(step_path)
    if not (math.isfinite(size_mm) and size_mm > 0):
        raise GeometryError(f"size_mm must be positive, got {size_mm}")

    ranges = shell_x_ranges_mm(step_path)
    coverages = {name: _sanitize_elset(name) for name in ranges}

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.model.add("step_mesh")
        gmsh.model.occ.importShapes(str(step_path.resolve()))
        gmsh.model.occ.synchronize()
        surfs = gmsh.model.getEntities(2)
        if len(surfs) < 1:
            raise GeometryError(f"no surfaces in {step_path}")

        # Imprint overlapping shells: creates vertices along ply-drop curves.
        if len(surfs) > 1:
            gmsh.model.occ.fragment([surfs[0]], list(surfs[1:]))
            gmsh.model.occ.synchronize()

        tiles = gmsh.model.getEntities(2)
        if not tiles:
            raise GeometryError("fragment left no surfaces to mesh")

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
            etags, flat = gmsh.model.mesh.getElementsByType(_GMSH_QUAD8, tag)
            if len(etags) == 0:
                # Allow empty tile only if some other type snuck in — refuse tris.
                types, type_tags, _ = gmsh.model.mesh.getElements(dim, tag)
                bad = [
                    int(t)
                    for t, tags in zip(types, type_tags, strict=True)
                    if int(t) != _GMSH_QUAD8 and len(tags)
                ]
                if bad:
                    raise GeometryError(
                        f"surface {tag} meshed as non-quad8 types {bad}; "
                        "refuse triangles for S8R"
                    )
                continue
            conn = np.array(flat, dtype=np.int64).reshape(len(etags), 8)
            tile_elems[tag] = [tuple(int(n) for n in row) for row in conn]

        if not tile_elems:
            raise GeometryError(
                f"no quad8 elements in {step_path}; try a smaller size_mm"
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
                # Element centroid (corners only) vs each shell x-span.
                cx = sum(nodes[n][0] for n in mapped[:4]) / 4.0
                names = {
                    coverages[name]
                    for name, (x0, x1) in ranges.items()
                    if x0 - coverage_tol_mm <= cx <= x1 + coverage_tol_mm
                }
                elements[next_eid] = mapped
                elem_coverage[next_eid] = set(names)
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
                f"coverage ELSETs empty for {missing}; shell x-ranges were "
                f"{ {k: (round(v[0],3), round(v[1],3)) for k,v in ranges.items()} }"
            )

        # Tip = free tip edge at xmax. Clamp = all nodes under clamp_coverage
        # (default HEAL), else the xmin edge.
        xs = [p[0] for p in nodes.values()]
        xmin, xmax = min(xs), max(xs)
        tip_tol = max(1e-6, 1e-4 * (xmax - xmin))

        def on_plane(x: float, target: float) -> bool:
            return abs(x - target) <= tip_tol

        tip_nodes = tuple(
            sorted(nid for nid, (x, _, _) in nodes.items() if on_plane(x, xmax))
        )
        if not tip_nodes:
            raise GeometryError(
                f"no tip-edge nodes at xmax={xmax:g}"
            )

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
                sorted(nid for nid, (x, _, _) in nodes.items() if on_plane(x, xmin))
            )
        if not root_nodes:
            raise GeometryError("clamp NSET (fixed_end) is empty")

        default_heading = (
            f"STEP import {step_path.name}: {len(elements)} S8R, "
            f"coverages {sorted(coverages.values())}; "
            f"clamp={'mask '+clamp_name if clamp_name else 'xmin edge'}, "
            "tip=xmax edge; nested ACP masks via OCC fragment imprint"
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

    return _orient_normals_outward(mesh)


def _flip_quad8(conn: tuple[int, ...]) -> tuple[int, ...]:
    """Reverse a serendipity quad so the normal flips."""
    a, b, c, d, e, f, g, h = conn
    return (a, d, c, b, h, g, f, e)


def _orient_normals_outward(mesh: Mesh) -> Mesh:
    """Force shell normals toward +z on the flat skin (first ply = -z face).

    Onshape surfaces in this export face -z; HEAL is bent so its normal is
    mostly in-plane. Flip any element with a clear -z normal component.
    """
    flipped = dict(mesh.elements)
    n_flip = 0
    for eid, normal in mesh.element_normals().items():
        if normal[2] < -0.1:
            flipped[eid] = _flip_quad8(mesh.elements[eid])
            n_flip += 1
    if n_flip == 0:
        return mesh
    return Mesh(
        nodes=mesh.nodes,
        elements=flipped,
        nsets=mesh.nsets,
        elsets=mesh.elsets,
        heading=mesh.heading,
    )
