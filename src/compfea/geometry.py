"""Planform outline -> quad8 shell mesh -> CalculiX .inp fragment.

Emits *NODE / *ELEMENT / *NSET / *ELSET only. Materials, sections and steps come
from `layup.py` and `deck.py`; `deck.assemble` takes what `Mesh.to_inp` returns.

Units: mm, N, tonne, MPa, s. Coordinates are mm.

## Axes

x is the chord direction, y is the span (the long axis), z is the shell normal.
That matches every case deck in the repo. It transposes CLAUDE.md's "+x along
the long axis", which is why `layup.py` takes an explicit `long_axis` -- a mesh
from here wants `long_axis="y"`.

## Three things here are load-bearing

- **Element type 16, not 10.** Type 16 is the 8-node serendipity quad that maps
  to S8R; type 10 is the 9-node quad, and it is what
  `gmsh.model.mesh.getElementType("Quadrangle", 2)` hands you. The 8-node form
  comes from `Mesh.SecondOrderIncomplete = 1`. Anything that is not type 16 is
  refused: only quad8 elements are read back, so a triangle would be silently
  *dropped*, leaving a hole in the part and nodes attached to nothing -- and a
  deck with a hole in it is still a deck ccx will solve.

- **The curve-loop direction sets the element normal, and the normal decides
  which ply is at -z.** A clockwise loop gives every element a -z normal, which
  turns the section card's first ply line -- documented as the -z ply -- into the
  +z ply, inverting every unsymmetric stack in the model. Nothing downstream
  would complain. The outline is therefore forced counter-clockwise by signed
  area before meshing, and every element's normal is checked afterwards.

- **gmsh's quad8 node order already is CalculiX's**: four corners, then mid-side
  node k between corners k and k+1. Verified against coordinates rather than
  assumed, in `tests/test_geometry.py`. Connectivity passes straight through.

gmsh is a process-global singleton, so this module initialises and finalises
around each call and never leaves state behind. `sweep.py` runs jobs in separate
processes, which is what makes that safe.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import gmsh
import numpy as np

Point = tuple[float, float]

EDGE_ROLES = ("root", "leading", "tip", "trailing")

# Coordinates are snapped to this many decimals. Transfinite meshing returns
# 3.1249999999877 for 3.125, which is physically irrelevant at 1e-9 mm but makes
# output non-reproducible and node sorting unstable.
_SNAP_DECIMALS = 9

# CalculiX quad8 shell.
_ELEMENT_TYPE = "S8R"
# gmsh's tag for the same element.
_GMSH_QUAD8 = 16


class GeometryError(ValueError):
    """The outline, the mesh request, or what gmsh returned is unusable."""


# --------------------------------------------------------------------------
# Curvature
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Camber:
    """Single-curvature lift of a flat planform, preserving arc length.

    The blade is meshed flat and then bent, and the bend is an **isometry**: a
    point at developed spanwise distance ``s`` moves to ``(x, Y(s), Z(s))`` where
    ``(Y, Z)`` traces a plane curve parameterised by its own arc length,

        theta(s) = integral_0^s kappa,   Y = integral cos theta,  Z = integral sin theta

    so distances along the span are unchanged and ``x`` is untouched. That is the
    right model for a laminate laid up flat and put in a mould: the fibres do not
    stretch, and the developed length is what the layup sees.

    Setting ``z = f(y)`` instead would stretch the span by the ratio of arc length
    to chord -- 4.6% at a 60 degree bend -- and quietly change the stiffness as
    the cube of the length. Do not "simplify" it back to that.

    Only single curvature is supported. A doubly curved surface is not
    developable, so a flat laminate cannot take that shape without in-plane
    strain, and the residual stress that implies is modelled nowhere in this repo.
    """

    def curvature(self, s: np.ndarray, span: float) -> np.ndarray:
        """Curvature in 1/mm at developed distance ``s``. Subclasses define this."""
        raise NotImplementedError

    def _grid(self, span: float) -> np.ndarray:
        return np.linspace(0.0, span, 4001)

    def tangent_angle(self, s: np.ndarray, span: float) -> np.ndarray:
        """theta(s), the integral of curvature. Always integrated on the same
        dense grid and then sampled: integrating over the node positions
        directly would give two nodes at equal s different answers, because they
        would not share a quadrature."""
        grid = self._grid(span)
        theta = _cumulative_trapezoid(self.curvature(grid, span), grid)
        return np.interp(s, grid, theta)

    def map_span(self, s: np.ndarray, span: float) -> tuple[np.ndarray, np.ndarray]:
        """Developed distance ``s`` -> ``(y, z)`` on the curved mid-surface."""
        grid = self._grid(span)
        theta = self.tangent_angle(grid, span)
        y = _cumulative_trapezoid(np.cos(theta), grid)
        z = _cumulative_trapezoid(np.sin(theta), grid)
        return np.interp(s, grid, y), np.interp(s, grid, z)


def _cumulative_trapezoid(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    steps = np.diff(x) * (values[:-1] + values[1:]) / 2.0
    return np.concatenate([[0.0], np.cumsum(steps)])


@dataclass(frozen=True)
class CircularCamber(Camber):
    """Constant curvature: the span becomes a circular arc of radius ``radius``.

    Positive radius bends toward +z. The arc subtends ``span / radius`` radians,
    so the tip lands at ``(R sin(L/R), R (1 - cos(L/R)))``.
    """

    radius: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.radius) or self.radius == 0.0:
            raise GeometryError(f"radius must be finite and nonzero, got {self.radius}")

    def curvature(self, s: np.ndarray, span: float) -> np.ndarray:
        return np.full_like(np.asarray(s, dtype=float), 1.0 / self.radius)


@dataclass(frozen=True)
class CurvatureProfile(Camber):
    """Piecewise-linear curvature against fraction of developed span.

    ``points`` is ``((s_fraction, kappa), ...)`` with fractions in [0, 1] and
    kappa in 1/mm. Curvature rather than a z profile on purpose: any kappa(s)
    integrates to a valid arc-length-preserving shape, while an arbitrary z(s)
    generally does not correspond to one.

    Curvature is linear between the given fractions and **held constant outside
    them**, not tapered to zero: a profile starting at 0.3 curves the inboard 30%
    of the blade by its first kappa. Give an explicit point at 0.0 and 1.0 if you
    want the ends flat.
    """

    points: tuple[tuple[float, float], ...]

    def __post_init__(self) -> None:
        if len(self.points) < 2:
            raise GeometryError("a curvature profile needs at least two points")
        fractions = [f for f, _ in self.points]
        if list(fractions) != sorted(fractions):
            raise GeometryError("curvature profile fractions must be increasing")
        if not all(math.isfinite(v) for pair in self.points for v in pair):
            raise GeometryError(f"curvature profile must be finite, got {self.points}")
        if fractions[0] < 0.0 or fractions[-1] > 1.0:
            raise GeometryError("curvature profile fractions must lie in [0, 1]")

    def curvature(self, s: np.ndarray, span: float) -> np.ndarray:
        fractions = np.array([f for f, _ in self.points])
        kappas = np.array([k for _, k in self.points])
        return np.interp(np.asarray(s, dtype=float) / span, fractions, kappas)


# --------------------------------------------------------------------------
# Outline
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Outline:
    """Planform boundary, as four edges in loop order.

    Each edge is a sequence of ``(x, y)`` points and may be any shape; the four
    *roles* are fixed because they are what the node sets are built from. The
    edges join end to end: ``root`` ends where ``leading`` starts, and
    ``trailing`` ends where ``root`` starts.

    Naming an edge in ``splines`` fits a Catmull-Rom spline through its points
    instead of joining them with straight segments.

    Orientation is normalised to counter-clockwise seen from +z on construction,
    so the mesh normal is +z whichever way the caller listed the points. That is
    not cosmetic -- see the module docstring.
    """

    root: tuple[Point, ...]
    leading: tuple[Point, ...]
    tip: tuple[Point, ...]
    trailing: tuple[Point, ...]
    splines: frozenset[str] = field(default_factory=frozenset)

    def __post_init__(self) -> None:
        for role in EDGE_ROLES:
            edge = getattr(self, role)
            if len(edge) < 2:
                raise GeometryError(f"edge {role!r} needs at least two points")
            for x, y in edge:
                if not (math.isfinite(x) and math.isfinite(y)):
                    raise GeometryError(f"edge {role!r} has a non-finite point")
        unknown = set(self.splines) - set(EDGE_ROLES)
        if unknown:
            raise GeometryError(f"unknown spline edges: {sorted(unknown)}")
        for role, nxt in zip(EDGE_ROLES, EDGE_ROLES[1:] + EDGE_ROLES[:1], strict=True):
            end = getattr(self, role)[-1]
            start = getattr(self, nxt)[0]
            if not _close(end, start):
                raise GeometryError(
                    f"edge {role!r} ends at {end} but {nxt!r} starts at {start}; "
                    "the four edges must join end to end in loop order"
                )
        if self.signed_area() == 0.0:
            raise GeometryError("outline encloses no area")
        self._check_spanwise_convention()
        if self.signed_area() < 0.0:
            # Walk the same loop the other way, so the normal comes out +z. Each
            # edge reverses, and the two side edges swap loop position: the new
            # order is root, trailing, tip, leading. root and tip keep their
            # identity -- they are what fixed_end and far_face are built from --
            # while the side labels are interchangeable and may end up swapped.
            reversed_edges = {r: tuple(reversed(getattr(self, r))) for r in EDGE_ROLES}
            moved = {"root": "root", "leading": "trailing", "tip": "tip",
                     "trailing": "leading"}
            for role, source in moved.items():
                object.__setattr__(self, role, reversed_edges[source])
            # `splines` names roles, so it has to move with them. Leaving it
            # alone splines whichever edge lands in the named slot, which on a
            # curved-leading-edge blade listed clockwise means the curve is
            # chorded and the straight edge is bulged -- a different part, with
            # no error.
            object.__setattr__(
                self,
                "splines",
                frozenset(
                    role for role, source in moved.items() if source in self.splines
                ),
            )

    def _check_spanwise_convention(self) -> None:
        """The root must be the inboard end of the blade and the tip the outboard.

        Everything spanwise -- zone fractions, the camber map -- measures from the
        smallest y in the outline. If the root edge is not the edge that sits
        there, a zone boundary asked for at 40% of span lands somewhere else and
        the camber anchors its zero curvature at the wrong end, both silently.
        The edges themselves may be any shape; only which end they occupy is
        fixed.
        """
        ys = [y for _, y in self.polygon()]
        low, high = min(ys), max(ys)
        root_y = [y for _, y in self.root]
        tip_y = [y for _, y in self.tip]
        if abs(min(root_y) - low) > 1e-9:
            raise GeometryError(
                f"the root edge must reach the inboard end of the blade "
                f"(y = {low:g}); its lowest point is y = {min(root_y):g}. "
                "Spanwise position is measured from there."
            )
        if abs(max(tip_y) - high) > 1e-9:
            raise GeometryError(
                f"the tip edge must reach the outboard end of the blade "
                f"(y = {high:g}); its highest point is y = {max(tip_y):g}."
            )
        if max(root_y) >= min(tip_y):
            raise GeometryError(
                "the root and tip edges overlap in y; the blade has no "
                "unambiguous spanwise direction"
            )

    def root_station(self) -> float:
        """The y the span is measured from."""
        return min(y for _, y in self.polygon())

    def polygon(self) -> tuple[Point, ...]:
        """The closed boundary as one point sequence, corners not repeated."""
        points: list[Point] = []
        for role in EDGE_ROLES:
            points.extend(getattr(self, role)[:-1])
        return tuple(points)

    def signed_area(self) -> float:
        """Shoelace area. Positive means counter-clockwise seen from +z."""
        pts = self.polygon()
        return 0.5 * sum(
            pts[i][0] * pts[(i + 1) % len(pts)][1]
            - pts[(i + 1) % len(pts)][0] * pts[i][1]
            for i in range(len(pts))
        )

    def span(self) -> float:
        ys = [y for _, y in self.polygon()]
        return max(ys) - min(ys)

    @classmethod
    def rectangle(cls, chord: float, span: float) -> Outline:
        """Constant-chord strip: x from 0 to ``chord``, y from 0 to ``span``."""
        return cls(
            root=((0.0, 0.0), (chord, 0.0)),
            leading=((chord, 0.0), (chord, span)),
            tip=((chord, span), (0.0, span)),
            trailing=((0.0, span), (0.0, 0.0)),
        )

    @classmethod
    def tapered(
        cls,
        root_chord: float,
        tip_chord: float,
        span: float,
        *,
        sweep: float = 0.0,
    ) -> Outline:
        """Straight-tapered blade. ``sweep`` offsets the tip in +x."""
        return cls(
            root=((0.0, 0.0), (root_chord, 0.0)),
            leading=((root_chord, 0.0), (sweep + tip_chord, span)),
            tip=((sweep + tip_chord, span), (sweep, span)),
            trailing=((sweep, span), (0.0, 0.0)),
        )


def _close(a: Point, b: Point, tol: float = 1e-9) -> bool:
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol


# --------------------------------------------------------------------------
# Mesh
# --------------------------------------------------------------------------


# Shell element by connectivity length. Serendipity quad8 -> S8R, and the
# 6-node triangle -> S6, which ccx accepts in a *SHELL SECTION, COMPOSITE and
# in the same ELSET as S8R (verified against ccx, not assumed). Triangles are
# read rather than dropped: deleting one interior element of 128 moved the
# reported force 2.5%, far more than its 0.8% of area, because it severs load
# path. A slightly stiff S6 is a much smaller error than a hole.
_SHELL_TYPE_BY_NODES = {8: "S8R", 6: "S6"}


def _corners(conn: tuple[int, ...]) -> tuple[int, ...]:
    """Corner nodes of a shell element: 4 for a quad8, 3 for a tri6."""
    n = 4 if len(conn) == 8 else 3
    return conn[:n]


def shell_type(conn: tuple[int, ...]) -> str:
    try:
        return _SHELL_TYPE_BY_NODES[len(conn)]
    except KeyError:
        raise GeometryError(
            f"element has {len(conn)} nodes; only quad8 (S8R) and tri6 (S6) "
            "shells are supported"
        ) from None


@dataclass(frozen=True)
class Mesh:
    """A quad8 shell mesh, ready to hand to `deck.assemble` as ``mesh_inp``."""

    nodes: dict[int, tuple[float, float, float]]
    elements: dict[int, tuple[int, ...]]
    nsets: dict[str, tuple[int, ...]]
    elsets: dict[str, tuple[int, ...]]
    heading: str = ""

    def element_normals(self) -> dict[int, tuple[float, float, float]]:
        """Unit normal per element, from the corner ordering. Should be ~+z."""
        out = {}
        for eid, conn in self.elements.items():
            c = [np.array(self.nodes[n]) for n in _corners(conn)]
            # Last corner either way: c[3] on a quad, c[2] on a triangle. Both
            # give the outward normal for counter-clockwise ordering.
            normal = np.cross(c[1] - c[0], c[-1] - c[0])
            length = float(np.linalg.norm(normal))
            if length == 0.0:
                raise GeometryError(
                    f"element {eid} is degenerate: its corners are collinear or "
                    "coincident, so it has no normal and no area"
                )
            out[eid] = tuple(normal / length)
        return out

    def element_areas(self) -> dict[int, float]:
        """Planform area per element, from its corners."""
        out = {}
        for eid, conn in self.elements.items():
            c = [np.array(self.nodes[n]) for n in _corners(conn)]
            area = 0.5 * float(np.linalg.norm(np.cross(c[1] - c[0], c[2] - c[0])))
            if len(c) == 4:
                area += 0.5 * float(
                    np.linalg.norm(np.cross(c[2] - c[0], c[3] - c[0]))
                )
            out[eid] = area
        return out

    def to_inp(self) -> str:
        lines: list[str] = []
        for line in self.heading.strip().splitlines():
            lines.append(line if line.startswith("**") else f"** {line}")
        lines.append("*NODE, NSET=all_nodes")
        for nid in sorted(self.nodes):
            x, y, z = self.nodes[nid]
            lines.append(f"{nid}, {x:.8f}, {y:.8f}, {z:.8f}")
        by_type: dict[str, list[int]] = {}
        for eid, conn in self.elements.items():
            by_type.setdefault(shell_type(conn), []).append(eid)
        # S8R first so an all-quad mesh emits byte-identical text to before.
        for etype in sorted(by_type, key=lambda t: t != "S8R"):
            lines.append(f"*ELEMENT, TYPE={etype}, ELSET=blade")
            for eid in sorted(by_type[etype]):
                lines.append(
                    f"{eid}, " + ", ".join(str(n) for n in self.elements[eid])
                )
        for name, members in self.nsets.items():
            lines.append(f"*NSET, NSET={name}")
            lines.extend(_wrap(members))
        for name, members in self.elsets.items():
            if name == "blade":
                continue  # already declared by the *ELEMENT card
            lines.append(f"*ELSET, ELSET={name}")
            lines.extend(_wrap(members))
        return "\n".join(lines) + "\n"


def _wrap(values: tuple[int, ...], per_line: int = 8) -> list[str]:
    return [
        ", ".join(str(v) for v in values[i : i + per_line])
        for i in range(0, len(values), per_line)
    ]


def mesh_outline(
    outline: Outline,
    *,
    n_chord: int | None = None,
    n_span: int | None = None,
    size: float | None = None,
    camber: Camber | None = None,
    zones: tuple[float, ...] = (),
    heading: str = "",
) -> Mesh:
    """Mesh ``outline`` with quad8 shells.

    Give ``n_chord`` and ``n_span`` for a structured transfinite grid -- the
    element counts across and along the blade -- or ``size`` for an unstructured
    mesh with that target edge length. Structured is preferred wherever the shape
    allows it: the grid is predictable and reproducible, which matters for node
    sets and for a sweep that caches on the design vector.

    ``zones`` are spanwise boundaries as fractions of the span measured from the
    root station, so ``(0.4, 0.7)`` yields ``zone_1``, ``zone_2``, ``zone_3`` from
    root to tip. Each gets its own ELSET, which is what `layup.py` hangs a
    ``*SHELL SECTION, COMPOSITE`` on.

    A boundary lands on whichever element row it falls in, so the achieved ply
    drop sits at a multiple of the element length, not exactly at the fraction
    asked for: 0.47 on a 10-row mesh puts it at 0.5. Refine the span if a drop
    has to land somewhere precise. A zone thinner than one row is an error, not
    an empty ELSET.
    """
    if (n_chord is None) != (n_span is None):
        raise GeometryError("give both n_chord and n_span, or neither")
    if n_chord is None and size is None:
        raise GeometryError("give either n_chord/n_span or size")
    if n_chord is not None and (n_chord < 1 or n_span < 1):
        raise GeometryError("n_chord and n_span must be >= 1")
    if size is not None and not (math.isfinite(size) and size > 0.0):
        raise GeometryError(f"size must be finite and positive, got {size}")

    span = outline.span()
    root_y = outline.root_station()
    # Everything that depends on position along the blade is decided on the flat
    # mesh, where the spanwise coordinate *is* the developed distance. Bending is
    # the last step and moves nodes only.
    flat_nodes, elements, edge_nodes = _generate(outline, n_chord, n_span, size)
    zone_sets = _zone_elsets(elements, flat_nodes, zones, span, root_y)

    nodes = (
        flat_nodes
        if camber is None
        else _apply_camber(flat_nodes, camber, span, root_y)
    )
    mesh = Mesh(
        nodes=nodes,
        elements=elements,
        nsets={"fixed_end": edge_nodes["root"], "far_face": edge_nodes["tip"]},
        elsets={"blade": tuple(sorted(elements)), **zone_sets},
        heading=heading,
    )
    _check_normals(mesh, flat_nodes, camber, span, root_y)
    check_watertight(mesh, what="meshed outline")
    return mesh


def _generate(
    outline: Outline,
    n_chord: int | None,
    n_span: int | None,
    size: float | None,
) -> tuple[
    dict[int, tuple[float, float, float]],
    dict[int, tuple[int, ...]],
    dict[str, tuple[int, ...]],
]:
    """Run gmsh; return snapped, renumbered nodes, elements and edge node sets."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("compfea_planform")
        geo = gmsh.model.geo

        curves = []
        # One gmsh point per distinct coordinate: edges share their corners, and
        # adding a second point at the same place leaves the loop open.
        registry: dict[tuple[float, float], int] = {}

        def point(x: float, y: float) -> int:
            key = (round(x, _SNAP_DECIMALS), round(y, _SNAP_DECIMALS))
            if key not in registry:
                registry[key] = geo.addPoint(x, y, 0.0, size or 1.0)
            return registry[key]

        for role in EDGE_ROLES:
            pts = getattr(outline, role)
            tags = [point(x, y) for x, y in pts]
            if role in outline.splines and len(tags) > 2:
                curves.append(geo.addSpline(tags))
            elif len(tags) == 2:
                curves.append(geo.addLine(tags[0], tags[1]))
            else:
                curves.append(geo.addPolyline(tags))
        surface = geo.addPlaneSurface([geo.addCurveLoop(curves)])

        if n_chord is not None:
            # Opposite sides must carry the same node count: root/tip are the
            # chordwise pair, leading/trailing the spanwise pair.
            counts = {
                "root": n_chord + 1,
                "tip": n_chord + 1,
                "leading": n_span + 1,
                "trailing": n_span + 1,
            }
            for role, curve in zip(EDGE_ROLES, curves, strict=True):
                geo.mesh.setTransfiniteCurve(curve, counts[role])
            geo.mesh.setTransfiniteSurface(surface)
        geo.mesh.setRecombine(2, surface)
        geo.synchronize()

        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.Algorithm", 8)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)
        gmsh.option.setNumber("Mesh.ElementOrder", 2)
        # The whole point: 8-node serendipity quads, not 9-node.
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
        try:
            gmsh.model.mesh.generate(2)
        except Exception as exc:  # gmsh raises bare Exception
            raise GeometryError(
                f"gmsh could not mesh this outline: {exc}. A `size` larger than "
                "the shortest boundary edge is the usual cause -- a narrow tip "
                "chord cannot be divided by it. Use a smaller size, or give "
                "n_chord/n_span for the structured path, which has no such limit."
            ) from exc

        present = [int(t) for t in gmsh.model.mesh.getElements(2)[0]]
        if present != [_GMSH_QUAD8]:
            raise GeometryError(
                f"gmsh returned 2D element types {present}, wanted only "
                f"{_GMSH_QUAD8} (8-node quad). Only quad8 elements are read back, "
                f"so anything else would be dropped, leaving holes in the part "
                f"and nodes attached to nothing. Refine the outline or set an "
                f"explicit size."
            )

        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        raw_xyz = {
            int(t): _snap(coords[3 * i : 3 * i + 3]) for i, t in enumerate(node_tags)
        }
        element_tags, connectivity = gmsh.model.mesh.getElementsByType(_GMSH_QUAD8)
        conn = np.array(connectivity, dtype=np.int64).reshape(len(element_tags), 8)
        # Edge membership comes from the curve entities, not from a distance test
        # against the control points: a splined edge leaves its control polyline,
        # and the nodes that bulge away from it are exactly the ones a distance
        # test drops -- silently, since the control points themselves always
        # match. far_face is the driven set in every case deck, so losing part of
        # it turns an edge drive into a few-point drive.
        raw_edges = {
            role: set(int(t) for t in gmsh.model.mesh.getNodes(1, curve, True)[0])
            for role, curve in zip(EDGE_ROLES, curves, strict=True)
        }
    finally:
        gmsh.finalize()

    # Renumber: nodes up the span then across the chord, elements the same way,
    # so the file is identical between runs and readable next to the geometry.
    order = sorted(raw_xyz, key=lambda t: (raw_xyz[t][1], raw_xyz[t][0], t))
    renumber = {old: new for new, old in enumerate(order, start=1)}
    nodes = {renumber[t]: raw_xyz[t] for t in order}

    mapped = [tuple(renumber[int(n)] for n in row) for row in conn]
    centroid = [
        (
            sum(nodes[n][1] for n in row[:4]) / 4.0,
            sum(nodes[n][0] for n in row[:4]) / 4.0,
        )
        for row in mapped
    ]
    element_order = sorted(range(len(mapped)), key=lambda i: (centroid[i], mapped[i]))
    elements = {new: mapped[i] for new, i in enumerate(element_order, start=1)}

    edge_nodes = {}
    for role in ("root", "tip"):
        members = tuple(sorted(renumber[t] for t in raw_edges[role] if t in renumber))
        if not members:
            raise GeometryError(f"no nodes landed on the {role!r} edge")
        edge_nodes[role] = members
    return nodes, elements, edge_nodes


def _snap(values) -> tuple[float, float, float]:
    return tuple(round(float(v), _SNAP_DECIMALS) + 0.0 for v in values)


def _apply_camber(
    nodes: dict[int, tuple[float, float, float]],
    camber: Camber,
    span: float,
    root_y: float,
) -> dict[int, tuple[float, float, float]]:
    """Bend the flat mesh, preserving spanwise arc length.

    ``s`` is measured from the root station, not from y = 0: an outline that does
    not start at the origin would otherwise run off the end of the integration
    grid, where np.interp clamps and collapses every outboard node onto the tip.
    """
    ids = sorted(nodes)
    s = np.array([nodes[n][1] - root_y for n in ids])
    y, z = camber.map_span(s, span)
    return {
        n: (
            nodes[n][0],
            round(float(yy), _SNAP_DECIMALS),
            round(float(zz), _SNAP_DECIMALS),
        )
        for n, yy, zz in zip(ids, y, z, strict=True)
    }


def quad_fraction(mesh: Mesh) -> float:
    """Share of elements that are quad8. 1.0 for an all-quad mesh."""
    if not mesh.elements:
        return 0.0
    quads = sum(1 for conn in mesh.elements.values() if len(conn) == 8)
    return quads / len(mesh.elements)


def check_quad_fraction(
    mesh: Mesh, minimum: float = 0.98, *, what: str = "mesh"
) -> None:
    """Refuse a mesh that has degenerated into a triangle-dominated one.

    Triangles are read, not dropped, so a few of them cost a little local
    bending accuracy rather than a hole. That tolerance is not open-ended: S6 is
    a stiffer, less accurate bending element than S8R, and a mesh that is mostly
    triangles is a different model than the one this repo validates. Recombine
    better, or change the size, rather than raising this number.
    """
    fraction = quad_fraction(mesh)
    if fraction < minimum:
        tris = sum(1 for conn in mesh.elements.values() if len(conn) == 6)
        raise GeometryError(
            f"{what} is {fraction:.1%} quad8 ({tris} triangles of "
            f"{len(mesh.elements)} elements), below the {minimum:.0%} floor. S6 "
            "is a stiffer bending element than S8R; this many of them changes "
            "the model. Adjust the mesh size or recombination instead"
        )


def check_watertight(mesh: Mesh, *, what: str = "mesh") -> None:
    """Refuse a shell mesh with a hole in it.

    Checks the property rather than any particular cause. A quad mesh of one
    connected sheet has exactly one loop of free edges -- its outline. Every
    interior edge is shared by exactly two elements. A second loop is a hole.

    This is deliberately not a check on element types. Holes have arrived here
    from a tile that recombined to quad8 *plus* triangles (only quad8 is read
    back, so the triangles vanish), and can equally arrive from a tile that
    meshes to nothing at all, or from a seam where two tiles fail to share
    nodes. Enumerating those causes is how the first one survived: it was
    checked only when a tile produced zero quad8, so the mixed case -- the one
    that actually happened -- was never tested. ccx solves a holed deck without
    complaint, no node is orphaned by an interior hole, and gmsh's own display
    is complete because gmsh still has the elements that the deck lost. Nothing
    downstream notices. So assert the property, once, on the way out.

    An edge shared by three or more elements is non-manifold and is refused for
    the same reason: it is not a sheet.
    """
    from collections import deque

    shared: dict[frozenset[int], int] = {}
    for conn in mesh.elements.values():
        corners = _corners(conn)
        for i in range(len(corners)):
            edge = frozenset((corners[i], corners[(i + 1) % len(corners)]))
            shared[edge] = shared.get(edge, 0) + 1

    over = [e for e, n in shared.items() if n > 2]
    if over:
        raise GeometryError(
            f"{what} has {len(over)} edge(s) shared by more than two elements; "
            "that is not a single shell sheet"
        )

    adjacent: dict[int, set[int]] = {}
    for edge, count in shared.items():
        if count == 1:
            a, b = tuple(edge)
            adjacent.setdefault(a, set()).add(b)
            adjacent.setdefault(b, set()).add(a)
    if not adjacent:
        raise GeometryError(f"{what} has no free edges at all; it has no outline")

    seen: set[int] = set()
    loops: list[list[int]] = []
    for start in adjacent:
        if start in seen:
            continue
        queue = deque([start])
        seen.add(start)
        loop = []
        while queue:
            node = queue.popleft()
            loop.append(node)
            for other in adjacent[node]:
                if other not in seen:
                    seen.add(other)
                    queue.append(other)
        loops.append(loop)

    if len(loops) > 1:
        interior = sorted(loops, key=len)[:-1]
        where = []
        for loop in interior[:3]:
            xs = [mesh.nodes[n][0] for n in loop]
            ys = [mesh.nodes[n][1] for n in loop]
            where.append(
                f"{len(loop)} nodes near x={sum(xs) / len(xs):.1f}, "
                f"y={sum(ys) / len(ys):.1f}"
            )
        raise GeometryError(
            f"{what} has {len(loops) - 1} hole(s): free edges form {len(loops)} "
            f"loops, not one outline ({'; '.join(where)}). Elements are missing "
            "from the sheet -- ccx would solve the holed part without complaint"
        )


def _zone_elsets(
    elements: dict[int, tuple[int, ...]],
    flat_nodes: dict[int, tuple[float, float, float]],
    zones: tuple[float, ...],
    span: float,
    root_y: float,
) -> dict[str, tuple[int, ...]]:
    """Split elements into spanwise bands at the requested span fractions.

    Called on the flat mesh, where y is the developed distance along the blade,
    so a band means the same length of material whether the blade is bent or not.
    """
    if not zones:
        return {}
    if len(set(zones)) != len(zones) or list(zones) != sorted(zones):
        raise GeometryError(f"zone fractions must be strictly increasing: {zones}")
    if not all(0.0 < f < 1.0 for f in zones):
        raise GeometryError(f"zone fractions must lie in (0, 1): {zones}")

    position = {
        eid: (sum(flat_nodes[n][1] for n in conn[:4]) / 4.0 - root_y) / span
        for eid, conn in elements.items()
    }
    bounds = (0.0,) + tuple(zones) + (1.0,)
    out: dict[str, tuple[int, ...]] = {}
    assigned: set[int] = set()
    for index in range(len(bounds) - 1):
        low, high = bounds[index], bounds[index + 1]
        last = index == len(bounds) - 2
        members = tuple(
            sorted(
                eid
                for eid, frac in position.items()
                if low <= frac < high or (last and frac == high)
            )
        )
        if not members:
            # An empty ELSET is written as a bare *ELSET card with no data lines.
            # ccx accepts that, layup.py hangs a section on it, and the deck
            # solves with the ply drop simply missing -- a converged number that
            # is not the laminate anyone asked for.
            raise GeometryError(
                f"zone_{index + 1} spans {low:g}-{high:g} of the span and caught "
                f"no elements; it is thinner than one element row. Move the "
                f"boundary or refine the mesh."
            )
        out[f"zone_{index + 1}"] = members
        assigned |= set(members)
    if assigned != set(elements):
        raise GeometryError(
            f"zones covered {len(assigned)} of {len(elements)} elements; "
            "every element must land in exactly one zone"
        )
    return out


def _check_normals(
    mesh: Mesh,
    flat_nodes: dict[int, tuple[float, float, float]],
    camber: Camber | None,
    span: float,
    root_y: float,
) -> None:
    """Every element normal must point out of the +z face of the laminate.

    This is the check that stands between a clockwise outline and a model whose
    every unsymmetric stack is upside down. On a bent blade the target is the
    local surface normal ``(0, -sin theta, cos theta)``, taken at the element's
    developed position, which is its flat y.
    """
    for eid, normal in mesh.element_normals().items():
        if camber is None:
            reference = np.array([0.0, 0.0, 1.0])
        else:
            conn = mesh.elements[eid]
            s = sum(flat_nodes[n][1] for n in conn[:4]) / 4.0 - root_y
            theta = float(camber.tangent_angle(np.array([s]), span)[0])
            reference = np.array([0.0, -math.sin(theta), math.cos(theta)])
        # `not (... > 0)` rather than `<= 0`: a degenerate element has a
        # zero-length normal, which comes back NaN, and every comparison against
        # NaN is False. Written the other way this guard waves those through.
        if not float(np.dot(normal, reference)) > 0.0:
            raise GeometryError(
                f"element {eid} has normal {normal}, which faces away from the "
                "laminate's +z side; the first ply line would land on the wrong "
                "face and every unsymmetric stack in the model would be inverted"
            )
