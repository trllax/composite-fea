# step_fin

Nested ACP-style STEP import via `compfea.step_mesh.mesh_step`.

## Vertices at ply drops

Overlapping named shells alone do **not** put vertices on the mesh. The
importer runs OpenCASCADE ``fragment`` so coverage boundaries are imprinted
and the mesh gets edges (and vertices) along those curves. No separate CAD
imprint step is required if the shells already meet at the intended drops.

## Mesh size: 15 mm, and why

`mesh_step` refuses a tile that recombines to quad8 **plus** triangles. Only
quad8 is read back, so those triangles would be dropped, leaving a hole in the
part that ccx solves without complaint.

On `test_fin_2.step` this is not hypothetical. Measured, deterministic across
repeats:

| `--size-mm` | 15 | 16 | 16.5 | 17 | 17.25 | 17.5 | 18 | 20 | 40 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| elements | 1922 | 1786 | 1762 | 1566 | — | — | — | — | — |
| result | clean | clean | clean | clean | tris | tris | tris | tris | tris |

**17.0 mm is the coarsest that comes out all quad8**, and the failing band
starts immediately above it at 17.25. The default is **16.0**: still clean,
7% lighter than 15, and far enough from the cliff that a small change in the
STEP or gmsh version does not land on it. The earlier runs at `--size-mm 40`
produced a mesh
with a one-element hole in it -- its free edges form two loops, the outline and
a 5-node interior loop -- so `results/fin_ubend_90_saved` (F_90 = 0.3719 N) and
the 180-degree run were both solved on a holed blade. Those numbers need
redoing at 16 mm before they mean anything.

The check that catches this used to be gated behind "this tile produced zero
quad8", so the mixed case -- the one that actually occurs -- was never tested.
Note that no node is orphaned by the dropped elements, because an interior
triangle shares every node with its neighbours; counting orphan nodes does not
detect this. Free-edge loop counting does.

Nor does inspecting the mesh in gmsh: gmsh keeps the triangles, so its display
is complete and correct (`results/test_fin_2.msh` holds 335 quad8 **plus** the
2 tri6). The hole appears only in the `.inp`, because `mesh_step` reads back
`getElementsByType(quad8)` alone. Compare element counts between the two files,
or open `deck.inp` in gmsh rather than the `.msh`.

## The guard is on the property, not the cause

`mesh_step` refuses a tile carrying non-quad8 elements, and separately calls
`geometry.check_watertight`, which requires the free edges of the quad mesh to
form exactly one loop -- the outline. A second loop is a hole; an edge shared by
more than two elements is not a sheet.

Both checks exist because the type check alone is not enough. It enumerates one
cause. A tile that meshes to **zero** elements passes it and still drops a hole,
and so would a seam whose two tiles fail to share nodes. `check_watertight` does
not care which of those happened.

Note that counting orphan nodes does **not** detect this, and neither does
looking at the mesh in gmsh -- an interior element shares every node with its
neighbours, so deleting it orphans nothing, and gmsh still holds the elements
the deck lost. Both are pinned in `tests/test_geometry.py`.

## Coverage

Each Onshape product name (`FULL`, `HALF`, `TIP`, …) becomes an ELSET mask.
An element may belong to several masks. `3_4ths` is sanitized to `z_3_4ths`
for CalculiX.

Drop the STEP at the repo root (or pass a path):

```python
from compfea.step_mesh import mesh_step
mesh = mesh_step("test_fin_2.step", size_mm=40)
```

Note: in `test_fin_2.step`, `HALF` and `QUARTER` share the same x-span, so
their ELSETs match until the CAD spans differ.

## Covered plies -> COMPOSITE

```python
from compfea.step_mesh import mesh_step
from compfea.layup import Ply, coverages_from_mesh, layup_from_coverage, mesh_elsets_for_stacks
from compfea.geometry import Mesh
from compfea.deck import assemble

mesh = mesh_step("test_fin_2.step", size_mm=40)
plies = [
    Ply(0.15, 0.0, coverage="FULL"),
    Ply(0.15, 90.0, coverage="FULL"),
    Ply(0.10, 0.0, coverage="z_3_4ths"),
    Ply(0.20, 0.0, coverage="TIP"),
]
cov = coverages_from_mesh(mesh.elsets)
layup, stacks = layup_from_coverage(plies, cov, long_axis="x")  # span is +x on this STEP
mesh = Mesh(
    nodes=mesh.nodes, elements=mesh.elements, nsets=mesh.nsets,
    elsets=mesh_elsets_for_stacks(stacks, all_elements=mesh.elsets["blade"]),
    heading=mesh.heading,
)
deck = assemble(mesh_inp=mesh.to_inp(), layup=layup, initial_bc="*BOUNDARY\nfixed_end, 1, 6")
```

Use the sanitized mask names from the mesh (`z_3_4ths`, not `3_4ths`).

## Triangles are kept, not dropped

Recombination does not always give all quads on a fragmented tile. Those strays
are now emitted as **S6** alongside the S8R, in the same ELSET under one
`*SHELL SECTION, COMPOSITE` -- verified against ccx, which accepts both mixed.
`mesh_step` refuses anything that is neither, and refuses a mesh below
`quad_floor` (default 98%), because S6 is a stiffer bending element than S8R and
a triangle-dominated mesh is a different model.

Dropping them was never cheap. On the strip, deleting one interior element of
128 -- 0.8% of the area -- moved the reported force **2.5%**, because a hole
severs load path rather than just removing material.

| `--size-mm` | elements | tri6 | quads |
| --- | --- | --- | --- |
| 40 | 337 | 2 | 99.4% |
| 30 | 540 | 2 | 99.6% |
| 20 | 1174 | 4 | 99.7% |
| 16 | 1786 | 0 | 100% |

So 40 mm is usable again: 337 elements against 1786 at 16 mm. Pick the size on
mesh convergence, not on triangle avoidance -- that study has not been done.

## Solver increments: the strip's calibration does not transfer

The sweep's `DEFAULT_STATIC_LINE` (max 0.5, and 0.25 when this was measured)
**diverges on the fin**. Measured
here at 40 mm with 1-degree steps, minimum increment held at `1.E-10`:

| max increment | result |
| --- | --- |
| 0.25 | diverges at ~88% of step 1 |
| 0.1 | converges, 21.5 s, 44 increments (**adopted**) |
| 0.01 | converges, 51.8 s, 211 increments (the old value) |

2.4x here, not the strip's 5.4x. The minimum increment matters as much as the
maximum: it is what lets ccx cut back through the first step, where a flat blade
takes its first bend.

**It is not the triangles.** None of the 46 diverging nodes belongs to an S6
element (0 of 296 mentions), and the all-quad 16 mm mesh diverges at the same
point. Calibrate per model; do not carry a solver setting across geometries.

## Tip U-bend runs

```bash
# HEAL mask clamped, tip edge driven; span along +x
python cases/step_fin/run_ubend.py --start-deg 1 --end-deg 90 --step-deg 1 --threads 4
```

Saved 90° solve (do not overwrite): `results/fin_ubend_90_saved/`
(`results/fin_ubend_90/` is the same run). Post reads energy from `ccx/job.dat`.

Both were meshed at 40 mm and are therefore holed -- keep them as a record of
the run, not as a result.

## How a named shell becomes an ELSET (measured, 2026-09-07)

Zone membership used to be a **1-D bounding box**: an element belonged to a
named shell if its centroid's `x` fell inside that shell's x-range. Names came
from a regex scan of the STEP text, geometry came from OpenCASCADE, and the two
were never joined by entity identity. On this very STEP that produced wrong
zones for three of six shells:

| shell | old x-span mask | true extent | true area (mm^2) |
| --- | --- | --- | --- |
| `FULL` | 0 .. 1174.92 | 0 .. 1174.92 | 336223.808 |
| `3_4ths` | 0 .. 700.00 | 0 .. 700.00 | 193746.626 |
| `HALF` | 0 .. **674.92** | 0 .. **450.00** | 118746.626 |
| `QUARTER` | 0 .. **674.92** | 0 .. **300.00** | 73746.626 |
| `TIP` | **674.92** .. 1174.92 | **974.92** .. 1174.92 | 60000.000 |
| `HEAL` | 0 .. 172.66 | 0 .. 172.66 | 37502.639 |

`HALF` and `QUARTER` came out with **identical** ELSETs, and the `TIP` mask was
2.5x its true area. The cause is `_collect_points`: it walks the entity graph
collecting `CARTESIAN_POINT`s, which for a B-spline face are **control points**,
not the trimmed face. Their hull overhangs the real surface, so the derived
x-range is an over-estimate — silently, and by 225 mm here.

### What replaced it

The STEP is regular: each `SHELL_BASED_SURFACE_MODEL` names an `OPEN_SHELL`
holding an ordered list of `ADVANCED_FACE`s.

```
FULL    #483 -> 7 faces  #462..#468      QUARTER #486 -> 3 faces  #475..#477
3_4ths  #484 -> 5 faces  #469..#473      HALF    #487 -> 4 faces  #478..#481
HEAL    #485 -> 1 face   #474            TIP     #488 -> 1 face   #482
```

21 faces, and OCC imports exactly 21 surfaces as tags 1..21 **in declaration
order**. `occ.fragment` returns `outDimTagsMap`, parallel to `objects + tools`,
whose entry *i* lists the tiles input face *i* became. Composing the two gives
name -> tiles -> elements exactly, with no geometric tolerance anywhere.

### Three things probed, and what each was worth

- **`Geometry.OCCImportLabels = 1` does not carry the names.** Only 2 of the 6
  arrive (`Shapes/HEAL`, `Shapes/TIP`), and OCC attaches them to *every*
  geometrically coincident face — 5 tags are labelled `HEAL`, including faces
  that belong to `FULL`, `HALF`, `QUARTER` and `3_4ths`. Useless as the primary
  mapping. Useful as a check: the tiles of a labelled tag must sit inside the
  tiles assigned to that name.
- **Area conservation across `fragment` has no discriminating power.** It looks
  like a good invariant and it is a tautology: `fragment` conserves area, so
  *any* block partition of the inputs conserves it. Swapping `HEAL` and `TIP`
  passes with `rel = 0.000e+00`. Do not re-add this check believing it tests
  something.
- **Control-point hull containment is the label-free check that works.** By the
  convex-hull property a B-spline face lies inside its control points, so the
  bbox of the tiles assigned to a shell must sit inside that shell's
  control-point bbox. It is a superset test, so it cannot false-alarm; measured
  overhang on the correct assignment is exactly `0.0000` for all six shells.
  It catches a `HEAL`/`TIP` swap (1002 mm overhang), a one-step rotation
  (675 mm), and a `FULL`/`3_4ths` block swap (275 mm).

`shell_x_ranges_mm` is kept — it is still how the names are read — but it is a
**hull**, not a mask, and nothing decides membership from it any more.

### What the old mask did to this fin's laminate

`default_plies()` resolved through `layup_from_coverage` at `--size-mm 16`,
before and after. Same 1786 elements, same six shells, same ply list:

| | zones | stacks |
| --- | --- | --- |
| **old** (x-range mask) | 5 | `[0/90/0]` 754 el · `[0/90/0/0]` **40 el** · `[0/90/0/90]` 784 el · `[0/90/0/90/0]` **3 el** · `[0/90/0/90/45]` 205 el |
| **new** (fragment map) | 5 | `[0/90]` 423 el · `[0/90/0]` 375 el · `[0/90/0/90]` 452 el · `[0/90/0/90/45]` 205 el · `[0/90/0]` 331 el (TIP) |

The two tiny zones in the old column — 40 elements and **3 elements** — are not
design intent. They are the ±2 mm centroid band where the `HALF` and `TIP`
hulls overlapped, and they put an extra 0-degree ply on a 3-element sliver.
The `TIP` reinforcement also covered about 2.5x the area it was drawn on, and
`HALF` and `QUARTER` were the same set, so a ply on one was a ply on both.

Any force previously computed from this STEP was computed on that laminate.
