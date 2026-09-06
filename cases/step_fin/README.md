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

The sweep's `DEFAULT_STATIC_LINE` (max 0.25) **diverges on the fin**. Measured
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
