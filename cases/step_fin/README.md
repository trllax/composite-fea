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

## Tip U-bend runs

```bash
# HEAL mask clamped, tip edge driven; span along +x
python cases/step_fin/run_ubend.py --start-deg 1 --end-deg 90 --step-deg 1 --threads 4
```

Saved 90° solve (do not overwrite): `results/fin_ubend_90_saved/`
(`results/fin_ubend_90/` is the same run). Post reads energy from `ccx/job.dat`.

Both were meshed at 40 mm and are therefore holed -- keep them as a record of
the run, not as a result.
