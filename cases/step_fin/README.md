# step_fin

Nested ACP-style STEP import via `compfea.step_mesh.mesh_step`.

## Vertices at ply drops

Overlapping named shells alone do **not** put vertices on the mesh. The
importer runs OpenCASCADE ``fragment`` so coverage boundaries are imprinted
and the mesh gets edges (and vertices) along those curves. No separate CAD
imprint step is required if the shells already meet at the intended drops.

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
python cases/step_fin/run_ubend.py --start-deg 1 --end-deg 90 --step-deg 1 --threads 5
```

Saved 90° solve (do not overwrite): `results/fin_ubend_90_saved/`
(`results/fin_ubend_90/` is the same run). Post reads energy from `ccx/job.dat`.
