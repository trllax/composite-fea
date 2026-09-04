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
