# tip_clamp_u_drive

Experiment on branch `exp/tip-clamp-u-drive`: tip-edge clamp via **prescribed U
only** (no shell UR), circular-arc tip-tangent poses, read tip RF TOTALS.

## Results (ccx 2.23, strip from `cantilever_89deg`, 2026-09-03)

| drive | outcome |
| --- | --- |
| single shot θ = 90° | **converged** (58 incs, ~1.4 s). \|F_tip\| ≈ 141.5 N (fy≈−114.6, fz≈−83.0) |
| single shot θ = 180° | **failed** early (t ≈ 0.021) — residual divergence |
| two-step 90° → 180° | step 1 OK; died in step 2 at tot time ≈ 1.36 / 2 |
| single shots 120–170° | all failed early |

So: **U-clamp clears the 90° UR wall for a 90° pose**, but a straight
circular-arc jump (or coarse continuation) to 180° is not yet stable on this
mesh/pose. Next levers: finer mesh, smaller angle steps after 90°, path closer
to the bench (force-dominated clamp motion), or contact if the strip
self-intersects.

Forces above are for the circular-arc *pose*, not fin test correlation.
Target later: \|F_FEA − F_test\| ≤ 2 N on the real fin case.

## Run

```sh
mamba activate composite-fea
python cases/tip_clamp_u_drive/build_and_run.py 90
python cases/tip_clamp_u_drive/build_and_run.py 180
```
