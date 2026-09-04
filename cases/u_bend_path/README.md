# u_bend_path

Stay in **ccx**. 180° tip tangent = **U-shape** (no self-contact).

## Metric (fixed)

Bench force is **always perpendicular to the tip**, scaled by a moment arm:

```
M = 2U / θ          (from *EL PRINT ELSE totals)
F = M / arm         (arm = undeformed length L by default)
```

**Do not use tip |RF|** as F_90 / F_180 — that lab-frame magnitude does not
double with θ. Helper: `metrics.py`.

## Recipe

Tip-edge **prescribed U** on circular-arc poses (no shell UR), S8R ≥32×2,
multi-step ≤5° from 5°→180°, `*EL PRINT, TOTALS=ONLY` / `ELSE`.

## 0.1 mm plies / 0.4 mm total (placeholder CFRP [0/90]s)

| θ | F = M/L (N) | \|RF\| (wrong) |
| --- | --- | --- |
| 90° | **1.75** | 9.46 |
| 180° | **3.43** | 11.73 |
| **ratio** | **1.96 ≈ 2** | 1.24 |

Deck: `u_bend_t0.1mm_energy.inp`. CSV: `results_t0.1mm_Fspring.csv`.
