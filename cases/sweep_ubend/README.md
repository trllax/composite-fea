# sweep_ubend

Ply-count sweep on the strip U-bend path, driven by `compfea.sweep`.

## Metric

```
M = 2U / θ      # U from *EL PRINT ELSE
F = M / arm     # arm = undeformed length L
```

Do not use tip |RF|. See `cases/u_bend_path/README.md`.

## Design

Symmetric cross-ply `[0/90]_ns` at fixed ply thickness (default 0.1 mm):

| `--n-pairs` | stack | plies |
| --- | --- | --- |
| 1 | `[0/90]s` | 4 |
| 2 | `[0/90]_2s` | 8 |
| 3 | `[0/90]_3s` | 12 |

Mesh: 100×20 mm strip, 32×2 S8R, `long_axis=y`. Path: 5° steps to 180°.

## Run

From the repo root, with the `composite-fea` env active (`ccx` on PATH):

```bash
# detaches; prints run id
compfea-sweep --n-pairs 1 2 3 --ply-mm 0.1

# foreground (tests / short demos)
compfea-sweep --sync --n-pairs 1 --ply-mm 0.1

# decks only (no solve)
compfea-sweep --deck-only --n-pairs 1 2
```

Poll a detached run:

```bash
jq . results/<run_id>/status.json
tail -n 40 results/<run_id>/sweep.log
```

Results: `results/<run_id>/results.csv` with `f_90`, `f_180`, `f_ratio_180_90`.
