# composite-fea

Parametric large-deflection FEA of laminated composite parts.
gmsh (mesh) -> CalculiX `ccx` (solve) -> pandas (post). CLI only.

Thin laminated shells, loaded well past small-deflection theory. Sweep layups,
get force-displacement curves, pick candidates for layup.

See `CLAUDE.md` for units, modeling conventions, solver rules, what counts as
a valid result, and the working loop. Read the units section before touching a
material card.

## Setup

Python and the solver both come from conda-forge:

```sh
mamba env create -f environment.yml
mamba activate composite-fea
pip install -e . --no-deps
ccx -v                            # must work before anything else
```

## Cases

| case | what it is |
| --- | --- |
| `cases/smoke_cantilever` | 8 elements, <5 s, checked against a closed-form large-deflection answer |
| `cases/fin_20n` | freediving fin, pinned against a physical measurement |
| `cases/cantilever_ansys` | cross-check against ANSYS; ply-stress moment probe |
| `cases/u_bend_path` | the validated tip-U clamp path and its energy metric |
| `cases/tip_clamp_u_drive` | tip-clamp drive kinematics |
| `cases/tip_force_disp` | tip force-displacement |
| `cases/step_fin` | STEP-imported fin U-bend runner |
| `cases/sweep_ubend` | the laminate sweep on the strip U-bend |

Both are regression cases. `CLAUDE.md` explains why neither tolerance moves.

## Modules

| module | what it does |
| --- | --- |
| `geometry.py` | planform outline -> quad8 shell mesh + zone ELSETs -> `.inp` |
| `step_mesh.py` | STEP import, OCC imprint, ACP-style ply coverage masks |
| `layup.py` | design vector -> `*SHELL SECTION, COMPOSITE` |
| `deck.py` | assemble a complete `.inp` |
| `ubend.py` | tip-U clamp path -> multi-step NLGEOM deck |
| `run.py` | run `ccx`, validate convergence, parse `.dat` |
| `metrics.py` | ELSE energy -> secant/tangent moment, linearity deviation |
| `sweep.py` | parameter grid -> process pool -> `results.parquet` |
| `post.py` | one solve dir -> F(θ) CSV + SVG |
| `sweep_post.py` | a sweep's `results.parquet` -> ranked CSV + figures |
| `frd.py` | streaming `.frd` reader; displacements only |
| `shapes.py` | planform, target-arc and deformed-shape plots |

## Sweep (U-bend laminate)

`compfea-sweep` varies the laminate on the fixed 100×20 strip U-bend and
reports tip-normal spring forces **F_90** and **F_180** from ELSE energy
(`F = M/L`, not tip |RF|). The grid is the Cartesian product of ply count per
spanwise zone, ply angles, ply thickness, and fibre type (UD or woven). See
`cases/sweep_ubend/README.md`.

```bash
compfea-sweep --n-pairs 1 2 3 --ply-mm 0.1

compfea-sweep --zones 0.5 --zone-pairs 2,1 3,1 \
              --fiber ud woven --ply-mm 0.1 0.15
```

## Post-processing

Tip U-bend metric is tip-normal spring force from ELSE energy (not tip |RF|):

```
compfea-post results/fin_ubend_90_saved --arm-mm 1002.26
```

Writes `post_F_theta.csv`, `post_F_theta.svg`, and `post_meta.json` in the run dir.

`M = 2U/θ` is a secant identity that is exact only for a linear spring, so the
output also carries the Castigliano tangent `dU/dθ` and `max_linearity_dev`,
which says how far apart the two ran. Read `F` as a force only where that
deviation is small.

### Whole sweeps

`compfea-post` handles one solve directory. `compfea-sweep-post` reads a
sweep's `results.parquet` and ranks the grid:

```bash
compfea-sweep-post <run_id>              # ranked CSV + figures in results/<run_id>/
compfea-sweep-post <run_id> --shapes 3   # also plot the top 3 deformed shapes
```

Writes `sweep_post_rank.csv`, `sweep_post_rank.svg`, `sweep_post.json`, and a
`sweep_post_ratio_check.svg`. Ranking is on `f_90` / `f_180`; `f_ratio_180_90`
is ~2 by construction and is refused as a sort key. Each force carries its
`linearity_dev` as a whisker in the same units, so two designs whose whiskers
cover each other's point are visibly not separable.

A run whose design columns are all constant (only `static_line` varying) is
reported as a **solver calibration** rather than a design comparison, and no
ranking figure is written.

`--shapes N` reads the cached `.frd` for the top N designs and plots the real
deformed mid-surface against the circular arc the tip was driven along. Those
files run 15-130 MB each, so it is off by default. Shapes exist only from the
first angle the deck asked `*NODE FILE` for — typically 90° — and an angle with
no DISP is skipped rather than substituted.
