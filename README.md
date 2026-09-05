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

Both are regression cases. `CLAUDE.md` explains why neither tolerance moves.

## Status

Skeleton only. `src/compfea/` holds the package marker and nothing else;
modules land in this order, with a review after each:

1. `run.py`
2. `layup.py`, `deck.py`
3. `geometry.py`
4. `sweep.py`, `post.py`

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
