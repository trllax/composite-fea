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
