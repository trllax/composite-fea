# smoke_cantilever

32 elements, 0.9 s of solver time, checked against two closed-form references.

This is the fast gate: it must pass after any change to `deck.py`, `layup.py`,
or `geometry.py`. Run it with `pytest tests/test_smoke_cantilever.py`.

The deck is **built by `layup.py` and `deck.py`**, not checked in. That is the
whole point: a hand-written `.inp` would keep passing after a refactor scrambled
the ply order.

## Model

100 x 20 mm flat strip, 32 S8R elements along the length, 1 across the width
(built by `geometry.py`, not checked in). Axes: x = width, **y = the long axis**,
z = normal, so `long_axis = "y"` and a 0-degree ply runs along the beam.

Material is `PLACEHOLDER_CFRP` from `layup.py`: E1 135000, E2 = E3 9000,
nu12 = nu13 0.30, nu23 0.45, G12 = G13 4500, G23 3000 MPa, density 1.6e-9
tonne/mm^3. Four 0.25 mm plies, 1.0 mm total, in four stacks.

## Load case

Displacement control (see `CLAUDE.md`).

| generic name | this part |
| --- | --- |
| driven reference node | `far_face`, the three nodes on the y = 100 edge |
| driven DOF | 3 (transverse, +z), to 1 mm and to 30 mm |
| free at the tip | DOF 1 and 2, so the strip draws in as it bends |
| reaction node set | `fixed_end`, the y = 0 edge, all 6 DOF fixed |
| compared quantities | `fz` total on `fixed_end`; tip `u_y` for the draw-in |

The tip is free axially on purpose: the reaction is then a pure transverse force
and the classic tip-loaded elastica is the right reference. Constraining DOF 2
would make it a different problem with no closed form.

## Which CLPT reduction — read this before touching a tolerance

There are two ways to turn `D11` into a beam `EI`, and they differ by 0.247% on
this strip:

- **wide plate**, `EI = D11 * b`, which assumes transverse curvature is
  suppressed;
- **narrow strip with free edges**, `EI = b (D11 - D12^2/D22)`, which lets the
  strip curl anticlastically.

The long edges here are free, so the second is correct. It matters: measured on a
converged mesh, ccx sits **+0.004%** from the narrow-strip value and **-0.243%**
from `D11 * b`. That -0.24% is the reduction being left out of the hand
calculation, **not** a property of ccx's shell. `CLAUDE.md` previously recorded it
as a formulation offset; it is not one, and the note there has been corrected.

## Expected values

Four stacks, bottom (-z) ply first, each four 0.25 mm plies:

| stack | D11 (N.mm^2/mm) | EI narrow-strip (N.mm^2) | b11/d11 (mm) |
| --- | --- | --- | --- |
| `[0/90]s` | 9997.4849 | 1.994558e5 | 0 |
| `[90/0]s` | 2074.9497 | 4.139649e4 | 0 |
| `[0/0/90/90]` | 6036.2173 | 1.205546e5 | +0.21875 |
| `[90/90/0/0]` | 6036.2173 | 1.205546e5 | -0.21875 |

`clpt.py` computes the full 6x6 ABD. The first two are symmetric (`B = 0`) and
4.8x apart in stiffness; the last two are each other's reverse, so they have
**identical D and opposite B**.

`elastica.py` is the closed-form large-deflection cantilever under a transverse
tip load, by elliptic integrals: at `delta/L = 0.30` it gives a 26.28 deg tip
angle and `P = 19.8079 N`, 10.34% above the small-deflection 17.9510 N.

### Measured

ccx 2.23, conda-forge build, SPOOLES, `OMP_NUM_THREADS=1`, 2026-09-04. Five
solves, 0.90 s of solver time.

| check | reference | ccx | error |
| --- | --- | --- | --- |
| `[0/90]s`, delta = 1 mm | 0.598367 N (CLPT) | 0.598328 N | **-0.007%** |
| `[90/0]s`, delta = 1 mm | 0.124189 N (CLPT) | 0.124489 N | **+0.242%** |
| stiffness ratio | 4.8182 (CLPT) | 4.8063 | **-0.248%** |
| `[0/90]s`, delta = 30 mm | 19.8079 N (elastica) | 19.8012 N | **-0.034%** |
| reversed pair, force | equal (D is reversal-invariant) | agree to 0.009% | — |
| reversed pair, draw-in | -6.5625e-3 mm (CLPT coupling) | -6.5607e-3 mm | **-0.03%** |

## Mesh

Do **not** transplant the convergence table from `cases/cantilever_ansys`. That
one is a pure end moment, where curvature is constant and a quadratic element is
near-exact at 8 elements. Under a tip load the curvature varies and convergence
is slower. Measured here, error against each reduction at `delta` = 1 mm:

| mesh | solver time | `[0/90]s` vs `D11 b` | vs narrow | `[90/0]s` vs `D11 b` | vs narrow | elastica |
| --- | --- | --- | --- | --- | --- | --- |
| 1x8 | 0.23 s | -0.091% | +0.157% | +0.723% | +0.972% | +0.144% |
| 1x16 | 0.36 s | -0.199% | +0.048% | +0.235% | +0.483% | +0.024% |
| **1x32** | **0.64 s** | **-0.254%** | **-0.007%** | **-0.006%** | **+0.242%** | **-0.034%** |
| 2x32 | 1.28 s | -0.218% | +0.029% | +0.009% | +0.256% | +0.061% |
| 2x64 | 2.55 s | -0.243% | +0.004% | -0.109% | +0.138% | +0.034% |

Two things to read off it. The error against `D11 * b` is monotone and settles at
-0.24% -- the reduction -- while the error against the narrow-strip value settles
at zero, which is what identifies the -0.24% as a reference artifact. And
refining across the width (1x32 -> 2x32) moves the answer 0.02%, so one element
across is enough: the quad8 mid-side node carries the anticlastic curvature.

32 along the length is the shipped mesh: converged, and still under a second.

The mesh used to be a checked-in `mesh_strip_32el.inp`. It now comes from
`geometry.py` (`Outline.rectangle` + `mesh_outline`, `n_chord=1`, `n_span=32`),
which is what makes this case the gate for the mesher that `CLAUDE.md` already
says it is. Swapping the source moved nothing: every measured error in the table
above is identical to four decimal places on the hand-built mesh and the
generated one, and the draw-in check moved from -0.030% to -0.028%.

## Tolerances

| check | tolerance | measured | why |
| --- | --- | --- | --- |
| CLPT force, both stacks, and the ratio | 0.5% | -0.007%, +0.242%, -0.248% | residual is discretization, not solver bias -- it falls monotonically with refinement |
| elastica | 0.3% | -0.034% | `cases/cantilever_ansys` records 0.15-0.22% for this comparison on a finer mesh |
| reversed pair, force | 0.1% | +0.009% | `D` invariance is exact in theory; the band is solver repeatability |
| reversed pair, draw-in | 1% | -0.03% | small-deflection beam theory, so it carries an idealisation error |

That 0.5% is `CLPT_TOL` in the test, 2x the largest measured error. Do not assert a direction: ccx is above the reference on one stack and
below it on another, and which way it sits depends on the mesh.

All of these are 20x or more inside the effect the case is looking for. A
ply-order or long-axis fault moves the answer by 380%, not by 1%.

## What each check catches, and what nothing here catches

| check | catches |
| --- | --- |
| CLPT, both stacks | wrong material card, wrong ply thickness, wrong units, a section that never reached the elements |
| stiffness ratio | any bug that loses which ply sits at which z; a flipped `long_axis` |
| reversed pair, force | the two decks not actually being the same laminate |
| reversed pair, draw-in | **the first ply line not being the -z ply** |
| elastica at delta/L = 0.3 | a wrong large-deflection path |
| NLGEOM check | a `*STEP` that silently lost `NLGEOM` -- ccx accepts it and the answer stays plausible |

Verified by sabotage, each reverted afterwards:

| sabotage | `test_smoke_cantilever` | `test_layup_deck` |
| --- | --- | --- |
| plies emitted top-to-bottom | 1 of 8 fails | 1 of 18 |
| plies sorted by angle | 6 of 8 fail | 2 of 18 |
| `long_axis = "x"` | 6 of 8 fail | — |
| `NLGEOM` dropped | 3 of 8 fail | — |
| orientation mirrored (angle negation) | **none** | 3 of 18 |

### What this case does *not* cover in the mesher

The case builds its mesh with `geometry.py`, so a mesher bug that changes this
strip's stiffness fails here. But it exercises exactly one path: a rectangular
outline, transfinite 1x32, counter-clockwise as constructed, no camber, no zones,
no splines. Everything outside that -- the unstructured path, spline edges, zone
ELSETs, curvature, and the outline-reorientation and normal checks -- is pinned
only by `tests/test_geometry.py`. A green run here is not evidence about any of
them. Both files have to pass before a mesh change is safe.

### The blind spot

**No check in this case can see the angle sign convention.** The mirrored form of
`orientation_card` is exactly angle negation: it emits every +theta ply as -theta.
A stack of 0 and 90 degree plies is completely immune (+90 and -90 are the same
fibre direction), and for any other stack negation changes only the *sign* of
bend-twist coupling, never the magnitude of a reaction. Restoring the mirror
leaves all 8 tests here green.

It is pinned in `tests/test_layup_deck.py` instead, on the direction cosines
themselves, by the invariant that `long_axis="y"` at theta must equal
`long_axis="x"` at theta + 90. Closing it *here* would need an unbalanced
off-axis stack and a signed tip twist, and a reference for that sign that does
not come from ccx -- the ANSYS cross-check in `cases/cantilever_ansys` is the
natural place.

Second, smaller blind spot: reversing a **symmetric** stack is a no-op, correctly
-- it is the same list of plies.

## Files

| file | what it is |
| --- | --- |
| — | the mesh is generated by `compfea.geometry` in `build_deck.mesh_inp()`, so this case gates the mesher too |
| `build_deck.py` | assembles the deck through `layup.py` and `deck.py` |
| `clpt.py` | lamination theory reference: ABD, both reductions, the coupling ratio |
| `elastica.py` | closed-form large-deflection reference |

Both references are deliberately outside `src/compfea/`. A reference that shares
a module with the code under test can be wrong in the same direction and agree.
