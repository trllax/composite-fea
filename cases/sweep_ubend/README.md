# sweep_ubend

Ply-count sweep on the strip U-bend path, driven by `compfea.sweep`.

## Metric

```
M = 2U / θ      # U from *EL PRINT ELSE
F = M / arm     # arm = undeformed length L
```

Do not use tip |RF|. See `cases/u_bend_path/README.md`.

### `M = 2U/θ` is an assumption, and the results say how good it was

That formula is the *secant* moment. It is exact for a linear torsional spring
(`U = Mθ/2`) and approximate otherwise. The moment that is always right is
Castigliano's tangent `M = dU/dθ`, so every run now reports both:

| column | meaning |
| --- | --- |
| `max_linearity_dev` | worst `(M_sec - M_tan)/M_tan` over the path |
| `linearity_dev_f_90`, `linearity_dev_f_180` | the same at the reported angles |
| `linearity_dev_theta_f_90`, `..._f_180` | the angle each of those was taken at |

The largest angle is always the last row, and there is no central difference at
an endpoint, so its deviation falls back to the nearest interior angle. On the
default 5-degree path that is 175, not 180 -- read
`linearity_dev_theta_f_180` rather than trusting the column name.

Small means `F` really is the force at that angle. Large means `F` is a secant
average across a curved `M(θ)` and should not be read as a force. On the saved
90-degree fin solve the deviation at 90 degrees is 0.14%, while the worst point
on the whole path is 1.1% and sits at θ=2°, where a 1-degree step is a large
*relative* change in θ. Judge a run by the value at the reported angle.

`F_N` stays the secant value so numbers already on disk keep their meaning.

### `f_ratio_180_90` is a diagnostic, not an objective

Under that same linear-spring assumption `M ∝ θ`, so the ratio is **2 by
construction** for every layup -- `tests/test_post.py` asserts exactly that on
synthetic data. The information is in the *deviation* from 2 (the validated
strip run gives 1.961). Do not rank layups on it; rank on `f_90` / `f_180`.

## Design

The grid is the Cartesian product of four axes. Geometry is deliberately not
one of them: this sweep varies the laminate, not the part.

| flag | axis | example |
| --- | --- | --- |
| `--zone-pairs` | repeats of the angle unit per spanwise zone, root first | `3,2,1` |
| `--zones` | span fractions splitting those zones | `0.5` |
| `--angles` | the repeating ply unit | `0/90`, `45/-45` |
| `--ply-mm` | ply thickness | `0.1 0.15` |
| `--fiber` | fibre architecture | `ud woven` |

Each zone is mirrored, so `n` repeats of `0/90` is `[0/90]_ns` -- `4n` plies.
`--n-pairs` is the old single-zone spelling of `--zone-pairs` and still works;
giving both is refused rather than silently resolved.

Zone boundaries must land on element rows, and this is enforced: with 32
spanwise elements a fraction must be a multiple of 1/32, and anything else is
refused with the nearest valid value. A boundary between rows would snap to one
silently, so `--zones 0.47` would deliver the ply drop at 0.5 while the results
row still said 0.47 -- and 0.47/0.48/0.49 would mint three cache keys for one
identical deck. `geometry.mesh_outline` separately raises if a zone catches no
elements at all.

**Fibre types.** `woven` is not a second set of invented constants: it is
derived from the UD card by membrane (A-matrix) equivalence to a `[0/90]` pair
-- see `layup.woven_from_ud`. At equal total thickness the two match in `A` but
not in `D`, because the symmetric UD stack clusters like plies at the outer
fibre (66% higher `D11` at n=2, 33% at n=3). That gap is what a UD-vs-woven
point measures.

There is no `n_plies` column: on a zoned design it would have to pick a zone.
`n_plies_root` / `n_plies_tip` and `thickness_root_mm` / `thickness_tip_mm`
are reported instead.

Mesh: 100×20 mm strip, 32×2 S8R, `long_axis=y`. Path: 5° steps to 180°.

## Caching

`Design.cache_key()` hashes everything that reaches the deck: geometry, mesh,
zones, stacks, angles, ply thickness, the full `*STATIC` line, the reported
angles, **and all nine engineering constants plus density of the material**.

Two ways this went wrong before, both the same shape -- an input that changes
the deck but not the key, so a stale row is served for a deck never solved:

- it hashed the literal string `mat=placeholder_cfrp`, so editing a modulus
  reused the old answer;
- `report_deg` was a function argument rather than part of the design, and it
  sets how many `*STEP` blocks the deck carries, so a cheap 90-degree
  evaluation would cache under the key a 180-degree one looks up.

Everything that reaches the deck now lives on `Design` and is hashed. `deck_for`
is the single deck-building path, shared by `--deck-only` and the solve, so what
you are shown is what gets solved.

## Solver increment: calibrated

`*STATIC` caps the increment size, and that cap sets a *floor* on increments per
angle step. The old `0.02` forced 53 per step on a solve that converges in two
iterations essentially every increment.

`results/incr-calib` solved the same `[0/90]s` design at five caps:

| max increment | wall s | increments | per step | f_90 | f_180 | vs base | speedup |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 25.2 | 223 | 6.2 | 1.759634 | 3.451304 | 0.00e+00 | 5.43x |
| 0.25 | 28.2 | 288 | 8.0 | 1.759634 | 3.451304 | 0.00e+00 | 4.86x |
| 0.1 | 42.3 | 506 | 14.1 | 1.759634 | 3.451304 | 0.00e+00 | 3.24x |
| 0.05 | 66.3 | 865 | 24.0 | 1.759634 | 3.451304 | 0.00e+00 | 2.07x |
| 0.02 | 137.0 | 1908 | 53.0 | 1.759634 | 3.451304 | — | 1.00x |

Bit-identical forces, because a converged elastic equilibrium does not depend on
the increment path. The default is now `0.05, 1.0, 1.E-8, 0.25`: at that cap ccx
uses 8.0 increments per step against a floor of 4, so the cap is no longer
binding and the solver is choosing its own step, with a ceiling still in reserve
for stiffer designs.

Re-calibrate if the design space moves a long way from this one -- the check is
cheap and `--static-line` is a sweep axis, so it runs through the normal machinery.

## Run

From the repo root, with the `composite-fea` env active (`ccx` on PATH):

```bash
# detaches; prints run id
compfea-sweep --n-pairs 1 2 3 --ply-mm 0.1

# a real grid: 2 stacks x 2 thicknesses x 2 fibre types = 8 points
compfea-sweep --zones 0.5 --zone-pairs 2,1 3,1 \
              --fiber ud woven --ply-mm 0.1 0.15

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

Results: `results/<run_id>/results.csv` (and `.parquet`) with `f_90`, `f_180`,
`f_ratio_180_90`, `max_linearity_dev`, and the per-design identity columns.

## Reading the results

`sweep.py` writes the table; `compfea-sweep-post` turns it into a ranking.

```bash
compfea-sweep-post <run_id>              # ranked CSV + figures in results/<run_id>/
compfea-sweep-post <run_id> --shapes 3   # also plot the top 3 deformed shapes
```

| file | what it is |
| --- | --- |
| `sweep_post_rank.csv` | every solved design, sorted by `f_180`, with the linearity columns and a `suspect_f_<deg>` flag beside each force |
| `sweep_post_rank.svg` | one panel per reported angle; a dot at `F` with a whisker of `|linearity_dev| * F` |
| `sweep_post_ratio_check.svg` | `f_180` against `f_90` with the `y = 2x` the model predicts — diagnostic only |
| `sweep_post.json` | mode, counts, varying axes, missing columns, best/worst, and every failed design's error string |
| `sweep_post_shape_<cache_key>.svg` | solved mid-surface against the arc the tip was driven along |

The whisker is the secant-vs-tangent disagreement expressed in the units of the
objective, so it can be read against the gaps between designs: two designs whose
whiskers cover each other's dot are not separable by this sweep. The panel title
names the angle the deviation was actually measured at — on the default
5-degree path the `f_180` value comes from 175, for the endpoint reason above.

### Modes

The tool refuses to present a run as something it is not.

- **calibration** — every design column is constant and only `static_line`
  varies. That is a solver-increment study, not a design comparison, so no
  ranking figure is written. `results/incr-calib` is one: five settings, forces
  identical to the last digit, and the coarsest is 5.4x faster than the finest.
- **single** — one solved design, nothing to compare.
- **design** — everything else.

### Shapes come from the `.frd`, and only from 90 degrees up

`--shapes N` reads DISP out of each design's cached `job.frd`. `build_deck`
asks for `*NODE FILE` at 90 and 180 only, and ccx carries that request forward
once made, so **there is no DISP below 90 degrees at all** and an angle without
one is skipped rather than filled in from a neighbour. These files run 15-130 MB
each, which is why `--shapes` is off by default.

On a spanwise-uniform layup the designs bend into the *same* shape and differ
only in the force needed — stiffness scales `F`, not the deflected curve. The
shape comparison earns its keep on zoned designs, where a ply drop moves
curvature outboard.
