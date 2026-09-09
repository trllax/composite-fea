# Designing a fin blade to a target flex and kick point

This is the operator recipe. The physics, the bench, and the validation live in
`README.md`; the repo-wide conventions live in `../../CLAUDE.md`. Read those once.
This file is the loop you actually run.

## What you get out

Two numbers per candidate layup, from one solve:

| number | what it is | range |
| --- | --- | --- |
| **flex** | `W`, the dead weight (N) hung from the tip that bends the tip tangent to 90 deg on the buckling bench | ~5 N supersoft → ~15–20 N hard |
| **kick point** | where the blade curves most: `f`, the fraction of free span, bucketed **heel** (`f < 1/3`) / **mid** / **tip** (`f > 2/3`) | a location, plus a migration with load |

Both come from `results/<run>/flex_kick.json`. Two SVGs go with it:
`kick_kappa.svg` (curvature vs span) and `kick_shape.svg` (the bent blade in
side view, coloured by curvature, kick point ringed).

## What you provide

Three files, same as `compfea-build` everywhere:

1. **A STEP** with the named zone shells. `FIN_TEST_3.step` has
   `FULL / 3_4ths / MID / QUARTER / HEAL / TIP`. Long axis is **+y** (heal near
   `+y`, tip near `−y`), chord `±x`, camber `z`. Mesh at 16 mm.
2. **A ply book CSV** — one row per ply in global stacking order, each naming
   the zone it covers (`plybook_shop.csv`). Ply 1 is the `−z` (showing) face. An
   element's stack is the subsequence of plies whose zone covers it, so a ply
   drop is a subsequence and continuity across a zone boundary is automatic.
3. **A materials CSV** — lamina cards in the deck's units (MPa, tonne/mm³),
   produced from `materials/shop_inventory.csv`:

   ```sh
   python -c "from compfea.shop_inventory import load_shop_inventory, \
     write_materials_csv as w; \
     w(load_shop_inventory('materials/shop_inventory.csv'), 'materials/from_shop.csv')"
   ```

## The loop

```sh
# 1. edit the ply book (or the STEP zones, or the materials)

# 2. solve + measure  (~2 min; does NOT self-detach -- launch it yourself)
nohup python cases/fin_test_3/run_tipweight.py \
        --plybook cases/fin_test_3/plybook_shop.csv \
        --materials materials/from_shop.csv \
        --uy-frac 0.60 --bow-tip-mm 0.5 --kick-bands 41 \
        --run-dir results/<name> > results/<name>.log 2>&1 &

# 3. poll
tail -f results/<name>.log
# ... look for:  flex = W(90 deg) = X.X N   kick = MID  s_kick = ... (f = 0.NN)

# 4. read results/<name>/flex_kick.json and the two SVGs

# 5. adjust the layup, go to 1
```

## Sweeping many layups at once

Instead of hand-editing one ply book, describe a family of them as a
**`TaperDesign`** and let `sweep_layups.py` solve the grid. A `TaperDesign` is
the same knobs as below, made explicit:

- **`skins`** -- plies per face on the through-going zone (`FULL`), the
  continuous outer surface;
- **`core`** -- a half-stack on `FULL`, mirrored about the mid-plane (the flex
  knob);
- **`pads`** -- a half-stack per inboard/tip zone, mirrored (the kick knob).
  `pad_order` is `-z` → mid-plane order, i.e. through-thickness position.

Everything is **symmetric by construction** (`B ≈ 0`, no warp off the mould) and
every ply thickness is the stocked value -- you pick material and angle, not
thickness. What this cannot express is an odd ply count in a pad zone or an
arbitrary unsymmetric order; hand-write a ply book for those.

```sh
# designs_example.json: a "base" design + "axes" whose Cartesian product is the
# grid. Axis keys are top-level fields (core, skins, ...) or "pads.<ZONE>".
python cases/fin_test_3/sweep_layups.py cases/fin_test_3/designs_example.json \
       --uy-frac 0.60 --bow-tip-mm 0.5 --kick-bands 41 --jobs 6
# prints a run id and detaches (like compfea.sweep). Poll:
tail -f results/<run-id>/sweep.log
jq . results/<run-id>/status.json          # state: running -> done (or error)
```

Each design is one `run_tipweight.py` solve (single core; `--jobs` at once).
Results land in `results/<run-id>/results.parquet` -- one row per **distinct
laminate** (designs that resolve identically are collapsed), with the flex/kick
numbers and the design parameters flattened. The cache dir is keyed on the
resolved laminate **and** the solve settings (`--uy-frac`, `--bow-tip-mm`,
`--kick-bands`, `--size-mm`, the materials file, the STEP), so re-running with an
added grid point is cheap and changing a solve setting correctly re-solves.

Then rank against a target:

```sh
python cases/fin_test_3/rank_layups.py results/<run-id>/results.parquet \
       --target-flex 12 --target-kick tip           # or --target-kick-frac 0.8
```

`dist = hypot((flex_n - target_flex)/flex_tol, kick_s_frac - target_frac)` with
`--flex-tol` (default 2 N) setting the trade rate. Rows with a kick-point
warning, an error, or no 90° are dropped and **counted** in the printout.
Writes `ranked.csv` + a scatter SVG (flex vs kick fraction, target starred,
top-N labelled). `--require-bucket heel|mid|tip` is a **hard filter** -- drop
designs not in that bucket before scoring -- for when the bucket is a
requirement, not a preference (the distance metric alone lets a close flex
match outrank a whole bucket). There is **no change-cost term yet** -- a later
version will rank "move a zone boundary < change a ply count < swap a fabric".

### Which zone a pad tiers -- FIN_TEST_3

The kick point is where the stiff root region **steps down**, so it lands just
outboard of the last padded zone. On `FIN_TEST_3.step`, measured from the clamp
as a fraction of the 649 mm free span:

| zone | covers free span | a pad here puts the kick near |
| --- | --- | --- |
| `HEAL` | inside the fixture (`fixed_end`) | nothing -- those plies are clamped |
| `QUARTER` | 0 – 15 % | f ≈ 0.21 (**heel**) |
| `MID` | 0 – 39 % | f ≈ 0.45 (**mid**) |
| `z_3_4ths` | 0 – 68 % | further out; needs a thin `FULL` core to read **tip** |
| `TIP` | 84 – 100 % | holds the last sixth straight (keeps the kick inboard) |
| `FULL` | whole blade | the flex knob; does not move the kick |

Measured on an 18-design sweep (`designs_example.json`): `MID` unpadded →
**heel** at f ≈ 0.21 for every flex level; `MID` padded → **mid** at
f ≈ 0.45 – 0.46, dead stable from 12 N to 56 N of flex; `MID` + `z_3_4ths`
padded **with a thin (`[UD]`) `FULL` core** → **tip** at f ≈ 0.71. Flex and
kick are close to separable once you pad the right zone: pick the bucket with
the pad reach, then trim flex with the `FULL` core count.

### The knobs

You are always changing the ply book (or, less often, the zone geometry). Two
independent effects:

- **Flex** scales with laminate bending stiffness, roughly as **thickness³** in
  the open blade. To make the whole fin stiffer/softer, add or drop plies (or
  swap a fabric for a stiffer/softer one) **uniformly** across `FULL` / the
  through-going zones. A single `FULL` ply of ~0.2 mm is a large flex move.
- **Kick point** is set by the **stiffness taper** — how far out the padded
  root region reaches (see the zone table above). The blade folds just outboard
  of the last padded zone.
  - pad out to `MID` → kick at mid; add `z_3_4ths` (light) over a thin `FULL`
    core → kick at tip; pad only `QUARTER` → kick stays at the heel.
  - a `TIP` pad holds the last sixth straight, so it pulls the kick **inboard**
    off the tip.
  - a pad on `HEAL` does nothing here — it is inside the fixture.
  - moving a **zone boundary** in the STEP is the cheapest kick move; changing
    a **ply count** in a zone is next; **swapping a fabric** (UD ↔ woven, or to
    one the shop does not stock) is the expensive one. (A future cost function
    will rank these; for now, prefer them in that order by hand.)

Changing the taper changes flex a little and changing thickness changes the kick
a little — they are not fully separable, so expect 2–3 iterations to hit both.

## Reading the outputs

`flex_kick.json`:

```jsonc
{
  "flex_n": 9.00,                 // the number of record; W at 90 deg
  "free_len_mm": 649.1,
  "headline": {                   // evaluated at theta = 90 deg
    "bucket": "mid",
    "s_kick_mm": 298, "s_kick_frac": 0.46,
    "s_rotmed_frac": 0.53,        // context only, not an agreement check
    "plateau_frac": 0.17,         // width of the >=80%-of-peak band / span
    "theta_deg": 90.0,            // band-centreline tip tangent (transparency)
    "warning": null               // <-- if not null, see below
  },
  "migration": { ... },           // same, at theta ~ 15 deg (light load)
  "migration_delta_mm": -2.4      // headline s_kick minus migration s_kick
}
```

- **`bucket` + `s_kick_frac`** is the answer. `s_rotmed_frac` is a different
  statistic (where half the tip rotation has accumulated); it legitimately sits
  outboard of the curvature peak and is there for context, not to match.
- **`migration_delta_mm`** near zero means the kick point is stable under load
  (good — a predictable fin). A large value means it walks up the span as you
  load it.
- **`warning` not null** means the curvature peak is a broad plateau or is
  bimodal — "kick point" is not a single location for that blade. Open
  `kick_kappa.svg` and look at the actual `|kappa(s)|` curve; the layup probably
  needs a sharper taper somewhere.

`kick_kappa.svg` — `|kappa(s)|` (curvature, 1/mm) from clamp to tip at both load
states, with the kick point and rotation median marked. This is the diagnostic:
a clean single hump = a well-defined kick; a flat top or two humps = ambiguous.

`kick_shape.svg` — the deformed mid-chord centreline in side view (span vs rise),
the line coloured by local `|kappa|`, the kick point ringed at both load states.
This is the human picture: you can see the fold sit in the open blade.

## Tip twist stiffness (`--twist-probe`)

Two blades can land on the same flex and kick point and still differ in how the
loaded blade resists a **twist** — the thing biax / ±45 / transverse plies buy.
`run_tipweight.py --twist-probe` appends a small `+F / −F` force couple on the
tip chord edges of the buckled blade and writes `twist.json` with
`k_twist_nmm_per_rad = couple / twist_angle` (a moment about global z, over the
deformed chord).

- **Comparative only.** No closed form; the absolute value is not meaningful.
  Larger is stiffer is better. Compare `k_twist` only between layups run the same
  way (same `--twist-cload-n`, same `--size-mm`) **and that reach a similar
  `theta_at_probe_deg`** — the z-couple is a pure twist only near 90°, so a wide
  θ spread across the compared set makes the ranking meaningless. Check the
  `theta_at_probe_deg` column before using `--twist-tiebreak`.
- It rewards `D66` / `A66`: swapping a 0 or 90 ply for a ±45 pair should raise
  `k_twist` sharply. If it does not, look hard at the layup before trusting it.
- A **force** couple, not a prescribed rotation, because an unsymmetric layup has
  rolled the tip edge by ~90° so its nodes sit at unknown axial positions. Tune
  `--twist-cload-n` so the reported `twist_angle` lands in ~0.1..2°.
- `--twist-probe` does **not** change the drive `--uy-frac`; the probe is read
  wherever the drive step ends (θ recorded).
- **Warnings** in `twist.json` (`warning` field): `theta_at_probe` outside
  45..135° (`k_twist` not usable); twist angle below 0.03° (then
  `k_twist_nmm_per_rad` is **`null`** — raise `--twist-cload-n`) or above 3°
  (lower it); leading/trailing edge motion very asymmetric (strong bend-twist
  coupling / roll — the comparison is soft); large heel RF swing. A warned or
  `null`-`k_twist` row is **kept** by the ranker but does **not** compete on
  twist in `--twist-tiebreak` — it keeps its `dist` order.
- Non-convergence of the twist step means the section has no positive torsional
  stiffness there. It fails loudly (no `twist.json`); do not work around it.

In the sweep: `sweep_layups.py --twist-probe` adds `k_twist_nmm_per_rad`,
`twist_angle_deg`, `theta_at_probe_deg`, `twist_le_te_asymmetry`,
`twist_warning` to `results.parquet` (roughly 2× the solve cost; adding
`--twist-probe` also changes `solve_phash`, so the first twist run re-solves
every design). `rank_layups.py --twist-tiebreak` then re-orders rows in the same
`--tie-eps`-wide `dist` bin by descending `k_twist` — stiffer in torsion wins
when flex and kick are a wash.

## What kappa is (one paragraph)

`kappa` is curvature — 1 / (local bend radius), in 1/mm. It is built per span
station: the deformed centreline's local tangent angle `phi` (measured off that
segment's *undeformed* direction, so initial camber reads zero), differentiated
along the developed span, `kappa = dphi/ds`. `phi` climbs from ~0 at the clamp
to the tip angle; `kappa` is how fast it climbs. The kick point is the peak of a
lightly smoothed `|kappa(s)|`, refined by a parabola fit, evaluated at the
increment where the canonical clip→`tip_band` angle is 90 deg.

## Gotchas

- **`--kick-bands`.** 21 is enough for the bucket; 41 pins `s_kick_frac` to
  about ±0.02. The deck and the post read the same value, so don't re-analyse an
  old `.dat` with a different count (it now raises rather than silently
  truncating).
- **`--uy-frac`** must drive the tip past a 90 deg tangent or `flex_n` comes
  back `null` (90 deg not bracketed). ~0.39 L reaches 90 deg for the shop layup;
  0.60 is safe headroom. A much stiffer layup may need more.
- **`--bow-tip-mm 0.5`** seeds the buckle. An unsymmetric ply book usually
  starts on its own; a symmetric or very stiff one needs the bow. Keep it
  sub-millimetre — it is a `z = f(s)` lift, not the arc-length-preserving camber
  map, so it stretches the span as the square of the amplitude.
- **`abd.py` first.** Run `abd` on the deck before spending a solve — it prints
  A/B/D per zone and is the only check that catches a mirrored angle convention
  (`D16`). A ply-book typo that halves a zone's thickness produces a clean,
  wrong flex number.
- **`zone_report()`** — check each zone's meshed area against the CAD before the
  first solve on a new STEP.
- **What is not modelled:** the blade's own weight (pass `--gravity 9.81` to
  add it; on this thin blade it barely moves `W`), double curvature, twist about
  the long axis, and any residual cure stress. The STEP dimensions here are not
  final and `EI ~ thickness³`, so absolute numbers will move; trends between
  layups are the usable output today.

## The gates — do not weaken these

- `pytest tests/test_compressed_elastica.py tests/test_kick.py` must pass after
  any change to `deck.axial_drive_body`, `src/compfea/kick.py`, or
  `compressed_elastica.py`. The first checks the buckling reaction against the
  elliptic-integral closed form on a flat strip; the second checks the flex +
  kick reduction on top of it (uniform column → heel; synthetic hinge at
  `0.5 L` → mid; at `0.78 L` → tip). There is deliberately no skip on a missing
  `ccx`.
- `pytest tests/test_twist.py` gates the tip twist metric: the `twist_couple_body`
  cards, the chordwise `twist_le` / `twist_te` split, `k_twist` (and the `None`
  case) from a hand-written `.dat`, and end to end that a `[±45]s` strip is more
  than 1.5× the tip twist stiffness of a `[0]4` strip of equal thickness
  (observed ~2×), repeatable under `n_span` + `n_chord` refinement. Also no skip
  on a missing `ccx`.
- If a change breaks one, fix the change — do not widen a tolerance or move an
  expected value.
