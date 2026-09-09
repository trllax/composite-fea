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

### The knobs

You are always changing the ply book (or, less often, the zone geometry). Two
independent effects:

- **Flex** scales with laminate bending stiffness, roughly as **thickness³** in
  the open blade. To make the whole fin stiffer/softer, add or drop plies (or
  swap a fabric for a stiffer/softer one) **uniformly** across `FULL` / the
  through-going zones. A single `FULL` ply of ~0.2 mm is a large flex move.
- **Kick point** is set by the **stiffness taper** — where the stack steps down
  along the span. The blade bends where it is most compliant:
  - thick, long inboard zones (`HEAL`, `QUARTER`) + thin open blade → kick
    moves **outboard** (toward mid / tip).
  - a stiff `TIP` pad (paired UD) holds the last ~15 % straight → kick moves
    **inboard** off the tip.
  - to move the kick **toward the heel**, shorten the inboard zones or thin
    them so the blade starts bending sooner.
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
- If a change breaks one, fix the change — do not widen a tolerance or move an
  expected value.
