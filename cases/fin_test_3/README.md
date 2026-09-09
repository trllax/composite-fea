# fin_test_3

Realistic geometry STEP (`FIN_TEST_3.step`) without final dimensions.
Named shells: `FULL`, `3_4ths`, `MID`, `QUARTER`, `HEAL`, `TIP`.

**Long axis is +y** (heal near +y, tip near −y). Chord is ±x; camber in z.
Use `--long-axis y`. Mesh at 16 mm for quads (40 mm falls under the 98% floor).

This file is the physics and the validation. **`DESIGNING.md` is the operator
recipe** -- the loop you run to hit a target flex and kick point, by hand
(`run_tipweight.py`) or over a grid (`sweep_layups.py` + `rank_layups.py`, fed a
`TaperDesign` JSON like `designs_example.json`).

## Shop layup (initial)

Fabrics: `twill_3k_198` skins + `hexcel_uni_231` UD only. No `3_4ths` / no `MID`
plies (zones still mesh; coverage comes from nested FULL / QUARTER / HEAL).

| region | thickness | stack idea |
| --- | --- | --- |
| open blade | ~0.63 mm | `[woven / UD / woven]` ≈ 0.7 mm apex target |
| tip | ~1.09 mm | same + paired tip UD pads (high kick) |
| heal | ~1.55 mm | + quarter + heal UD pairs (1.5–1.7 mm) |

```sh
python -c "from compfea.shop_inventory import load_shop_inventory, write_materials_csv as w; w(load_shop_inventory('materials/shop_inventory.csv'), 'materials/from_shop.csv')"

compfea-build --step FIN_TEST_3.step \
              --plybook cases/fin_test_3/plybook_shop.csv \
              --materials materials/from_shop.csv \
              --long-axis y --size-mm 16 \
              --clamp-coverage HEAL \
              --out results/fin_test_3_shop
```

Expected tip force ~5 N at 90° is a bench target for later ubend solves, not
asserted by the smoke build.

## Clamp / tip NSETs

CAD origin is **not** at the heal edge. That is fine:

- `fixed_end` = every node under the `HEAL` mask
- `far_face` = tip edge at the span end **opposite** HEAL along `--long-axis y`
  (ymin here), not max-x

`mesh_step(..., long_axis=...)` and `tip_length_mm` both follow that sense.

## Tip-weight bench (buckling cantilever)

`run_tipweight.py` models the actual bench apparatus, which is **not** a fold and
not the circular-arc `ubend.py` path:

- the blade is clamped at the HEAL and held so its axis is **vertical, tip up**;
- a dead weight hangs from the tip and pulls **straight down** -- *along* the
  undeformed axis, in compression;
- the operator nudges the tip so it buckles the intended (flatwise) way;
- weight is added until the tip **tangent is horizontal** (90 deg) and holds
  there stably. That weight is the bench number.

This is a large-deflection **clamped-free buckling column**. Its small-angle
limit is the Euler load `pi^2 EI / 4L^2`; `W` then rises monotonically with tip
angle (no limit point below 180 deg), so a dead weight can hold any angle up to
there. It is a different deformed shape from both `smoke_cantilever` (transverse
load, curvature at the root) and `ubend.py` (uniform-curvature arc).

### What the model does

`axial_drive_body` (in `deck.py`) ramps the **axial** displacement of a small
tip **clip patch** (a few mid-chord nodes, where a string would tie) under
displacement control, leaving it free to rotate and swing -- free transverse
translation is what keeps this clamped-*free* rather than clamped-pinned (~8x
stiffer). The equivalent hung weight is read directly as the axial reaction `RF`
at the HEAL: no moment read-back, no `M = 2U/theta`, no arm, and `theta` is an
*output*, so ccx's 90 deg prescribed-rotation wall never comes up.

`theta` is the angle between the deformed and undeformed **clip -> `tip_band`**
segments -- two `*NODE PRINT` U stations ~6-45 mm apart, read from the `.dat`
(deck nodes). It is deliberately **not** taken from the `.frd` midsurface:
`deformed_midsurface` returns only a handful of stations over a 650 mm blade and
fans each one into many wherever the shell normal tilts, so a slope fit there is
a chord across half the blade, not a tip tangent (it inflated an earlier
W(90 deg) from ~9 to ~13 N). `run_tipweight` prints a `tangent_ratio` -- deformed
over undeformed segment length -- and warns if it strays from 1, and warns if
`W(theta)` is not monotonic.

The buckling bifurcation is broken by an imperfection. The unsymmetric shop ply
book couples the axial drive into bending on its own, but near Euler that is not
enough for ccx -- pass `--bow-tip-mm` for a small first-mode `z` pre-bow (zero at
the clamp face, clipped to the free length), applied before meshing the section.
Keep it sub-millimetre: it is a `z = f(s)` lift, not the arc-length-preserving
map real camber uses, so it stretches the span as the square of the amplitude.

**Not modelled:** the blade's own weight. On the bench it hangs vertically, so
~2-4 N of blade mass sits in the compression path and, once the blade is over,
becomes a distributed transverse load. Against a total of ~10 N that is a
15-40% systematic offset -- material in a case whose whole point is a ~2x
question. Record the blade mass with the bench number.

Output: `W_vs_theta.csv` (with `tangent_ratio`), `W_vs_theta.svg`, `W`
interpolated at `theta = 90 deg`, and a `dU/d(delta)` Castigliano cross-check.

### Kick point

`--kick-bands N` (default 21) adds a row of `N` mid-chord `*NODE PRINT` U
stations along the span (`compfea.kick.span_band_nsets`, clamp/tip/clip nodes
excluded). After the solve, `kick.flex_kick` builds the deformed centreline from
those, takes `kappa(s)` as the rate of change of the local material tangent
angle with developed span, and reports the arc-length location of the peak of a
lightly smoothed `|kappa(s)|` -- bucketed **heel** (`< 1/3`) / **mid** /
**tip** (`> 2/3`) -- at the increment where the canonical clip->`tip_band`
`theta = 90 deg`, and again at `~15 deg` so migration with load is visible. A
`rotation median` station (where the tip rotation is half done) is reported
alongside as context, not an agreement check. It **warns** when the peak is a
broad plateau or two comparable maxima far apart -- then "kick point" is
ill-posed for that blade. Outputs are `flex_kick.json`, `kick_kappa.svg`
(`|kappa|` vs span) and `kick_shape.svg` (the deformed mid-chord centreline in
side view, coloured by `|kappa|`, kick point ringed at both load states).

The gate is `curvature_profile` / `rotation_median_s` in `compressed_elastica.py`
plus `tests/test_kick.py`: a uniform clamped-free column's curvature is maximal
at the clamp, so its kick point is the **heel**; a synthetic hinge at `0.5 L`
reads **mid**, at `0.78 L` reads **tip**. The metric is band-count converged to
within one bucket from ~11 to ~41 bands on the real fin.

First run, shop plybook + 41 bands: `flex = W(90 deg) = 9.00 N`, `kick = MID`
at `s ~ 298 mm` (`f ~ 0.46`), no warning, and the kick point barely moves
(`~2 mm`) between `theta = 15 deg` and `90 deg`. The thick heel and the tip UD
pads push the peak curvature into the open blade, just inboard of mid-span --
which is the "high kick" the shop layup was aiming for.

### The analytic gate

`compressed_elastica.py` is the closed form (elliptic integrals, pure numpy,
shares no module with `compfea`). `tests/test_compressed_elastica.py` drives a
flat prismatic `[0/90]s` strip to the shortening the closed form predicts for a
target tip angle and checks the reaction is `W(phi_L)`, at 45 / 65 / 85 deg.
`tests/test_kick.py` reuses that strip for the flex + kick-point reduction.

The FE strip carries a deliberate `CircularCamber` bow of radius `1000 L` (tip
offset ~L/2000). A bowed column is genuinely a little soft: the gate lands
**-0.15% to -0.17%** below the perfect-column closed form, **mesh-independent**
from 60 to 120 elements along the span, and scaling linearly with bow amplitude
(L/1500 -> -0.19%, L/1000 -> -0.24%, L/500 -> -0.34%). Tolerance is 0.4%, the
measured offset plus margin -- not what passes. Re-measure before tightening.

### Bench comparison -- unpinned

**Measured (2026-09-08, eyeballed angle):** a 3 lb mass (13.345 N) hung from the
tip curled the blade to ~110 deg. `run_tipweight.py --bench-n 13.345
--bench-deg 110` plots the point and prints the ratio. Not a pinned tolerance --
record the fin serial, the layup, the fixture grip length, the blade mass, and a
*measured* (photo/protractor) tip angle before this becomes a regression case.

```sh
# real solve (~2 min) -- launch detached and poll the run dir
python cases/fin_test_3/run_tipweight.py --uy-frac 0.60 --bow-tip-mm 0.5 \
       --bench-n 13.345 --bench-deg 110 --kick-bands 41
```

First run, shop plybook + `from_shop.csv` IM cards: the clip fixture force `W`
rises monotonically from ~5.8 N near `theta = 0` (the blade's own Euler load)
through **~9.0 N at 90 deg** to **~10.9 N at 110 deg**, `tangent_ratio` =
1.000000 and `dU/d(delta)` within 0.2% of `RF`. The plateau-then-rise shape is
the compressed-column signature; the 90 deg / Euler ratio (~1.55) is just above
the uniform-column closed form's 1.39.

Against the bench: model **10.9 N vs 13.3 N at 110 deg, ratio ~0.82**. That is
the whole story of the "2-3x too weak" this case was opened for -- it was the
*boundary condition* (the circular-arc `ubend.py` path had no compressive
geometric term), not the material or the mesh. What is left:

- **eyeballed 110 deg.** `W` is ~0.4 N/deg here, so +-5 deg is +-2 N -- enough on
  its own to close or open the gap. A measured angle is the first thing to fix.
- **`FIN_TEST_3.step` dimensions are not final** (`EI ~ thickness^3`).
- **self-weight is small here and does not help.** `--gravity 9.81` models it:
  the thin FIN_TEST_3 blade weighs only ~1.8 N, and because that substitutes for
  fixture force in the buckling relation, `W(110 deg)` barely moves (~10.6 N).
  With `--gravity` the `dU/d(delta)` check loosens (gravity does work the clip RF
  does not see) -- expected, not a fault.
- the IM cards (`im_lamina.csv`) are ANSYS-class, not datasheet-fitted.

(An earlier draft read ~13 N at 90 deg from a `.frd`-midsurface slope fit; that
was a chord across 60% of the blade, not a tip tangent -- see above.)

