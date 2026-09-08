# fin_test_3

Realistic geometry STEP (`FIN_TEST_3.step`) without final dimensions.
Named shells: `FULL`, `3_4ths`, `MID`, `QUARTER`, `HEAL`, `TIP`.

**Long axis is +y** (heal near +y, tip near −y). Chord is ±x; camber in z.
Use `--long-axis y`. Mesh at 16 mm for quads (40 mm falls under the 98% floor).

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

### The analytic gate

`compressed_elastica.py` is the closed form (elliptic integrals, pure numpy,
shares no module with `compfea`). `tests/test_compressed_elastica.py` drives a
flat prismatic `[0/90]s` strip to the shortening the closed form predicts for a
target tip angle and checks the reaction is `W(phi_L)`, at 45 / 65 / 85 deg.

The FE strip carries a deliberate `CircularCamber` bow of radius `1000 L` (tip
offset ~L/2000). A bowed column is genuinely a little soft: the gate lands
**-0.15% to -0.17%** below the perfect-column closed form, **mesh-independent**
from 60 to 120 elements along the span, and scaling linearly with bow amplitude
(L/1500 -> -0.19%, L/1000 -> -0.24%, L/500 -> -0.34%). Tolerance is 0.4%, the
measured offset plus margin -- not what passes. Re-measure before tightening.

### Bench comparison -- unpinned

Only a weight and a rough (eyeballed) 90 deg angle have been recorded, so this is
a comparison, not a pinned tolerance. `run_tipweight.py --bench-n <N>` draws the
line on the plot and prints the ratio. The IM shop cards (`from_shop.csv`)
already closed most of the earlier gap; record the weight, the fixture grip
length, the fin serial and the layup here, and the tip angle if it can be
measured, before this becomes a regression case.

```sh
# real solve (~2 min, 65 increments) -- launch detached and poll the run dir
python cases/fin_test_3/run_tipweight.py --bow-tip-mm 0.5 --bench-n <N>
```

First run, shop plybook + `from_shop.csv` IM cards, `--bow-tip-mm 0.5`, 65
increments: `W` rises monotonically from ~5.8 N near `theta = 0` (the blade's own
Euler load) to **~9.0 N at `theta = 90 deg`**, `tangent_ratio` = 1.000000
throughout and `dU/d(delta)` within **0.2%** of `RF`. The plateau-then-rise shape
is the compressed-column signature; the 90 deg / Euler ratio (~1.55) is just
above the uniform-column closed form's 1.39 -- the zoned heel stiffens the root a
little, not by 2x. Not pinned -- see above. (An earlier draft read ~13 N here
from a `.frd`-midsurface slope fit; that was a chord across 60% of the blade, not
a tip tangent.)

