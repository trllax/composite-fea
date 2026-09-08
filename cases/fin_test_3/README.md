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

