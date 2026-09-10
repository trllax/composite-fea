# materials/

Lamina cards and shop stock for `compfea-build` / ply books.

| file | role |
| --- | --- |
| `generic.csv` | Textbook UD + derived woven cards already in the repo |
| `ematerials_all.xml` | Raw ANSYS Workbench Engineering Data offload; source for `ansys_composites.csv` |
| `ansys_composites.csv` | The ANSYS Epoxy Carbon library flattened to lamina cards (`# units: stress=Pa density=kg/m3`), plus two `im7_*_276` IM-class rows adapted from it. Paste more ANSYS rows here. |
| `shop_inventory.csv` | What you actually stock: product, architecture, tow, gsm and/or oz/yd², thickness, moduli |

## Units

Every library file that carries elastic constants **must** start with a
`# units: stress=… density=…` line. `compfea.materials` converts into
**MPa** and **tonne/mm³** and range-checks. No units line is refused.

`shop_inventory.csv` uses shop-friendly **GPa** and **g/cm³**. A future
loader can emit a `materials.csv` row in those units (or convert first);
until that lands, copy a filled shop row into a materials file by hand
and keep the `# units:` line honest.

## Shop → ply book

Design parameters are discrete:

- **zones / coverage** and **ply count / angle** in the ply book
- **thickness** only from `shop_inventory.thickness_mm` (or, when blank,
  `areal_weight_gsm / 1000 * CURED_PLY_FACTOR`)
- **material** name must match a materials library `name`

Outputs that matter for the fin (F_90 and kick locus) are post metrics,
not columns here.


### Weight and tow

- **`tow`**: filament count label only (`3k`, `6k`, `12k`, …). Provenance for
  the roll; not a deck input.
- **`width_in`**: roll or panel width in inches. Cutting / nesting only;
  not a deck input.
- **Areal weight**: fill `areal_weight_gsm` and/or `areal_weight_oz_yd2`
  (oz per square yard). Conversion used by a loader:
  `gsm = oz_yd2 * 33.9057`.
- **Thickness**: the as-laid **cured** ply thickness. When `thickness_mm` is
  blank the loader derives `areal_weight_gsm / 1000 * CURED_PLY_FACTOR`
  (`shop_inventory.CURED_PLY_FACTOR`, currently **1.1** — the built FIN_TEST_3
  prototype measured ~1.2× the dry `gsm/1000` rule on 2026-09-10; 1.1 is the
  adopted compromise, part real ply loft and part resin). The card moduli and
  density are NOT scaled with it, so since flex ~ `E·t³` this factor moves
  stiffness directly — keep it small. An explicit `thickness_mm` is trusted
  as-is (no factor); the shipped rows carry the ×1.1 cured value. Prefer a
  per-fabric micrometer + Vf re-derivation over the uniform factor.

### Architecture → ply angle

| architecture | plybook `kind` | meaning of `angle_deg` |
| --- | --- | --- |
| `ud` | `ud` | fibre direction |
| `0_90` | `woven` | warp along this angle, weft at +90 |
| `biax` | `woven` | fibres at ±45 to this angle; use `0` for bisector on the long axis |

`biax` is one fabric SKU in inventory, not two UD plies in the book.

Fibre tensile modulus / strength on a roll label is **not** a lamina card.
Put those in `fibre_e_msi` / `fibre_xt_ksi` (never copied into `e1`). Leave
lamina `e1`… blank to use the default UD / woven card, or fill measured
lamina engineering constants when you have them.

## Default lamina when E is blank

`compfea.shop_inventory.load_shop_inventory` fills missing lamina elastic
fields from ``UD_CFRP_GENERIC`` (``ud``) or ``woven_from_ud`` (``0_90`` /
``biax``). That is enough to run a deck on scarce stock; it is **not** a
datasheet. ``fibre_e_msi`` / ``fibre_xt_ksi`` stay out of ``e1``.

```python
from compfea.shop_inventory import load_shop_inventory, write_materials_csv
skus = load_shop_inventory("materials/shop_inventory.csv")
write_materials_csv(skus, "materials/from_shop.csv")  # for compfea-build
```
