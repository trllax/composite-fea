# fin_zoned

The design front door, on a real part: a STEP of named zones, a ply book, and a
materials library. Nothing about the laminate is on the command line.

```sh
compfea-build --step test_fin_2.step \
              --plybook cases/fin_zoned/plybook.csv \
              --materials materials/generic.csv \
              --long-axis x --size-mm 16 \
              --out results/fin_zoned
```


## Shop-stock demo

Scarce inventory only (`hexcel_uni_231` + `twill_3k_198`; no 379),
blank-E filled from default lamina cards via `materials/from_shop.csv`:

```sh
python -c "from compfea.shop_inventory import load_shop_inventory, write_materials_csv as w; w(load_shop_inventory('materials/shop_inventory.csv'), 'materials/from_shop.csv')"

compfea-build --step test_fin_2.step \
              --plybook cases/fin_zoned/plybook_shop.csv \
              --materials materials/from_shop.csv \
              --long-axis x --size-mm 40 \
              --out results/fin_shop_demo
```

Or `python cases/fin_zoned/run_ubend.py` (defaults to that plybook / materials;
add `--end-deg 5` for a cheap tip-force smoke before a full 90° hunt).


## Long axis: +x

A 0-degree ply runs along the span, root to tip, which on this STEP is **+x**.
`layup.py` takes `long_axis` with no default and `compfea-build` requires the
flag, because this repo's cases disagree — `cases/cantilever_ansys` and the
strip sweep use +y — and a default is silently right for half of them.

## The zones

Six named `SHELL_BASED_SURFACE_MODEL`s, nested root to tip, plus one outboard:

| zone | x range (mm) | area (mm²) | elements @16 mm |
| --- | --- | --- | --- |
| `HEAL` | 0 – 172.66 | 37502.6 | 205 |
| `QUARTER` | 0 – 300 | 73733.4 | 415 |
| `HALF` | 0 – 450 | 118733.4 | 657 |
| `z_3_4ths` | 0 – 700 | 193733.4 | 1032 |
| `FULL` | 0 – 1174.92 | 336210.5 | 1786 |
| `TIP` | 974.92 – 1174.92 | 60000.0 | 331 |

`3_4ths` starts with a digit, so it is sanitized to `z_3_4ths` for ccx; write
either in the ply book, the loader sanitizes what you write.

`HEAL` is also the clamp: `fixed_end` is every node under it, not just the
root edge, because the fin is bonded into a foot pocket rather than knife-edge
clamped. `--clamp-coverage ''` falls back to the min-x edge.

These areas come from the OCC fragment map and match the CAD faces to 0.007%.
Read them out of `zone_report.csv` and check them against the CAD **before**
spending a solve — see `cases/step_fin/README.md` for what the old bounding-box
mask did to this same geometry.

## The ply book

14 plies. A ply covers its zone, so an inboard zone carries every ply naming it
or any zone containing it: the root sees 13 of the 14 (ply 12 is TIP-only), the outboard blade sees 4.

**Ply 1 is the −z ply** — the first line of the `*SHELL SECTION, COMPOSITE`
card. Reversing the file is not a no-op: it leaves `D` unchanged and flips `B`,
so the blade bends identically and draws in the other way. That is the only
place ply order at z shows up in a displacement-controlled bend.

Three properties the file is built to have, and `tests/test_build.py` asserts
rather than trusting:

- **Every zone symmetric about its own mid-plane**, so no zone warps flat.
- **Every UD off-axis ply paired with its negation in the same zone.** The
  woven skins need no partner *at 45°*: for a balanced weave `Q11 == Q22`, which
  reduces `Qbar16` to `(Q11 - Q12 - 2 Q66)·cs(c² - s²)` — zero at 0°, 45° and
  90° and nowhere else. A weave at 22.5° would genuinely be unbalanced, its tows
  lying at 22.5° and 112.5°. `abd.is_balanced` reads this off `A16`/`A26`, which
  is ACP's definition and needs no material metadata.
- **Drops buried.** Plies 1, 2, 13, 14 run everywhere, so the outer surface is
  continuous and every drop sits under the skin.

## What comes out

| file | what it is |
| --- | --- |
| `zone_report.csv` | elements, area, extent per zone — check against CAD |
| `stack_table.csv` | resolved stack per zone, ply by ply, with z stations |
| `laminate_abd.csv` | A, B, D, flexural moduli, symmetric/balanced/D16 flags |

| `mesh.inp`, `layup.inp` | the deck pieces |
| `build.json` | inputs, counts, layup fingerprint, warnings |

Resolved stacks at 16 mm, root to outboard:

| section | t (mm) | D11 (N·mm) | D16 (N·mm) | stack |
| --- | --- | --- | --- | --- |
| `cov_5` (HEAL) | 2.300 | 81758.2 | 1426.1 | `[45/0/0/45/-45/0/0/0/-45/45/0/0/45]` |
| `cov_4` (QUARTER) | 2.000 | 52735.0 | 998.2 | `[45/0/0/45/-45/0/0/-45/45/0/0/45]` |
| `cov_3` (HALF) | 1.600 | 26172.9 | 427.8 | `[45/0/0/45/-45/-45/45/0/0/45]` |
| `cov_2` (3_4ths) | 1.000 | 5193.6 | 0.0 | `[45/0/0/0/0/45]` |
| `cov_1` (FULL) | 0.700 | 1413.6 | 0.0 | `[45/0/0/45]` |
| `cov_6` (TIP) | 0.800 | 2295.2 | 0.0 | `[45/0/0/0/45]` |

Flexural moduli come out as `ef_1_mpa` / `ef_2_mpa`, **not** `ef_x` / `ef_y`:
direction 1 is the laminate's 0° direction, which is whichever global axis
`long_axis` names. Here that is x, but for a `long_axis="y"` case — half this
repo — direction 1 is global y, and a column called `ef_x` would be read as the
wrong axis. Every row carries `long_axis` so the mapping travels with the CSV.

| section | ef_1 (MPa) | ef_2 (MPa) |
| --- | --- | --- |
| `cov_5` (HEAL) | 65754 | 23187 |
| `cov_4` (QUARTER) | 63784 | 23360 |
| `cov_3` (HALF) | 60672 | 23615 |
| `cov_2` (3_4ths) | 42199 | 23660 |
| `cov_1` (FULL) | 25602 | 20436 |
| `cov_6` (TIP) | 31200 | 22010 |

`D16` is not zero on the zones carrying ±45 UD pairs, and that is correct. A
balanced symmetric angle-ply laminate still has bend–twist coupling: the +45
and −45 plies sit at different z, so they cancel in membrane and not in
bending. It is 1.7% of `D11` at the root here. Reported, not gated.

## Validating against ANSYS

`laminate_abd.csv` is the comparison to make **first**, because it needs no
solver on either side. Every number in it is one ACP reports for a laminate, so
it establishes that ccx and ANSYS were handed the same laminate before anyone
argues about a force. Compare per zone: total thickness, `A`, `B`, `D`, and the
symmetric/balanced flags.

`D16`'s **sign** is worth checking specifically. `layup.py` notes that writing
the mirrored `*ORIENTATION` form flips bend–twist coupling and that no force
check in this repo can catch it. `D16` can, and it is the only quantity here
that can.

Only then compare solved quantities — `F_90`/`F_180` from ELSE energy, and the
per-ply peaks from `compfea.stress`.

## Does the deck solve?

Yes. A smoke run — 5 steps of 1 degree, `0.001, 1.0, 1.E-10, 0.1`, one core —
reaches TOT TIME 5.0 and exits clean, so both of `CLAUDE.md`'s validity
conditions hold. **This is the first deck in this repo to carry two
`*MATERIAL` cards**, and ccx takes six `*SHELL SECTION, COMPOSITE` blocks over
one `blade` ELSET without complaint.

| theta | U (N·mm) | M = 2U/theta (N·mm) | F = M/arm (N) |
| --- | --- | --- | --- |
| 1° | 3.462e-1 | 39.68 | 0.0396 |
| 2° | 1.472e+0 | 84.33 | 0.0841 |
| 3° | 3.390e+0 | 129.50 | 0.1292 |
| 4° | 6.110e+0 | 175.05 | 0.1747 |
| 5° | 9.634e+0 | 220.79 | 0.2203 |

Arm 1002.26 mm (tip minus the outboard `HEAL` station). `M` increments by
44.7, 45.2, 45.6, 45.7 — near-linear with about 11% geometric stiffening from
1° to 5°, which is what early large-deflection behaviour looks like.

**These are not validated forces.** 5 degrees is the start of a 180-degree
path, the material is `materials/generic.csv` — textbook constants, explicitly
not a datasheet — and nothing here is pinned to a measurement. The number the
run establishes is that the deck is well formed and converges, not what the fin
does.

## Not done here

`compfea-build` does not solve; the U-bend path lives in
`cases/step_fin/run_ubend.py` and this case is not wired to it. The full
180-degree run also needs the fin's own increments (`FIN_STATIC_LINE`) — the
strip's calibrated cap diverges on this geometry — and it is much heavier than
the smoke run above, because ccx expands each S8R into one solid per ply and
the root zone here is 13 plies deep.
