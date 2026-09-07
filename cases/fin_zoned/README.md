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
or any zone containing it: the root sees all 14, the outboard blade sees 4.

**Ply 1 is the −z ply** — the first line of the `*SHELL SECTION, COMPOSITE`
card. Reversing the file is not a no-op: it leaves `D` unchanged and flips `B`,
so the blade bends identically and draws in the other way. That is the only
place ply order at z shows up in a displacement-controlled bend.

Three properties the file is built to have, and `tests/test_build.py` asserts
rather than trusting:

- **Every zone symmetric about its own mid-plane**, so no zone warps flat.
- **Every UD off-axis ply paired with its negation in the same zone.** The
  woven skins need no partner: a woven ply at 45° carries tows at 45° and 135°
  and is balanced alone. `abd.is_balanced` needs the material kinds to know
  that; without them it reports a woven-skinned laminate as unbalanced.
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

## Not done here

`compfea-build` does not solve. The U-bend path exists in
`cases/step_fin/run_ubend.py`, and this case has not been run through it, so
there is no force on this page — a number here would be one nobody checked.
The fin also needs its own solver increments (`FIN_STATIC_LINE`); the strip's
calibrated cap diverges on this geometry.
