# Next-gen fin blade — layup study (2026-09-10)

Branch `exp/nextgen-fin`. Full report (tables, ply-stack diagrams, deformed-blade
shapes, cut sheets): the published artifact
<https://claude.ai/code/artifact/ee343b2d-6441-4609-8475-388929fd5ec6>.
The design playbook this produced is the repo skill
`.claude/skills/fin-layup-design/`.

## Goal

One fin on the unchanged `FIN_TEST_3.step` planform: ~2 N softer than the
current prototype, same mid-span kick, more tip-twist stiffness, ideally no new
fabric purchase.

## Two model changes (committed here)

1. **Cured ply thickness ×1.2.** The built prototype laminate measured ~20 %
   over the `gsm/1000` rule. `shop_inventory.CURED_PLY_FACTOR = 1.2` applies to
   the derived thickness; explicit `thickness_mm` in a CSV is trusted as-is, and
   the shipped `shop_inventory.csv` rows were updated to the ×1.2 cured value.
   Flex ~ thickness³, so the modelled prototype flex rose from ~9 N to 15.5 N.
2. **In-stock UD is intermediate-modulus.** `hexcel_uni_231` / `_379` now carry
   the vetted IM lamina constants (E1 121 → 153 GPa, reusing
   `hexcel_im2_uni_193`; the roll is an ~43 Msi IM fibre). The 3k twill stays
   standard-modulus AS4C.

The 53 solving regression tests still pass.

## Baseline

`plybook_shop_measured.csv` (the prototype ply book, every thickness ×1.2):
**flex 15.5 N, kick MID (f 0.46), migration −2.2 mm, k_twist 16 250 N·mm/rad at
θ_at_probe 118.3°.** Target: flex ~13.5 N, MID kick, k_twist > 16 250 at a
matched θ.

## Result — the ladder

All no-biax, all ±45 `twill_3k_198` skins over a UD (or all-twill) core,
single-ply MID + QUARTER steps, **interleaved so each zone stack is a palindrome
in angle and material → B ≈ 0** (cures flat). Ply books in `ax/`.

| design | flex N | kick | k_twist (matched θ) | B/√(AD) | fabrics |
| --- | --- | --- | --- | --- | --- |
| `ax0`  | 7.1  | MID f0.45 | ~1.3× baseline | ~0 | twill only (cheapest) |
| `ax5r` | 11.2 | MID f0.45 | 1.50× | 0.19† | twill + IM UD |
| `ax3r` | 13.3 | MID f0.45 | 1.56× | ~0 | twill + IM UD |
| **`ax2r`** | **14.1** | **MID f0.45** | **1.63×** | **~0** | **twill + IM UD — recommendation** |
| `ax4`  | 21.2 | MID f0.45 | ~2.2× | 0.32‡ | twill + IM UD (stiff end) |

† `ax5r` keeps a `twill@0` softening ply; the 0°-material mismatch with the UD
is an irreducible B floor. `ax3r` drops it (all-UD) and reaches B = 0 at the
same flex band. ‡ `ax4` and the first-pass `ax1/ax3/ax5(+s)` books (flex
9.9–11.5 N) are **superseded** — they carry their un-interleaved ply order.

`ax2r` stack (ply 1 = −z / mould face):
`FULL tw±45 · QUARTER UD0 · MID UD0 · QUARTER tw±45 · FULL UD0 · FULL UD0 ·
TIP tw±45 · FULL tw±45` → open blade `[45/0/0/45]`, root `[45/0/0/45/0/0/45]`.

## Findings

- The **MID pad places the kick at mid-span** (no MID pad → the fold sits at the
  heel) and is also the **largest flex lever**: a mirrored UD pair on MID ≈ +10 N.
  A symmetric `TaperDesign` sweep can only add mirrored pairs, so it cannot seat
  a MID kick at the soft target — the asymmetric single-ply taper is what fills
  the gap.
- **Full ±45 twill skins are the twist lever** — ~2.4× the k_twist of a
  0/90-skinned blade at the same flex and a matched θ.
- **k_twist only compares at a matched `theta_at_probe`** (~±5°). Re-probe a
  finalist by nudging `--uy-frac` before quoting a ratio.
- **Biax earns nothing** at the current library: `hexcel_himax_biax_100` is a
  stand-in card identical in stiffness to the twill, only thinner. `ax2rb` (its
  biax analog) lands 25 N / 1.31× — worse on every axis. A *measured* HiMax ±45
  lamina card (authored through `materials.py`) would be needed to revisit.
- **Interleaving is not always a ~0.3 N flex move.** `ax2 → ax2r` moved 0.3 N;
  `ax3 → ax3r` moved 1.8 N (its ±45 plies travelled further to the faces).
  Re-solve after reordering.

## Caveats

`FIN_TEST_3.step` is not dimensionally final; flex ~ thickness³ ~ length³, so
the absolute newtons will move — the **trends and ratios** are the deliverable.
Blade self-weight, residual cure stress and double curvature are not modelled;
the B ≈ 0 result is the only guard against mould warp.
