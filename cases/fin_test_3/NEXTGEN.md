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

1. **Cured ply thickness ×1.1.** The built prototype laminate measured ~1.2× the
   dry `gsm/1000` rule; **1.1 is the adopted compromise** — part real ply loft,
   part resin. `shop_inventory.CURED_PLY_FACTOR = 1.1` applies to the derived
   thickness only (an explicit `thickness_mm` is trusted as-is); the shipped
   `shop_inventory.csv` rows carry the ×1.1 value. **The card moduli and density
   are not scaled with it**, so this factor moves bending stiffness directly
   (`D ~ E·t³`) — a full resin-pickup model would also drop `E1` ~1.2× and give
   a smaller net stiffening. Flag for a future per-fabric measurement + Vf
   re-derivation.
2. **In-stock UD is intermediate-modulus.** `hexcel_uni_231` / `_379` now carry
   the vetted IM lamina constants (E1 121 → 153 GPa, reusing
   `hexcel_im2_uni_193`; the roll is an ~43 Msi IM fibre). The 3k twill stays
   standard-modulus AS4C.

The solving regression gates still pass; the baseline is unaffected by the UD
change (its ply book uses IM2/IM7/biax, not `hexcel_uni_231`).

## Baseline

`plybook_shop_measured.csv` (the prototype ply book, every thickness ×1.1),
`--uy-frac 0.60`: **flex 11.96 N, kick MID (f 0.46), migration −2.3 mm,
k_twist 12 774 N·mm/rad at θ_at_probe 118.3°.**
Target: **flex ~10 N** (baseline − 2), MID kick, k_twist > 12 774 at a matched θ.

## Result — the buildable ladder

All no-biax, all ±45 `twill_3k_198` skins over a UD (or all-twill) core,
single-ply MID + QUARTER steps, **interleaved so each zone stack is a palindrome
in angle and material → B ≈ 0** (cures flat). k_twist quoted at a drive
re-matched to θ ≈ 118° (baseline's), the only condition it compares under.
Ply books in `ax/`.

| design | flex N | kick | k_twist (θ≈118) | ×base | B/√(AD) | fabrics / note |
| --- | --- | --- | --- | --- | --- | --- |
| `ax2r1` | 3.56 | MID f0.45 | ~9 960 | 0.78× | ~0† | ultra-soft, ax2r with a single UD core ply — **below the twist floor** |
| `ax0`   | 5.49 | MID f0.45 | ~16 600 | 1.30× | ~0 | super-soft, **twill only** (cheapest, single fabric) |
| `ax5r`  | 8.67 | MID f0.45 | 19 050 (θ117) | 1.49× | 0.19‡ | twill + IM UD |
| `ax3r`  | 10.29 | MID f0.45 | ~20 300 | 1.59× | ~0 | twill + IM UD — closest to target |
| **`ax2r`** | **10.85** | **MID f0.45** | **~20 640** | **1.62×** | **~0** | **twill + IM UD — recommendation** |
| `ax4`   | 16.33 | MID f0.45 | ~30 750 | 2.41× | 0.32§ | twill + IM UD (stiff end) |

† `ax2r1` root (QUARTER∩MID) carries a 0.068 residual; open blade / MID / QUARTER
regions are B = 0. ‡ `ax5r` keeps a `twill@0` softening ply — the 0°-material
mismatch with the UD is an irreducible B floor; `ax3r` drops it (all-UD) for
B = 0 at the same flex. § `ax4` and the first-pass `ax1/ax3/ax5(+s)` books are
**superseded** — un-interleaved ply order, B ≈ 0.3.

`ax2r` stack (ply 1 = −z / mould face):
`FULL tw±45 · QUARTER UD0 · MID UD0 · QUARTER tw±45 · FULL UD0 · FULL UD0 ·
TIP tw±45 · FULL tw±45` → open blade `[45/0/0/45]`, root `[45/0/0/45/0/0/45]`.
`ax2r1` is the same minus one `FULL UD0` ply.

Migration is −1.4 to −2.6 mm on every candidate: the kick point barely moves
between light load and 90°.

## Findings

- The **MID pad places the kick at mid-span** (no MID pad → the fold sits at the
  heel) and is also the **largest flex lever**: a mirrored UD pair on MID ≈ +8 N;
  a single UD ply on `FULL` moved `ax2r → ax2r1` by ~7 N. A symmetric
  `TaperDesign` sweep can only add mirrored pairs, so it cannot seat a MID kick
  at the soft target — the asymmetric single-ply taper is what fills the gap.
- **Full ±45 twill skins are the twist lever** — ~2.4× the k_twist of a
  0/90-skinned blade at the same flex and a matched θ.
- **k_twist only compares at a matched `theta_at_probe`** (~±5°). Re-probe a
  finalist by nudging `--uy-frac` before quoting a ratio. The 1.1× / 1.2× factor
  change moved every absolute flex ~−23% but the **k_twist ratios were
  unchanged** (baseline and candidate scale together).
- **Biax earns nothing** at the current library: `hexcel_himax_biax_100` is a
  stand-in card identical in stiffness to the twill, only thinner. `ax2rb` (its
  biax analog) lands 19.5 N / 1.34× — worse on every axis. A *measured* HiMax
  ±45 lamina card (authored through `materials.py`) would be needed to revisit.
- **Interleaving is not always a ~0.3 N flex move.** `ax2 → ax2r` moved 0.3 N;
  `ax3 → ax3r` moved ~1.5 N (its ±45 plies travelled further to the faces).
  Re-solve after reordering.

## Caveats

`FIN_TEST_3.step` is not dimensionally final; flex ~ thickness³ ~ length³, so
the absolute newtons will move — the **trends and ratios** are the deliverable.
The ×1.1 cured factor is a modulus-free compromise (see model change 1); a
per-fabric micrometer + Vf pass would firm the absolute numbers. Blade
self-weight, residual cure stress and double curvature are not modelled; the
B ≈ 0 result is the only guard against mould warp.
