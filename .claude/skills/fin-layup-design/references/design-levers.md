# Design levers — measured on FIN_TEST_3

All numbers from the 2026-09-10 study: `FIN_TEST_3.step` planform, tip-weight
buckling bench, `--uy-frac 0.60`, cured ply thickness = `gsm/1000 * 1.1`,
materials `hexcel_uni_231` (IM UD, E1 153 GPa, 0.277 mm) + `twill_3k_198`
(AS4C 3k 0/90 woven, E1 59 GPa, 0.238 mm). **Absolute newtons move with the
planform; use these as ratios and directions.**

Zones, free-span reach from the clamp (heel):
`QUARTER` 0–15 % · `MID` 0–39 % · `z_3_4ths` 0–68 % · `TIP` 84–100 % ·
`FULL` whole blade · `HEAL` inside the fixture (its plies do nothing).

## Cured ply thickness

`shop_inventory.CURED_PLY_FACTOR = 1.1` — the built prototype laminate measured
~1.2× the dry `gsm/1000` rule; **1.1 is the adopted compromise** (part real ply
loft, part resin). Applied to the *derived* thickness only; an explicit
`thickness_mm` in the CSV is trusted as-is. **The card moduli and density are
NOT scaled with it**, so since flex ~ `E·t³` this factor moves stiffness
directly — a full resin-pickup model would also drop `E1` ~1.2× (A₁₁ ≈ E·t
conserved) and give a smaller net stiffening. Going 1.0 → 1.1 lifted the
modelled prototype flex ~9 N → 12 N; 1.2 gave ~15.5 N. **Re-measure and
re-derive the card (Vf) per fabric when you have micrometer data.**

## Flex (`W` at 90° tip tangent)

- Scales roughly as **through-section thickness³** in the open blade, so it is
  brutally sensitive near the target: dropping **one `FULL` UD ply** took
  `ax2r → ax2r1` from 10.9 N to 3.6 N on a ~1 mm blade.
- **One `hexcel_uni_231` (IM UD) ply on `FULL` ≈ +5 to +8 N** at this thickness.
- **One `twill_3k_198 @ 0` ply ≈ +2–4 N** (softer per ply — lower E1, thinner).
- Dropping the **`QUARTER` UD stiffener** (0–15 % span, part-clamped) ≈ **−1 N**,
  kick unmoved — the fine adjustment.
- The symmetric `TaperDesign` sweep can only add plies in **mirrored pairs**, so
  its flex quantum is ~2 UD plies. Single-ply asymmetric taper halves that step.

## Kick point

- **Set by the stiffness taper**, not the overall thickness. The blade folds
  just outboard of the last padded zone.
- **No `MID` pad → the kick sits at the HEEL** (f ≈ 0.21). A `MID` pad of any
  weight moves it to **MID** (f ≈ 0.45–0.46), and it stays there dead stable
  from ~12 N to ~55 N of flex.
- The `MID` pad is also the **biggest single flex lever**: a mirrored UD pair on
  `MID` ≈ **+10 N**; a single UD ply ≈ +4–6 N; a single `twill@45` or
  `himax_biax@45` ≈ +2–4 N. Use the lightest pad that still buckets MID.
- A `TIP` pad holds the last sixth straight → pulls the kick **inboard** off the
  tip. Pad only `QUARTER` → kick stays at the heel.
- Migration with load (`migration_delta_mm`) is small (−2 to −3 mm) for every
  design that reached a clean MID bucket — a predictable fin.

## Tip-twist stiffness (`k_twist`)

- **Comparative only, and only across layups reaching a similar
  `theta_at_probe_deg`.** Re-probe the finalist to θ ≈ the baseline's (nudge
  `--uy-frac`) before quoting a ratio.
- Rewards `D66` / `A66` → ±45 and transverse plies.
- **Full ±45 skins are the lever.** `twill_3k_198 @ 45` skins front+back over a
  UD bending core give **k_twist ~2.4× a 0/90-skinned blade of the same flex**,
  at a matched θ, for ~+1.5 N of flex vs 0/90 skins. A single ±45 twill ply is
  already balanced (0/90 woven rotated 45°), so it also keeps `D16 = 0`.
- **Biax earns nothing today.** `hexcel_himax_biax_100` in the library is a
  *stand-in* — its stiffness card is byte-identical to `twill_3k_198`, only
  thinner (0.12 vs 0.24 mm) and with no ±45 baked in (angle comes from the ply
  book). So biax = a thinner twill: the thinner ±45 shell carries *less*
  torsion, needs more UD to hold flex, and ends up heavier and less twist-stiff
  than the no-biax `twill@45`-skin design. A **measured HiMax ±45 lamina card**
  (a real ±45 architecture has ~5–10× the in-plane shear of a rotated 0/90
  woven) would be needed to revisit — author it through `materials.py` with a
  declared `# units:` line.

## The winning architecture (2026-09-10)

**±45 `twill_3k_198` skins · IM UD spanwise core · single-ply MID + QUARTER
steps · interleaved for B ≈ 0.** Family in `cases/fin_test_3/ax/`:

Numbers at `CURED_PLY_FACTOR = 1.1`; baseline (prototype ×1.1) is flex 11.96 N,
kick MID, k_twist 12 774 at θ 118.3. k_twist ratios are θ-matched to ~118°.

| design | flex N | kick | ×base k_twist | B/√AD | note |
| --- | --- | --- | --- | --- | --- |
| ax2r1 | 3.56  | mid | 0.78× | ~0 (root .07) | ultra-soft, ax2r − one FULL UD ply; below the twist floor |
| ax0   | 5.49  | mid | 1.30× | ~0 | super-soft, single-fabric (twill only) |
| ax5r  | 8.67  | mid | 1.49× | 0.19 | keeps a twill@0 softener (B floor) |
| ax3r  | 10.29 | mid | 1.59× | ~0 | all-UD 0°, interleaved — closest to a −2 N target |
| **ax2r** | **10.85** | **mid** | **1.62×** | **~0** | **the recommendation** |
| ax4   | 16.33 | mid | 2.41× | 0.32 | stiff end |

`ax2rb` (biax analog of ax2r) lands 19.5 N / 1.34× — the biax overshoot
described above. The 1.1×→1.2× factor change moves every flex ~+30 % but leaves
the k_twist ratios unchanged.

## Solver / workflow traps

- `ccx` exits `3221225781` without the full conda PATH (missing DLL).
- `run_tipweight.py` default `--uy-frac` is 0.45; sweeps use 0.60 — match them.
- No `--step` / `--long-axis` on either script (hard-wired to FIN_TEST_3, +y).
- `jq` absent — parse `status.json` in python.
- A thin soft open blade under a heavy stiff MID pad → `ccx` "too many
  cutbacks". That taper is genuinely bad; do not chase it.
- A run is valid only if `ccx` exits 0 **and** the `.sta` reached total time
  1.0 (`run.py` enforces this and raises rather than returning a partial).
