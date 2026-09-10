---
name: fin-layup-design
description: >-
  Design or retune a laminated fin-blade layup to a target flex, kick point and
  tip-twist stiffness on the cases/fin_test_3 tip-weight buckling bench. Use when
  asked to hit a flex/kick/twist target, sweep layups, build an asymmetric
  single-ply taper, cut a fin's cure warp (B coupling), or produce a fin-design
  report with cut sheets. Assumes the fin method in cases/fin_test_3 already
  exists (DESIGNING.md, README.md, run_tipweight.py, sweep_layups.py).
---

# Fin layup design loop

The repo-wide conventions are in `../../CLAUDE.md`; the bench physics and the
manual loop are in `cases/fin_test_3/README.md` and `DESIGNING.md`. **Read those
once.** This skill is the design *playbook* on top of them — the levers, the
traps, and the finishing steps that a green solve does not cover — kept out of
`CLAUDE.md` so it only loads when you are actually designing a fin.

## Environment

`ccx` needs the **full conda env PATH** or it exits `3221225781` (missing DLL):

```sh
CF=$(conda info --base)/envs/composite-fea      # or the env you built
export PATH="$CF:$CF/Library/mingw-w64/bin:$CF/Library/usr/bin:$CF/Library/bin:$CF/Scripts:$CF/bin:$PATH"
export OMP_NUM_THREADS=1
PY="$CF/python"
```

`jq` is not installed — read `status.json` with `python -c`. `sweep_layups.py`
and `run_tipweight.py` have **no `--step` / `--long-axis`** (hard-wired to
`FIN_TEST_3.step`, long axis +y). `run_tipweight.py` default `--uy-frac` is
**0.45**; the sweep uses **0.60** — match them or the flex/twist numbers are not
comparable.

## The loop

1. **Baseline.** Re-solve the reference layup **at the current cured ply
   thickness** with the twist probe:
   ```sh
   $PY cases/fin_test_3/run_tipweight.py --plybook <ref>.csv \
       --materials materials/from_shop.csv --uy-frac 0.60 --bow-tip-mm 0.5 \
       --kick-bands 41 --twist-probe --run-dir results/baseline
   ```
   Record `flex_n`, `kick_bucket` + `kick_s_frac`, `migration_delta_mm`,
   `k_twist_nmm_per_rad`, and **`theta_at_probe_deg`**. Cured thickness =
   `gsm/1000 * shop_inventory.CURED_PLY_FACTOR` — if a physical layup was
   measured, that factor is the calibration (see `references/design-levers.md`).

2. **Targets.** Flex = baseline ± your delta (`--flex-tol` sets the trade rate,
   default 2 N). Kick bucket. **k_twist > baseline, compared only at a matched
   `theta_at_probe`** (± ~5°) — the twist couple is a pure twist only near a
   90° tip tangent, so a wide θ spread makes `k_twist` meaningless.

3. **Coarse search — `sweep_layups.py` (`TaperDesign` grids).** Fast, but
   symmetric-only and every pad is a **mirrored pair**. Consequence: the
   smallest MID step it can build is 2 plies and that pad is a **~+10 N flex
   lever**, so a MID-kick design at a soft target may not exist in the grid.
   Add a light `QUARTER` pad + thin the through-section to force a clean MID
   kick everywhere, then read where flex lands.

4. **Fill the gap — asymmetric single-ply taper.** When the target sits between
   symmetric grid points (the flex quantum is ~one UD ply, 6–11 N), hand-write
   a ply-book CSV like `cases/fin_test_3/plybook_shop.csv`: one row per ply in
   `-z -> +z` order, each naming its zone. `cases/fin_test_3/ax/*.csv` are the
   worked examples — ±45 twill skins, UD spanwise core, **single-ply** MID /
   QUARTER steps. Loop `run_tipweight.py --twist-probe` over the family.

5. **Interleave for B ≈ 0** (cure warp). See `references/b-coupling.md`. Reorder
   the plies so **each zone's through-thickness stack is a palindrome in both
   angle and material**. Verify: `compfea-build --plybook … --out /tmp/x` then
   read `/tmp/x/laminate_abd.csv`; `B/sqrt(A11*D11)` should be < ~0.08 in every
   load-bearing zone (the tip zone often keeps a small residual — low stress,
   acceptable). Reordering shifts flex ~0.3 N and does not move the kick.

6. **Matched-θ twist re-probe.** If the finalist's `theta_at_probe_deg` is more
   than ~5° off the baseline's, nudge `--uy-frac` (up = more tip rotation) and
   re-run `--twist-probe` until θ matches, then quote `k_twist` from *that*
   run. `k_twist` grows with θ, so an unmatched comparison flatters the softer
   blade.

7. **Confirm.** `--kick-bands 41`, check `flex_kick.json headline.warning` is
   null, `|migration_delta_mm|` is small (a stable kick), `twist.json warning`
   is null.

8. **Adversarial review + commit.** Fresh-context **opus** subagent on the diff
   vs the stated target — units, ply order / `-z` face, `B`, `D16` convention,
   convergence handling, whether the `k_twist` comparison is θ-matched. Fix real
   findings, then commit.

## Reporting

`scripts/report.py <results-dir> ...` builds the HTML candidate report:
the candidate table, each finalist's **deformed-blade side view**
(`kick_shape.svg`, produced by `run_tipweight.py --kick-bands N`), its ply
stack, and its cut sheet. `scripts/cutsheet.py <plybook.csv>` prints one cut
sheet (per-ply cut list + per-fabric area / roll yield / mass) on its own.

## Do not

- Sweep load and solve for deflection — displacement control only, always.
- Trust a `k_twist` ratio across an unmatched `theta_at_probe`.
- Work around a "too many cutbacks" non-convergence — a thin, soft open blade
  under a heavy stiff MID pad is a genuinely bad taper, not a solver problem.
- Add a `twill@0` ply and a `UD@0` ply to the same stack if you need B ≈ 0 —
  the 0° material mismatch is irreducible; pick one (all-UD, or single-fabric).
- Report absolute newtons as final — `FIN_TEST_3.step` is not dimensionally
  fixed; flex ~ thickness³ ~ length³. Trends and ratios are the deliverable.
