# B — the extension–bending coupling matrix, and why an asymmetric laminate warps

## What B is

Classical laminate theory writes the plate response as

```
| N |   | A  B | | eps0 |          N  = in-plane force / width  (N/mm)
|   | = |      | |      |          M  = moment / width          (N)
| M |   | B  D | | kappa|          eps0 = mid-plane strain, kappa = curvature
```

- **A** (extensional) and **D** (bending) are always there.
- **B** couples them: a pure in-plane load `N` produces curvature, and a pure
  moment produces in-plane strain. `B_ij = (1/2) Σ_k Qbar_ij(k) (z_k² − z_{k−1}²)`
  — each ply weighted by the **signed** difference of the squares of its
  through-thickness bounds.

`B = 0` **iff the stack is a mirror image of itself about the mid-plane** — same
ply *stiffness* (material **and** angle) at `+z` and `−z` for every depth. Any
mismatch leaves a residual `B`.

## Why it matters for a fin

`B ≠ 0` means:

1. **The part cures warped.** Cure shrinkage is an in-plane `N`; through `B` it
   becomes curvature. The blade comes off the single-sided mould with a twist
   or a cup. You can force it flat in the fixture, but that locks in residual
   stress — **which this model captures nowhere** (no thermal step, no cure
   strain). So `B` is the *only* warning you get that the real part will not
   match the solved shape.
2. **Extension–bending coupling under load.** The axial drive on the buckling
   bench feeds bending directly. A little of this is useful (it seeds the
   buckle), a lot of it muddies the flex reading.

`D16` / `D26` (bend–twist coupling) is separate: it comes from an **unbalanced**
±θ distribution, not from asymmetry. A balanced ±45 (a 0/90 woven rotated 45°
counts as balanced on its own) keeps `D16 = 0`. `D16` is the one quantity that
catches a mirrored `*ORIENTATION` sign convention — see `CLAUDE.md`.

## Reading it

`compfea-build` writes `laminate_abd.csv` with `a11 … d66` per zone. The useful
dimensionless warp index is

```
B / sqrt(A11 * D11)
```

- **< 0.03** — negligible, cures flat.
- **~0.05–0.15** — noticeable but tolerable; the shipped prototype
  (`plybook_shop.csv`) sits here in its thick zones and is accepted.
- **> 0.2** — the part will visibly warp; fix the ply order.

Compute it per zone; the number inflates on thin stacks (small `A11·D11`), so a
0.19 on a 3-ply open blade is milder than 0.19 on the 8-ply root.

## The interleave rule

Given a fixed set of plies for a zone, `B` depends only on their **order**.
Minimise it by making the stack a **palindrome in angle and material**:

1. **±45 skins on both faces** (ply 1 and ply N). A ±45 woven ply is balanced
   within itself, so a mirrored pair of them contributes no `B`.
2. **0° plies in even count, single material per stack.** `[45/0/0/45]` with
   both `0` the same fabric → `B = 0`. `[45/0/0/45]` with one `0` a `twill@0`
   (0.238 mm) and the other a `UD@0` (0.277 mm) → **`B ≠ 0`, irreducible**: the
   two "0" positions differ in stiffness and thickness, no reordering fixes it.
   Choose all-UD, or single-fabric.
3. **Odd centre ply on the mid-plane.** One extra 0° ply → put it dead centre
   (`core_center` in a `TaperDesign`, or the middle row of a hand book).
4. **Zone steps interleave, they don't stack at the end.** A ply book that
   lists *all* FULL rows, then *all* MID, then *all* QUARTER clusters the
   inboard-zone plies on one face → `B`. Write the rows in true through-
   thickness order so a QUARTER∩MID element's stack comes out palindromic,
   e.g. root `[45/0/45/0/45/0/45]` for an all-UD family.

Worked example: `cases/fin_test_3/ax/ax2.csv` (first-pass order, root
`[45/0/0/0/0/45/45]`, `B/√AD ≈ 0.20`) → `ax2r.csv` (interleaved, root
`[45/0/0/45/0/0/45]`, `B ≈ 0`), same plies, flex 13.7 → 14.1 N, kick unchanged.
`ax3.csv` → `ax3r.csv` reaches `B = 0` in every load-bearing zone because it has
no `twill@0` ply; `ax5.csv` → `ax5r.csv` cannot (it keeps the `twill@0`
softening ply) and holds a ~0.19 residual in the open blade.

## The tip-zone residual

The last-sixth `TIP` zone stacks as a symmetric FULL core plus **one** ±45 pad
ply → 5 plies, odd, always a small `B` (~0.1–0.26 on that thin stack). It is a
low-stress region; either accept it, or drop the `TIP` pad and let the tip run
the open-blade stack.
