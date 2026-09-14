# cases/fin_test_3 — the tip-weight buckling-cantilever bench

This is the method for characterising a fin blade. It models the actual shop
bench: the blade is clamped at the heel and held so its axis lies along a hung
weight's line of action; the weight compresses it **past its Euler load** and it
lies over until the tip tangent has turned some angle from the load line.

- Drive a small tip **clip patch** (a few mid-chord nodes, where a string ties)
  in **axial translation** under `*BOUNDARY`, toward the clamp. Leave every other
  DOF on that patch free — free transverse translation is what makes this a
  clamped-*free* column and not a clamped-pinned one (~8x stiffer).
- Read `W`, the equivalent hung weight, as the **axial reaction** — the clip's
  own `RF` along the drive axis (the fixture force, no body load), cross-checked
  against the heel `RF`. No moment read-back, no `M = 2U/theta`, no moment arm.
- `theta`, the tip tangent, is an **output**, read from two or more `*NODE PRINT`
  U stations in the `.dat` (deck nodes — never the `.frd`, whose expanded solid
  mesh fans one planform station into many and turns a slope fit into a chord
  across curvature). So ccx's 90-degree prescribed-rotation wall never comes up.
- The buckling bifurcation is broken by an imperfection: an unsymmetric layup
  couples the axial drive into bending on its own; near Euler, add a
  sub-millimetre first-mode bow (`run_tipweight.py --bow-tip-mm`).
- **Flex** = `W` at `theta = 90 deg`, in newtons (~5 N supersoft, ~15–20 N
  hard). **Kick point** = where the blade curves most: `compfea.kick` builds the
  deformed centreline from a row of mid-chord U stations (`--kick-bands`),
  takes the peak of a smoothed `|kappa(s)|`, and buckets it heel / mid / tip. It
  warns when the peak is a broad plateau or bimodal — then "kick point" is
  ill-posed for that blade. `flex_kick` also reports how the kick point migrates
  up the span between a light load and 90 deg.

```
*STEP, NLGEOM, INC=40000
*STATIC
0.002, 1.0, 1.E-10, 0.05
*BOUNDARY
clip, 2, 2, -180.0          ** axial drive toward the clamp; dof = long axis
*NODE PRINT, NSET=fixed_end, TOTALS=YES
RF
*NODE PRINT, NSET=clip, TOTALS=YES
RF
*NODE PRINT, NSET=clip
U
*NODE PRINT, NSET=tip_band
U
*EL PRINT, ELSET=blade, TOTALS=ONLY
ELSE
*END STEP
```

`deck.axial_drive_body` writes this (plus an optional `*DLOAD, GRAV` settle
step for self-weight). The analytic gate is
`cases/fin_test_3/compressed_elastica.py` — the clamped-free buckling elastica
in elliptic-integral form, checked on a flat prismatic strip by
`tests/test_compressed_elastica.py` and `tests/test_kick.py`. Per-part node set
names and the drive DOF live in the case directory, not here.

- **Tip twist stiffness** (`run_tipweight.py --twist-probe`, `compfea.twist`,
  `deck.twist_couple_body`) is a second/third `*STATIC` step on the buckled
  state: a `+F`/`-F` **force couple** (`*CLOAD`, second step `OP=NEW`) on the
  chord-split halves of the tip edge. A *force*, not a prescribed rotation: an
  unsymmetric layup has rolled the tip edge by the time it is 90 deg over, so its
  nodes sit at unknown axial positions and a prescribed total `*BOUNDARY`
  displacement there cannot be formed without first solving. `k_twist = couple /
  twist_angle`, all on the **deformed** chord width `s`: `couple = |F|·n·s`,
  `twist_angle` the half-swing of the `twist_le`/`twist_te` `uy` difference over
  `s`. Read from `*NODE PRINT` `U` — never `RM`. A moment about **global z**,
  which is the twist axis only near `theta = 90 deg` — so `k_twist` is
  **comparative only** *and only across layups reaching a similar
  `theta_at_probe`* (recorded). `--twist-probe` does **not** change `--uy-frac`
  (the sweep needs 0.60 to bracket 90 deg or `flex_n` is null); the probe reads
  wherever the drive ends. `k_twist_nmm_per_rad` is `None` when the swing is
  below ~0.03 deg; `twist_stiffness` warns on that, on θ outside 45..135, on a
  lopsided le/te swing, and on a large heel RF swing. `rank_layups
  --twist-tiebreak` re-orders same-`dist`-bin rows by `k_twist` but skips warned
  / `None` rows. Gated by `tests/test_twist.py` (do not weaken). `collect()`
  bounds `W_vs_theta.csv` and the kick map to the drive step, so flex/kick are
  untouched when the flag is off.

## Legacy / superseded

Kept for the strip cases and for reading old results, **not** the path for new
fin work:

- `ubend.py` — the circular-arc tip-U clamp path. It prescribes the tip onto an
  undeformed-length arc and reads `F = M / arm` with `M = 2U/theta`. It carries
  **no compressive geometric-stiffness term**, so on a fin it read 2–3x too
  soft; that was a boundary-condition error, not material or mesh. The
  `M = 2U/theta` identity is also only a secant spring.
- `metrics.py` — the secant/tangent moment and linearity-deviation machinery for
  that `F(theta)`.
- `post.py`, `shapes.py` — the `F(theta)` CSV/SVG and the circular-arc
  target-arc poses that go with the U-bend path.
- `cases/fin_20n` — its binding drives DOF 6 (rotation) to 1.5708 rad at
  `tip_ref`. The case stays (it is the only hardware-pinned one) but new fin
  characterisation uses the buckling bench above.
