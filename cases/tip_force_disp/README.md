# tip_force_disp

Compare drive methods on the `cantilever_89deg` strip (ccx 2.23).

## Perpendicular tip force → measure displacement

`*CLOAD` +Fz on `far_face`, tip free to draw in.

| Fz (N) | tip θ (deg) | tip uz (mm) | notes |
| --- | --- | --- | --- |
| 0.5–2 | 0.7–2.8 | matches FL³/3EI to ~0.1% | linear OK |
| 5 | 6.9 | 8.28 | mild NL |
| ≥10 | — | — | fails (too-slow Newton; increment → min) |

**Verdict:** fine for small-deflection stiffness checks. **Cannot** reach 90°/180°.
Vertical tip-force elastica also asymptotes near 90° tip slope, so load
control is the wrong primary path for the bend-test outputs.

## Tip Uz displacement → measure force (and θ)

Prescribe tip-edge Uz only; Uy free.

| Uz (mm) | δ/L | tip θ (deg) | \|F_tip\| (N) |
| --- | --- | --- | --- |
| 5 | 0.05 | 4.3 | 3.0 |
| 20 | 0.20 | 17.3 | 12.5 |
| 40 | 0.40 | 35.5 | 28.8 |
| 60 | 0.60 | 56.0 | 59.0 |
| ≥66 | ≥0.66 | — | fails before end of step |

**Verdict:** much better than force control; reaches ~56° here (your earlier
elastica notes got ~63°). Still short of 90°, far from 180°.

## Implications for F_90 / F_180

- Bench outputs need **kinematic control to a tip-tangent target**, then tip RF.
- Pure force + read δ inverts the test and stalls early.
- Pure tip Uz helps but plateaus before 90° on this mesh/path.
- Pure moment / tip UR works to 89° (see `cantilever_ansys`) but hits the 90°
  UR wall and does not give tip *force* as the natural metric.
- Full circular U-clamp: 90° OK, 180° not yet (`tip_clamp_u_drive`).

Next levers for 90°+: finer mesh, continuation in Uz with tighter STATIC
controls, or clamp kinematics that match the bench (not circular-arc assumption).
