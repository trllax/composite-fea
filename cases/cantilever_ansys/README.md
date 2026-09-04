# cantilever_ansys

Flat laminated strip in large deflection, cross-checked against ANSYS for UD
and woven CF. This is a code-to-code validation case: a third reference class
alongside `smoke_cantilever` (closed-form) and `fin_20n` (hardware).

Status: **awaiting ANSYS values and material data.** The model definition below
is fixed; the numbers are not yet filled in.

## Why this shape

A strip is the only cantilever where we have an independent answer that owes
nothing to either code. ccx already reproduces the closed-form elastica to
0.15-0.22% from 6 to 63 degrees of tip rotation on a 100 x 20 x 1 mm [0/90]s
strip, so a disagreement with ANSYS points at model setup, not at the solver.

## Load case

Prescribed rotation of the far edge, translations free. This is a **pure end
moment**: curvature is constant, so `M = EI * theta / L` is exact at any
rotation, with no elastica integral. That gives the moment reaction probe a
closed-form target.

| generic name | this part |
| --- | --- |
| driven node set | `far_face`, all nodes at y = L |
| driven DOF | 4 (rotation about x), to **1.5533430343 rad = 89 deg** |
| free at the far edge | all translations |
| reaction node set | `fixed_end`, all nodes at y = 0, all 6 DOF fixed |
| compared quantity | bending moment, sampled **mid-span** |

Axes: x short (width), y long (beam axis), z normal, plies stack in z. Note
this transposes CLAUDE.md's "+x along the long axis" convention, so 0 deg here
means fibre along **y**. Confirm ANSYS uses the same reference direction or the
two models are 90 deg apart.

### Do not prescribe 90 degrees

ccx has a hard wall at exactly 90.000 deg of prescribed shell rotation.
Verified: 91 deg dies at t=0.989, 100 deg at t=0.900, 179 deg at t=0.503 --
all of them at 90.00 deg of accumulated rotation. Independent of mesh size,
increment size and solver. 89.999 deg converges every time. M is exactly linear
in theta for a pure end moment, so drive 89 deg and scale.

Do **not** add an `x` translation constraint on the far edge. It makes ccx die
with exit 255 and no diagnostic, and it buys nothing: with x free the
anticlastic displacement is +/-2.9e-4 mm.

### An earlier load case, rejected

A tip-displacement elastica case was tried first and works (ccx matches the
closed-form elastica to 0.15-0.22% from 6 to 63 deg). It was dropped because a
pure vertical end load cannot reach 90 deg at all -- it asymptotes, needing
delta/L ~ 0.92 for 89.9 deg. These were also tried and rejected:

- Driving DOF 6 is the drilling rotation for a strip in the x-y plane. It
  converges and returns a root reaction of exactly 0.000000E+00.
- Driving the true bending rotation (DOF 5) at a single node diverges at
  t = 0.027.
- A prescribed rotation with the tip free to translate transmits a pure moment,
  and ccx 2.23 has no `RM` label for `*NODE PRINT` -- there is no moment to read
  back. A pure vertical end load also cannot reach 90 degrees of tip rotation at
  all: it asymptotes, needing delta/L ~ 0.92 for 89.9 degrees.

So: drive the whole tip edge transversely with the axial direction free, which
is the classic elastica problem and gives an unambiguous force to compare.

| generic name | this part |
| --- | --- |
| driven node set | `tip_edge`, all nodes at x = L |
| driven DOF | 3 (transverse), to delta/L = 0.66 |
| free at the tip | DOF 1, so the strip draws in as it bends |
| reaction node set | `fixed_end`, all nodes at x = 0, all 6 DOF fixed |

Compare the **force-displacement curve**, not one endpoint. Both codes produce
the full history for free and a curve catches a stiffness error that a single
point can hide.

## Reading the moment out of ccx

`*SECTION PRINT` does not work with `*SHELL SECTION, COMPOSITE`. Verified on
two decks identical but for the section card, 1 N tip load: a plain shell
section gives 99.99933 N.mm against an exact 100, a composite section gives
0.000000. ccx finds the surface and computes its area, centre of gravity and
normal correctly in both cases, but never populates the stress array it
integrates for composite sections. No surface definition works around it.
`*NODE PRINT` has no `RM` label, and `*RIGID BODY` is rejected on shell nodes.

So the moment comes from ply stresses, cross-checked against energy:

```
*EL PRINT, ELSET=blade
S
*EL PRINT, ELSET=blade, TOTALS=ONLY
ELSE
```

Two things to get right:

- ccx expands each S8R into one C3D20R **per ply**, so a 4-ply shell prints 32
  integration points per element. Thickness is the slowest index. The **first
  ply listed in the section card is the -z ply**.
- `*EL PRINT S` writes stress in the `*ORIENTATION` frame **as defined at t=0**;
  it does not co-rotate. At 85 deg of local rotation the axial stress has moved
  from `sxx` into `szz`. Use the rotation-invariant `tr(sigma) - sigma_ww`.

Sample **mid-span**, not the root: the root element reads +1.26% high and does
not improve with refinement -- a real clamped-edge boundary layer where fixing
the rotational DOFs suppresses anticlastic curvature.

Cross-check every moment against `2U/theta` from the energy output. The two
routes agree to 0.006% when ply ordering and z-stations are right, and diverge
loudly when they are not.

## Mesh

Converged by 8 elements along the beam: 4->8 moves the answer 0.2%, 8->64 moves
it 0.05%. Keep more than one element across the width if you refine further.

| n | wall | M_mid (N.mm) | err vs exact |
| --- | --- | --- | --- |
| 4 | 0.41 s | 3092.5 | -0.43% |
| 8 | 0.85 s | 3099.1 | -0.22% |
| 16 | 2.72 s | 3098.5 | -0.24% |
| 64 | 9.86 s | 3097.7 | -0.27% |

Reference values above use the placeholder CFRP, not the real material:
CLPT gives EI = 1.9995e5 N.mm^2 and M_exact(89 deg) = 3105.9 N.mm.

## Geometry

Proposed, and yours to change -- the only hard requirement is that both codes
use identical numbers. Keep L/t >= 50 so shell and solid formulations agree.

| quantity | value |
| --- | --- |
| length L | 100 mm |
| width w | 20 mm |
| ply thickness | 0.25 mm |
| stack | TODO: ply count and sequence for each material |

## Materials

Units are mm, N, tonne, MPa, s. Moduli in MPa, density in tonne/mm^3.
Woven is modeled as one homogenized balanced orthotropic ply per fabric layer.

| constant | UD CF | woven CF |
| --- | --- | --- |
| E1 | TODO | TODO |
| E2 | TODO | TODO |
| E3 | TODO | TODO |
| nu12 | TODO | TODO |
| nu13 | TODO | TODO |
| nu23 | TODO | TODO |
| G12 | TODO | TODO |
| G13 | TODO | TODO |
| G23 | TODO | TODO |
| density | TODO | TODO |

These must be the exact constants used in the ANSYS run. A material card that
differs in the fourth digit makes the comparison meaningless.

## What to export from ANSYS

1. Reaction force at the fixed end vs prescribed tip displacement -- the whole
   history, not the endpoint.
2. Tip axial draw-in at the final substep, which pins the deformed shape rather
   than just the stiffness.
3. Element type (SHELL181 / SHELL281 / SOLSH190), and whether large deflection
   was on.
4. Mesh density and how the root was clamped (which DOFs, on nodes or an edge).
5. Whether the laminate was defined by ACP or by a section stack.

## Expected disagreement

ccx authors shells as S8R but expands them into stacked solid elements with
knot MPCs before solving; ANSYS SHELL181/281 are true shear-deformable shells.
On a slender strip these agree closely. Near the clamped edge, or on a thick
laminate, they will not agree exactly, and that difference is real rather than
a defect in either code. Set the tolerance after seeing the first comparison --
do not pick a number now and then rationalise toward it.

## Pinned value

TODO: the ANSYS reaction at delta/L = 0.66 for each material, the tolerance,
the ANSYS version, and who ran it.
