"""Classical lamination theory for a laminated strip, by hand.

The reference the linear gate is checked against. It imports nothing from
compfea -- layup.py writes cards and computes no stiffness at all, so this is a
genuinely independent number rather than the same arithmetic run twice.

Only 0 and 90 degree plies are supported, which is all this case uses; that
keeps Q16 = Q26 = 0 and the transformation to one swap.

Units: mm, N, MPa. D is N.mm^2 per mm of width, A is N/mm, B is N.

## Which reduction to a beam

There are two, and they differ by 0.25% on this strip -- enough to matter at the
tolerances here, and enough to look like a solver offset if you pick the wrong
one.

- **Wide plate / cylindrical bending**, ``EI = D11 * b``. Assumes the transverse
  curvature is suppressed. Right for a plate held flat across its width.
- **Narrow strip with free edges**, ``EI = b / d11_compliance = b (D11 -
  D12^2/D22)``. The long edges of this strip carry no moment, so the transverse
  curvature is free and the strip curls anticlastically. This is the correct one
  here, and ccx agrees with it to 0.01% on a converged mesh while sitting 0.24%
  below ``D11 * b`` -- see README.md. That 0.24% is the reduction, not a
  property of ccx's shell.
"""

from __future__ import annotations

import numpy as np

ANGLES = (0.0, 90.0)


def q_ply(angle_deg: float, *, e1: float, e2: float, nu12: float, g12: float):
    """Reduced stiffness (Q11, Q12, Q22, Q66) with fibres at 0 or 90 degrees."""
    if angle_deg not in ANGLES:
        raise ValueError(f"only {ANGLES} degree plies are supported, got {angle_deg}")
    nu21 = nu12 * e2 / e1
    denom = 1.0 - nu12 * nu21
    q11, q22 = e1 / denom, e2 / denom
    q12 = nu12 * e2 / denom
    if angle_deg == 90.0:
        q11, q22 = q22, q11  # a 90 degree ply is the same ply, swapped
    return q11, q12, q22, g12


def abd(
    angles_deg: tuple[float, ...] | list[float],
    ply_thickness: float,
    *,
    e1: float,
    e2: float,
    nu12: float,
    g12: float,
) -> np.ndarray:
    """The 6x6 [[A, B], [B, D]] laminate stiffness. Stack given bottom (-z) up."""
    total = ply_thickness * len(angles_deg)
    z = [-0.5 * total + k * ply_thickness for k in range(len(angles_deg) + 1)]
    a = np.zeros((3, 3))
    b = np.zeros((3, 3))
    d = np.zeros((3, 3))
    for k, angle in enumerate(angles_deg):
        q11, q12, q22, q66 = q_ply(angle, e1=e1, e2=e2, nu12=nu12, g12=g12)
        q = np.array([[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, q66]])
        a += q * (z[k + 1] - z[k])
        b += q * (z[k + 1] ** 2 - z[k] ** 2) / 2.0
        d += q * (z[k + 1] ** 3 - z[k] ** 3) / 3.0
    return np.block([[a, b], [b, d]])


def d_matrix(angles_deg, ply_thickness, **material) -> np.ndarray:
    """The 3x3 bending block."""
    return abd(angles_deg, ply_thickness, **material)[3:, 3:]


def d11(angles_deg, ply_thickness, **material) -> float:
    """D11 in N.mm^2 per mm of width.

    Exactly invariant under reversing the stack: ply k spans [z_{k-1}, z_k] and
    its mirror spans [-z_k, -z_{k-1}], and the two contribute the same
    (z^3 - z^3)/3. So no bending-stiffness check can ever see a reversal -- see
    `axial_coupling_ratio`, which can.
    """
    return float(d_matrix(angles_deg, ply_thickness, **material)[0, 0])


def ei_wide_plate(angles_deg, ply_thickness, *, width: float, **material) -> float:
    """D11 * b. Transverse curvature suppressed; not this strip. See module docs."""
    return d11(angles_deg, ply_thickness, **material) * width


def ei_narrow_strip(angles_deg, ply_thickness, *, width: float, **material) -> float:
    """b (D11 - D12^2/D22). Free long edges, transverse curvature free."""
    d = d_matrix(angles_deg, ply_thickness, **material)
    return width * float(d[0, 0] - d[0, 1] ** 2 / d[1, 1])


def axial_coupling_ratio(angles_deg, ply_thickness, **material) -> float:
    """b11/d11 from the inverted ABD, in mm. Zero for a symmetric stack.

    Under pure bending of a narrow strip every stress resultant vanishes except
    M1, so the inverted laminate compliance gives mid-plane axial strain
    eps1 = b11 M1 and curvature kappa1 = d11 M1, hence eps1 = (b11/d11) kappa1.
    Integrating along the beam, the mid-plane stretching adds

        du_y = (b11/d11) * theta_tip

    to the tip's axial position. B changes sign when the stack is reversed and D
    does not, so a stack and its reverse bend identically and draw in
    *differently*, by 2 (b11/d11) theta_tip. That is the only place ply order at
    z shows up in a displacement-controlled bend test, and it is what pins "the
    first ply line is the -z ply" through the solver rather than in a string.
    """
    compliance = np.linalg.inv(abd(angles_deg, ply_thickness, **material))
    return float(compliance[0, 3] / compliance[3, 3])


def tip_load_small_deflection(delta: float, *, ei: float, length: float) -> float:
    """P = 3 EI delta / L^3, the small-deflection cantilever result."""
    return 3.0 * ei * delta / length**3


def tip_angle_small_deflection(delta: float, *, length: float) -> float:
    """theta = 3 delta / (2 L) radians, for the same tip-loaded cantilever."""
    return 1.5 * delta / length
