"""Bend metrics for tip-normal spring force.

Bench force is perpendicular to the tip and scaled by a moment arm::

    M = 2U / theta
    F = M / arm

Do not use tip |RF| magnitude as F_90 / F_180.

``M = 2U/theta`` is a *secant* identity. It is exact for a linear torsional
spring (``U = M*theta/2``) and only approximate otherwise, so it is an
assumption about the response, not a definition. The generalized moment that is
always right is Castigliano's ``M = dU/dtheta``. ``moment_tangent`` computes
that, and ``linearity_deviation`` reports how far apart the two are, so the
assumption is a measured number in every run rather than an unexamined one.

Keep reporting the secant value as F: it is what the existing results mean, and
a redefinition would silently move every number already on disk. The deviation
is the flag.
"""

from __future__ import annotations

from collections.abc import Sequence


def moment_from_energy(u_nmm: float, theta_rad: float) -> float:
    if theta_rad <= 0:
        raise ValueError("theta must be > 0")
    return 2.0 * u_nmm / theta_rad


def tip_normal_force(m_nmm: float, arm_mm: float) -> float:
    if arm_mm == 0:
        raise ValueError("arm must be nonzero")
    return m_nmm / arm_mm


def f_spring(
    u_nmm: float, theta_rad: float, arm_mm: float
) -> tuple[float, float]:
    m = moment_from_energy(u_nmm, theta_rad)
    return m, tip_normal_force(m, arm_mm)


def moment_tangent(
    u_nmm: Sequence[float], theta_rad: Sequence[float]
) -> list[float]:
    """``dU/dtheta`` along the path -- the generalized moment, per Castigliano.

    Central difference, which on a **uniformly spaced** theta is exact for the
    quadratic ``U(theta)`` of a linear spring and so contributes no error of
    its own to the linearity check. ``theta_grid_deg`` pins the last angle to
    ``end_deg``, so a span that is not a whole number of steps leaves one short
    final interval; the second-to-last point then carries a first-order error
    of its own, and that is exactly the point the reported deviation falls back
    to for the largest angle. The two endpoints come back ``nan`` on purpose: there is no central
    difference to take there, and substituting a one-sided slope would inject
    its own first-order error -- on a 1-degree grid starting at 1 degree that
    error is 33% at the first point, which would swamp the real nonlinearity
    this function exists to expose. Callers must skip the ends rather than
    treat a worse estimator as the same quantity.
    """
    u = [float(v) for v in u_nmm]
    th = [float(v) for v in theta_rad]
    if len(u) != len(th):
        raise ValueError(f"u has {len(u)} points, theta has {len(th)}")
    if len(u) < 3:
        raise ValueError(
            f"moment_tangent needs at least three points for a central "
            f"difference, got {len(u)}"
        )
    out = [float("nan")] * len(u)
    for i in range(1, len(u) - 1):
        dth = th[i + 1] - th[i - 1]
        if dth == 0.0:
            raise ValueError(f"repeated theta around index {i}: {th[i]}")
        out[i] = (u[i + 1] - u[i - 1]) / dth
    return out


def linearity_deviation(
    m_secant: Sequence[float], m_tangent: Sequence[float]
) -> list[float]:
    """``(M_secant - M_tangent) / M_tangent``, the error in the 2U/theta model.

    Near zero means the response really is the linear spring that ``M = 2U/theta``
    assumes. A large value means the reported F is a secant average over a
    curved M(theta) and should not be read as the force at that angle.
    """
    sec = [float(v) for v in m_secant]
    tan = [float(v) for v in m_tangent]
    if len(sec) != len(tan):
        raise ValueError(f"{len(sec)} secant values, {len(tan)} tangent values")
    return [
        float("nan") if t == 0.0 else (s - t) / t for s, t in zip(sec, tan, strict=True)
    ]
