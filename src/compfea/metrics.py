"""Bend metrics for tip-normal spring force.

Bench force is perpendicular to the tip and scaled by a moment arm::

    M = 2U / theta
    F = M / arm

Do not use tip |RF| magnitude as F_90 / F_180.
"""

from __future__ import annotations


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
