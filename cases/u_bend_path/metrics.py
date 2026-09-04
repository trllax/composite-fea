"""Bend metrics: tip-normal spring force from energy.

Bench force is perpendicular to the tip and scaled by a moment arm:
    M = 2U / theta      (pure-moment / linear M-theta path)
    F = M / arm        (arm default = undeformed length L)

Do not use tip |RF| as F_90/F_180.
"""
from __future__ import annotations
import math

def moment_from_energy(U_Nmm: float, theta_rad: float) -> float:
    if theta_rad <= 0:
        raise ValueError("theta must be > 0")
    return 2.0 * U_Nmm / theta_rad

def tip_normal_force(M_Nmm: float, arm_mm: float) -> float:
    if arm_mm == 0:
        raise ValueError("arm must be nonzero")
    return M_Nmm / arm_mm

def f_spring(U_Nmm: float, theta_rad: float, arm_mm: float) -> tuple[float, float]:
    M = moment_from_energy(U_Nmm, theta_rad)
    return M, tip_normal_force(M, arm_mm)
