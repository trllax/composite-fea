"""post.py force curve from synthetic ELSE energy."""

from __future__ import annotations

import math

import pandas as pd

from compfea.post import force_curve, plot_force_curve, summary_at
from compfea.ubend import theta_grid_deg


def test_force_curve_spring_ratio(tmp_path):
    angles = theta_grid_deg(step_deg=45.0, start_deg=45.0, end_deg=180.0)
    arm = 100.0
    # U = 0.5 * k * θ^2 => M = kθ, F = kθ/arm  (linear in θ)
    k = 2.0
    rows = []
    for i, deg in enumerate(angles, start=1):
        th = math.radians(deg)
        u = 0.5 * k * th * th
        rows.append({"increment": i, "time": float(i), "elset": "blade", "energy": u})
    energy = pd.DataFrame(rows)
    curve = force_curve(energy, angles, arm)
    assert list(curve["theta_deg"]) == angles
    # F should ~double 90→180
    f90 = float(curve.loc[curve.theta_deg == 90.0, "F_N"].iloc[0])
    f180 = float(curve.loc[curve.theta_deg == 180.0, "F_N"].iloc[0])
    assert abs(f180 / f90 - 2.0) < 1e-9
    summ = summary_at(curve)
    assert abs(summ["F_180_over_F_90"] - 2.0) < 1e-9
    svg = plot_force_curve(curve, tmp_path / "f.svg")
    assert svg.is_file() and svg.stat().st_size > 100
