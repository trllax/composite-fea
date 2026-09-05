"""post.py force curve from synthetic ELSE energy."""

from __future__ import annotations

import math

import pandas as pd

from compfea.post import (
    force_curve,
    linearity_dev_near,
    max_linearity_dev,
    plot_force_curve,
    summary_at,
)
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


def _curve(u_of_theta, *, step_deg=1.0, end_deg=90.0, arm=100.0):
    angles = theta_grid_deg(step_deg=step_deg, start_deg=1.0, end_deg=end_deg)
    rows = [
        {
            "increment": i,
            "time": float(i),
            "elset": "blade",
            "energy": u_of_theta(math.radians(deg)),
        }
        for i, deg in enumerate(angles, start=1)
    ]
    return force_curve(pd.DataFrame(rows), angles, arm)


def test_a_linear_spring_reports_essentially_zero_deviation():
    """M = 2U/theta is exact here, so the check must not invent a deviation."""
    curve = _curve(lambda th: 0.5 * 2.0 * th * th)
    assert max_linearity_dev(curve) < 1e-12


def test_a_stiffening_response_is_flagged():
    """A curved M(theta) must show up, or the check is decorative."""
    curve = _curve(lambda th: 0.5 * 2.0 * th * th * (1.0 + 0.3 * th))
    assert max_linearity_dev(curve) > 0.05


def test_endpoints_are_nan_not_a_one_sided_slope():
    """A one-sided slope at the ends is a different, worse estimator.

    Mixing it in reads as 33% nonlinearity on an exactly linear spring, which
    would bury the real signal. The ends must abstain instead.
    """
    curve = _curve(lambda th: 0.5 * 2.0 * th * th)
    assert math.isnan(curve["linearity_dev"].iloc[0])
    assert math.isnan(curve["linearity_dev"].iloc[-1])
    assert curve["linearity_dev"].iloc[1:-1].notna().all()


def test_the_reported_angle_falls_back_to_the_nearest_interior_point():
    """F_180 is always the last row, so it never has a central difference."""
    curve = _curve(lambda th: 0.5 * 2.0 * th * th, step_deg=1.0, end_deg=180.0)
    near = linearity_dev_near(curve)
    assert "F_90" in near and "F_180" in near
    assert not math.isnan(near["F_180"])
    assert abs(near["F_180"]) < 1e-12


def test_an_unsolved_angle_is_not_reported():
    curve = _curve(lambda th: 0.5 * 2.0 * th * th, end_deg=90.0)
    near = linearity_dev_near(curve)
    assert "F_90" in near
    assert "F_180" not in near
