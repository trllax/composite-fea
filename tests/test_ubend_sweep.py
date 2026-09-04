"""Unit tests for U-bend path + energy parse + sweep deck build (no ccx)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from compfea.layup import cross_ply_symmetric
from compfea.metrics import f_spring
from compfea.run import SolveError, parse_dat_energy
from compfea.ubend import (
    build_deck,
    force_at_theta,
    step_index_for,
    tip_displacements,
    tip_length_mm,
    theta_grid_deg,
)


SAMPLE_ENERGY_DAT = """
 total force (fx,fy,fz) for set FAR_FACE and time  0.1000000E+01

       -1.000000E+00  0.000000E+00 -2.000000E+00

 total internal energy for set BLADE and time  0.1000000E+01

        1.372843E+02

 total force (fx,fy,fz) for set FAR_FACE and time  0.2000000E+01

       -1.000000E+00  0.000000E+00 -3.000000E+00

 total internal energy for set BLADE and time  0.2000000E+01

        5.391728E+02
"""


def test_parse_dat_energy_reads_else_totals(tmp_path: Path):
    path = tmp_path / "job.dat"
    path.write_text(SAMPLE_ENERGY_DAT)
    frame = parse_dat_energy(path)
    assert list(frame.columns) == ["increment", "time", "elset", "energy"]
    assert len(frame) == 2
    assert frame.iloc[0]["elset"] == "blade"
    assert frame.iloc[0]["time"] == pytest.approx(1.0)
    assert frame.iloc[0]["energy"] == pytest.approx(137.2843)
    assert frame.iloc[1]["energy"] == pytest.approx(539.1728)


def test_parse_dat_energy_rejects_empty(tmp_path: Path):
    path = tmp_path / "job.dat"
    path.write_text("no energy here\n")
    with pytest.raises(SolveError, match="internal-energy"):
        parse_dat_energy(path)


def test_theta_grid_includes_90_and_180():
    grid = theta_grid_deg(step_deg=5.0, start_deg=5.0, end_deg=180.0)
    assert grid[0] == 5.0
    assert grid[-1] == 180.0
    assert 90.0 in grid
    assert len(grid) == 36


def test_tip_displacements_90_deg_on_fake_mesh():
    # Minimal stand-in: root at y=0, tip at y=100, width nodes at x=0 and 20.
    from compfea.geometry import Mesh

    mesh = Mesh(
        nodes={
            1: (0.0, 0.0, 0.0),
            2: (20.0, 0.0, 0.0),
            10: (0.0, 100.0, 0.0),
            11: (20.0, 100.0, 0.0),
        },
        elements={1: (1, 2, 11, 10, 1, 2, 11, 10)},
        nsets={"fixed_end": (1, 2), "far_face": (10, 11)},
        elsets={"blade": (1,)},
    )
    assert tip_length_mm(mesh) == pytest.approx(100.0)
    u = tip_displacements(mesh, math.radians(90.0))
    # Circular arc: R = L/θ, y_t = R sinθ, z_t = R (1-cosθ).
    # At 90°, y_t = z_t = 2L/π ≈ 63.662.
    y_t = 2.0 * 100.0 / math.pi
    assert u[10][0] == pytest.approx(0.0)
    assert u[10][1] == pytest.approx(y_t - 100.0)
    assert u[10][2] == pytest.approx(y_t)
    assert u[11][2] == pytest.approx(y_t)


def test_force_at_theta_matches_known_u_bend_reference(tmp_path: Path):
    path = tmp_path / "job.dat"
    path.write_text(SAMPLE_ENERGY_DAT)
    energy = parse_dat_energy(path)
    angles = [90.0, 180.0]
    u90, m90, f90 = force_at_theta(
        energy, theta_deg=90.0, step_index=step_index_for(angles, 90.0), arm_mm=100.0
    )
    u180, m180, f180 = force_at_theta(
        energy, theta_deg=180.0, step_index=step_index_for(angles, 180.0), arm_mm=100.0
    )
    assert f90 == pytest.approx(1.75, rel=2e-3)
    assert f180 == pytest.approx(3.43, rel=1e-2)
    assert f180 / f90 == pytest.approx(2.0, rel=2e-2)
    # Cross-check helper directly.
    assert f_spring(u90, math.pi / 2, 100.0)[1] == pytest.approx(f90)


def test_build_deck_emits_multi_step_tip_u_and_else():
    from compfea.geometry import Mesh

    mesh = Mesh(
        nodes={
            1: (0.0, 0.0, 0.0),
            2: (20.0, 0.0, 0.0),
            10: (0.0, 100.0, 0.0),
            11: (20.0, 100.0, 0.0),
        },
        elements={1: (1, 2, 11, 10, 1, 2, 11, 10)},
        nsets={"fixed_end": (1, 2), "far_face": (10, 11)},
        elsets={"blade": (1,)},
    )
    layup = cross_ply_symmetric(0.1, long_axis="y")
    deck = build_deck(mesh, layup, [5.0, 10.0])
    assert deck.count("*STEP") == 2
    assert "*BOUNDARY\nfixed_end, 1, 6" in deck or "fixed_end, 1, 6" in deck
    assert "*EL PRINT, ELSET=blade, TOTALS=ONLY" in deck
    assert "ELSE" in deck
    assert "10, 3, 3," in deck


def test_sweep_layup_ply_counts():
    from compfea.sweep import Design, layup_for

    d1 = Design(n_pairs=1, ply_mm=0.1)
    d2 = Design(n_pairs=2, ply_mm=0.1)
    assert d1.n_plies == 4
    assert d2.n_plies == 8
    assert d1.stack_label == "[0/90]s"
    plies = layup_for(d2).zones[0].plies
    assert len(plies) == 8
    assert [p.angle_deg for p in plies] == [0, 90, 0, 90, 90, 0, 90, 0]
