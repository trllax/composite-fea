"""Gate for the flex + kick-point metric (``compfea.kick``).

The kick point is where a bent blade curves most: the peak of a smoothed
``|kappa(s)|`` on the deformed centreline, bucketed heel / mid / tip. The
reference is ``cases/fin_test_3/compressed_elastica.py`` -- a clamped-free
buckling column's curvature is largest at the clamp, so a *uniform* blade kicks
at the heel; a stiffness taper that thins the tip moves the peak outboard.

Most of this file is pure unit work on synthetic centrelines sampled from that
closed form. The last test drives a real flat strip with ccx and has, like
``tests/test_compressed_elastica.py``, **no skipif on a missing solver**: a gate
that quietly skips is not a gate.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from compfea import kick
from compfea.deck import StaticStep, assemble, axial_drive_body
from compfea.geometry import CircularCamber, Outline, mesh_outline
from compfea.layup import PLACEHOLDER_CFRP, Layup, Ply
from compfea.run import solve

ROOT = Path(__file__).resolve().parents[1]


def _module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


clpt = _module(ROOT / "cases/smoke_cantilever/clpt.py", "kick_gate_clpt")
ce = _module(ROOT / "cases/fin_test_3/compressed_elastica.py", "kick_gate_ce")

LENGTH_MM = 120.0
EI_NMM2 = 1.0e6


# --------------------------------------------------------------------------
# Synthetic centreline builders
# --------------------------------------------------------------------------


def _integrate_shape(s0: np.ndarray, kappa: np.ndarray) -> pd.DataFrame:
    """A flat-along-y blade bent by ``kappa(s0)`` -> a kick centreline frame.

    Undeformed centreline is the straight line ``y = s0`` (x = z = 0). The
    deformed centreline integrates the prescribed curvature: ``phi = int kappa``,
    then ``y = int cos phi``, ``z = int sin phi``. Columns are exactly what
    ``kick.tangent_and_curvature`` consumes.
    """
    def _cumtrap(values: np.ndarray) -> np.ndarray:
        step = 0.5 * (values[:-1] + values[1:]) * np.diff(s0)
        return np.concatenate([[0.0], np.cumsum(step)])

    phi = _cumtrap(kappa)
    return pd.DataFrame(
        {
            "s0_mm": s0,
            "x0_mm": 0.0,
            "y0_mm": s0,
            "z0_mm": 0.0,
            "x_mm": 0.0,
            "y_mm": _cumtrap(np.cos(phi)),
            "z_mm": _cumtrap(np.sin(phi)),
        }
    )


def _uniform_column(phi_l_deg: float, *, n_bands: int) -> pd.DataFrame:
    """Centreline of a uniform clamped-free buckling column at tip angle phi_l."""
    phi_l = math.radians(phi_l_deg)
    s_fine, kap_fine = ce.curvature_profile(
        phi_l, ei=EI_NMM2, length=LENGTH_MM, n=8001
    )
    fine = _integrate_shape(s_fine, kap_fine)
    s0 = np.linspace(6.0, LENGTH_MM - 6.0, n_bands)
    return pd.DataFrame(
        {
            "s0_mm": s0,
            "x0_mm": 0.0,
            "y0_mm": s0,
            "z0_mm": 0.0,
            "x_mm": 0.0,
            "y_mm": np.interp(s0, s_fine, fine["y_mm"]),
            "z_mm": np.interp(s0, s_fine, fine["z_mm"]),
        }
    )


def _gaussian_hinge(
    center_frac: float, *, n_bands: int = 41, sigma_mm: float = 10.0
) -> pd.DataFrame:
    s0 = np.linspace(3.0, LENGTH_MM - 3.0, n_bands)
    kappa = 0.03 * np.exp(-(((s0 - center_frac * LENGTH_MM) / sigma_mm) ** 2))
    return _integrate_shape(s0, kappa)


# --------------------------------------------------------------------------
# The closed form the gate leans on
# --------------------------------------------------------------------------


@pytest.mark.parametrize("deg", [30.0, 60.0, 90.0])
def test_curvature_profile_is_consistent_with_its_own_scalars(deg):
    phi_l = math.radians(deg)
    s, kappa = ce.curvature_profile(phi_l, ei=EI_NMM2, length=LENGTH_MM, n=6001)
    assert s[0] == pytest.approx(0.0, abs=1e-6)
    assert s[-1] == pytest.approx(LENGTH_MM, rel=1e-4)
    # Curvature is monotone decreasing clamp -> tip, maximal at the clamp.
    assert np.all(np.diff(kappa) <= 1e-9)
    assert kappa[0] == pytest.approx(
        ce.curvature_root(phi_l, ei=EI_NMM2, length=LENGTH_MM), rel=1e-3
    )
    # int kappa ds over the span is the tip rotation.
    assert np.trapezoid(kappa, s) == pytest.approx(phi_l, rel=1e-4)


# --------------------------------------------------------------------------
# tangent_and_curvature
# --------------------------------------------------------------------------


@pytest.mark.parametrize("deg", [20.0, 45.0, 90.0])
def test_tangent_angle_recovers_the_tip_angle(deg):
    cur = kick.tangent_and_curvature(_uniform_column(deg, n_bands=25), long_axis="y")
    tip_deg = math.degrees(cur["phi_rad"].to_numpy()[-1])
    assert tip_deg == pytest.approx(deg, rel=0.03)


def test_tangent_angle_is_not_a_running_sum_of_slopes():
    """phi is the local material rotation, not a cumsum of absolute slopes.

    A cumsum bug inflates the tip angle several-fold; guard the regression.
    """
    cur = kick.tangent_and_curvature(_uniform_column(90.0, n_bands=25), long_axis="y")
    assert math.degrees(cur["phi_rad"].abs().max()) < 120.0


# --------------------------------------------------------------------------
# kick_point
# --------------------------------------------------------------------------


@pytest.mark.parametrize("deg", [20.0, 45.0, 90.0])
@pytest.mark.parametrize("n_bands", [15, 21, 31])
def test_uniform_column_kicks_at_the_heel(deg, n_bands):
    cl = _uniform_column(deg, n_bands=n_bands)
    cur = kick.tangent_and_curvature(cl, long_axis="y")
    kp = kick.kick_point(cur, free_len_mm=LENGTH_MM)
    assert kp.bucket == "heel"
    assert kp.s_kick_frac < 0.25
    assert kp.warning is None
    # rotation median of a uniform column is ~0.32 L (closed form), and well
    # outboard of the curvature peak -- the two are different statistics.
    expected = ce.rotation_median_s(math.radians(deg), length=LENGTH_MM)
    assert kp.s_rotmed_mm == pytest.approx(expected, rel=0.06)
    assert kp.s_rotmed_mm > kp.s_kick_mm


def test_mid_span_hinge_reads_as_mid():
    kp = kick.kick_point(
        kick.tangent_and_curvature(_gaussian_hinge(0.5), long_axis="y"),
        free_len_mm=LENGTH_MM,
    )
    assert kp.bucket == "mid"
    assert kp.s_kick_frac == pytest.approx(0.5, abs=0.08)
    assert kp.warning is None


def test_outboard_hinge_reads_as_tip():
    kp = kick.kick_point(
        kick.tangent_and_curvature(_gaussian_hinge(0.78), long_axis="y"),
        free_len_mm=LENGTH_MM,
    )
    assert kp.bucket == "tip"
    assert 0.65 < kp.s_kick_frac < 0.9
    assert kp.warning is None


def test_two_far_apart_maxima_warn():
    s0 = np.linspace(3.0, LENGTH_MM - 3.0, 41)
    kappa = 0.02 * np.exp(-(((s0 - 0.25 * LENGTH_MM) / 8.0) ** 2)) + 0.02 * np.exp(
        -(((s0 - 0.75 * LENGTH_MM) / 8.0) ** 2)
    )
    kp = kick.kick_point(
        kick.tangent_and_curvature(_integrate_shape(s0, kappa), long_axis="y"),
        free_len_mm=LENGTH_MM,
    )
    assert kp.warning is not None


def test_broad_plateau_warns():
    s0 = np.linspace(3.0, LENGTH_MM - 3.0, 61)
    kappa = np.where((s0 > 0.2 * LENGTH_MM) & (s0 < 0.8 * LENGTH_MM), 0.015, 0.001)
    kp = kick.kick_point(
        kick.tangent_and_curvature(_integrate_shape(s0, kappa), long_axis="y"),
        free_len_mm=LENGTH_MM,
    )
    assert kp.warning is not None
    assert kp.plateau_frac > 0.4


# --------------------------------------------------------------------------
# span_band_nsets
# --------------------------------------------------------------------------


def test_span_band_nsets_are_ordered_and_disjoint_along_the_span():
    mesh = mesh_outline(
        Outline.rectangle(chord=20.0, span=LENGTH_MM), n_chord=2, n_span=40
    )
    clamped = set(mesh.nsets["fixed_end"])
    tip = set(mesh.nsets["far_face"])
    bands = kick.span_band_nsets(
        mesh, long_axis="y", n_bands=15, chord_frac=1.0, exclude=clamped | tip
    )
    assert list(bands) == sorted(bands)
    s0_clamp, sign = kick.clamp_face(mesh, "y")
    centres = [
        sign * (np.mean([mesh.nodes[n][1] for n in ids]) - s0_clamp)
        for ids in bands.values()
    ]
    assert centres == sorted(centres)
    assert centres[0] > 0.0 and centres[-1] < LENGTH_MM
    # No node in two bands, and none clamped or on the tip edge.
    seen: set[int] = set()
    for ids in bands.values():
        assert seen.isdisjoint(ids), "a node is shared between bands"
        seen.update(ids)
    assert seen.isdisjoint(clamped) and seen.isdisjoint(tip)


def test_kick_is_independent_of_which_way_the_blade_buckled():
    """A centreline mirrored in z gives the same kick point (sign-normalised)."""
    cl = _uniform_column(90.0, n_bands=21)
    up = kick.kick_point(
        kick.tangent_and_curvature(cl, long_axis="y"), free_len_mm=LENGTH_MM
    )
    cl_down = cl.copy()
    cl_down["z_mm"] = -cl_down["z_mm"]
    cl_down["z0_mm"] = -cl_down["z0_mm"]
    down = kick.kick_point(
        kick.tangent_and_curvature(cl_down, long_axis="y"), free_len_mm=LENGTH_MM
    )
    assert down.bucket == up.bucket
    assert down.s_kick_mm == pytest.approx(up.s_kick_mm, abs=1e-6)
    assert down.theta_deg == pytest.approx(up.theta_deg, abs=1e-6)


# --------------------------------------------------------------------------
# End to end: a real flat strip, solved.
# --------------------------------------------------------------------------

_STRIP_STACK = (0.0, 90.0, 90.0, 0.0)
_PLY_MM = 0.25
_STRIP_WIDTH = 20.0


@pytest.fixture(scope="module")
def strip_run(tmp_path_factory):
    """A flat [0/90]s strip driven as a buckling column past a 90 deg tip tangent.

    Mirrors ``tests/test_compressed_elastica.py`` (same bow, same drive) but adds
    a row of ``span_band_nsets`` U prints so ``flex_kick`` has a centreline.
    """
    run_dir = tmp_path_factory.mktemp("kick_strip")
    ei = clpt.ei_narrow_strip(
        _STRIP_STACK,
        _PLY_MM,
        width=_STRIP_WIDTH,
        e1=PLACEHOLDER_CFRP.e1,
        e2=PLACEHOLDER_CFRP.e2,
        nu12=PLACEHOLDER_CFRP.nu12,
        g12=PLACEHOLDER_CFRP.g12,
    )
    # Drive to a tip tangent past 90 deg so flex_kick can interpolate W(90).
    delta = ce.shortening_ratio(math.radians(100.0)) * LENGTH_MM
    mesh = mesh_outline(
        Outline.rectangle(chord=_STRIP_WIDTH, span=LENGTH_MM),
        n_chord=2,
        n_span=60,
        camber=CircularCamber(radius=1000.0 * LENGTH_MM),
    )
    bands = kick.span_band_nsets(
        mesh,
        long_axis="y",
        n_bands=19,
        chord_frac=1.0,
        exclude=set(mesh.nsets["fixed_end"]) | set(mesh.nsets["far_face"]),
    )
    import dataclasses

    mesh = dataclasses.replace(mesh, nsets={**mesh.nsets, **bands})
    layup = Layup.uniform(
        [Ply(_PLY_MM, a) for a in _STRIP_STACK],
        long_axis="y",
        elset="blade",
        material=PLACEHOLDER_CFRP,
    )
    body = axial_drive_body(
        "far_face", 2, -abs(delta), tangent_nsets=tuple(bands)
    )
    deck = assemble(
        mesh_inp=mesh.to_inp(),
        layup=layup,
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=[StaticStep(body, inc=6000, static_line="0.01, 1.0, 1.E-9, 0.05")],
    )
    (run_dir / "deck.inp").write_text(deck)
    result = solve(run_dir / "deck.inp", run_dir, final_time=1.0)
    return mesh, run_dir, ei, result


def test_flat_strip_kicks_at_the_heel_end_to_end(strip_run):
    from compfea.run import parse_dat_totals

    mesh, run_dir, ei, _ = strip_run
    dat = run_dir / "job.dat"

    # Locate the load states the way run_tipweight does: a theta(time) curve and
    # a W(time) curve, then interpolate the target angles.
    theta_t = kick.tip_tangent_series(dat, mesh, long_axis="y")
    rf = parse_dat_totals(dat)
    w_by_time = (
        rf[rf["nset"] == "fixed_end"].set_index("time")["fy"].abs()
    )
    theta_t = theta_t.assign(W_N=theta_t["time"].map(w_by_time)).dropna()
    assert theta_t["theta_deg"].max() > 90.0, "strip did not reach a 90 deg tip tangent"
    t_h = float(np.interp(90.0, theta_t["theta_deg"], theta_t["time"]))
    t_m = float(np.interp(25.0, theta_t["theta_deg"], theta_t["time"]))
    w90 = float(np.interp(90.0, theta_t["theta_deg"], theta_t["W_N"]))

    fk = kick.flex_kick(
        dat,
        mesh,
        time_headline=t_h,
        time_migration=t_m,
        flex_n=w90,
        long_axis="y",
        out_dir=run_dir,
    )
    assert fk.headline.bucket == "heel"
    assert fk.headline.s_kick_frac < 0.3
    assert fk.headline.warning is None
    # Both load states resolve to a heel kick for a uniform strip; migration is
    # small (the kick point of a uniform column does not walk up the span).
    assert fk.migration.bucket == "heel"
    assert abs(fk.migration_delta_mm) < 0.2 * LENGTH_MM
    # flex is W at a 90 deg tip tangent; the closed form for this EI and length.
    expected = ce.axial_load_at_angle(math.radians(90.0), ei=ei, length=LENGTH_MM)
    assert fk.flex_n == pytest.approx(expected, rel=0.05)
    for name in ("kick_kappa.csv", "kick_kappa.svg", "kick_shape.svg"):
        assert (run_dir / name).is_file(), f"{name} not written"


def test_flat_strip_gate_is_fast(strip_run):
    _, _, _, result = strip_run
    assert result.wall_time_s < 30.0, (
        f"kick strip solve took {result.wall_time_s:.1f} s"
    )
