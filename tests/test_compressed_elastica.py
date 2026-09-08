"""Gate for the axial (buckling) tip-weight load case, solved with ccx.

``cases/fin_test_3/run_tipweight.py`` models a blade held along a hung weight's
line of action: the weight compresses it past its Euler load and it lies over,
tip tangent rotating toward horizontal. This checks that ccx reproduces the
closed-form buckling elastica on a flat prismatic strip, so the only unknowns
left in the fin run are geometry and material -- not the boundary condition.

Like ``tests/test_smoke_cantilever.py`` there is no skipif on a missing ccx: a
gate that quietly skips is not a gate. The reference is
``cases/fin_test_3/compressed_elastica.py`` (elliptic-integral form, shares no
code with ``compfea``); ``EI`` is the free-edge narrow-strip reduction from
``cases/smoke_cantilever/clpt.py``, the same one that case pins.

The bifurcation is broken with a small developable ``CircularCamber`` bow
(~L/2000 at the tip). That is a real imperfect column and lands ~0.16% below the
perfect-column closed form, independent of mesh from 60 to 120 elements along
the span and scaling with the bow -- see the study in
``cases/fin_test_3/README.md``. The tolerance is that measured offset plus
margin, not what passes.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

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


clpt = _module(ROOT / "cases/smoke_cantilever/clpt.py", "ce_gate_clpt")
ce = _module(ROOT / "cases/fin_test_3/compressed_elastica.py", "compressed_elastica")
rtw = _module(ROOT / "cases/fin_test_3/run_tipweight.py", "ce_gate_run_tipweight")

LENGTH_MM = 120.0
WIDTH_MM = 20.0
PLY_MM = 0.25
N_SPAN = 60
N_CHORD = 1
# Symmetric, so B = 0 and the strip is a pure column: no bend-extension coupling
# to muddy the comparison with a coupling-free closed form.
STACK = (0.0, 90.0, 90.0, 0.0)
LONG_AXIS = "y"

# Imperfection that breaks the buckling bifurcation: a constant-curvature bow of
# radius BOW_RADIUS_MULT * L, i.e. a tip offset near L/2000. Small enough that
# the column is ~0.16% soft against the perfect closed form and large enough that
# ccx never has to pick a branch and diverge.
BOW_RADIUS_MULT = 1000.0

MATERIAL = dict(
    e1=PLACEHOLDER_CFRP.e1,
    e2=PLACEHOLDER_CFRP.e2,
    nu12=PLACEHOLDER_CFRP.nu12,
    g12=PLACEHOLDER_CFRP.g12,
)
EI = clpt.ei_narrow_strip(STACK, PLY_MM, width=WIDTH_MM, **MATERIAL)

# Firmly post-buckling, where the bench sits (90 deg). Below ~35 deg the imperfect
# response is near-Euler and sensitive; that is not the regime being validated.
TARGET_ANGLES_DEG = (45.0, 65.0, 85.0)

# Measured -0.15% .. -0.17% over these angles, mesh-independent 60->120 elements.
# Tighten only against a fresh study, and never toward what happens to pass.
ELASTICA_TOL = 0.004


def _mesh_inp(radius_mult: float = BOW_RADIUS_MULT) -> str:
    return mesh_outline(
        Outline.rectangle(chord=WIDTH_MM, span=LENGTH_MM),
        n_chord=N_CHORD,
        n_span=N_SPAN,
        camber=CircularCamber(radius=radius_mult * LENGTH_MM),
        heading=(
            f"compressed_elastica gate: {LENGTH_MM:g} x {WIDTH_MM:g} mm strip, "
            f"{N_SPAN} S8R along the length, bowed R = {radius_mult:g} L. "
            "Built by compfea.geometry."
        ),
    ).to_inp()


def _deck(delta_mm: float, radius_mult: float = BOW_RADIUS_MULT) -> str:
    """Strip driven axially to ``delta_mm`` of shortening, tip free to swing.

    ``Outline.rectangle`` puts the root at y = 0 and the tip edge at y = span, so
    compression drives ``far_face`` in -y. The tip keeps free translation in x
    and z and all three rotations; ``fixed_end`` is clamped in all six and is
    where the reaction is read.
    """
    layup = Layup.uniform(
        [Ply(PLY_MM, angle) for angle in STACK],
        long_axis=LONG_AXIS,
        elset="blade",
        material=PLACEHOLDER_CFRP,
    )
    body = axial_drive_body("far_face", 2, -abs(delta_mm))
    step = StaticStep(body, inc=6000, static_line="0.01, 1.0, 1.E-9, 0.05")
    return assemble(
        mesh_inp=_mesh_inp(radius_mult),
        layup=layup,
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=[step],
        heading=(
            f"compressed_elastica gate: [0/90]s strip, axial shortening "
            f"{delta_mm:g} mm (delta/L = {delta_mm / LENGTH_MM:g}). "
            "Built by tests/test_compressed_elastica.py."
        ),
    )


def _solve_to_angle(
    run_dir: Path, phi_l_deg: float, radius_mult: float = BOW_RADIUS_MULT
):
    delta = ce.shortening_ratio(np.radians(phi_l_deg)) * LENGTH_MM
    run_dir.mkdir(parents=True, exist_ok=True)
    deck = run_dir / "deck.inp"
    deck.write_text(_deck(delta, radius_mult))
    return solve(deck, run_dir, final_time=1.0)


def _weight(result) -> float:
    """Equivalent hung weight: the axial reaction at the clamp, N, sign dropped."""
    return abs(result.total_force("fixed_end", "fy"))


@pytest.fixture(scope="module")
def solved(tmp_path_factory):
    """One solve per target tip angle. Well under half a minute in total."""
    root = tmp_path_factory.mktemp("compressed_elastica")
    return {
        deg: _solve_to_angle(root / f"phi_{deg:g}", deg)
        for deg in TARGET_ANGLES_DEG
    }


@pytest.mark.parametrize("deg", TARGET_ANGLES_DEG)
def test_buckling_load_matches_closed_form(solved, deg):
    """The axial reaction at each driven shortening is the elastica's W(phi_L)."""
    expected = ce.axial_load_at_angle(np.radians(deg), ei=EI, length=LENGTH_MM)
    assert _weight(solved[deg]) == pytest.approx(expected, rel=ELASTICA_TOL)


def test_it_is_post_buckling_not_a_stiffer_mode(solved):
    """W at 45 deg sits just above the clamped-free Euler load, and above it.

    A model that is not actually post-buckling in the clamped-free mode -- NLGEOM
    lost, or the tip transversely pinned -- misses this badly: clamped-pinned
    would be ~8x pi^2 EI / 4L^2, small-deflection theory would sit below it.
    """
    euler = ce.euler_load(ei=EI, length=LENGTH_MM)
    ratio = _weight(solved[45.0]) / euler
    assert 1.03 < ratio < 1.15
    assert "NLGEOM" in _deck(5.0)


def test_imperfection_amplitude_barely_moves_the_load(tmp_path):
    """Post-buckling W is set by the column, not by the bow that seeds it.

    phi_L = 65 deg at three bow radii over a 3x range; the reaction must barely
    move. This is what lets the fin run treat the operator's hand-bias as
    immaterial to the reported force.
    """
    weights = [
        _weight(_solve_to_angle(tmp_path / f"bow_{i}", 65.0, radius_mult=rm))
        for i, rm in enumerate((1500.0, 1000.0, 500.0))
    ]
    assert max(weights) / min(weights) - 1.0 < 0.005


def test_tip_angle_from_two_stations_matches_the_closed_form(tmp_path):
    """`run_tipweight` reads theta from two .dat U stations, not the .frd.

    Drive the strip to the shortening the closed form predicts for a 65 deg tip
    angle and check the two-station tangent recovers ~65 deg -- with the segment
    length within 2% of inextensible, i.e. a tangent and not a chord across
    curvature. This is the check the frd-midsurface path failed on the fin.
    """
    phi_deg = 65.0
    delta = ce.shortening_ratio(np.radians(phi_deg)) * LENGTH_MM
    mesh = mesh_outline(
        Outline.rectangle(chord=WIDTH_MM, span=LENGTH_MM),
        n_chord=N_CHORD,
        n_span=N_SPAN,
        camber=CircularCamber(radius=BOW_RADIUS_MULT * LENGTH_MM),
    )
    clip = rtw.tip_clip(mesh, long_axis=LONG_AXIS, n=3)
    band = rtw.tip_band(
        mesh, long_axis=LONG_AXIS, near_mm=4.0, far_mm=18.0
    )
    mesh = dataclasses.replace(
        mesh, nsets={**mesh.nsets, "clip": clip, "tip_band": band}
    )
    layup = Layup.uniform(
        [Ply(PLY_MM, angle) for angle in STACK],
        long_axis=LONG_AXIS,
        elset="blade",
        material=PLACEHOLDER_CFRP,
    )
    body = axial_drive_body(
        "clip", 2, -delta, tangent_nsets=("tip_band",)
    )
    deck = assemble(
        mesh_inp=mesh.to_inp(),
        layup=layup,
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=[StaticStep(body, inc=6000, static_line="0.01, 1.0, 1.E-9, 0.05")],
    )
    run_dir = tmp_path / "angle"
    run_dir.mkdir()
    (run_dir / "deck.inp").write_text(deck)
    solve(run_dir / "deck.inp", run_dir, final_time=1.0)

    dat = run_dir / "job.dat"
    clip_u = rtw.read_set_disp(dat, "clip")
    band_u = rtw.read_set_disp(dat, "tip_band")
    t = max(clip_u)
    clip0 = np.mean([mesh.nodes[n] for n in clip], axis=0)
    band0 = np.mean([mesh.nodes[n] for n in band], axis=0)
    theta, ratio = rtw.station_tangent_deg(
        tuple(clip0),
        tuple(clip0 + np.asarray(clip_u[t])),
        tuple(band0),
        tuple(band0 + np.asarray(band_u[t])),
    )
    assert ratio == pytest.approx(1.0, abs=0.02)
    assert theta == pytest.approx(phi_deg, rel=0.05)


def test_the_gate_is_fast(solved):
    total = sum(r.wall_time_s for r in solved.values())
    assert total < 30.0, f"compressed_elastica gate took {total:.1f} s of solver"
