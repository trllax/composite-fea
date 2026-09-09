"""Tip twist-stiffness metric: deck cards, chordwise node split, and the
comparative k_twist read-back.

The string tests prove the twist step says what it was asked to say and that the
regression-gated ``axial_drive_body`` text did not move. The end-to-end tests
solve two flat strips and check the physics that motivates the metric: a +/-45
laminate resists a tip twist far better than a cross-ply of equal thickness.
There is deliberately no skipif on a missing ccx.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from compfea import twist
from compfea.deck import StaticStep, assemble, axial_drive_body, twist_couple_body
from compfea.geometry import CircularCamber, Mesh, Outline, mesh_outline
from compfea.layup import PLACEHOLDER_CFRP, Layup, Ply
from compfea.run import solve

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "cases" / "fin_test_3"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ce = _load("twist_gate_compressed_elastica", CASE / "compressed_elastica.py")


def _replace_nsets(mesh: Mesh, extra: dict) -> Mesh:
    import dataclasses

    return dataclasses.replace(mesh, nsets={**mesh.nsets, **extra})


# --------------------------------------------------------------------------
# 1. twist_couple_body card text
# --------------------------------------------------------------------------


def test_twist_couple_body_cards():
    body = twist_couple_body(le_nset="twist_le", te_nset="twist_te", force=5.0)
    lines = body.splitlines()
    assert lines[0] == "*CLOAD"
    assert lines[1] == "twist_le, 2, 5.0000000000"
    assert lines[2] == "twist_te, 2, -5.0000000000"
    # RF on le/te (the free-DOF check) and on fixed_end (load path)
    assert "*NODE PRINT, NSET=twist_le, TOTALS=YES\nRF" in body
    assert "*NODE PRINT, NSET=twist_te, TOTALS=YES\nRF" in body
    assert "*NODE PRINT, NSET=fixed_end, TOTALS=YES\nRF" in body
    assert "*NODE PRINT, NSET=twist_le\nU" in body
    assert "*NODE PRINT, NSET=twist_te\nU" in body
    assert body.rstrip().endswith("*EL PRINT, ELSET=blade, TOTALS=ONLY\nELSE")

    # it perturbs the held buckled state -- no restated drive/clamp, no *BOUNDARY
    assert "*BOUNDARY" not in body
    assert "clip," not in body


def test_twist_couple_body_op_new_and_dof():
    body = twist_couple_body(
        le_nset="le", te_nset="te", force=-3.0, dof=1, op_new=True,
        read_nsets=("root", "mid"), energy_elset="strip",
    )
    assert body.splitlines()[0] == "*CLOAD, OP=NEW"
    assert "le, 1, -3.0000000000" in body
    assert "te, 1, 3.0000000000" in body
    assert "*NODE PRINT, NSET=root, TOTALS=YES\nRF" in body
    assert "*NODE PRINT, NSET=mid, TOTALS=YES\nRF" in body
    assert "*EL PRINT, ELSET=strip, TOTALS=ONLY\nELSE" in body


def test_twist_couple_body_rejects_bad_input():
    with pytest.raises(ValueError):
        twist_couple_body(le_nset="a", te_nset="b", force=0.0)
    with pytest.raises(ValueError):
        twist_couple_body(le_nset="a", te_nset="b", force=1.0, dof=7)


# --------------------------------------------------------------------------
# 2. the regression-gated axial_drive_body text is unchanged
# --------------------------------------------------------------------------


def test_axial_drive_body_text_did_not_move():
    """A frozen snapshot so a reviewer sees the gated body is byte-for-byte as it
    was. If this fails, tests/test_compressed_elastica.py and test_kick.py are
    the authority -- do not "fix" this by editing the expected string."""
    body = axial_drive_body("clip", 2, -180.0, tangent_nsets=("tip_band",))
    assert body == (
        "*BOUNDARY\n"
        "clip, 2, 2, -180.0000000000\n"
        "*NODE PRINT, NSET=fixed_end, TOTALS=YES\n"
        "RF\n"
        "*NODE PRINT, NSET=clip, TOTALS=YES\n"
        "RF\n"
        "*NODE PRINT, NSET=clip\n"
        "U\n"
        "*NODE PRINT, NSET=tip_band\n"
        "U\n"
        "*EL PRINT, ELSET=blade, TOTALS=ONLY\n"
        "ELSE"
    )


# --------------------------------------------------------------------------
# 3-5. chordwise node split, arms, couple (no solver)
# --------------------------------------------------------------------------

_STRIP_CHORD = 20.0
_STRIP_SPAN = 120.0


def _flat_strip_mesh(*, n_chord: int, n_span: int) -> Mesh:
    return mesh_outline(
        Outline.rectangle(chord=_STRIP_CHORD, span=_STRIP_SPAN),
        n_chord=n_chord,
        n_span=n_span,
    )


def test_tip_chord_nsets_split_disjoint_equal_arms():
    mesh = _flat_strip_mesh(n_chord=6, n_span=40)
    exclude = set(mesh.nsets["fixed_end"])
    sets = twist.tip_chord_nsets(mesh, long_axis="y", k=4, exclude=exclude)
    le, te = sets["twist_le"], sets["twist_te"]

    assert len(le) >= 2 and len(te) >= 2
    assert len(le) == len(te)  # equal counts -> a +F/-F couple has zero net force
    assert set(le).isdisjoint(te)
    assert set(le).isdisjoint(exclude) and set(te).isdisjoint(exclude)
    assert set(le) <= set(mesh.nsets["far_face"])
    assert set(te) <= set(mesh.nsets["far_face"])

    mesh = _replace_nsets(mesh, sets)
    x_mid, arm_le, arm_te = twist.section_arms(mesh, "twist_le", "twist_te", long_axis="y")
    assert all(mesh.nodes[n][0] > x_mid for n in le)
    assert all(mesh.nodes[n][0] < x_mid for n in te)
    assert arm_le == pytest.approx(arm_te, rel=0.2)
    assert arm_le > 0.3 * (_STRIP_CHORD / 2)
    assert all(abs(mesh.nodes[n][1] - _STRIP_SPAN) <= 1e-6 for n in le + te)


def test_tip_chord_nsets_rejects_a_degenerate_pick():
    mesh = _flat_strip_mesh(n_chord=6, n_span=40)
    with pytest.raises(ValueError):
        twist.tip_chord_nsets(mesh, long_axis="y", chord_frac_min=1.1)
    with pytest.raises(ValueError):
        twist.tip_chord_nsets(
            mesh, long_axis="y", exclude=set(mesh.nsets["far_face"])
        )


def test_couple_nmm_is_force_times_arm_sum():
    mesh = Mesh(
        nodes={
            1: (10.0, 100.0, 0.0), 2: (10.0, 101.0, 0.0),
            3: (-10.0, 100.0, 0.0), 4: (-10.0, 101.0, 0.0),
        },
        elements={},
        nsets={"twist_le": (1, 2), "twist_te": (3, 4)},
        elsets={},
        heading="",
    )
    # arms |x - x_mid| = 10 for all four nodes, force 5 -> 5 * (10*4) = 200
    assert twist.couple_nmm(
        mesh, "twist_le", "twist_te", force=5.0, long_axis="y"
    ) == pytest.approx(200.0)


# --------------------------------------------------------------------------
# 6. twist_stiffness from a hand-written .dat
# --------------------------------------------------------------------------


def _synthetic_dat() -> str:
    """le at x = +10, te at x = -10, span_def 20 mm (ux = 0), force 5 N/node,
    2 nodes per side.

    couple(+F) = 5 * 2 * 20 = 200 N.mm. Baseline uy = 50 on both edges; +F moves
    le +0.1 / te -0.1 (rotation +0.01 rad), -F the reverse. phi_half = 0.01 rad
    -> k_twist = 200 / 0.01 = 20000 N.mm/rad. le/te RF and heel RF swing = 0.
    """
    out = []

    def disp(nset, t, vy):
        out.append(
            f"\n displacements (vx,vy,vz) for set {nset} and time  {t:.6f}\n\n"
            f"      1  0.0000000E+00  {vy: .7E}  0.0000000E+00\n"
            f"      2  0.0000000E+00  {vy: .7E}  0.0000000E+00\n"
        )

    def total(nset, t, fy):
        out.append(
            f"\n total force (fx,fy,fz) for set {nset} and time  {t:.6f}\n\n"
            f"  0.0000000E+00  {fy: .7E}  0.0000000E+00\n"
        )

    disp("TWIST_LE", 1.0, 50.0)
    disp("TWIST_TE", 1.0, 50.0)
    disp("TWIST_LE", 2.0, 50.1)
    disp("TWIST_TE", 2.0, 49.9)
    disp("TWIST_LE", 3.0, 49.9)
    disp("TWIST_TE", 3.0, 50.1)
    for t in (2.0, 3.0):
        total("TWIST_LE", t, 0.0)
        total("TWIST_TE", t, 0.0)
        total("FIXED_END", t, 0.0)
    return "".join(out)


def _synthetic_mesh() -> Mesh:
    return Mesh(
        nodes={
            1: (10.0, 100.0, 0.0), 2: (10.0, 101.0, 0.0),
            3: (-10.0, 100.0, 0.0), 4: (-10.0, 101.0, 0.0),
        },
        elements={},
        nsets={"twist_le": (1, 2), "twist_te": (3, 4)},
        elsets={},
        heading="",
    )


def test_twist_stiffness_from_synthetic_dat(tmp_path):
    dat = tmp_path / "job.dat"
    dat.write_text(_synthetic_dat())

    res = twist.twist_stiffness(
        dat, _synthetic_mesh(),
        drive_end_time=1.0, plus_time=2.0, minus_time=3.0,
        force_n=5.0, theta_at_probe_deg=90.0,
    )
    assert res.couple_nmm == pytest.approx(200.0, rel=1e-9)
    assert res.chord_span_def_mm == pytest.approx(20.0, rel=1e-9)
    assert res.k_twist_nmm_per_rad == pytest.approx(20000.0, rel=1e-6)
    assert res.phi_deg == pytest.approx(math.degrees(0.01), rel=1e-6)
    assert res.le_te_asymmetry == pytest.approx(0.0, abs=1e-9)
    assert res.rf_le_n == pytest.approx(0.0, abs=1e-9)
    assert res.heel_rf_swing_n == pytest.approx(0.0, abs=1e-9)
    assert res.warning is None


def test_twist_stiffness_warns_theta_out_of_band(tmp_path):
    dat = tmp_path / "job.dat"
    dat.write_text(_synthetic_dat())
    res = twist.twist_stiffness(
        dat, _synthetic_mesh(),
        drive_end_time=1.0, plus_time=2.0, minus_time=3.0,
        force_n=5.0, theta_at_probe_deg=40.0,
    )
    assert res.warning is not None and "outside" in res.warning


def test_twist_stiffness_none_k_when_swing_unresolved(tmp_path):
    # a swing far below phi_lo_deg -> k_twist is None, not a huge number
    dat = tmp_path / "job.dat"
    out = _synthetic_dat().replace(
        " 5.0100000E+01", " 5.0000100E+01"  # le +F: +1e-4 instead of +0.1
    ).replace(" 4.9900000E+01", " 4.9999900E+01")  # te +F: -1e-4
    dat.write_text(out)
    res = twist.twist_stiffness(
        dat, _synthetic_mesh(),
        drive_end_time=1.0, plus_time=2.0, minus_time=3.0,
        force_n=5.0, theta_at_probe_deg=90.0,
    )
    assert res.k_twist_nmm_per_rad is None
    assert res.warning is not None and "None" in res.warning


# --------------------------------------------------------------------------
# 7-9. end to end: +/-45 is torsionally stiffer than 0/90 at equal thickness
# --------------------------------------------------------------------------

_E2E_CHORD = 24.0
_E2E_SPAN = 84.0
_E2E_PLY = 0.25
_E2E_FORCE_N = 8.0
_E2E_DRIVE_DEG = 88.0  # near the metric's design point; strip is prismatic so it holds


def _twist_run(run_dir: Path, stack, *, n_span: int, n_chord: int = 6):
    import dataclasses

    run_dir.mkdir(parents=True, exist_ok=True)
    length = _E2E_SPAN
    mesh = mesh_outline(
        Outline.rectangle(chord=_E2E_CHORD, span=length),
        n_chord=n_chord,
        n_span=n_span,
        camber=CircularCamber(radius=1000.0 * length),
    )
    sets = twist.tip_chord_nsets(
        mesh, long_axis="y", exclude=set(mesh.nsets["fixed_end"])
    )
    # Drive the tip edge minus its two chord extremes: stable like the buckling
    # gate, but leaves twist_le / twist_te free for the *CLOAD couple to twist
    # (a *CLOAD on a *BOUNDARY-constrained DOF does nothing).
    edges = set(sets["twist_le"]) | set(sets["twist_te"])
    drive_set = tuple(sorted(n for n in mesh.nsets["far_face"] if n not in edges))
    mesh = dataclasses.replace(
        mesh, nsets={**mesh.nsets, "drive_set": drive_set, **sets}
    )

    layup = Layup.uniform(
        [Ply(_E2E_PLY, a) for a in stack],
        long_axis="y",
        elset="blade",
        material=PLACEHOLDER_CFRP,
    )
    delta = ce.shortening_ratio(math.radians(_E2E_DRIVE_DEG)) * length

    steps = [
        StaticStep(
            axial_drive_body(
                "drive_set", 2, -abs(delta),
                tangent_nsets=("twist_le", "twist_te"),
            ),
            inc=6000,
            static_line="0.01, 1.0, 1.E-9, 0.05",
        )
    ]
    for sgn in (1.0, -1.0):
        steps.append(
            StaticStep(
                twist_couple_body(
                    le_nset="twist_le", te_nset="twist_te",
                    force=sgn * _E2E_FORCE_N, op_new=sgn < 0.0,
                ),
                inc=2000,
                static_line="0.25, 1.0, 1.E-6, 1.0",
            )
        )

    deck = assemble(
        mesh_inp=mesh.to_inp(),
        layup=layup,
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=steps,
        heading=f"twist gate {stack} n_span={n_span}",
    )
    (run_dir / "deck.inp").write_text(deck)
    result = solve(run_dir / "deck.inp", run_dir, final_time=3.0)
    res = twist.twist_stiffness(
        run_dir / "job.dat",
        mesh,
        drive_end_time=1.0,
        plus_time=2.0,
        minus_time=3.0,
        force_n=_E2E_FORCE_N,
        theta_at_probe_deg=float(_E2E_DRIVE_DEG),
    )
    return result, res


@pytest.fixture(scope="module")
def e2e(tmp_path_factory):
    root = tmp_path_factory.mktemp("twist_e2e")
    return {
        "pm45": _twist_run(root / "pm45", (45.0, -45.0, -45.0, 45.0), n_span=18),
        "uni0": _twist_run(root / "uni0", (0.0, 0.0, 0.0, 0.0), n_span=18),
        # refine both directions: n_span AND n_chord (the latter moves the
        # load-introduction patch, which n_span alone cannot).
        "uni0_fine": _twist_run(
            root / "uni0_fine", (0.0, 0.0, 0.0, 0.0), n_span=26, n_chord=8
        ),
    }


def test_pm45_is_torsionally_stiffer_than_fiber_aligned_end_to_end(e2e):
    _, pm45 = e2e["pm45"]
    _, uni0 = e2e["uni0"]
    assert pm45.k_twist_nmm_per_rad is not None and pm45.k_twist_nmm_per_rad > 0.0
    assert uni0.k_twist_nmm_per_rad is not None and uni0.k_twist_nmm_per_rad > 0.0
    ratio = pm45.k_twist_nmm_per_rad / uni0.k_twist_nmm_per_rad
    # +/-45 feeds D66 with fibres; a 0-only stack leaves it to the matrix. The
    # tip-twist metric dilutes the laminate D66 ratio (~5-10x) but must still
    # see it clearly. Observed on this coarse strip: ratio ~2.0 (record it here
    # so a regression from 2x to 1.1x is visible, not silently green).
    assert ratio > 1.5, f"pm45/uni0 k_twist ratio {ratio:.2f} (expected ~2)"
    assert uni0.le_te_asymmetry < 0.15
    assert pm45.le_te_asymmetry < 0.4
    assert pm45.warning is None
    assert uni0.warning is None


def test_twist_metric_is_repeatable_under_refinement(e2e):
    _, coarse = e2e["uni0"]
    _, fine = e2e["uni0_fine"]
    rel = abs(fine.k_twist_nmm_per_rad - coarse.k_twist_nmm_per_rad) / coarse.k_twist_nmm_per_rad
    assert rel < 0.20, f"k_twist moved {rel * 100:.0f}% from 18x6 -> 26x8"


def test_twist_gate_is_fast(e2e):
    # three buckling solves plus two couple steps each; a coarse strip. The
    # single-solve gates budget 30 s; this one gets a little more headroom.
    total = sum(result.wall_time_s for result, _ in e2e.values())
    assert total < 80.0, f"twist gate took {total:.1f} s of solver"
