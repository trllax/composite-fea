"""stress.py: parsing, ply arithmetic, and the material-frame correction.

The correction is the point of this module. ccx prints Cauchy stress in a basis
frozen at t=0, so at large rotation the axial stress leaks into the
through-thickness component -- 392 MPa of it at 89 degrees on a strip that
carries none. Every test that touches rotation is really testing that the leak
is undone.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from compfea.run import SolveError
from compfea.stress import (
    IPS_PER_PLY,
    in_plane_to_fibre,
    ips_per_ply_from,
    material_frame,
    parse_dat_stress,
    phi_from_stress,
    ply_of_ip,
    stress_summary,
    unrotate,
)

FIX = Path(__file__).parent / "fixtures"


# --- parsing ---------------------------------------------------------------

def test_local_rows_carry_an_orientation_name_and_global_rows_do_not():
    """Field count is the discriminator: 9 fields local, 8 global."""
    local = parse_dat_stress(FIX / "stress_block.dat")
    assert set(local["frame"]) == {"local"}
    assert local["ori"].iloc[0] == "ORI_P00_shell_000000"

    glob = parse_dat_stress(FIX / "stress_block_global.dat")
    assert set(glob["frame"]) == {"global"}
    assert glob["ori"].isna().all()


def test_the_blank_line_after_the_header_does_not_end_the_block():
    """ccx writes a blank line between the header and the first record."""
    df = parse_dat_stress(FIX / "stress_block.dat")
    first = df[df.time == 0.5]
    assert len(first) == 3  # all three rows, not zero


def test_a_following_energy_block_is_not_swallowed_as_stress():
    """The fixture puts an ELSE block between two stress blocks."""
    df = parse_dat_stress(FIX / "stress_block.dat")
    assert sorted(df["time"].unique()) == [0.5, 1.0]
    assert len(df) == 5  # 3 + 2, and none of the energy lines
    assert 123.4 not in df["sxx"].values


def test_increment_is_derived_from_distinct_times():
    df = parse_dat_stress(FIX / "stress_block.dat")
    assert sorted(df["increment"].unique()) == [1, 2]


def test_a_file_with_no_stress_blocks_names_the_card_it_needs(tmp_path):
    empty = tmp_path / "job.dat"
    empty.write_text(" total internal energy for set BLADE and time  0.1E+01\n\n 1.0\n")
    with pytest.raises(SolveError, match="GLOBAL=YES"):
        parse_dat_stress(empty)


def test_a_missing_file_raises_rather_than_returning_empty(tmp_path):
    with pytest.raises(SolveError, match="no results file"):
        parse_dat_stress(tmp_path / "absent.dat")


# --- ply arithmetic --------------------------------------------------------

def test_ply_zero_is_the_minus_z_ply_and_thickness_is_the_slowest_index():
    """Load-bearing: reversing this mirrors every unsymmetric stack silently."""
    assert ply_of_ip(1, 4) == 0                 # first IP -> first section line
    assert ply_of_ip(IPS_PER_PLY, 4) == 0       # last IP of ply 0
    assert ply_of_ip(IPS_PER_PLY + 1, 4) == 1   # first IP of ply 1
    assert ply_of_ip(32, 4) == 3                # last IP -> +z ply


def test_integration_points_per_ply_are_measured_not_assumed():
    """S8R gives 8 per ply, S6 gives 6 -- both live in the same ELSET.

    Assuming 8 on a triangle is worse than wrong: a 4-ply S6 tops out at ip 24,
    so ``ply_of_ip``'s bounds check never fires, plies are silently mis-assigned
    and the top ply is never reached.
    """
    assert ips_per_ply_from(32, 4) == IPS_PER_PLY   # S8R, 4 plies
    assert ips_per_ply_from(24, 4) == 6             # S6, 4 plies
    assert ips_per_ply_from(48, 6) == IPS_PER_PLY   # S8R, 6 plies


def test_an_indivisible_integration_point_count_is_refused():
    """Silently flooring would mis-assign every ply above the first."""
    with pytest.raises(ValueError, match="do not divide"):
        ips_per_ply_from(30, 4)
    with pytest.raises(ValueError, match="at least one ply"):
        ips_per_ply_from(32, 0)


def test_a_six_point_element_maps_plies_correctly():
    """The S6 case the default would get wrong."""
    assert ply_of_ip(1, 4, 6) == 0
    assert ply_of_ip(6, 4, 6) == 0
    assert ply_of_ip(7, 4, 6) == 1      # the default of 8 would say ply 0
    assert ply_of_ip(24, 4, 6) == 3     # the default would never reach ply 3


def test_an_ip_past_the_stack_is_refused():
    with pytest.raises(ValueError, match="ips_per_ply"):
        ply_of_ip(33, 4)


def test_integration_points_are_one_based():
    with pytest.raises(ValueError, match="1-based"):
        ply_of_ip(0, 4)


# --- the rotation ----------------------------------------------------------

def _bending_tensor(sigma_axial: float) -> np.ndarray:
    """Pure axial stress in the material frame, long_axis='y' (axial = 1)."""
    t = np.zeros((3, 3))
    t[1, 1] = sigma_axial
    return t


@pytest.mark.parametrize("deg", [0.0, 30.0, 89.0, 90.0, 180.0])
def test_a_rotated_tensor_round_trips(deg):
    """Rotate a pure axial state by phi, then undo it."""
    mat = _bending_tensor(950.0)
    phi = math.radians(deg)
    # forward-rotate into the frozen frame the way the solver's output is posed
    fixed = unrotate(mat, -phi, width_axis=0)
    back = unrotate(fixed, phi, width_axis=0)
    assert back[1, 1] == pytest.approx(950.0, abs=1e-9)
    assert back[2, 2] == pytest.approx(0.0, abs=1e-9)


def test_the_through_thickness_leak_is_real_before_correction():
    """At 89 degrees a pure axial state looks like large szz if left alone."""
    fixed = unrotate(_bending_tensor(950.0), -math.radians(89.0), width_axis=0)
    assert abs(fixed[2, 2]) > 0.4 * 950.0     # the artefact this module removes
    assert fixed[1, 1] + fixed[2, 2] == pytest.approx(950.0, abs=1e-9)  # invariant


def test_rotating_the_wrong_way_leaves_a_large_residual():
    """The sign guard: undoing +phi with -phi must not look correct."""
    phi = math.radians(60.0)
    fixed = unrotate(_bending_tensor(950.0), -phi, width_axis=0)
    right = unrotate(fixed, phi, width_axis=0)
    wrong = unrotate(fixed, -phi, width_axis=0)
    assert abs(right[2, 2]) < 1e-9
    assert abs(wrong[2, 2]) > 0.5 * 950.0


def test_phi_is_recoverable_from_the_stress_itself():
    """Independent of geometry, so it can cross-check the geometric angle."""
    for deg in (10.0, 40.0, 75.0):
        fixed = unrotate(_bending_tensor(950.0), -math.radians(deg), width_axis=0)
        row = type("R", (), dict(
            sxx=fixed[0, 0], syy=fixed[1, 1], szz=fixed[2, 2],
            sxy=fixed[0, 1], sxz=fixed[0, 2], syz=fixed[1, 2],
        ))
        assert math.degrees(phi_from_stress(row, "y")) == pytest.approx(deg, abs=1e-6)


# --- fibre transformation --------------------------------------------------

def test_a_zero_degree_ply_puts_axial_stress_into_s11():
    s11, s22, t12 = in_plane_to_fibre(100.0, 10.0, 0.0, 0.0)
    assert (s11, s22, t12) == pytest.approx((100.0, 10.0, 0.0))


def test_a_ninety_degree_ply_swaps_the_two_normals():
    s11, s22, t12 = in_plane_to_fibre(100.0, 10.0, 0.0, 90.0)
    assert (s11, s22) == pytest.approx((10.0, 100.0))
    assert t12 == pytest.approx(0.0, abs=1e-12)


def test_a_forty_five_degree_ply_turns_axial_stress_into_shear():
    """Signed, not abs(): ``abs()`` here cannot see a mirrored convention."""
    s11, s22, t12 = in_plane_to_fibre(100.0, 0.0, 0.0, 45.0)
    assert s11 == pytest.approx(50.0)
    assert s22 == pytest.approx(50.0)
    assert t12 == pytest.approx(-50.0)


@pytest.mark.parametrize("long_axis", ["x", "y"])
@pytest.mark.parametrize("angle", [0.0, 30.0, 45.0, -45.0, 90.0])
def test_the_fibre_transform_matches_the_orientation_card_ccx_is_given(
    long_axis, angle
):
    """Handedness guard, against the very card `layup.py` writes into the deck.

    The transform must agree with sigma_11 = e1.sigma.e1 built from the
    ``*ORIENTATION`` direction cosines, or the post-processing is describing a
    different laminate than the solver was given. A cross-ply cannot detect a
    mirror -- tau_12 is ~0 at 0 and 90 degrees -- so this sweeps off-axis angles
    and a state with real shear.
    """
    from compfea.layup import orientation_card

    nums = [
        float(v)
        for v in orientation_card(angle, long_axis=long_axis).splitlines()[1].split(",")
    ]
    e1, e2 = np.array(nums[:3]), np.array(nums[3:6])

    w, a, n, _ = __import__("compfea.stress", fromlist=["_axes"])._axes(long_axis)
    sigma = np.zeros((3, 3))
    sigma[a, a], sigma[w, w] = 900.0, 60.0
    sigma[a, w] = sigma[w, a] = 25.0          # real shear: the mirror shows here

    from compfea.stress import _axes

    _, _, _, w_sign = _axes(long_axis)
    got = in_plane_to_fibre(sigma[a, a], sigma[w, w], w_sign * sigma[w, a], angle)
    want = (e1 @ sigma @ e1, e2 @ sigma @ e2, e1 @ sigma @ e2)
    # 1e-5 is loose only against the card's 8-decimal direction cosines; a
    # mirrored convention misses by ~50 MPa here, not by rounding.
    assert got == pytest.approx(want, abs=1e-5)


# --- against the real probe solve, when it is on disk ----------------------

PROBE = Path("results/stress_probe/probe_b.dat")
probe_only = pytest.mark.skipif(not PROBE.is_file(), reason=f"{PROBE} not on disk")


@probe_only
def test_the_probe_solve_reproduces_beam_theory_and_passes_the_shell_check():
    """End to end on real ccx output: cantilever_89deg, [0/90/90/0] at 0.25 mm.

    Independent check, not a regression against our own output: a beam of
    L=100 mm at 89 degrees has r = L/theta, and the outermost integration point
    sits at the 2-point Gauss station inside the outer ply, not at the surface.
    E1*z/r there is 937.7 MPa.
    """
    df = parse_dat_stress(PROBE)
    last = df[df.time == df.time.max()]
    assert len(last) == 512                       # 16 elements x 32 IPs
    assert set(last["elem"]) == set(range(1, 17))  # shell ids, not expanded ids

    phi = {e: phi_from_stress(g.iloc[0], "y") for e, g in last.groupby("elem")}
    mat = material_frame(last, phi, [0.0, 90.0, 90.0, 0.0], long_axis="y")
    summary = stress_summary(mat)

    # tension on one face, compression on the other, both in the 0-degree plies
    assert summary["s11_max_t"] == pytest.approx(937.7, rel=0.03)
    assert summary["s11_max_c"] == pytest.approx(-937.7, rel=0.03)
    assert summary["s11_max_t_angle"] == 0.0
    assert summary["s11_max_c_angle"] == 0.0
    assert summary["s11_max_t_ply"] != summary["s11_max_c_ply"]

    # the shell cannot carry through-thickness normal stress
    assert summary["s33_residual_max"] < 0.05


@probe_only
def test_the_recovered_rotation_ramps_along_the_span():
    """Uniform curvature under a prescribed end rotation: phi grows with span."""
    df = parse_dat_stress(PROBE)
    last = df[df.time == df.time.max()]
    phi = [
        math.degrees(phi_from_stress(g.iloc[0], "y"))
        for _, g in sorted(last.groupby("elem"), key=lambda kv: kv[0])
    ]
    assert phi == sorted(phi)                  # monotone root to tip
    assert phi[0] < 5.0 and phi[-1] > 80.0     # spans nearly the full 89 degrees


# --- against a real stress sweep, when one is on disk ----------------------

SWEEP = Path("results/stressrun/cache/ffd8059ba2187fcd")
sweep_only = pytest.mark.skipif(
    not (SWEEP / "ccx" / "job.dat").is_file(), reason=f"{SWEEP} not on disk"
)


@sweep_only
@pytest.mark.parametrize(("step", "deg"), [(18, 90), (36, 180)])
def test_geometry_and_stress_agree_on_the_rotation_angle(step, deg):
    """The two sources of phi are independent, so agreement is real evidence.

    One comes from the deformed mid-surface in the .frd, the other from the
    principal direction of the stress tensor. Nothing links them but the
    physics, so if the un-rotation were wrong -- wrong axis, wrong sign, wrong
    element-to-span mapping -- they would diverge.

    The principal direction is pi-periodic while the geometric angle runs past
    90 degrees, so the comparison has to be wrapped; without that the two look
    ~90 degrees apart for no reason.
    """
    from compfea.frd import (
        deformed_midsurface,
        disp_at_step,
        index_disp_blocks,
        read_nodes,
    )
    from compfea.stress import bend_angle_by_span, element_spans

    frd = SWEEP / "ccx" / "job.frd"
    blocks = index_disp_blocks(frd)
    _, disp = disp_at_step(frd, step, blocks=blocks)
    phi_at = bend_angle_by_span(
        deformed_midsurface(read_nodes(frd), disp, long_axis="y")
    )
    spans = element_spans(SWEEP / "deck.inp", long_axis="y")

    df = parse_dat_stress(SWEEP / "ccx" / "job.dat")
    block = df[(df["time"] - float(step)).abs() <= 1e-6]
    assert not block.empty, f"no stress printed at step {step}"

    diffs = []
    for elem, group in block.groupby("elem"):
        if elem not in spans:
            continue
        geo = math.degrees(phi_at(spans[elem]))
        sig = math.degrees(phi_from_stress(group.iloc[0], "y"))
        diffs.append(abs((geo - sig + 90.0) % 180.0 - 90.0))
    assert len(diffs) == 64
    assert max(diffs) < 5.0, f"phi sources diverge by up to {max(diffs):.1f} deg"


@sweep_only
def test_the_peak_stress_matches_local_beam_theory():
    """Independent magnitude check, at the curvature the laminate actually took.

    The strip is driven onto a circular target but does not adopt it, so the
    mean curvature L/theta under-predicts by ~20%. Against the local maximum
    curvature of the solved shape, E1*z*kappa lands within a couple of percent.
    """
    from compfea.frd import (
        deformed_midsurface,
        disp_at_step,
        read_nodes,
    )
    from compfea.stress import bend_angle_by_span

    frd = SWEEP / "ccx" / "job.frd"
    _, disp = disp_at_step(frd, 36)
    mid = deformed_midsurface(read_nodes(frd), disp, long_axis="y")
    phi_at = bend_angle_by_span(mid)
    span = mid["span0_mm"].to_numpy()
    kappa = np.gradient(np.unwrap([phi_at(v) for v in span]), span)

    e1 = 135000.0
    z = 0.15 + 0.05 / math.sqrt(3)      # outer 0-deg ply's Gauss station, t=0.4
    predicted = e1 * z * kappa.max()
    assert predicted == pytest.approx(963.8, rel=0.05)


PROBE_LOCAL = Path("results/stress_probe/probe_a.dat")
two_frame_only = pytest.mark.skipif(
    not (PROBE_LOCAL.is_file() and PROBE.is_file()),
    reason="both probe frames needed",
)


@two_frame_only
def test_the_two_frames_of_one_solve_transform_into_each_other():
    """probe_a and probe_b are the SAME solve printed in two frames.

    ccx's local frame *is* the ply's ``*ORIENTATION`` frame, so rotating the
    global block by the ply angle must reproduce the local block exactly --
    element by element, point by point, to machine precision. Nothing else in
    the suite pins the sign of the shear term this hard: before the width-basis
    handedness was fixed, the normals agreed to 1e-16 while tau_12 was out by
    exactly twice sigma_xy for every row.
    """
    local = parse_dat_stress(PROBE_LOCAL)
    glob = parse_dat_stress(PROBE)
    local = local[(local["time"] - 1.0).abs() <= 1e-9].reset_index(drop=True)
    glob = glob[(glob["time"] - 1.0).abs() <= 1e-9].reset_index(drop=True)
    assert len(local) == len(glob) == 512

    by_name = {"ORI_P00_shell_000000": 0.0, "ORI_P90_shell_000000": 90.0}
    from compfea.stress import _axes

    _, _, _, w_sign = _axes("y")
    worst = [0.0, 0.0, 0.0]
    for lo, gl in zip(local.itertuples(), glob.itertuples(), strict=True):
        assert (lo.elem, lo.ip) == (gl.elem, gl.ip)
        got = in_plane_to_fibre(gl.syy, gl.sxx, w_sign * gl.sxy, by_name[lo.ori])
        for i, (g, ref) in enumerate(zip(got, (lo.sxx, lo.syy, lo.sxy), strict=True)):
            worst[i] = max(worst[i], abs(g - ref))
    assert max(worst) < 1e-9, f"frames disagree: {worst}"
