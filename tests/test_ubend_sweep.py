"""Unit tests for U-bend path + energy parse + sweep deck build (no ccx)."""

from __future__ import annotations

import math
from pathlib import Path

from dataclasses import replace

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

    d1 = Design(zone_pairs=(1,), ply_mm=0.1)
    d2 = Design(zone_pairs=(2,), ply_mm=0.1)
    assert d1.n_plies_by_zone == (4,)
    assert d2.n_plies_by_zone == (8,)
    assert d1.stack_label == "[0/90]s"
    assert d1.zone_elsets == ("blade",)
    plies = layup_for(d2).zones[0].plies
    assert len(plies) == 8
    assert [p.angle_deg for p in plies] == [0, 90, 0, 90, 90, 0, 90, 0]


def test_zoned_design_stacks_are_symmetric_and_root_first():
    from compfea.sweep import Design, layup_for

    d = Design(zone_pairs=(3, 1), zones=(0.5,), ply_mm=0.1)
    assert d.zone_elsets == ("zone_1", "zone_2")
    assert d.n_plies_by_zone == (12, 4)
    assert d.stack_label == "[0/90]_3s|[0/90]s"
    layup = layup_for(d)
    for zone in layup.zones:
        angles = [p.angle_deg for p in zone.plies]
        assert angles == angles[::-1], f"{zone.elset} stack is not symmetric"
    assert layup.zones[0].elset == "zone_1"
    assert layup.zones[0].thickness > layup.zones[1].thickness


def test_zone_boundary_count_is_enforced():
    from compfea.sweep import Design

    with pytest.raises(ValueError, match="boundaries"):
        Design(zone_pairs=(2, 1), zones=())
    with pytest.raises(ValueError, match="boundaries"):
        Design(zone_pairs=(2,), zones=(0.5,))


def test_the_cache_key_moves_with_everything_that_reaches_the_deck():
    """A key that misses an input silently returns a stale force.

    The old key hashed the literal string "mat=placeholder_cfrp", so a material
    edit reused the previous answer. Each of these must produce a new key.
    """
    from compfea.sweep import Design

    base = Design(zone_pairs=(2,), ply_mm=0.1, angles=(0.0, 90.0), fiber="ud")
    variants = [
        Design(zone_pairs=(3,), ply_mm=0.1, angles=(0.0, 90.0), fiber="ud"),
        Design(zone_pairs=(2,), ply_mm=0.15, angles=(0.0, 90.0), fiber="ud"),
        Design(zone_pairs=(2,), ply_mm=0.1, angles=(45.0, -45.0), fiber="ud"),
        Design(zone_pairs=(2,), ply_mm=0.1, angles=(0.0, 90.0), fiber="woven"),
        Design(zone_pairs=(2, 2), zones=(0.5,), ply_mm=0.1, fiber="ud"),
        # The solver setting changes the answer, so it must change the key.
        Design(zone_pairs=(2,), ply_mm=0.1, static_line="0.02, 1.0, 1.E-8, 0.1"),
    ]
    keys = {base.cache_key()} | {v.cache_key() for v in variants}
    assert len(keys) == len(variants) + 1
    assert base.cache_key() == Design(zone_pairs=(2,), ply_mm=0.1).cache_key()


def test_a_malformed_static_line_is_refused():
    from compfea.sweep import Design

    with pytest.raises(ValueError, match="initial, period, min, max"):
        Design(zone_pairs=(1,), static_line="0.02, 1.0, 0.1")


def test_material_fingerprint_covers_all_nine_constants():
    from compfea.layup import UD_CFRP_GENERIC
    from compfea.sweep import _material_fingerprint

    base = _material_fingerprint(UD_CFRP_GENERIC)
    for field in (
        "e1", "e2", "e3", "nu12", "nu13", "nu23", "g12", "g13", "g23", "density",
    ):
        bumped = replace(
            UD_CFRP_GENERIC, **{field: getattr(UD_CFRP_GENERIC, field) * 1.01}
        )
        assert _material_fingerprint(bumped) != base, f"{field} not in fingerprint"


def test_the_grid_is_the_cartesian_product():
    from compfea.sweep import build_parser, designs_from_args

    args = build_parser().parse_args(
        ["--zone-pairs", "2,1", "3,1", "--zones", "0.5",
         "--fiber", "ud", "woven", "--ply-mm", "0.1", "0.15"]
    )
    designs = designs_from_args(args)
    assert len(designs) == 2 * 2 * 2
    assert len({d.cache_key() for d in designs}) == 8


def test_zone_pairs_and_n_pairs_are_not_both_accepted():
    from compfea.sweep import build_parser, designs_from_args

    args = build_parser().parse_args(["--zone-pairs", "2", "--n-pairs", "3"])
    with pytest.raises(SystemExit, match="not both"):
        designs_from_args(args)


def test_n_pairs_still_means_what_it_used_to():
    from compfea.sweep import Design, build_parser, designs_from_args

    args = build_parser().parse_args(["--n-pairs", "1", "2", "3"])
    designs = designs_from_args(args)
    assert [d.zone_pairs for d in designs] == [(1,), (2,), (3,)]
    assert designs[0].cache_key() == Design(zone_pairs=(1,)).cache_key()


def test_report_deg_is_in_the_cache_key():
    """It sets how many *STEP blocks the deck has, so it must move the key.

    Otherwise a cheap 90-degree evaluation caches under the key a 180-degree
    one looks up, and the 180 request is served a row with no f_180 in it.
    """
    from compfea.sweep import Design

    assert (
        Design(zone_pairs=(1,), report_deg=(90.0,)).cache_key()
        != Design(zone_pairs=(1,), report_deg=(90.0, 180.0)).cache_key()
    )


def test_a_zone_boundary_off_the_element_rows_is_refused():
    """A fraction between rows snaps, so requested != delivered ply drop."""
    from compfea.sweep import Design

    with pytest.raises(ValueError, match="not a multiple"):
        Design(zone_pairs=(2, 1), zones=(0.47,))
    Design(zone_pairs=(2, 1), zones=(0.5,))  # on a row: fine


def test_build_deck_is_told_the_long_axis(monkeypatch):
    """It must be passed explicitly, not left to its own default of 'y'.

    layup_for and tip_length_mm both use LONG_AXIS. If build_deck is left to
    default while LONG_AXIS says 'x', the fibre direction and the moment arm
    use x while the prescribed tip path still drives y -- a deck that converges
    and reports a force for a load 90 degrees off the axis the cache key
    records. The strip spans y, so flipping the constant cannot show this; what
    matters is that the argument is handed over at all.
    """
    import compfea.sweep as sweep

    seen = {}
    real = sweep.build_deck

    def spy(mesh, layup, angles, **kwargs):
        seen.update(kwargs)
        return real(mesh, layup, angles, **kwargs)

    monkeypatch.setattr(sweep, "build_deck", spy)
    sweep.deck_for(sweep.Design(zone_pairs=(1,), report_deg=(10.0,)))
    assert seen.get("long_axis") == sweep.LONG_AXIS


def test_deck_only_and_the_solve_build_the_same_deck():
    """They were separate paths and had already drifted on end_deg."""
    from compfea.sweep import Design, angles_for, deck_for

    design = Design(zone_pairs=(1,), report_deg=(90.0,))
    _deck, angles, _arm = deck_for(design)
    assert angles == angles_for(design)
    assert max(angles) == 90.0
