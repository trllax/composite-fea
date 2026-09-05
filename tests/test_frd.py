"""frd.py: fixed-width parsing, block indexing, mid-surface reduction.

The fixture is ``tiny_disp.frd.txt``, not ``.frd``: ``.gitignore`` line 8 is
``*.frd``, so a fixture with the real extension passes locally and vanishes on
a fresh clone. The frd format has no comment syntax, so that note lives here.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from compfea.frd import (
    deformed_midsurface,
    disp_at_step,
    index_disp_blocks,
    parse_record,
    read_disp,
    read_nodes,
)

FIXTURE = Path(__file__).parent / "fixtures" / "tiny_disp.frd.txt"

# A verified record from a real ccx run. The third component runs together with
# the second: there is no space before the minus sign.
REAL_LINE = " -1       272 0.00000E+00 0.00000E+00-3.73579E-06"


def test_a_record_is_fixed_width_not_whitespace_separated():
    """The load-bearing test. A negative value abuts the field before it.

    Asserting the split count too, so the test says *why* ``.split()`` is
    wrong rather than only that this implementation happens to disagree.
    """
    assert len(REAL_LINE.split()) == 4  # four fields where there are five
    node, (d1, d2, d3) = parse_record(REAL_LINE)
    assert node == 272
    assert d1 == 0.0 and d2 == 0.0
    assert d3 == pytest.approx(-3.73579e-06, rel=1e-12)


def test_a_node_id_that_abuts_its_first_value_still_parses():
    """In the fixture the id itself runs into the value: ``110-1.00000E-05``."""
    line = " -1       110-1.00000E-05 0.00000E+00 4.00000E-01"
    # split fuses the id onto the value: four fields, and the second is junk
    assert line.split()[1] == "110-1.00000E-05"
    node, (d1, _, d3) = parse_record(line)
    assert node == 110
    assert d1 == pytest.approx(-1.0e-05)
    assert d3 == pytest.approx(0.4)


def test_the_short_format_uses_a_narrower_id_field():
    """NOT verified against real output -- this build of ccx never writes it.

    ``strings`` on the conda-forge ccx shows one record writer,
    ``%3s%10d%12.5E%12.5E%12.5E``, i.e. format 1 only. The format-0 slices
    follow the documented cgx convention (5-wide id, 12-wide values) and are
    here so a short-format file is read correctly rather than mis-sliced in
    silence, but this test only checks the implementation against itself. Treat
    it as a guard against accidental edits, not as evidence the layout is right.
    """
    line = " -1  272 0.00000E+00 0.00000E+00-3.73579E-06"
    node, (_, _, d3) = parse_record(line, fmt=0)
    assert node == 272
    assert d3 == pytest.approx(-3.73579e-06)


def test_a_binary_frd_is_refused_rather_than_mis_sliced():
    with pytest.raises(ValueError, match="binary"):
        parse_record(REAL_LINE, fmt=2)


def test_the_element_block_is_not_read_as_nodal_coordinates():
    """``-1`` is overloaded; the 3C block uses it with a different layout."""
    nodes = read_nodes(FIXTURE)
    assert len(nodes) == 18
    assert min(nodes) == 101 and max(nodes) == 118
    # element ids 318/319 live in the 3C block and must not appear as nodes
    assert 318 not in nodes and 319 not in nodes


def test_index_finds_every_disp_block_with_its_step_and_time():
    blocks = index_disp_blocks(FIXTURE)
    assert [(b.step, b.time) for b in blocks] == [(2, 1.4), (2, 1.8), (2, 2.0)]
    assert all(b.n_nodes == 18 for b in blocks)
    # step 1 carries a STRESS block only, and must not be indexed as DISP
    assert all(b.step != 1 for b in blocks)


def test_disp_at_step_takes_the_end_of_step_not_an_intermediate_increment():
    block, disp = disp_at_step(FIXTURE, 2)
    assert block.time == pytest.approx(2.0)
    assert len(disp) == 18
    assert disp[110][2] == pytest.approx(1.0)  # the t=2.0 lift, not 0.4 or 0.8


def test_a_step_that_never_reached_its_end_raises_instead_of_substituting(
    tmp_path,
):
    """The forbidden fallback: handing back the last converged increment.

    A solve killed part way through a step leaves increments but no end-of-step
    block. Returning the nearest one gives the caller a ~177-degree pose under
    a 180-degree label -- a curve that looks entirely reasonable and is wrong.
    """
    lines = FIXTURE.read_text().splitlines()
    # drop the t=2.0 block, leaving step 2 with only its 1.4 and 1.8 increments
    out, drop = [], False
    for line in lines:
        if line.startswith("  100CL"):
            drop = abs(float(line[12:24]) - 2.0) < 1e-9
        if not drop:
            out.append(line)
        if line.startswith(" -3"):
            drop = False
    partial = tmp_path / "partial.frd.txt"
    partial.write_text("\n".join(out) + "\n")

    blocks = index_disp_blocks(partial)
    assert [b.time for b in blocks] == [1.4, 1.8]  # increments, no end of step
    with pytest.raises(LookupError, match="never reached its end"):
        disp_at_step(partial, 2)


def test_a_step_without_disp_raises_and_names_the_steps_that_have_it():
    with pytest.raises(LookupError, match="steps with DISP: 2"):
        disp_at_step(FIXTURE, 1)


def test_a_file_with_no_disp_blocks_says_to_re_solve(tmp_path):
    stripped = tmp_path / "nodisp.frd.txt"
    keep, mode = [], True
    for line in FIXTURE.read_text().splitlines():
        if line.startswith(" -4"):
            mode = "DISP" not in line
        if mode:
            keep.append(line)
        if line.startswith(" -3"):
            mode = True
    stripped.write_text("\n".join(keep) + "\n")
    assert index_disp_blocks(stripped) == []
    with pytest.raises(LookupError, match="predates \\*NODE FILE"):
        disp_at_step(stripped, 2)


def test_read_disp_seeks_to_one_block_without_rescanning():
    blocks = index_disp_blocks(FIXTURE)
    first = read_disp(FIXTURE, blocks[0])
    last = read_disp(FIXTURE, blocks[-1])
    assert first[110][2] == pytest.approx(0.4)
    assert last[110][2] == pytest.approx(1.0)


def test_mid_surface_uses_the_extreme_pair_not_the_mean_of_the_column():
    """The fixture's mid-thickness node is deliberately displaced 4 mm off.

    Both reductions agree on an evenly weighted symmetric column, so the
    fixture breaks that symmetry: the extreme pair gives 1.0 and a mean over
    the column would give (0.9 + 5.0 + 1.1)/3 = 2.333.
    """
    nodes = read_nodes(FIXTURE)
    _, disp = disp_at_step(FIXTURE, 2)
    mid = deformed_midsurface(nodes, disp, long_axis="y")
    tip = mid.iloc[-1]
    assert tip.z_mm == pytest.approx(1.0)
    assert tip.z_mm != pytest.approx(7.0 / 3.0, rel=1e-3)


def test_a_repeated_extreme_level_is_averaged_not_picked_arbitrarily():
    """Defensive branch: neither the fixture nor the real mesh repeats an extreme.

    ccx gives each station a single node at +t/2 and one at -t/2, so the
    averaging path never runs on real output. Exercise it directly, because an
    untested branch that silently picks one node of a pair would tilt the
    mid-surface without changing anything visible.
    """
    nodes = {
        1: (0.0, 0.0, -0.1), 2: (0.0, 0.0, -0.1),   # repeated bottom
        3: (0.0, 0.0, 0.1), 4: (0.0, 0.0, 0.1),     # repeated top
    }
    # bottom pair averages to z=1.0, top pair to z=3.0 -> mid-surface 2.0
    disp = {
        1: (0.0, 0.0, 1.5), 2: (0.0, 0.0, 0.7),
        3: (0.0, 0.0, 3.1), 4: (0.0, 0.0, 2.7),
    }
    mid = deformed_midsurface(nodes, disp, long_axis="y")
    assert len(mid) == 1
    assert mid.iloc[0].z_mm == pytest.approx(2.0)


def test_mid_surface_keeps_only_the_chord_centre_column():
    nodes = read_nodes(FIXTURE)
    _, disp = disp_at_step(FIXTURE, 2)
    mid = deformed_midsurface(nodes, disp, long_axis="y")
    assert list(mid.span0_mm) == [0.0, 50.0]  # one row per span station
    assert mid.transverse_mm.iloc[0] == pytest.approx(10.0)  # x centre of 0..20


def test_mid_surface_refuses_a_disp_block_from_a_different_mesh():
    nodes = read_nodes(FIXTURE)
    _, disp = disp_at_step(FIXTURE, 2)
    disp.pop(110)
    with pytest.raises(ValueError, match="no displacement"):
        deformed_midsurface(nodes, disp, long_axis="y")


# --- against real solver output, when a sweep run is on disk -----------------
# results/ is gitignored, so these skip on a fresh clone.

REAL = Path("results/incr-calib/cache/6418231dde201796/ccx/job.frd")
real_only = pytest.mark.skipif(not REAL.is_file(), reason=f"{REAL} not on disk")


@real_only
def test_frd_node_ids_are_the_expanded_solid_mesh_not_the_shell_mesh():
    """Guards the 'never join frd onto deck.inp' rule with actual numbers."""
    nodes = read_nodes(REAL)
    assert len(nodes) == 2484  # 261 planform stations, expanded per ply
    assert min(nodes) == 265  # ids do not start at 1: not the deck's numbering
    assert max(nodes) == 4176


@real_only
@pytest.mark.parametrize(("step", "theta_deg"), [(18, 90.0), (36, 180.0)])
def test_the_deformed_tip_matches_the_circular_arc_it_was_driven_to(
    step, theta_deg
):
    """The tip is displacement-driven onto r=L/θ, so the shape must land there.

    This is the end-to-end check that the fixed-width parse, the expanded-mesh
    handling and the mid-surface reduction are all right at once: any one of
    them wrong moves the tip off the arc.
    """
    nodes = read_nodes(REAL)
    _, disp = disp_at_step(REAL, step)
    mid = deformed_midsurface(nodes, disp, long_axis="y")

    theta = math.radians(theta_deg)
    r = 100.0 / theta
    assert len(mid) == 65  # 32 spanwise elements, quad8 -> 65 stations
    root = mid.iloc[0]
    assert root.span_mm == pytest.approx(0.0, abs=1e-6)
    assert root.z_mm == pytest.approx(0.0, abs=1e-6)
    tip = mid.iloc[-1]
    assert tip.span_mm == pytest.approx(r * math.sin(theta), abs=1e-3)
    assert tip.z_mm == pytest.approx(r * (1.0 - math.cos(theta)), abs=1e-3)
