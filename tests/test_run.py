"""Tests for compfea.run.

Sample solver output is inlined below rather than committed as files:
.gitignore drops *.sta and *.cvg on purpose, so a checked-in sample would
silently disappear from the repo.

There is deliberately no skipif on a missing ccx. The solver is a hard
dependency of this repo, and a suite that goes green by skipping itself is the
same failure mode as a solve that goes green by stopping early.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from compfea.run import (
    CcxFailed,
    NotConverged,
    SolveError,
    parse_dat_totals,
    parse_sta,
    solve,
)

FIXTURES = Path(__file__).parent / "fixtures"
STRIP = FIXTURES / "strip_4el.inp"
STRIP_SHORT = FIXTURES / "strip_4el_short.inp"

# The converged strip fixture. Hand-checked: the laminate gives EI = 1.999e5
# N.mm^2, so 3EI/L^3 = 0.5998 N/mm, and the solver returns 0.6002 N/mm in the
# small-deflection limit. At 20 mm the response has stiffened to -12.53 N.
# A value the code must keep reproducing, not a physical validation.
STRIP_FZ_N = -12.525
STRIP_INCREMENTS = 22
STRIP_LINEAR_STIFFNESS_N_PER_MM = 0.5998
STRIP_DRIVEN_MM = 20.0
# What the truncated fixture leaves in its .dat. Nothing may ever return this.
SHORT_PARTIAL_FZ_N = -1.381295

SAMPLE_STA = """SUMMARY OF JOB INFORMATION
  STEP      INC     ATT  ITRS     TOT TIME     STEP TIME      INC TIME
     1          1     1     2  0.200000E-01  0.200000E-01  0.200000E-01
     1          2     1     3  0.400000E-01  0.400000E-01  0.200000E-01
     1          3     1     2  0.100000E+01  0.100000E+01  0.350000E-01
"""

SAMPLE_DAT = """

                        S T E P       1


                                INCREMENT     1


 forces (fx,fy,fz) for set FIXED_END and time  0.2000000E-01

         1  1.505959E-03  4.561112E-04  2.837116E-02
        10 -3.011077E-03 -4.069159E-11 -2.968419E-01

 total force (fx,fy,fz) for set FIXED_END and time  0.2000000E-01

        8.414649E-07  6.254762E-16 -2.400996E-01

                                INCREMENT     2


 forces (fx,fy,fz) for set FIXED_END and time  0.1000000E+01

         1  3.820788E+00  1.157285E+00  1.550800E+00

 total force (fx,fy,fz) for set FIXED_END and time  0.1000000E+01

       -2.009884E-06 -9.170442E-14 -1.252544E+01
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path


def _sta(*rows: tuple[int, str, float]) -> str:
    """A .sta file from (increment, attempt, total time) rows."""
    head = (
        "SUMMARY OF JOB INFORMATION\n"
        "  STEP      INC     ATT  ITRS     TOT TIME     STEP TIME      INC TIME\n"
    )
    body = "".join(
        f"     1     {inc:6d}  {att:>4}     2  {t:.6E}  {t:.6E}  0.200000E-01\n"
        for inc, att, t in rows
    )
    return head + body


def _dat(*times: float, nset: str = "FIXED_END", fz: float = -12.52544) -> str:
    """A .dat holding one total-force block per given time."""
    return "".join(
        f"\n total force (fx,fy,fz) for set {nset} and time  {t:.7E}\n"
        f"\n       -2.009884E-06 -9.170442E-14 {fz:.6E}\n"
        for t in times
    )


def _stub_ccx(
    tmp_path: Path, *, sta: str = "", dat: str = "", exit_code: int = 0
) -> str:
    """A fake solver that writes canned output. Invoked as `ccx -i <job>`.

    This is the only way to reach the code paths where ccx exits zero on a solve
    that did not finish -- the real solver exits nonzero for those, which would
    let the exit-status branch satisfy tests meant to pin the .sta check.
    """
    script = tmp_path / "stub_ccx"
    lines = ["#!/bin/sh"]
    if sta:
        lines.append(f'cat > "$2.sta" <<\'STA\'\n{sta}STA')
    if dat:
        lines.append(f'cat > "$2.dat" <<\'DAT\'\n{dat}DAT')
    lines.append(f"exit {exit_code}")
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o755)
    return str(script)


# --- parsers, no solver needed -------------------------------------------


def test_parse_sta_reads_every_increment(tmp_path):
    status = parse_sta(_write(tmp_path, "job.sta", SAMPLE_STA))
    assert list(status.columns) == [
        "step", "inc", "att", "itrs", "tot_time", "step_time", "inc_time",
    ]
    assert len(status) == 3
    assert status["itrs"].tolist() == [2, 3, 2]
    assert status["tot_time"].iloc[-1] == pytest.approx(1.0)


def test_parse_sta_drops_unconverged_attempts(tmp_path):
    """A U-suffixed row is an attempt, not an increment."""
    text = _sta((1, "1", 0.5), (2, "1U", 1.0))
    status = parse_sta(_write(tmp_path, "job.sta", text))
    assert len(status) == 1
    assert status["tot_time"].iloc[-1] == pytest.approx(0.5)


def test_parse_sta_rejects_a_file_with_no_increments(tmp_path):
    header = "\n".join(SAMPLE_STA.splitlines()[:2])
    with pytest.raises(SolveError, match="no increments"):
        parse_sta(_write(tmp_path, "job.sta", header))


def test_parse_dat_totals_reads_totals_not_per_node_lines(tmp_path):
    reactions = parse_dat_totals(_write(tmp_path, "job.dat", SAMPLE_DAT))
    assert list(reactions.columns) == ["increment", "time", "nset", "fx", "fy", "fz"]
    assert len(reactions) == 2  # two totals blocks, not the per-node lines
    assert reactions["increment"].tolist() == [1, 2]
    assert reactions["nset"].tolist() == ["fixed_end", "fixed_end"]  # lower-cased
    assert reactions["time"].tolist() == pytest.approx([0.02, 1.0])
    assert reactions["fz"].iloc[-1] == pytest.approx(-12.52544)


def test_parse_dat_totals_numbers_increments_consistently_across_nsets(tmp_path):
    text = _dat(0.5, 1.0, nset="FIXED_END") + _dat(0.5, 1.0, nset="FOOT_POCKET")
    reactions = parse_dat_totals(_write(tmp_path, "job.dat", text))
    at_end = reactions[reactions["time"] == 1.0]
    assert len(at_end) == 2
    assert at_end["increment"].nunique() == 1  # same increment for both sets


def test_parse_dat_totals_rejects_an_empty_dat(tmp_path):
    with pytest.raises(SolveError, match="no reaction totals"):
        parse_dat_totals(_write(tmp_path, "job.dat", ""))


# --- the convergence guard, pinned with a stub solver ---------------------
# The real ccx exits nonzero when it stops early, so these paths are otherwise
# unreachable and the .sta check could be deleted with the suite still green.


def test_exit_zero_with_a_short_sta_is_rejected(tmp_path):
    """Condition 2 standing alone: clean exit, converged-looking .dat, short step."""
    stub = _stub_ccx(tmp_path, sta=_sta((1, "1", 0.5)), dat=_dat(0.5))
    with pytest.raises(NotConverged, match="stopped at total time 0.5, not 1"):
        solve(STRIP, tmp_path / "short_sta", ccx=stub)


def test_exit_zero_with_an_unconverged_attempt_at_the_final_time_is_rejected(tmp_path):
    """The last attempt sits at 1.0 but never converged, so it does not count."""
    stub = _stub_ccx(
        tmp_path, sta=_sta((1, "1", 0.5), (2, "1U", 1.0)), dat=_dat(0.5, 1.0)
    )
    with pytest.raises(NotConverged, match="stopped at total time 0.5"):
        solve(STRIP, tmp_path / "attempt", ccx=stub)


def test_reaching_the_final_time_without_printing_reactions_there_is_rejected(tmp_path):
    stub = _stub_ccx(tmp_path, sta=_sta((1, "1", 0.5), (2, "1", 1.0)), dat=_dat(0.5))
    with pytest.raises(NotConverged, match="printed no reactions"):
        solve(STRIP, tmp_path / "unprinted", ccx=stub)


def test_stale_outputs_are_cleared_before_the_solve(tmp_path):
    """A previous job's converged .sta and .dat must not be read as this run's."""
    run_dir = tmp_path / "stale"
    run_dir.mkdir()
    (run_dir / "job.sta").write_text(SAMPLE_STA)
    (run_dir / "job.dat").write_text(SAMPLE_DAT)
    stub = _stub_ccx(tmp_path)  # exits 0 having written nothing
    with pytest.raises(NotConverged, match="never ran"):
        solve(STRIP, run_dir, ccx=stub)


def test_a_deck_that_ends_at_another_total_time_can_say_so(tmp_path):
    """ccx accumulates TOT TIME across steps; 1.0 is a default, not a law."""
    stub = _stub_ccx(
        tmp_path, sta=_sta((1, "1", 1.0), (2, "1", 2.0)), dat=_dat(1.0, 2.0)
    )
    with pytest.raises(NotConverged):
        solve(STRIP, tmp_path / "two_step_default", ccx=stub)
    result = solve(STRIP, tmp_path / "two_step", ccx=stub, final_time=2.0)
    assert result.final["time"].iloc[0] == pytest.approx(2.0)
    assert result.total_force("fixed_end") == pytest.approx(-12.52544)


def test_a_nonzero_exit_is_rejected_even_with_a_complete_sta(tmp_path):
    """Condition 1 standing alone."""
    stub = _stub_ccx(tmp_path, sta=_sta((1, "1", 1.0)), dat=_dat(1.0), exit_code=201)
    with pytest.raises(CcxFailed, match="exited 201"):
        solve(STRIP, tmp_path / "crash", ccx=stub)


# --- the real solver ------------------------------------------------------


def test_solve_returns_the_converged_reaction(tmp_path):
    result = solve(STRIP, tmp_path / "converged")
    assert result.increments == STRIP_INCREMENTS
    assert result.total_force("fixed_end", "fz") == pytest.approx(STRIP_FZ_N, rel=0.01)
    assert result.wall_time_s > 0


def test_solve_keeps_the_whole_load_history(tmp_path):
    result = solve(STRIP, tmp_path / "history")
    assert result.reactions["time"].nunique() == STRIP_INCREMENTS
    assert result.reactions["time"].is_monotonic_increasing
    assert len(result.final) == 1


def test_the_response_stiffens_under_nlgeom(tmp_path):
    """NLGEOM is actually doing something: the strip is not linear at 20 mm.

    Monotonic |force| would hold for a linear response too, so compare the
    final secant stiffness against the closed-form small-deflection value.
    """
    result = solve(STRIP, tmp_path / "nlgeom")
    secant = abs(result.total_force("fixed_end")) / STRIP_DRIVEN_MM
    assert secant > STRIP_LINEAR_STIFFNESS_N_PER_MM * 1.02


def test_solve_raises_on_a_solve_that_stopped_early(tmp_path):
    # Running out of increments is also a nonzero exit, so this surfaces as
    # CcxFailed -- but the message has to say how far the step actually got.
    with pytest.raises(CcxFailed, match="total time 0.115"):
        solve(STRIP_SHORT, tmp_path / "short")


def test_the_partial_reaction_is_never_returned(tmp_path):
    """The truncated deck's .dat holds a plausible number. It must not escape."""
    run_dir = tmp_path / "short"
    with pytest.raises(SolveError):
        solve(STRIP_SHORT, run_dir)
    partial = parse_dat_totals(run_dir / "job.dat")
    assert partial["fz"].iloc[-1] == pytest.approx(SHORT_PARTIAL_FZ_N, rel=1e-4)


def test_solve_raises_when_ccx_exits_zero_without_analysing(tmp_path):
    """ccx skips keywords it does not know and reports a deck with no *STEP as
    "Job finished", exit 0. A clean exit status on its own proves nothing."""
    broken = _write(tmp_path, "broken.inp", "*NODE\n1, 0., 0., 0.\n*NONSENSE\n")
    with pytest.raises(NotConverged, match="never ran"):
        solve(broken, tmp_path / "broken")


def test_solve_is_single_threaded(tmp_path, monkeypatch):
    # ccx already defaults to one thread, so this only tests anything if the
    # ambient environment asks for more.
    monkeypatch.setenv("OMP_NUM_THREADS", "8")
    result = solve(STRIP, tmp_path / "threads")
    assert os.environ["OMP_NUM_THREADS"] == "8"
    stdout = (result.job_dir / f"{result.job_name}.stdout").read_text()
    assert "Using up to 1 cpu(s)" in stdout


def test_total_force_rejects_an_unprinted_nset(tmp_path):
    result = solve(STRIP, tmp_path / "nset")
    with pytest.raises(SolveError, match="not printed"):
        result.total_force("foot_pocket")


def test_solve_needs_a_deck(tmp_path):
    with pytest.raises(FileNotFoundError):
        solve(tmp_path / "absent.inp", tmp_path / "run")
