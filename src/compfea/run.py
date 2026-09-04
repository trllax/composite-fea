"""Run a CalculiX job and refuse to report anything but a converged result.

Units are CalculiX-consistent (mm, N, tonne, MPa, s), so reaction forces
parsed here are newtons.

A solve that stops early still writes a well-formed ``.dat``. Its last block
is a clean number that looks exactly like an answer and is meaningless, so
nothing in this module parses reactions before both validity conditions hold:
``ccx`` exited zero, and the last increment in the ``.sta`` file reached total
time 1.0. Failure raises. There is deliberately no path that returns the last
available increment.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# ccx writes these beside the deck. A stale one left by an earlier solve parses
# just as cleanly as a fresh one, so they are removed before every run.
_SOLVER_OUTPUTS = (
    ".dat", ".sta", ".cvg", ".frd", ".12d", ".equ", ".rout", ".rin", ".rfn",
)

# The load-case template in CLAUDE.md drives one step over a period of 1.0, so
# that is the default. It is a property of the deck, not a law: ccx accumulates
# TOT TIME across steps, so a two-step deck ends at 2.0 and must say so.
_DEFAULT_FINAL_TIME = 1.0
_TIME_TOL = 1e-6

_STA_COLUMNS = ("step", "inc", "att", "itrs", "tot_time", "step_time", "inc_time")
_REACTION_COLUMNS = ("increment", "time", "nset", "fx", "fy", "fz")

# " total force (fx,fy,fz) for set FIXED_END and time  0.1000000E+01"
_TOTAL_HEADER = re.compile(
    r"^\s*total force \(fx,fy,fz\) for set (?P<nset>\S+) and time\s+(?P<time>\S+)\s*$"
)


class SolveError(RuntimeError):
    """Base for every way a solve fails to produce a usable result."""


class CcxFailed(SolveError):
    """ccx exited nonzero."""


class NotConverged(SolveError):
    """ccx exited zero but the step never reached total time 1.0."""


class SolveTimeout(SolveError):
    """ccx exceeded the wall-clock budget and was killed."""


@dataclass(frozen=True, eq=False)
class SolveResult:
    """A converged solve. Constructing one is a claim that it converged."""

    job_dir: Path
    job_name: str
    reactions: pd.DataFrame
    increments: int
    wall_time_s: float

    @property
    def final(self) -> pd.DataFrame:
        """Reaction totals at the end of the analysis, one row per printed nset.

        Keyed off the last time actually present rather than a constant: solve()
        has already established that this is the deck's final time.
        """
        if self.reactions.empty:
            raise SolveError(f"no reaction totals in {self.job_dir}")
        last = self.reactions["time"].max()
        return self.reactions[(self.reactions["time"] - last).abs() <= _TIME_TOL]

    def total_force(self, nset: str, component: str = "fz") -> float:
        """Reaction force component on ``nset`` at the end of the step, in N."""
        if component not in ("fx", "fy", "fz"):
            raise ValueError(f"component must be fx, fy or fz, not {component!r}")
        rows = self.final
        rows = rows[rows["nset"] == nset.lower()]
        if rows.empty:
            printed = sorted(self.reactions["nset"].unique())
            raise SolveError(
                f"nset {nset!r} not printed in {self.job_name}.dat; has {printed}"
            )
        if len(rows) > 1:
            raise SolveError(
                f"nset {nset!r} has {len(rows)} totals blocks at the final time in "
                f"{self.job_name}.dat; the deck prints it more than once"
            )
        return float(rows[component].iloc[0])


def parse_sta(path: str | Path) -> pd.DataFrame:
    """Parse a CalculiX ``.sta`` increment summary."""
    path = Path(path)
    if not path.is_file():
        raise SolveError(f"ccx wrote no status file at {path}")
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) != len(_STA_COLUMNS):
            continue  # header lines
        # CalculiX flags an attempt that did not converge with a U suffix in the
        # ATT column ("1U"). Those rows are attempts, not increments, and are
        # dropped on purpose: an unconverged attempt sitting at the final time
        # must never satisfy the convergence check. Do not "repair" this by
        # stripping the suffix -- that turns condition 2 into a rubber stamp.
        if any(field.upper().endswith("U") for field in fields[:4]):
            continue
        try:
            rows.append(
                [int(v) for v in fields[:4]] + [float(v) for v in fields[4:]]
            )
        except ValueError:
            continue
    if not rows:
        raise SolveError(f"no increments recorded in {path}")
    return pd.DataFrame(rows, columns=list(_STA_COLUMNS))


def parse_dat_totals(path: str | Path) -> pd.DataFrame:
    """Parse the ``*NODE PRINT, TOTALS=YES`` reaction blocks from a ``.dat``.

    One row per increment per node set. Set names are lower-cased; ccx prints
    them upper-case regardless of how the deck spells them.
    """
    path = Path(path)
    if not path.is_file():
        raise SolveError(f"ccx wrote no results file at {path}")
    lines = path.read_text().splitlines()
    rows = []
    for i, line in enumerate(lines):
        header = _TOTAL_HEADER.match(line)
        if header is None:
            continue
        fx, fy, fz = _totals_after(lines, i + 1, path)
        rows.append(
            {
                "time": float(header["time"]),
                "nset": header["nset"].lower(),
                "fx": fx,
                "fy": fy,
                "fz": fz,
            }
        )
    if not rows:
        raise SolveError(
            f"no reaction totals in {path}; the deck needs "
            f"*NODE PRINT, NSET=..., TOTALS=YES with RF"
        )
    frame = pd.DataFrame(rows)
    # Keyed on time so the number means the same thing for every node set, even
    # when print blocks have different frequencies. A per-set running count
    # would silently disagree between sets.
    frame.insert(0, "increment", frame["time"].rank(method="dense").astype(int))
    return frame[list(_REACTION_COLUMNS)]


def _totals_after(
    lines: list[str], start: int, path: Path
) -> tuple[float, float, float]:
    """The three force components on the first non-blank line after a header."""
    for line in lines[start : start + 3]:
        fields = line.split()
        if not fields:
            continue
        if len(fields) == 3:
            try:
                fx, fy, fz = (float(v) for v in fields)
            except ValueError:
                break
            return fx, fy, fz
        break
    raise SolveError(f"malformed total-force block at line {start} of {path}")


def solve(
    inp_path: str | Path,
    run_dir: str | Path,
    *,
    job_name: str = "job",
    ccx: str = "ccx",
    timeout_s: float = 1800.0,
    final_time: float = _DEFAULT_FINAL_TIME,
) -> SolveResult:
    """Solve ``inp_path`` in its own directory and return the converged result.

    The deck is copied into ``run_dir`` as ``<job_name>.inp`` and solved there,
    so concurrent sweep jobs cannot collide on ccx's job-derived scratch file
    names. Each solve is single-threaded on purpose: parallelism belongs at the
    sweep level, not inside one solve.

    ``final_time`` is the total time the analysis must reach to count as
    converged. It defaults to 1.0, matching the single-step load-case template,
    but ccx accumulates TOT TIME across steps -- a deck with a pretension step
    plus a driven step, or a ``*STATIC`` period other than 1.0, must pass its
    own value or it will be rejected while perfectly converged.

    Raises ``CcxFailed``, ``NotConverged`` or ``SolveTimeout`` rather than
    returning a number from a solve that did not finish.
    """
    inp_path = Path(inp_path)
    run_dir = Path(run_dir)
    if not inp_path.is_file():
        raise FileNotFoundError(f"no deck at {inp_path}")
    executable = shutil.which(ccx)
    if executable is None:
        raise SolveError(f"solver {ccx!r} is not on PATH")

    run_dir.mkdir(parents=True, exist_ok=True)
    deck = run_dir / f"{job_name}.inp"
    stdout_path = run_dir / f"{job_name}.stdout"
    for suffix in _SOLVER_OUTPUTS:
        (run_dir / f"{job_name}{suffix}").unlink(missing_ok=True)
    stdout_path.unlink(missing_ok=True)
    if inp_path.resolve() != deck.resolve():
        shutil.copyfile(inp_path, deck)

    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [executable, "-i", job_name],
            cwd=run_dir,
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        output = _as_text(exc.stdout) + _as_text(exc.stderr)
        stdout_path.write_text(output)
        raise SolveTimeout(
            f"ccx exceeded {timeout_s:g} s on {deck}\n{_tail(output)}"
        ) from exc
    wall_time_s = time.monotonic() - started
    output = completed.stdout + completed.stderr
    stdout_path.write_text(output)

    if completed.returncode != 0:
        raise CcxFailed(
            f"ccx exited {completed.returncode} on {deck}"
            f"{_progress(run_dir, job_name, final_time)}; "
            f"see {stdout_path}\n{_tail(output)}"
        )
    # A clean exit proves nothing: ccx skips keywords it does not recognise, and
    # a deck that ends up with no *STEP is reported as "Job finished", exit 0.
    try:
        status = parse_sta(run_dir / f"{job_name}.sta")
    except SolveError as exc:
        raise NotConverged(
            f"ccx exited zero on {deck} but recorded no increments; the analysis "
            f"never ran. See {stdout_path}\n{_tail(output)}"
        ) from exc
    reached = float(status["tot_time"].iloc[-1])
    if abs(reached - final_time) > _TIME_TOL:
        raise NotConverged(
            f"{deck} stopped at total time {reached:g}, not {final_time:g}; "
            f"its reactions are meaningless. See {stdout_path}\n{_tail(output)}"
        )

    reactions = parse_dat_totals(run_dir / f"{job_name}.dat")
    if not (reactions["time"] - final_time).abs().le(_TIME_TOL).any():
        raise NotConverged(
            f"{deck} reached total time {reached:g} but printed no reactions "
            f"there; the last totals block is at {reactions['time'].max():g}"
        )
    return SolveResult(
        job_dir=run_dir,
        job_name=job_name,
        reactions=reactions,
        increments=int(len(status)),
        wall_time_s=wall_time_s,
    )


def _progress(run_dir: Path, job_name: str, final_time: float) -> str:
    """How far the analysis got, for error messages. Never raises."""
    try:
        status = parse_sta(run_dir / f"{job_name}.sta")
    except SolveError:
        return " having recorded no converged increment"
    reached = float(status["tot_time"].iloc[-1])
    if abs(reached - final_time) <= _TIME_TOL:
        return f" after reaching total time {reached:g}"
    return f", stopped at total time {reached:g} of {final_time:g}"


def _as_text(raw: str | bytes | None) -> str:
    if raw is None:
        return ""
    return raw if isinstance(raw, str) else raw.decode(errors="replace")


def _tail(text: str, lines: int = 20) -> str:
    return "\n".join(text.splitlines()[-lines:])
