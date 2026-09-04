"""Parameter sweep: strip U-bend ply count -> F_90 / F_180.

Default design: symmetric cross-ply stacks ``[0/90]_ns`` with fixed ply
thickness. Metric is tip-normal spring force from ELSE energy
(``compfea.metrics.f_spring``), not tip |RF|.

Long runs detach (``setsid nohup``) into ``results/<run_id>/`` and return the
run id immediately. Use ``--sync`` to stay in the foreground (tests, short
demos).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from compfea.geometry import Outline, mesh_outline
from compfea.layup import PLACEHOLDER_CFRP, Layup, Ply
from compfea.run import parse_dat_energy, solve
from compfea.ubend import (
    DEFAULT_END_DEG,
    DEFAULT_START_DEG,
    DEFAULT_STEP_DEG,
    build_deck,
    final_time_for,
    force_at_theta,
    step_index_for,
    tip_length_mm,
    theta_grid_deg,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = REPO_ROOT / "results"

# Strip that matched the validated U-bend path (32 along, 2 across).
LENGTH_MM = 100.0
WIDTH_MM = 20.0
N_SPAN = 32
N_CHORD = 2
LONG_AXIS = "y"
PLY_MM_DEFAULT = 0.1
REPORT_DEG = (90.0, 180.0)


@dataclass(frozen=True)
class Design:
    """One sweep point: n repeats of [0/90] before mirroring -> [0/90]_ns."""

    n_pairs: int
    ply_mm: float = PLY_MM_DEFAULT

    @property
    def n_plies(self) -> int:
        return 4 * self.n_pairs

    @property
    def stack_label(self) -> str:
        if self.n_pairs == 1:
            return "[0/90]s"
        return f"[0/90]_{self.n_pairs}s"

    def cache_key(self) -> str:
        payload = (
            f"ubend|L={LENGTH_MM}|W={WIDTH_MM}|ns={N_SPAN}|nc={N_CHORD}|"
            f"ply={self.ply_mm}|npairs={self.n_pairs}|step={DEFAULT_STEP_DEG}|"
            f"mat=placeholder_cfrp"
        )
        return hashlib.sha256(payload.encode()).hexdigest()[:16]


def layup_for(design: Design) -> Layup:
    half: list[Ply] = []
    for _ in range(design.n_pairs):
        half.append(Ply(design.ply_mm, 0.0))
        half.append(Ply(design.ply_mm, 90.0))
    plies = half + list(reversed(half))
    return Layup.uniform(
        plies,
        long_axis=LONG_AXIS,
        elset="blade",
        material=PLACEHOLDER_CFRP,
    )


def strip_mesh():
    return mesh_outline(
        Outline.rectangle(chord=WIDTH_MM, span=LENGTH_MM),
        n_chord=N_CHORD,
        n_span=N_SPAN,
        heading=(
            f"sweep U-bend strip {LENGTH_MM:g}x{WIDTH_MM:g} mm, "
            f"{N_SPAN}x{N_CHORD} S8R; y = long axis"
        ),
    )


def _default_jobs() -> int:
    return max(1, (os.cpu_count() or 2) - 2)


def _run_id() -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"ubend-plycount-{stamp}"


def evaluate_design(
    design: Design,
    *,
    run_dir: Path,
    timeout_s: float = 1800.0,
    report_deg: tuple[float, ...] = REPORT_DEG,
) -> dict:
    """Solve one design; return a result row dict (raises on solve failure)."""
    run_dir.mkdir(parents=True, exist_ok=True)
    cache_path = run_dir / "result.json"
    if cache_path.is_file():
        return json.loads(cache_path.read_text())

    mesh = strip_mesh()
    arm = tip_length_mm(mesh)
    angles = theta_grid_deg(
        step_deg=DEFAULT_STEP_DEG,
        start_deg=DEFAULT_START_DEG,
        end_deg=max(report_deg),
    )
    deck_text = build_deck(mesh, layup_for(design), angles)
    deck_path = run_dir / "deck.inp"
    deck_path.write_text(deck_text)

    started = time.monotonic()
    result = solve(
        deck_path,
        run_dir / "ccx",
        job_name="job",
        timeout_s=timeout_s,
        final_time=final_time_for(angles),
    )
    energy = parse_dat_energy(run_dir / "ccx" / "job.dat")

    row: dict = {
        "n_pairs": design.n_pairs,
        "n_plies": design.n_plies,
        "ply_mm": design.ply_mm,
        "stack": design.stack_label,
        "thickness_mm": design.n_plies * design.ply_mm,
        "cache_key": design.cache_key(),
        "arm_mm": arm,
        "wall_time_s": result.wall_time_s,
        "increments": result.increments,
        "status": "ok",
        "elapsed_s": time.monotonic() - started,
    }
    for deg in report_deg:
        idx = step_index_for(angles, deg)
        u, m, f = force_at_theta(
            energy, theta_deg=deg, step_index=idx, arm_mm=arm
        )
        tag = str(int(deg)) if float(deg).is_integer() else f"{deg:g}"
        row[f"u_{tag}"] = u
        row[f"m_{tag}"] = m
        row[f"f_{tag}"] = f
    if 90.0 in report_deg and 180.0 in report_deg:
        row["f_ratio_180_90"] = row["f_180"] / row["f_90"]
    cache_path.write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
    return row


def _worker(payload: dict) -> dict:
    design = Design(n_pairs=payload["n_pairs"], ply_mm=payload["ply_mm"])
    try:
        return evaluate_design(
            design,
            run_dir=Path(payload["run_dir"]),
            timeout_s=payload["timeout_s"],
        )
    except Exception as exc:  # noqa: BLE001 — sweep records failures as rows
        return {
            "n_pairs": design.n_pairs,
            "n_plies": design.n_plies,
            "ply_mm": design.ply_mm,
            "stack": design.stack_label,
            "thickness_mm": design.n_plies * design.ply_mm,
            "cache_key": design.cache_key(),
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }


def write_status(run_root: Path, **fields) -> None:
    path = run_root / "status.json"
    current = {}
    if path.is_file():
        current = json.loads(path.read_text())
    current.update(fields)
    current["updated_at"] = datetime.now(timezone.utc).isoformat()
    path.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n")


def run_sweep(
    designs: list[Design],
    *,
    run_root: Path,
    jobs: int,
    timeout_s: float,
) -> pd.DataFrame:
    run_root.mkdir(parents=True, exist_ok=True)
    log_path = run_root / "sweep.log"
    cache_dir = run_root / "cache"
    cache_dir.mkdir(exist_ok=True)

    write_status(
        run_root,
        state="running",
        n_designs=len(designs),
        jobs=jobs,
        done=0,
        ok=0,
        error=0,
    )
    payloads = []
    for design in designs:
        payloads.append(
            {
                "n_pairs": design.n_pairs,
                "ply_mm": design.ply_mm,
                "run_dir": str(cache_dir / design.cache_key()),
                "timeout_s": timeout_s,
            }
        )

    rows: list[dict] = []
    with log_path.open("a") as log:
        log.write(f"sweep start designs={len(designs)} jobs={jobs}\n")
        log.flush()
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(_worker, p): p for p in payloads}
            done = ok = err = 0
            for fut in as_completed(futures):
                row = fut.result()
                rows.append(row)
                done += 1
                if row.get("status") == "ok":
                    ok += 1
                    log.write(
                        f"ok {row['stack']} f_90={row.get('f_90', float('nan')):.4g} "
                        f"f_180={row.get('f_180', float('nan')):.4g}\n"
                    )
                else:
                    err += 1
                    log.write(f"error {row.get('stack')}: {row.get('error')}\n")
                log.flush()
                write_status(run_root, state="running", done=done, ok=ok, error=err)

    frame = pd.DataFrame(rows).sort_values(["ply_mm", "n_pairs"]).reset_index(drop=True)
    out_parquet = run_root / "results.parquet"
    out_csv = run_root / "results.csv"
    frame.to_parquet(out_parquet, index=False)
    frame.to_csv(out_csv, index=False)
    write_status(
        run_root,
        state="done",
        done=len(rows),
        ok=int((frame["status"] == "ok").sum()) if len(frame) else 0,
        error=int((frame["status"] != "ok").sum()) if len(frame) else 0,
        results_csv=str(out_csv),
        results_parquet=str(out_parquet),
    )
    return frame


def _detach_and_exit(argv: list[str], run_root: Path) -> int:
    """Re-exec under setsid/nohup; print run id and exit."""
    run_root.mkdir(parents=True, exist_ok=True)
    log_path = run_root / "sweep.log"
    worker_argv = [sys.executable, "-m", "compfea.sweep", "--sync", *argv]
    env = {**os.environ, "COMPFEA_SWEEP_WORKER": "1"}
    # Linux has setsid(1); macOS does not. Popen(start_new_session=True)
    # already starts a new process group on both, so nohup alone is enough.
    cmd = ["nohup", *worker_argv]
    write_status(run_root, state="starting", run_id=run_root.name, pid=None)
    with log_path.open("a") as log:
        log.write(f"detach: {' '.join(cmd)}\n")
        log.flush()
        proc = __import__("subprocess").Popen(
            cmd,
            cwd=str(REPO_ROOT),
            env=env,
            stdout=log,
            stderr=log,
            start_new_session=True,
        )
    write_status(run_root, state="detached", pid=proc.pid, run_id=run_root.name)
    print(run_root.name)
    print(f"status: {run_root / 'status.json'}")
    print(f"log:    {run_root / 'sweep.log'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="compfea-sweep",
        description=(
            "Sweep [0/90]_ns ply count on the strip U-bend; report F_90 and "
            "F_180 from ELSE energy (F=M/L)."
        ),
    )
    p.add_argument(
        "--n-pairs",
        type=int,
        nargs="+",
        default=[1, 2, 3],
        help="repeats of [0/90] before mirror (1 -> 4 plies). Default: 1 2 3",
    )
    p.add_argument("--ply-mm", type=float, default=PLY_MM_DEFAULT)
    p.add_argument("--jobs", type=int, default=None, help="default cpu_count()-2")
    p.add_argument("--timeout-s", type=float, default=1800.0)
    p.add_argument("--run-id", default=None)
    p.add_argument(
        "--sync",
        action="store_true",
        help="run in the foreground (default detaches with setsid nohup)",
    )
    p.add_argument(
        "--deck-only",
        action="store_true",
        help="write decks for each design and exit without solving",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    jobs = args.jobs if args.jobs is not None else _default_jobs()
    designs = [Design(n_pairs=n, ply_mm=args.ply_mm) for n in args.n_pairs]
    for d in designs:
        if d.n_pairs < 1:
            raise SystemExit("n-pairs must be >= 1")

    run_id = args.run_id or _run_id()
    run_root = RESULTS_ROOT / run_id

    if args.deck_only:
        run_root.mkdir(parents=True, exist_ok=True)
        mesh = strip_mesh()
        angles = theta_grid_deg(end_deg=DEFAULT_END_DEG)
        for design in designs:
            out = run_root / f"deck_{design.cache_key()}.inp"
            out.write_text(build_deck(mesh, layup_for(design), angles))
            print(f"wrote {out}")
        return 0

    # Detach unless --sync or we are already the worker.
    if not args.sync and os.environ.get("COMPFEA_SWEEP_WORKER") != "1":
        # Strip a possible existing --run-id; pin the id we created.
        forward = [a for a in argv if a != "--sync"]
        if "--run-id" not in forward:
            forward = ["--run-id", run_id, *forward]
        return _detach_and_exit(forward, run_root)

    frame = run_sweep(designs, run_root=run_root, jobs=jobs, timeout_s=args.timeout_s)
    print(frame.to_string(index=False))
    print(f"\nwrote {run_root / 'results.csv'}")
    return 0 if (frame["status"] == "ok").all() else 1


if __name__ == "__main__":
    raise SystemExit(main())
