#!/usr/bin/env python3
"""Sweep ``TaperDesign`` layups through the FIN_TEST_3 tip-weight bench.

One design -> a ply book (``compfea.plybook_gen``) -> one ``run_tipweight.py``
solve, unchanged -> ``flex_kick.json``. Every run is collected into
``results/<run_id>/results.parquet`` with the design's parameters flattened
alongside its flex and kick numbers.

Input is a JSON file, one of:

- a bare list ``[{design}, {design}, ...]``;
- a single ``{design}`` dict;
- ``{"base": {design}, "axes": {"core": [...], "pads.HEAL": [...], ...}}`` --
  the Cartesian product via ``plybook_gen.expand_grid``.

Each solve is single-core (``OMP_NUM_THREADS=1``, ``--threads 1``); parallelism
is N solves at once, ``--jobs`` of them (default ``cpu_count() - 2``), per
CLAUDE.md. Designs are cached on the resolved-laminate fingerprint, so a rerun
that only adds grid points is cheap.

Like ``compfea.sweep`` this detaches itself: launch it, it re-execs under
``nohup``, prints the run id, and writes ``status.json`` / ``sweep.log`` /
``results.parquet`` under ``results/<run_id>/``. Poll those. ``--sync`` stays in
the foreground (tests, short grids).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from compfea.plybook_gen import (
    TaperDesign,
    design_fingerprint,
    design_from_dict,
    design_rows,
    design_to_dict,
    expand_grid,
    write_plybook,
)
from compfea.shop_inventory import load_shop_inventory
from compfea.step_mesh import _sanitize_elset
from compfea.sweep import write_status

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).resolve()
RUN_TIPWEIGHT = SCRIPT.parent / "run_tipweight.py"
WORKER_ENV = "COMPFEA_LAYUP_SWEEP_WORKER"

DEFAULT_MATERIALS = ROOT / "materials" / "from_shop.csv"
DEFAULT_INVENTORY = ROOT / "materials" / "shop_inventory.csv"


# --------------------------------------------------------------------------
# design list


def load_designs(path: str | Path) -> list[TaperDesign]:
    """Read the designs JSON into a list of validated ``TaperDesign``."""
    obj = json.loads(Path(path).read_text())
    if isinstance(obj, list):
        return [design_from_dict(d) for d in obj]
    if isinstance(obj, dict) and "base" in obj:
        return expand_grid(obj["base"], obj.get("axes", {}))
    if isinstance(obj, dict):
        return [design_from_dict(obj)]
    raise SystemExit(
        "designs JSON must be a list, a single design dict, or {base, axes}"
    )


def _spec_label(specs: list) -> str:
    return "/".join(f"{material}@{angle:g}" for material, angle in specs)


def _flatten_design(design: TaperDesign) -> dict:
    """Design identity columns for the results table."""
    dd = design_to_dict(design)
    row = {
        "through_zone": dd["through_zone"],
        "skins": _spec_label(dd["skins"]),
        "core": _spec_label(dd["core"]),
        "core_center": (
            _spec_label([dd["core_center"]]) if dd["core_center"] else ""
        ),
        "pad_order": ",".join(dd["pad_order"]),
    }
    for zone, half in dd["pads"].items():
        row[f"pad_{zone}"] = _spec_label(half)
    return row


# --------------------------------------------------------------------------
# one solve


def _run_tipweight_solve(
    plybook_path: Path,
    run_dir: Path,
    *,
    materials: Path,
    size_mm: float,
    uy_frac: float,
    bow_tip_mm: float,
    kick_bands: int,
    timeout_s: float,
) -> None:
    """Invoke ``run_tipweight.py`` unchanged for one ply book. Raises on failure.

    Factored out so a test can substitute a stub that writes a canned
    ``flex_kick.json`` -- the sweep's dispatch / cache / collect logic is unit
    tested without a two-minute solve.
    """
    cmd = [
        sys.executable,
        str(RUN_TIPWEIGHT),
        "--plybook", str(plybook_path),
        "--materials", str(materials),
        "--size-mm", str(size_mm),
        "--uy-frac", str(uy_frac),
        "--bow-tip-mm", str(bow_tip_mm),
        "--kick-bands", str(kick_bands),
        "--threads", "1",
        "--timeout-s", str(timeout_s),
        "--run-dir", str(run_dir),
    ]
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    with (run_dir / "run.log").open("w") as log:
        subprocess.run(
            cmd,
            check=True,
            cwd=str(ROOT),
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=timeout_s + 120.0,
        )


def solve_phash(
    *,
    materials: Path,
    size_mm: float,
    uy_frac: float,
    bow_tip_mm: float,
    kick_bands: int,
) -> str:
    """8-hex hash of everything that changes the *solve* but not the laminate.

    ``design_fingerprint`` identifies the laminate; this identifies the bench
    run around it -- the drive fraction, the seed bow, the band count (which is
    baked into the ``*NODE PRINT`` stations, so a cached ``flex_kick.json`` from
    a different count is stale, not re-analysable), the mesh size, and the bytes
    of the materials file and the STEP. The cache dir is keyed on both, so
    ``--kick-bands 41`` then ``--kick-bands 81`` under one ``--run-id`` misses
    rather than serving the 41-band result.

    Everything else reaching ``run_tipweight`` is its default and constant across
    a sweep (``--clip-n``, ``--gravity``, ``--static-line``). If a flag for one
    of those is ever added here, add it to this hash too.
    """
    h = hashlib.sha256()
    h.update(f"{size_mm}|{uy_frac}|{bow_tip_mm}|{kick_bands}".encode())
    for path in (Path(materials), ROOT / "FIN_TEST_3.step"):
        h.update(path.read_bytes() if path.is_file() else b"<missing>")
    return h.hexdigest()[:8]


def _through_rows(design: TaperDesign, skus) -> tuple[list, list]:
    through = _sanitize_elset(design.through_zone)
    rows = list(design_rows(design, skus=skus))
    return rows, [r for r in rows if r.zone == through]


def solve_one(
    design: TaperDesign,
    *,
    cache_dir: Path,
    skus,
    materials: Path,
    size_mm: float,
    uy_frac: float,
    bow_tip_mm: float,
    kick_bands: int,
    timeout_s: float,
) -> dict:
    """Resolve, solve (or hit cache), read ``flex_kick.json`` -> a result row."""
    fp = design_fingerprint(design, skus=skus)
    phash = solve_phash(
        materials=materials, size_mm=size_mm, uy_frac=uy_frac,
        bow_tip_mm=bow_tip_mm, kick_bands=kick_bands,
    )
    run_dir = cache_dir / f"{fp}.{phash}"
    run_dir.mkdir(parents=True, exist_ok=True)
    fk_path = run_dir / "flex_kick.json"

    row: dict = {"fingerprint": fp, **_flatten_design(design)}
    row["cached"] = fk_path.is_file()
    if not fk_path.is_file():
        write_plybook(
            design, run_dir / "plybook.csv", skus=skus, note=f"sweep design {fp}"
        )
        _run_tipweight_solve(
            run_dir / "plybook.csv",
            run_dir,
            materials=materials,
            size_mm=size_mm,
            uy_frac=uy_frac,
            bow_tip_mm=bow_tip_mm,
            kick_bands=kick_bands,
            timeout_s=timeout_s,
        )
    if not fk_path.is_file():
        raise RuntimeError(f"{fk_path} not written; see {run_dir / 'run.log'}")

    fk = json.loads(fk_path.read_text())
    head = fk.get("headline") or {}
    all_rows, through_rows = _through_rows(design, skus)
    row.update(
        status="ok",
        flex_n=fk.get("flex_n"),
        free_len_mm=fk.get("free_len_mm"),
        kick_bucket=fk.get("kick_bucket"),
        kick_s_frac=fk.get("kick_s_frac"),
        kick_s_mm=fk.get("kick_s_mm"),
        migration_delta_mm=fk.get("migration_delta_mm"),
        plateau_frac=head.get("plateau_frac"),
        warning=head.get("warning"),
        n_plies_total=len(all_rows),
        through_n_plies=len(through_rows),
        through_thickness_mm=sum(r.thickness_mm for r in through_rows),
    )
    return row


def _worker(task: tuple[TaperDesign, dict]) -> dict:
    design, kw = task
    try:
        return solve_one(design, **kw)
    except Exception as exc:  # noqa: BLE001 -- a failed design is a row, not a crash
        return {
            "fingerprint": design_fingerprint(design, skus=kw["skus"]),
            **_flatten_design(design),
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }


# --------------------------------------------------------------------------
# the sweep


def run_sweep(
    designs: list[TaperDesign],
    *,
    run_root: Path,
    jobs: int,
    materials: Path,
    inventory: Path,
    size_mm: float,
    uy_frac: float,
    bow_tip_mm: float,
    kick_bands: int,
    timeout_s: float,
) -> pd.DataFrame:
    run_root.mkdir(parents=True, exist_ok=True)
    cache_dir = run_root / "cache"
    cache_dir.mkdir(exist_ok=True)
    log_path = run_root / "sweep.log"
    skus = load_shop_inventory(inventory)

    kw = {
        "cache_dir": cache_dir,
        "skus": skus,
        "materials": materials,
        "size_mm": size_mm,
        "uy_frac": uy_frac,
        "bow_tip_mm": bow_tip_mm,
        "kick_bands": kick_bands,
        "timeout_s": timeout_s,
    }
    # Collapse designs that resolve to the same laminate: they share a cache
    # dir, so solving both would race two ccx runs into one directory (job
    # scratch files collide). One solve, one row per distinct laminate.
    unique: dict[str, TaperDesign] = {}
    for d in designs:
        unique.setdefault(design_fingerprint(d, skus=skus), d)
    collapsed = len(designs) - len(unique)
    tasks = [(d, kw) for d in unique.values()]

    write_status(
        run_root, state="running", n_designs=len(designs),
        n_unique=len(unique), jobs=jobs, done=0, ok=0, error=0,
    )
    rows: list[dict] = []
    with log_path.open("a") as log:
        log.write(
            f"sweep start designs={len(designs)} unique={len(unique)} "
            f"collapsed={collapsed} jobs={jobs}\n"
        )
        log.flush()
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(_worker, t): t for t in tasks}
            done = ok = err = 0
            for fut in as_completed(futures):
                row = fut.result()
                rows.append(row)
                done += 1
                if row.get("status") == "ok":
                    ok += 1
                    log.write(
                        f"ok {row['fingerprint']} flex_n={row.get('flex_n')} "
                        f"kick={row.get('kick_bucket')} "
                        f"f={row.get('kick_s_frac')}"
                        f"{' (cached)' if row.get('cached') else ''}\n"
                    )
                else:
                    err += 1
                    log.write(
                        f"error {row['fingerprint']}: {row.get('error')}\n"
                    )
                log.flush()
                write_status(
                    run_root, state="running", done=done, ok=ok, error=err
                )

    frame = pd.DataFrame(rows)
    if "flex_n" in frame.columns:
        frame = frame.sort_values(
            "flex_n", na_position="last"
        ).reset_index(drop=True)
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
        results_parquet=str(out_parquet),
        results_csv=str(out_csv),
    )
    return frame


# --------------------------------------------------------------------------
# CLI + detach


def _run_id() -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    return f"layup-sweep-{stamp}"


def _detach_and_exit(argv: list[str], run_root: Path) -> int:
    run_root.mkdir(parents=True, exist_ok=True)
    log_path = run_root / "sweep.log"
    cmd = ["nohup", sys.executable, str(SCRIPT), "--sync", *argv]
    env = {**os.environ, WORKER_ENV: "1"}
    write_status(run_root, state="starting", run_id=run_root.name, pid=None)
    with log_path.open("a") as log:
        log.write(f"detach: {' '.join(cmd)}\n")
        log.flush()
        proc = subprocess.Popen(
            cmd,
            # Inherit the caller's cwd, not ROOT: a relative designs /
            # --materials / --inventory path the parent resolved must resolve
            # the same way in the detached worker.
            cwd=os.getcwd(),
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
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "designs", type=Path, help="designs JSON (list | dict | {base, axes})"
    )
    p.add_argument("--materials", type=Path, default=DEFAULT_MATERIALS)
    p.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    p.add_argument("--size-mm", type=float, default=16.0)
    p.add_argument("--uy-frac", type=float, default=0.60)
    p.add_argument("--bow-tip-mm", type=float, default=0.5)
    p.add_argument("--kick-bands", type=int, default=41)
    p.add_argument(
        "--jobs", type=int, default=None, help="default cpu_count() - 2"
    )
    p.add_argument(
        "--timeout-s", type=float, default=3600.0, help="per-solve timeout"
    )
    p.add_argument("--run-id", default=None)
    p.add_argument(
        "--sync", action="store_true", help="foreground (default detaches)"
    )
    return p


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    jobs = args.jobs if args.jobs is not None else max(1, (os.cpu_count() or 2) - 2)

    # Parse the designs before anything detaches, so a broken file fails here.
    designs = load_designs(args.designs)
    if not designs:
        raise SystemExit(f"{args.designs}: no designs")

    run_id = args.run_id or _run_id()
    run_root = ROOT / "results" / run_id

    if not args.sync and os.environ.get(WORKER_ENV) != "1":
        forward = [a for a in argv if a != "--sync"]
        if "--run-id" not in forward:
            forward = ["--run-id", run_id, *forward]
        return _detach_and_exit(forward, run_root)

    try:
        frame = run_sweep(
            designs,
            run_root=run_root,
            jobs=jobs,
            materials=args.materials,
            inventory=args.inventory,
            size_mm=args.size_mm,
            uy_frac=args.uy_frac,
            bow_tip_mm=args.bow_tip_mm,
            kick_bands=args.kick_bands,
            timeout_s=args.timeout_s,
        )
    except BaseException as exc:
        # A crash in the (usually detached) worker must not leave status.json
        # stuck on "detached" -- a poller would wait forever.
        write_status(
            run_root, state="error",
            error=f"{type(exc).__name__}: {exc}",
        )
        raise
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(frame.to_string(index=False))
    print(f"\nwrote {run_root / 'results.parquet'}")
    return 0 if len(frame) and (frame["status"] == "ok").all() else 1


if __name__ == "__main__":
    raise SystemExit(main())
