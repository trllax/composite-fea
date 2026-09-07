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
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import pandas as pd

from compfea.geometry import Outline, mesh_outline
from compfea.layup import (
    UD_CFRP_GENERIC,
    EngineeringConstants,
    Layup,
    Ply,
    ZoneLayup,
    canonical_angle,
    woven_from_ud,
)
from compfea.run import parse_dat_energy, solve
from compfea.ubend import (
    DEFAULT_END_DEG,
    DEFAULT_INC,
    DEFAULT_START_DEG,
    DEFAULT_STATIC_LINE,
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
# Geometry is deliberately fixed: this sweep varies the laminate, not the part.
LENGTH_MM = 100.0
WIDTH_MM = 20.0
N_SPAN = 32
N_CHORD = 2
LONG_AXIS = "y"
PLY_MM_DEFAULT = 0.1
REPORT_DEG = (90.0, 180.0)
DEFAULT_ANGLES = (0.0, 90.0)

# Fibre architectures. The woven card is derived from the UD one by membrane
# equivalence (see layup.woven_from_ud) -- no second set of invented constants.
# At equal total thickness the two differ in D, which is what a bending sweep
# is measuring.
FIBERS: dict[str, EngineeringConstants] = {
    "ud": UD_CFRP_GENERIC,
    "woven": woven_from_ud(UD_CFRP_GENERIC),
}


def _material_fingerprint(ec: EngineeringConstants) -> str:
    """Every number that reaches the deck, in the cache key.

    The old key hashed the literal string "mat=placeholder_cfrp", so editing a
    modulus reused the old answer under the old hash and the sweep reported a
    stale force for a material it never solved.
    """
    return (
        f"{ec.name}:{ec.e1!r},{ec.e2!r},{ec.e3!r},{ec.nu12!r},{ec.nu13!r},"
        f"{ec.nu23!r},{ec.g12!r},{ec.g13!r},{ec.g23!r},{ec.density!r}"
    )


def _angle_label(angles: Sequence[float]) -> str:
    return "/".join(f"{a:g}" for a in angles)


@dataclass(frozen=True)
class Design:
    """One sweep point.

    ``zone_pairs`` is repeats of the ``angles`` unit per spanwise zone, root
    first. Each zone is mirrored, so a zone with ``n`` repeats of ``(0, 90)``
    carries ``[0/90]_ns`` -- ``4n`` plies. ``zones`` are the span fractions
    between them and must number one fewer than ``zone_pairs``.
    """

    zone_pairs: tuple[int, ...] = (1,)
    ply_mm: float = PLY_MM_DEFAULT
    angles: tuple[float, ...] = DEFAULT_ANGLES
    fiber: str = "ud"
    zones: tuple[float, ...] = ()
    static_line: str = DEFAULT_STATIC_LINE
    report_deg: tuple[float, ...] = REPORT_DEG
    stress: bool = False

    def __post_init__(self) -> None:
        if not self.zone_pairs:
            raise ValueError("a design needs at least one zone")
        if any(n < 1 for n in self.zone_pairs):
            raise ValueError(f"every zone needs >= 1 pair, got {self.zone_pairs}")
        if len(self.zones) != len(self.zone_pairs) - 1:
            raise ValueError(
                f"{len(self.zone_pairs)} zones need {len(self.zone_pairs) - 1} "
                f"boundaries, got {len(self.zones)}: {self.zones}"
            )
        if not self.angles:
            raise ValueError("a design needs at least one ply angle")
        if not (self.ply_mm > 0.0 and math.isfinite(self.ply_mm)):
            raise ValueError(f"ply_mm must be finite and positive, got {self.ply_mm}")
        if self.fiber not in FIBERS:
            raise ValueError(
                f"fiber must be one of {sorted(FIBERS)}, not {self.fiber!r}"
            )
        if not self.report_deg:
            raise ValueError("report_deg needs at least one angle")
        # (3) A boundary that is not on an element row snaps to one, so the
        # requested fraction and the delivered ply drop differ. geometry.py only
        # objects when a zone catches nothing at all, so 0.47 and 0.50 would mint
        # two cache keys for one identical deck. Refuse the ambiguity.
        rows = 1.0 / N_SPAN
        for fraction in self.zones:
            if not 0.0 < fraction < 1.0:
                raise ValueError(
                    f"zone fractions must lie in (0, 1), got {fraction}"
                )
            if abs(fraction / rows - round(fraction / rows)) > 1e-9:
                raise ValueError(
                    f"zone fraction {fraction} is not a multiple of 1/{N_SPAN}, "
                    f"so it would snap to a different ply drop than requested; "
                    f"nearest valid: {round(fraction / rows) * rows:g}"
                )
        if len(self.static_line.split(",")) != 4:
            raise ValueError(
                f"static_line must be 'initial, period, min, max', got "
                f"{self.static_line!r}"
            )

    @property
    def material(self) -> EngineeringConstants:
        return FIBERS[self.fiber]

    @property
    def zone_elsets(self) -> tuple[str, ...]:
        """ELSETs the sections hang on, matching geometry.mesh_outline.

        A single-zone design uses the whole-blade ELSET; a zoned one uses the
        ``zone_N`` sets mesh_outline emits, root first.
        """
        if len(self.zone_pairs) == 1:
            return ("blade",)
        return tuple(f"zone_{i}" for i in range(1, len(self.zone_pairs) + 1))

    @property
    def n_plies_by_zone(self) -> tuple[int, ...]:
        return tuple(2 * len(self.angles) * n for n in self.zone_pairs)

    @property
    def thickness_by_zone_mm(self) -> tuple[float, ...]:
        return tuple(n * self.ply_mm for n in self.n_plies_by_zone)

    @property
    def stack_label(self) -> str:
        """Per-zone stacks, root first, joined by ``|``.

        Deliberately no single ``n_plies``: on a zoned design that number would
        have to silently pick one zone, which is exactly the kind of plausible
        wrong value this repo exists to avoid. Root and tip are reported as
        separate columns instead.
        """
        unit = _angle_label(self.angles)
        parts = [
            f"[{unit}]s" if n == 1 else f"[{unit}]_{n}s" for n in self.zone_pairs
        ]
        return "|".join(parts)

    def cache_key(self) -> str:
        payload = (
            f"ubend|L={LENGTH_MM}|W={WIDTH_MM}|ns={N_SPAN}|nc={N_CHORD}|"
            f"axis={LONG_AXIS}|ply={self.ply_mm!r}|"
            f"zone_pairs={self.zone_pairs!r}|zones={self.zones!r}|"
            f"angles={tuple(canonical_angle(a) for a in self.angles)!r}|"
            f"mat={_material_fingerprint(self.material)}|"
            f"step={DEFAULT_STEP_DEG}|start={DEFAULT_START_DEG}|"
            f"report={self.report_deg!r}|static={self.static_line}|inc={DEFAULT_INC}|"
            # Output requests belong in the key. Without this a run solved
            # without the stress card is served from cache to a request that
            # expects one, and the stress columns come back empty (or, once a
            # schema exists for them, stale) with nothing to say why.
            f"stress={self.stress!r}"
        )
        return hashlib.sha256(payload.encode()).hexdigest()[:16]


def layup_for(design: Design) -> Layup:
    """One symmetric stack per zone, root first. First ply is the -z ply."""
    material = design.material
    zone_layups = []
    for elset, n in zip(design.zone_elsets, design.zone_pairs, strict=True):
        half = [
            Ply(design.ply_mm, angle, material.name)
            for _ in range(n)
            for angle in design.angles
        ]
        zone_layups.append(ZoneLayup(elset, tuple(half + list(reversed(half)))))
    return Layup(
        materials=(material,),
        zones=tuple(zone_layups),
        long_axis=LONG_AXIS,
    )


def strip_mesh(zones: tuple[float, ...] = ()):
    """The fixed strip. ``zones`` are span fractions from the root.

    mesh_outline already refuses a zone thinner than one element row, so no
    second check here.
    """
    return mesh_outline(
        Outline.rectangle(chord=WIDTH_MM, span=LENGTH_MM),
        n_chord=N_CHORD,
        n_span=N_SPAN,
        zones=zones,
        heading=(
            f"sweep U-bend strip {LENGTH_MM:g}x{WIDTH_MM:g} mm, "
            f"{N_SPAN}x{N_CHORD} S8R; y = long axis"
        ),
    )


def angles_for(design: Design) -> list[float]:
    """The tip path this design is solved over."""
    return theta_grid_deg(
        step_deg=DEFAULT_STEP_DEG,
        start_deg=DEFAULT_START_DEG,
        end_deg=max(design.report_deg),
    )


def deck_for(design: Design) -> tuple[str, list[float], float]:
    """Deck text, angle path and moment arm for one design.

    One path, used by both ``--deck-only`` and the solve, so the thing you are
    shown is the thing that gets solved. They were separate and had already
    drifted on ``end_deg``.
    """
    mesh = strip_mesh(design.zones)
    angles = angles_for(design)
    deck = build_deck(
        mesh,
        layup_for(design),
        angles,
        long_axis=LONG_AXIS,
        static_line=design.static_line,
        stress_deg=design.report_deg if design.stress else None,
        # the un-rotation needs DISP wherever stress is printed
        file_deg=design.report_deg if design.stress else (90.0, 180.0),
    )
    return deck, angles, tip_length_mm(mesh, long_axis=LONG_AXIS)


def _default_jobs() -> int:
    return max(1, (os.cpu_count() or 2) - 2)


def _run_id() -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"ubend-plycount-{stamp}"


def _design_columns(design: Design) -> dict:
    """Design identity columns. Root and tip are separate on purpose.

    A zoned design has no single ply count, so there is no ``n_plies`` column
    to misread.
    """
    plies = design.n_plies_by_zone
    thick = design.thickness_by_zone_mm
    return {
        "stack": design.stack_label,
        "fiber": design.fiber,
        "ply_mm": design.ply_mm,
        "angles": _angle_label(design.angles),
        "n_zones": len(design.zone_pairs),
        "zone_bounds": ",".join(f"{z:g}" for z in design.zones),
        "zone_pairs": ",".join(str(n) for n in design.zone_pairs),
        "n_plies_root": plies[0],
        "n_plies_tip": plies[-1],
        "thickness_root_mm": thick[0],
        "thickness_tip_mm": thick[-1],
        "static_line": design.static_line,
        "cache_key": design.cache_key(),
    }


def _stress_columns(
    design: Design,
    run_dir: Path,
    angles: list[float],
    report_deg: Sequence[float],
) -> dict:
    """Peak per-ply stress at each reported angle, in the material frame.

    Imported here rather than at module scope for the same reason ``post`` is:
    a sweep should not fail to start because an optional analysis path has a
    bad import. Nothing here pulls in matplotlib.

    The rotation angle comes from the deformed shape in the ``.frd``, not from
    the stress tensor, so ``s33_residual`` stays an independent check rather
    than a tautology.
    """
    from compfea.frd import (
        deformed_midsurface,
        disp_at_step,
        index_disp_blocks,
        read_nodes,
    )
    from compfea.stress import (
        bend_angle_by_span,
        element_spans,
        material_frame,
        parse_dat_stress,
        stress_summary,
    )

    if len(design.zone_pairs) != 1:
        raise NotImplementedError(
            "stress output is single-zone only for now: a zoned design gives "
            "each zone a different ply count, so one integration-point index "
            "maps to different plies in different elements. Refusing rather "
            "than assigning plies that are plausibly wrong."
        )
    ply_angles = [p.angle_deg for p in layup_for(design).zones[0].plies]

    dat = run_dir / "ccx" / "job.dat"
    frd = run_dir / "ccx" / "job.frd"
    stress = parse_dat_stress(dat)
    spans = element_spans(run_dir / "deck.inp", long_axis=LONG_AXIS)
    blocks = index_disp_blocks(frd)
    nodes = read_nodes(frd)

    out: dict = {}
    for deg in report_deg:
        step = step_index_for(angles, float(deg))
        # No fallback. Each *STATIC step runs a unit period, so end-of-step N is
        # tot-time N; if that block is absent, summarising the newest block
        # instead would file 180-degree numbers under the 90-degree column names
        # with nothing on the row to say so.
        block = stress[(stress["time"] - float(step)).abs() <= 1e-6]
        if block.empty:
            have = ", ".join(f"{t:g}" for t in sorted(stress["time"].unique()))
            raise LookupError(
                f"no stress printed at step {step} (theta={deg:g}); "
                f"blocks exist at times: {have}"
            )
        _, disp = disp_at_step(frd, step, blocks=blocks)
        phi_at = bend_angle_by_span(
            deformed_midsurface(nodes, disp, long_axis=LONG_AXIS)
        )
        phi = {e: phi_at(s) for e, s in spans.items()}
        mat = material_frame(block, phi, ply_angles, long_axis=LONG_AXIS)
        tag = str(int(deg)) if float(deg).is_integer() else f"{deg:g}"
        for key, value in stress_summary(mat).items():
            out[f"{key}_{tag}"] = value
    return out


def evaluate_design(
    design: Design,
    *,
    run_dir: Path,
    timeout_s: float = 1800.0,
) -> dict:
    """Solve one design; return a result row dict (raises on solve failure).

    The reported angles come from the design, not from a parameter: they set
    how many ``*STEP`` blocks the deck carries, so a parameter here would
    change the deck without changing the cache key and a 90-degree run would be
    served from the cache to a 180-degree request.
    """
    report_deg = design.report_deg
    run_dir.mkdir(parents=True, exist_ok=True)
    cache_path = run_dir / "result.json"
    if cache_path.is_file():
        return json.loads(cache_path.read_text())

    deck_text, angles, arm = deck_for(design)
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
        **_design_columns(design),
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
    # F is the secant 2U/theta. Carry how well that held, so a design point
    # whose M(theta) actually curved is visible in the results table instead of
    # being ranked as if it were a spring.
    # Imported here, not at module scope: post pulls in matplotlib and seaborn,
    # and a sweep should not fail to start because the plotting stack does.
    from compfea.post import (
        force_curve,
        linearity_dev_near,
        linearity_dev_theta,
        max_linearity_dev,
    )

    blade = energy[energy["elset"] == "blade"]
    curve = force_curve(blade, angles, arm)
    row["max_linearity_dev"] = max_linearity_dev(curve)
    at_theta = linearity_dev_theta(curve, report_deg)
    for key, dev in linearity_dev_near(curve, report_deg).items():
        row[f"linearity_dev_{key.lower()}"] = dev
        # Which angle it came from: on a 5-degree path the value keyed f_180 is
        # taken at 175, and a column named for 180 would misreport that.
        row[f"linearity_dev_theta_{key.lower()}"] = at_theta.get(key)
    if design.stress:
        try:
            row.update(_stress_columns(design, run_dir, angles, report_deg))
        except Exception as exc:  # noqa: BLE001 - a stress failure is not a
            # solve failure: the forces are still good, so keep the row and say
            # why the stress columns are missing. Deliberately NOT cached: a
            # cached failure is served back forever and never retried, and the
            # run would keep reporting "ok" with no stress columns.
            row["stress_error"] = f"{type(exc).__name__}: {exc}"
            return row

    cache_path.write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
    return row


def _worker(payload: dict) -> dict:
    design: Design = payload["design"]
    try:
        return evaluate_design(
            design,
            run_dir=Path(payload["run_dir"]),
            timeout_s=payload["timeout_s"],
        )
    except Exception as exc:  # noqa: BLE001 — sweep records failures as rows
        return {
            **_design_columns(design),
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
                "design": design,
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
            done = ok = err = stress_failed = 0
            for fut in as_completed(futures):
                row = fut.result()
                rows.append(row)
                done += 1
                if row.get("status") == "ok":
                    ok += 1
                    note = ""
                    if row.get("stress_error"):
                        stress_failed += 1
                        note = f"  STRESS FAILED: {row['stress_error']}"
                    log.write(
                        f"ok {row['stack']} f_90={row.get('f_90', float('nan')):.4g} "
                        f"f_180={row.get('f_180', float('nan')):.4g}{note}\n"
                    )
                else:
                    err += 1
                    log.write(f"error {row.get('stack')}: {row.get('error')}\n")
                log.flush()
                write_status(
                    run_root, state="running", done=done, ok=ok, error=err,
                    stress_failed=stress_failed,
                )

    frame = (
        pd.DataFrame(rows)
        .sort_values(["fiber", "ply_mm", "thickness_root_mm", "stack"])
        .reset_index(drop=True)
    )
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


def _parse_int_list(text: str) -> tuple[int, ...]:
    try:
        return tuple(int(t) for t in text.split(",") if t.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"{text!r} must be comma-separated integers, e.g. 3,2,1"
        ) from exc


def _parse_angle_list(text: str) -> tuple[float, ...]:
    try:
        return tuple(float(t) for t in text.split("/") if t.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"{text!r} must be slash-separated angles, e.g. 0/90 or 45/-45"
        ) from exc


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="compfea-sweep",
        description=(
            "Sweep laminates on the fixed strip U-bend; report F_90 and F_180 "
            "from ELSE energy (F=M/L), not tip |RF|. The grid is the Cartesian "
            "product of --zone-pairs x --ply-mm x --angles x --fiber."
        ),
    )
    p.add_argument(
        "--zone-pairs",
        type=_parse_int_list,
        nargs="+",
        default=None,
        help=(
            "repeats of the angle unit per spanwise zone, root first, e.g. "
            "'3,2,1'. One entry per design. Default: 1 2 3 (single zone)"
        ),
    )
    p.add_argument(
        "--n-pairs",
        type=int,
        nargs="+",
        default=None,
        help="single-zone shorthand for --zone-pairs. Kept for older commands.",
    )
    p.add_argument(
        "--zones",
        type=float,
        nargs="*",
        default=[],
        help=(
            "span fractions from the root splitting the zones; needs one fewer "
            "than the entries in each --zone-pairs. With 32 spanwise elements "
            "these must be multiples of 1/32."
        ),
    )
    p.add_argument(
        "--angles",
        type=_parse_angle_list,
        nargs="+",
        default=[DEFAULT_ANGLES],
        help="repeating ply unit(s), slash separated, e.g. 0/90 45/-45",
    )
    p.add_argument(
        "--ply-mm",
        type=float,
        nargs="+",
        default=[PLY_MM_DEFAULT],
        help=f"ply thickness axis. Default: {PLY_MM_DEFAULT}",
    )
    p.add_argument(
        "--fiber",
        nargs="+",
        choices=sorted(FIBERS),
        default=["ud"],
        help="fibre architecture axis. Default: ud",
    )
    p.add_argument(
        "--static-line",
        nargs="+",
        default=[DEFAULT_STATIC_LINE],
        help=(
            "*STATIC line(s) 'initial, period, min, max'. The max increment "
            "sets a floor on increments per angle step, so it dominates "
            "runtime. Sweep it to calibrate; adopt only on matching forces."
        ),
    )
    p.add_argument(
        "--stress",
        action="store_true",
        help=(
            "also print per-ply stress at the reported angles (*EL PRINT S, "
            "GLOBAL=YES). Changes the cache key, so designs already solved "
            "without it are re-solved."
        ),
    )
    p.add_argument("--jobs", type=int, default=None, help="default cpu_count()-2")
    p.add_argument("--timeout-s", type=float, default=1800.0)
    p.add_argument("--run-id", default=None)
    p.add_argument(
        "--sync",
        action="store_true",
        help="run in the foreground (default detaches with nohup)",
    )
    p.add_argument(
        "--deck-only",
        action="store_true",
        help="write decks for each design and exit without solving",
    )
    return p


def designs_from_args(args) -> list[Design]:
    """The Cartesian product, in a stable order.

    ``--n-pairs`` is the old single-zone spelling; it and ``--zone-pairs`` mean
    the same thing, so taking both would leave the grid ambiguous.
    """
    if args.zone_pairs is not None and args.n_pairs is not None:
        raise SystemExit("give --zone-pairs or --n-pairs, not both")
    if args.zone_pairs is not None:
        zone_pairs = list(args.zone_pairs)
    elif args.n_pairs is not None:
        zone_pairs = [(n,) for n in args.n_pairs]
    else:
        zone_pairs = [(1,), (2,), (3,)]

    zones = tuple(args.zones)
    return [
        Design(
            zone_pairs=zp,
            ply_mm=ply,
            angles=ang,
            fiber=fiber,
            zones=zones,
            static_line=static,
            stress=bool(getattr(args, "stress", False)),
        )
        for fiber in args.fiber
        for ply in args.ply_mm
        for ang in args.angles
        for zp in zone_pairs
        for static in args.static_line
    ]


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    jobs = args.jobs if args.jobs is not None else _default_jobs()
    try:
        designs = designs_from_args(args)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    run_id = args.run_id or _run_id()
    run_root = RESULTS_ROOT / run_id

    if args.deck_only:
        run_root.mkdir(parents=True, exist_ok=True)
        for design in designs:
            out = run_root / f"deck_{design.cache_key()}.inp"
            out.write_text(deck_for(design)[0])
            print(f"wrote {out}  {design.fiber} {design.stack_label} "
                  f"ply={design.ply_mm:g}")
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
