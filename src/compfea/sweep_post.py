"""Sweep-level post: ``results/<run_id>/results.parquet`` -> ranking + figures.

``post.py``'s unit of work is one solve directory. This module's is a grid of
them: it reads the table ``sweep.run_sweep`` wrote and turns it into a ranked
CSV, comparison figures and a summary. Kept separate because the two have
different failure modes -- schema drift, error rows and degenerate grids are
problems only a sweep has -- and because ``sweep.py`` lazily imports ``post``
mid-solve to keep matplotlib off the solve path.

Two rules from ``cases/sweep_ubend/README.md`` are enforced here rather than
left to the reader:

- **Rank on ``f_90`` / ``f_180``, never on ``f_ratio_180_90``.** Under the
  linear-spring assumption behind ``M = 2U/theta`` that ratio is 2 by
  construction, so ranking on it sorts designs by their numerical noise.
  ``rank_designs`` refuses it as a sort key.
- **``F`` is a secant average**, exact only where ``M(theta)`` is straight. The
  ranking carries ``linearity_dev`` beside every force and the plot draws it as
  a whisker in the units of the objective, so two designs whose whiskers cover
  each other's dot are visibly not separable.

Nothing is written into ``cache/<cache_key>/``: that directory is the solve
record keyed on the design hash, and a derived file there would eventually be
read as cached truth.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from compfea.frd import (
    deformed_midsurface,
    disp_at_step,
    index_disp_blocks,
    read_nodes,
)
from compfea.shapes import plot_deformed_side, plot_shape_compare
from compfea.sweep import LENGTH_MM, LONG_AXIS, RESULTS_ROOT
from compfea.ubend import (
    DEFAULT_START_DEG,
    DEFAULT_STEP_DEG,
    step_index_for,
    theta_grid_deg,
)

# Columns that describe the design itself. A run whose ok rows are constant
# across all of these is not a design comparison.
DESIGN_AXES = (
    "stack",
    "fiber",
    "ply_mm",
    "angles",
    "n_zones",
    "zone_bounds",
    "zone_pairs",
    "n_plies_root",
    "n_plies_tip",
    "thickness_root_mm",
    "thickness_tip_mm",
    # older schema (results/ubend-demo-n1)
    "n_pairs",
    "n_plies",
    "thickness_mm",
)

# Solver settings and bookkeeping: varying these is a calibration, not a design.
SOLVER_COLS = (
    "static_line",
    "cache_key",
    "wall_time_s",
    "increments",
    "elapsed_s",
)

RATIO_COL = "f_ratio_180_90"
_F_COL = re.compile(r"^f_(\d+(?:\.\d+)?)$")

DEFAULT_LINEARITY_WARN = 0.05


def load_results(run_dir: str | Path) -> pd.DataFrame:
    """Read a sweep's ``results.parquet``; falls back to ``results.csv``."""
    run_dir = resolve_run_dir(run_dir)
    parquet = run_dir / "results.parquet"
    if parquet.is_file():
        return pd.read_parquet(parquet)
    csv = run_dir / "results.csv"
    if csv.is_file():
        return pd.read_csv(csv)
    raise FileNotFoundError(f"no results.parquet or results.csv in {run_dir}")


def resolve_run_dir(run: str | Path) -> Path:
    """Accept a path or a bare run id under ``results/``."""
    path = Path(run)
    if path.is_dir():
        return path
    candidate = RESULTS_ROOT / str(run)
    if candidate.is_dir():
        return candidate
    raise FileNotFoundError(f"no sweep run at {run!r} or {candidate}")


def ok_rows(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split into (solved, failed). A failed row has no trustworthy force."""
    if "status" not in frame:
        return frame.copy(), frame.iloc[0:0].copy()
    ok = frame[frame["status"] == "ok"].reset_index(drop=True)
    bad = frame[frame["status"] != "ok"].reset_index(drop=True)
    return ok, bad


def report_angles(frame: pd.DataFrame) -> list[float]:
    """Angles present as ``f_<deg>`` columns, ascending.

    Discovered rather than hardcoded: ``report_deg`` is a per-design property,
    so a sweep may carry angles other than 90 and 180. ``f_ratio_180_90`` is
    excluded by name -- it matches no angle, but only because of the underscore.
    """
    out = []
    for col in frame.columns:
        if col == RATIO_COL:
            continue
        m = _F_COL.match(str(col))
        if m:
            out.append(float(m.group(1)))
    return sorted(out)


def angle_tag(deg: float) -> str:
    """The suffix ``sweep.py`` builds its columns with."""
    return str(int(deg)) if float(deg).is_integer() else f"{deg:g}"


def column(frame: pd.DataFrame, name: str) -> pd.Series:
    """A column, or an all-NaN series if this run's schema lacks it."""
    if name in frame:
        return frame[name]
    return pd.Series([float("nan")] * len(frame), index=frame.index, name=name)


def varying_axes(frame: pd.DataFrame) -> list[str]:
    """Design axes that actually take more than one value in this run."""
    out = []
    for axis in DESIGN_AXES:
        if axis not in frame:
            continue
        if frame[axis].astype(str).nunique(dropna=False) > 1:
            out.append(axis)
    return out


# Preferred order for labelling: the axes a person reads a layup by first.
# The grid's axes are heavily correlated -- stack, zone_pairs, n_plies_root and
# thickness_root_mm all move together on a ply-count sweep -- so naming every
# varying one produces a label six times longer than it needs to be.
_LABEL_PRIORITY = (
    "stack",
    "fiber",
    "angles",
    "ply_mm",
    "zone_bounds",
    "zone_pairs",
    "n_plies_root",
    "n_plies_tip",
    "thickness_root_mm",
    "thickness_tip_mm",
    "n_pairs",
    "n_plies",
    "thickness_mm",
)


def _fmt(value) -> str:
    """Format a cell for a label: no 1.2000000000000002."""
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def label_axes(frame: pd.DataFrame) -> list[str]:
    """The fewest varying axes that still tell the designs apart.

    Greedy over ``_LABEL_PRIORITY``: take axes until every row has a distinct
    label, then stop. On a ply-count sweep this collapses six correlated
    columns to ``stack`` alone.
    """
    varying = varying_axes(frame)
    if not varying:
        return []
    chosen: list[str] = []
    ordered = [a for a in _LABEL_PRIORITY if a in varying]
    ordered += [a for a in varying if a not in ordered]
    for axis in ordered:
        chosen.append(axis)
        combined = frame[chosen].astype(str).agg("|".join, axis=1)
        if combined.nunique() == len(frame):
            return chosen
    return chosen


def design_label(row: pd.Series, axes: Sequence[str]) -> str:
    """Label a design by the given axes, falling back to its identity."""
    usable = [a for a in axes if a in row.index]
    if usable:
        if usable == ["stack"]:
            return _fmt(row["stack"])
        return "  ".join(f"{a.replace('_mm', '')}={_fmt(row[a])}" for a in usable)
    for fallback in ("stack", "cache_key"):
        if fallback in row.index and pd.notna(row[fallback]):
            return _fmt(row[fallback])
    return "design"


def rank_designs(frame: pd.DataFrame, *, by: str = "f_180") -> pd.DataFrame:
    """Solved designs sorted by one objective, strongest first.

    ``f_ratio_180_90`` is refused: it is ~2 for every design by construction,
    so sorting on it ranks numerical noise. See the module docstring.
    """
    if by == RATIO_COL:
        raise ValueError(
            f"{RATIO_COL} is a diagnostic, not an objective -- it is 2 by "
            "construction under the linear-spring assumption. Rank on an "
            "f_<angle> column instead."
        )
    if by not in frame:
        raise KeyError(f"{by!r} not in results (have: {report_angles(frame)})")
    return frame.sort_values(by, ascending=False).reset_index(drop=True)


def suspect_flags(
    frame: pd.DataFrame, angles: Sequence[float], warn: float
) -> pd.DataFrame:
    """Add ``suspect_f_<deg>`` where the secant/tangent gap exceeds ``warn``.

    Missing deviation data stays ``pd.NA``, not ``False``. ``NaN > warn`` is
    ``False``, so a plain comparison would stamp "not suspect" on a run that
    predates the linearity columns entirely -- reporting an unmeasured quantity
    as a clean bill of health.
    """
    out = frame.copy()
    for deg in angles:
        tag = angle_tag(deg)
        dev = column(out, f"linearity_dev_f_{tag}").abs()
        flag = (dev > warn).astype("boolean")
        out[f"suspect_f_{tag}"] = flag.mask(dev.isna(), pd.NA)
    return out


def classify(frame: pd.DataFrame) -> str:
    """``design``, ``calibration`` or ``single``.

    Calibration is checked before anything else: ``results/incr-calib`` varies
    only ``static_line`` and its forces are identical to the last digit, which
    would otherwise read as a design comparison of five identical bars.
    """
    if len(frame) <= 1:
        return "single"
    if not varying_axes(frame):
        return "calibration"
    return "design"


def calibration_table(frame: pd.DataFrame, angles: Sequence[float]) -> pd.DataFrame:
    """Solver cost against the finest-increment row, which is taken as truth."""
    # every solver column this run carries, except the design hash, which is
    # an identity rather than a setting anyone tuned
    cols = [c for c in SOLVER_COLS if c != "cache_key"]
    keep = [c for c in cols if c in frame] + [
        f"f_{angle_tag(d)}" for d in angles if f"f_{angle_tag(d)}" in frame
    ]
    table = frame[keep].copy()
    # The reference is the finest increment, which is only the last row once
    # the table is sorted by it. Without that column there is no ordering and
    # so no defensible reference to take deviations against.
    if "increments" in table and len(table):
        table = table.sort_values("increments").reset_index(drop=True)
        reference = table.iloc[-1]
    else:
        reference = None
    for deg in angles:
        col = f"f_{angle_tag(deg)}"
        if reference is not None and col in table and reference[col]:
            table[f"dev_{col}"] = (table[col] - reference[col]) / reference[col]
    if reference is not None and "wall_time_s" in table:
        table["speedup_vs_finest"] = reference["wall_time_s"] / table["wall_time_s"]
    return table


def deviation_note(frame: pd.DataFrame, deg: float) -> str:
    """``  (dev at 175°)`` when the deviation did not come from ``deg`` itself.

    There is no central difference at a path endpoint, so the largest angle's
    ``linearity_dev`` falls back to the nearest interior one -- 175 on the
    default 5-degree grid. A label reading "180" over a number measured at 175
    would overstate the plot's own precision, so the axis says which angle it
    is. Empty string when the value really is at ``deg``.
    """
    at = column(frame, f"linearity_dev_theta_f_{angle_tag(deg)}").dropna().unique()
    if not len(at) or all(abs(float(a) - float(deg)) <= 1e-9 for a in at):
        return ""
    return "  (dev at " + ", ".join(f"{float(a):g}" for a in at) + "°)"


def plot_rank(
    frame: pd.DataFrame,
    out_svg: str | Path,
    *,
    angles: Sequence[float],
    axes: Sequence[str],
    warn: float = DEFAULT_LINEARITY_WARN,
    title: str = "Tip-normal spring force by design",
) -> Path:
    """One panel per reported angle; dot at F with a linearity whisker.

    The whisker is ``|linearity_dev| * F`` -- the secant-vs-tangent
    disagreement expressed in the units of the objective, so it can be compared
    against the gaps between designs directly. Where the deviation was taken at
    a different angle than the column is named for (on a 5-degree path the
    ``f_180`` deviation comes from 175), the axis label says so.
    """
    out_svg = Path(out_svg)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")
    angles = list(angles)
    if not angles:
        raise ValueError(
            "no f_<angle> columns in these results, so there is nothing to "
            "rank; the sweep wrote no reported forces"
        )
    labels = [design_label(r, axes) for _, r in frame.iterrows()]
    # leave horizontal room for the longest label rather than letting it run
    # off the canvas, which is what a fixed figsize does on a zoned design
    label_in = min(0.09 * max((len(t) for t in labels), default=8), 3.6)
    fig, panels = plt.subplots(
        1, len(angles),
        figsize=(label_in + 3.2 * len(angles), 0.45 * len(frame) + 2.2),
        sharey=True, squeeze=False,
    )
    y = range(len(frame))
    for ax, deg in zip(panels[0], angles, strict=True):
        tag = angle_tag(deg)
        f = column(frame, f"f_{tag}").astype(float)
        dev = column(frame, f"linearity_dev_f_{tag}").astype(float).abs()
        # A zero-length whisker means "measured, and straight". Where the
        # deviation was never recorded there must be no whisker at all, or the
        # plot claims a check it never ran.
        measured = dev.notna()
        err = (dev * f).where(measured, 0.0)
        flagged = (dev > warn).fillna(False)
        ax.errorbar(
            f, list(y), xerr=err, fmt="o", ms=6, lw=0, elinewidth=1.4,
            capsize=3, color="#1f4e79", ecolor="#7f8c8d", zorder=2,
        )
        if flagged.any():
            ax.scatter(
                f[flagged], [i for i, s in zip(y, flagged, strict=True) if s],
                s=90, facecolors="none", edgecolors="#c0392b",
                linewidths=1.6, zorder=3,
                label=f"|linearity dev| > {warn:.0%}",
            )
            ax.legend(loc="lower right", fontsize=7)
        ax.set_xlabel("F (N)")
        note = deviation_note(frame, deg)
        if not measured.any():
            note = "  (no linearity data)"
        ax.set_title(f"F at θ={deg:g}°{note}", fontsize=10)
    panels[0][0].set_yticks(list(y))
    panels[0][0].set_yticklabels(labels, fontsize=8)
    panels[0][0].invert_yaxis()
    fig.suptitle(title)
    fig.tight_layout()
    fig.subplots_adjust(left=label_in / fig.get_figwidth())
    fig.savefig(out_svg, format="svg")
    plt.close(fig)
    return out_svg


def plot_ratio_check(
    frame: pd.DataFrame, out_svg: str | Path, *, low: float, high: float
) -> Path:
    """Diagnostic only: F_high vs F_low against the y=2x the model predicts."""
    out_svg = Path(out_svg)
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    x = column(frame, f"f_{angle_tag(low)}").astype(float)
    y = column(frame, f"f_{angle_tag(high)}").astype(float)
    ax.scatter(x, y, s=34, color="#1f4e79", zorder=3)
    if len(x.dropna()):
        lo, hi = float(x.min()), float(x.max())
        span = [lo, hi if hi > lo else lo + 1.0]
        ax.plot(span, [2.0 * v for v in span], ls="--", color="#c0392b",
                label="y = 2x (linear spring)")
        ax.legend(fontsize=8)
    ax.set_xlabel(f"F at θ={low:g}° (N)")
    ax.set_ylabel(f"F at θ={high:g}° (N)")
    ax.set_title("DIAGNOSTIC — deviation from 2, not an objective")
    fig.tight_layout()
    fig.savefig(out_svg, format="svg")
    plt.close(fig)
    return out_svg


def plot_calibration(
    table: pd.DataFrame, out_svg: str | Path, *, angles: Sequence[float]
) -> Path:
    """Wall time against increment count, with the force held up beside it."""
    out_svg = Path(out_svg)
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot(table["increments"], table["wall_time_s"], "o-", color="#1f4e79")
    for _, row in table.iterrows():
        ax.annotate(
            str(row.get("static_line", "")), (row["increments"], row["wall_time_s"]),
            textcoords="offset points", xytext=(6, -9), fontsize=7,
        )
    ax.set_xlabel("increments")
    ax.set_ylabel("wall time (s)")
    tags = [f"f_{angle_tag(d)}" for d in angles if f"f_{angle_tag(d)}" in table]
    spread = max(
        (float(table[t].max() - table[t].min()) for t in tags), default=0.0
    )
    ax.set_title(
        f"Solver increment calibration — force spread {spread:.3g} N over the grid"
    )
    fig.tight_layout()
    fig.savefig(out_svg, format="svg")
    plt.close(fig)
    return out_svg


def shape_curves(
    run_dir: Path,
    cache_key: str,
    angles_deg: Sequence[float],
    *,
    end_deg: float,
    step_deg: float = DEFAULT_STEP_DEG,
    start_deg: float = DEFAULT_START_DEG,
) -> dict[float, pd.DataFrame]:
    """Deformed mid-surfaces for one design, keyed by angle.

    ``end_deg`` must be the design's own deepest reported angle, because that
    is what ``sweep.angles_for`` passes to ``theta_grid_deg`` and
    ``theta_grid_deg`` pins its last entry to it. Rebuilding the path with a
    fixed 180 instead would shift every step index for a design whose deepest
    angle is off the regular grid, and the pose would come back under the wrong
    angle's name.

    An angle with no DISP is skipped rather than substituted: ``*NODE FILE`` is
    requested on selected steps only and ccx carries it forward, so a run
    typically has no shape at all below its first requested angle. Only that
    absence is caught -- a corrupt record or an angle that was never on the
    deck's path is a different problem and is raised, not quietly dropped.
    """
    frd_path = run_dir / "cache" / cache_key / "ccx" / "job.frd"
    if not frd_path.is_file():
        raise FileNotFoundError(f"no job.frd for design {cache_key} at {frd_path}")
    grid = theta_grid_deg(step_deg=step_deg, start_deg=start_deg, end_deg=end_deg)
    blocks = index_disp_blocks(frd_path)
    nodes = read_nodes(frd_path)
    out: dict[float, pd.DataFrame] = {}
    for deg in angles_deg:
        # Off-grid angles are the caller's error, not missing output: say so
        # rather than reporting it as "this run has no DISP there".
        step = step_index_for(grid, float(deg))
        try:
            _, disp = disp_at_step(frd_path, step, blocks=blocks)
        except LookupError:
            continue
        out[float(deg)] = deformed_midsurface(nodes, disp, long_axis=LONG_AXIS)
    return out


def write_shapes(
    run_dir: Path,
    ranked: pd.DataFrame,
    out: Path,
    *,
    n: int,
    angles_deg: Sequence[float],
    axes: Sequence[str],
    end_deg: float,
) -> dict[str, str]:
    """Per-design shape plots for the top ``n`` designs, plus a comparison."""
    files: dict[str, str] = {}
    compare: dict[str, pd.DataFrame] = {}
    # The comparison is at one angle for every design. Taking each design's own
    # deepest available angle instead would draw a 90-degree pose under a
    # 180-degree title as soon as one design stopped short.
    target = max(float(d) for d in angles_deg)
    for _, row in ranked.head(n).iterrows():
        key = str(row.get("cache_key", ""))
        if not key:
            continue
        try:
            curves = shape_curves(run_dir, key, angles_deg, end_deg=end_deg)
        except FileNotFoundError:
            continue
        if not curves:
            continue
        label = design_label(row, axes)
        # The target arc's radius is L/theta, so it must use this design's own
        # moment arm. LENGTH_MM is the same number today, but arm_mm is the tip
        # minus the outboard clamp edge: add a clamp patch and the two diverge,
        # and the overlay would be drawn at the wrong radius while the solved
        # curve stayed right -- which reads as physics.
        arm = float(row["arm_mm"]) if pd.notna(row.get("arm_mm")) else LENGTH_MM
        svg = plot_deformed_side(
            curves, out / f"sweep_post_shape_{key}.svg",
            length_mm=arm, long_axis=LONG_AXIS,
            title=f"{label} — solved shape vs target arc",
        )
        files[f"shape_{key}"] = str(svg)
        if target in curves:
            # Labels come from the varying axes only, so two designs can share
            # one -- in a calibration every design does. Keying the comparison
            # on the label alone would drop all but the last without saying so.
            name = label if label not in compare else f"{label} [{key[:8]}]"
            compare[name] = curves[target]
    if len(compare) >= 2:
        svg = plot_shape_compare(
            compare, out / f"sweep_post_shape_compare_{angle_tag(target)}.svg",
            theta_deg=target,
        )
        files["shape_compare"] = str(svg)
    return files


def _jsonable(value):
    """NaN is not JSON. Emit null so a poller sees "absent", not a bare NaN."""
    if value is None or pd.isna(value):
        return None
    return float(value)


def sweep_post(
    run: str | Path,
    *,
    top: int = 0,
    warn: float = DEFAULT_LINEARITY_WARN,
    out_dir: str | Path | None = None,
    plots: bool = True,
    shapes: int = 0,
    shape_deg: Sequence[float] = (90.0, 180.0),
) -> dict:
    """Post one sweep run. Returns the summary that is also written as json."""
    run_dir = resolve_run_dir(run)
    out = Path(out_dir) if out_dir else run_dir
    out.mkdir(parents=True, exist_ok=True)

    frame = load_results(run_dir)
    solved, failed = ok_rows(frame)
    angles = report_angles(frame)
    axes = varying_axes(solved)
    labels = label_axes(solved)
    mode = classify(solved)

    missing = [
        c
        for c in (
            [f"linearity_dev_f_{angle_tag(d)}" for d in angles]
            + ["max_linearity_dev"]
        )
        if c not in frame
    ]

    summary: dict = {
        "run_dir": str(run_dir),
        "mode": mode,
        "n_rows": int(len(frame)),
        "n_ok": int(len(solved)),
        "n_error": int(len(failed)),
        "report_angles": angles,
        "varying_axes": axes,
        "missing_columns": missing,
        "linearity_warn": warn,
        "files": {},
        "errors": [
            {
                "stack": str(r.get("stack", "")),
                "cache_key": str(r.get("cache_key", "")),
                "error": str(r.get("error", "")),
            }
            for _, r in failed.iterrows()
        ],
    }

    if solved.empty:
        summary["note"] = "every design failed; nothing to rank"
        (out / "sweep_post.json").write_text(json.dumps(summary, indent=2) + "\n")
        return summary

    rank_on = 180.0 if 180.0 in angles else (angles[-1] if angles else None)
    ranked = solved
    if rank_on is not None:
        ranked = rank_designs(solved, by=f"f_{angle_tag(rank_on)}")
    ranked = suspect_flags(ranked, angles, warn)
    summary["ranked_on"] = f"f_{angle_tag(rank_on)}" if rank_on is not None else None

    keep = (
        [a for a in DESIGN_AXES if a in ranked]
        + [c for c in ("cache_key", "increments", "wall_time_s") if c in ranked]
        + [f"f_{angle_tag(d)}" for d in angles if f"f_{angle_tag(d)}" in ranked]
        + [
            c
            for d in angles
            for c in (
                f"linearity_dev_f_{angle_tag(d)}",
                f"linearity_dev_theta_f_{angle_tag(d)}",
                f"suspect_f_{angle_tag(d)}",
            )
            if c in ranked
        ]
        + [c for c in ("max_linearity_dev", RATIO_COL) if c in ranked]
    )
    rank_csv = out / "sweep_post_rank.csv"
    ranked[keep].to_csv(rank_csv, index=False)
    summary["files"]["rank_csv"] = str(rank_csv)

    for deg in angles:
        tag = angle_tag(deg)
        col = f"f_{tag}"
        if col not in ranked:
            continue
        best = ranked.iloc[0]
        worst = ranked.iloc[-1]
        summary[f"best_f_{tag}"] = {
            "design": design_label(best, labels),
            "value": _jsonable(best[col]),
            "suspect": (
                None
                if pd.isna(best.get(f"suspect_f_{tag}", pd.NA))
                else bool(best[f"suspect_f_{tag}"])
            ),
        }
        summary[f"worst_f_{tag}"] = {
            "design": design_label(worst, labels),
            "value": _jsonable(worst[col]),
        }

    if mode == "calibration":
        table = calibration_table(ranked, angles)
        cal_csv = out / "sweep_post_calibration.csv"
        table.to_csv(cal_csv, index=False)
        summary["files"]["calibration_csv"] = str(cal_csv)
        summary["note"] = (
            "this run varies only the solver increment; it is a calibration, "
            "not a design comparison, so no ranking figure was written"
        )
        if plots and "increments" in table and len(table) > 1:
            svg = plot_calibration(table, out / "sweep_post_calibration.svg",
                                   angles=angles)
            summary["files"]["calibration_svg"] = str(svg)
    elif mode == "single":
        summary["note"] = "single solved design — nothing to compare"
    elif plots:
        svg = plot_rank(
            ranked, out / "sweep_post_rank.svg", angles=angles, axes=labels, warn=warn,
            title=f"{run_dir.name}: tip-normal spring force by design",
        )
        summary["files"]["rank_svg"] = str(svg)
        if len(angles) >= 2:
            ratio = plot_ratio_check(
                ranked, out / "sweep_post_ratio_check.svg",
                low=angles[0], high=angles[-1],
            )
            summary["files"]["ratio_svg"] = str(ratio)

    if top:
        summary["top"] = [
            {
                "design": design_label(r, labels),
                **{
                    f"f_{angle_tag(d)}": _jsonable(r[f"f_{angle_tag(d)}"])
                    for d in angles
                    if f"f_{angle_tag(d)}" in r.index
                    and pd.notna(r[f"f_{angle_tag(d)}"])
                },
            }
            for _, r in ranked.head(top).iterrows()
        ]

    if shapes:
        summary["files"].update(
            write_shapes(
                run_dir, ranked, out, n=shapes, angles_deg=shape_deg,
                axes=labels, end_deg=max(angles) if angles else 180.0,
            )
        )
        if not any(k.startswith("shape_") for k in summary["files"]):
            summary["shape_note"] = (
                "no deformed shapes written: the cached .frd files carry no DISP "
                "at these angles (the run predates *NODE FILE, or asked for it "
                "only at deeper angles)"
            )

    (out / "sweep_post.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    summary["files"]["summary_json"] = str(out / "sweep_post.json")
    return summary


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="compfea-sweep-post",
        description=(
            "Post a finished sweep: results.parquet -> ranked CSV and figures. "
            "Ranks on f_90/f_180; f_ratio_180_90 is a diagnostic only."
        ),
    )
    p.add_argument("run", help="results/<run_id> or a bare <run_id>")
    p.add_argument("--top", type=int, default=5,
                   help="designs to list in the summary (default 5)")
    p.add_argument("--linearity-warn", type=float, default=DEFAULT_LINEARITY_WARN,
                   help="flag a design above this |linearity dev| (default 0.05)")
    p.add_argument("--out-dir", default=None, help="default: the run dir")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument(
        "--shapes", type=int, default=0, metavar="N",
        help=(
            "read the cached .frd for the top N designs and plot the real "
            "deformed shape. Off by default: these files run 15-130 MB each."
        ),
    )
    p.add_argument(
        "--shape-deg", type=float, nargs="+", default=[90.0, 180.0],
        help="angles to draw shapes at (default 90 180)",
    )
    args = p.parse_args(argv)

    try:
        summary = sweep_post(
            args.run,
            top=args.top,
            warn=args.linearity_warn,
            out_dir=args.out_dir,
            plots=not args.no_plots,
            shapes=args.shapes,
            shape_deg=args.shape_deg,
        )
    except (FileNotFoundError, LookupError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"mode     {summary['mode']}")
    print(f"designs  {summary['n_ok']} ok, {summary['n_error']} error")
    if summary["varying_axes"]:
        print(f"varying  {', '.join(summary['varying_axes'])}")
    if summary.get("missing_columns"):
        print(f"missing  {', '.join(summary['missing_columns'])}")
    if summary.get("note"):
        print(f"note     {summary['note']}")
    if summary.get("shape_note"):
        print(f"shapes   {summary['shape_note']}")
    for key, path in summary["files"].items():
        print(f"  {key:16s} {path}")
    for err in summary["errors"]:
        print(f"  error {err['stack']}: {err['error']}")
    return 1 if summary["n_ok"] == 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
