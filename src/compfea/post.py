"""Post-process tip U-bend solves: ELSE energy -> tip-normal F(θ).

Reads a CalculiX ``.dat`` with ``*EL PRINT, TOTALS=ONLY`` / ELSE, maps each
end-of-step tot-time to the tip-clamp angle grid, and reports tip-normal
spring force::

    M = 2U / θ
    F = M / arm

Do not use tip |RF| as F_90 / F_180. CLI writes CSV + SVG for a run dir.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from compfea.run import parse_dat_energy
from compfea.ubend import force_at_theta, step_index_for, theta_grid_deg
from compfea.shapes import (
    parse_nodes_elements,
    parse_nset,
    plot_bend_side,
    plot_planform,
)

REPORT_DEG = (90.0, 180.0)


def force_curve(
    energy: pd.DataFrame,
    angles_deg: Sequence[float],
    arm_mm: float,
    *,
    elset: str = "blade",
) -> pd.DataFrame:
    """One row per angle: end-of-step U, M, F tip-normal."""
    rows = []
    for deg in angles_deg:
        idx = step_index_for(angles_deg, float(deg))
        u, m, f = force_at_theta(
            energy,
            theta_deg=float(deg),
            step_index=idx,
            arm_mm=arm_mm,
            elset=elset,
        )
        rows.append(
            {
                "theta_deg": float(deg),
                "step": idx,
                "time": float(idx),
                "U_nmm": u,
                "M_nmm": m,
                "F_N": f,
            }
        )
    return pd.DataFrame(rows)


def summary_at(
    curve: pd.DataFrame,
    report_deg: Sequence[float] = REPORT_DEG,
) -> dict[str, float]:
    """Named F_θ values plus F_180/F_90 when both exist."""
    out: dict[str, float] = {}
    by_theta = {float(r.theta_deg): float(r.F_N) for r in curve.itertuples()}
    for deg in report_deg:
        key = f"F_{int(deg) if float(deg).is_integer() else deg}"
        matches = [t for t in by_theta if abs(t - float(deg)) <= 1e-9]
        if matches:
            out[key] = by_theta[matches[0]]
    if "F_90" in out and "F_180" in out and out["F_90"] != 0:
        out["F_180_over_F_90"] = out["F_180"] / out["F_90"]
    return out


def plot_force_curve(
    curve: pd.DataFrame,
    out_svg: str | Path,
    *,
    title: str = "Tip-normal spring force",
    mark_deg: Sequence[float] = REPORT_DEG,
) -> Path:
    """Write F(θ) SVG; marks report angles when present on the curve."""
    out_svg = Path(out_svg)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(curve["theta_deg"], curve["F_N"], color="#1f4e79", lw=1.8)
    present = set(float(t) for t in curve["theta_deg"])
    for deg in mark_deg:
        if any(abs(t - float(deg)) <= 1e-9 for t in present):
            row = curve.loc[(curve["theta_deg"] - float(deg)).abs() <= 1e-9].iloc[0]
            ax.scatter([row.theta_deg], [row.F_N], zorder=3, color="#c0392b", s=36)
            ax.annotate(
                f"F_{int(deg) if float(deg).is_integer() else deg}="
                f"{row.F_N:.3g} N",
                (row.theta_deg, row.F_N),
                textcoords="offset points",
                xytext=(6, 8),
                fontsize=9,
            )
    ax.set_xlabel("tip tangent θ (deg)")
    ax.set_ylabel("tip-normal F (N)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_svg, format="svg")
    plt.close(fig)
    return out_svg


def post_run(
    run_dir: str | Path,
    *,
    arm_mm: float,
    start_deg: float = 1.0,
    step_deg: float = 1.0,
    end_deg: float | None = None,
    elset: str = "blade",
    dat_name: str = "job.dat",
) -> dict:
    """Post one solve directory: CSV + SVG + summary dict.

    ``end_deg`` defaults to the largest integer tot-time present in the ``.dat``
    (so a stopped-at-90 run still posts without claiming 180).
    """
    run_dir = Path(run_dir)
    dat = run_dir / "ccx" / dat_name
    if not dat.is_file():
        # allow run_dir itself to be the ccx folder
        alt = run_dir / dat_name
        if alt.is_file():
            dat = alt
            ccx_parent = run_dir
            out_dir = run_dir
        else:
            raise FileNotFoundError(f"no {dat_name} under {run_dir}")
    else:
        ccx_parent = run_dir / "ccx"
        out_dir = run_dir

    energy = parse_dat_energy(dat)
    blade = energy[energy["elset"] == elset.lower()]
    if blade.empty:
        raise LookupError(f"no energy rows for elset={elset!r} in {dat}")
    tmax = float(blade["time"].max())
    if end_deg is None:
        # end-of-step times are integers 1..N
        end_deg = float(math.floor(tmax + 1e-9))
        if end_deg < start_deg:
            raise LookupError(
                f"dat tmax={tmax:g} is before start_deg={start_deg:g}"
            )
    angles = theta_grid_deg(
        step_deg=step_deg, start_deg=start_deg, end_deg=end_deg
    )
    # drop angles whose end-of-step is past available time
    angles = [a for a in angles if float(step_index_for(angles, a)) <= tmax + 1e-6]
    if not angles:
        raise LookupError(f"no completed tip steps in {dat} (tmax={tmax:g})")

    curve = force_curve(energy, angles, arm_mm, elset=elset)
    csv_path = out_dir / "post_F_theta.csv"
    svg_path = out_dir / "post_F_theta.svg"
    meta_path = out_dir / "post_meta.json"
    curve.to_csv(csv_path, index=False)
    plot_force_curve(
        curve,
        svg_path,
        title=f"Tip-normal F(θ)  arm={arm_mm:g} mm",
    )
    summ = summary_at(curve)
    meta = {
        "dat": str(dat),
        "arm_mm": arm_mm,
        "start_deg": start_deg,
        "step_deg": step_deg,
        "end_deg": end_deg,
        "n_angles": len(angles),
        "tmax": tmax,
        "summary": summ,
        "csv": str(csv_path),
        "svg": str(svg_path),
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def shapes_run(
    run_dir: str | Path,
    *,
    arm_mm: float,
    long_axis: str = "x",
    angles_deg: Sequence[float] = (45.0, 90.0, 135.0, 180.0),
) -> dict:
    """Planform + circular-arc side poses from the run's ``deck.inp`` / ``job.inp``."""
    run_dir = Path(run_dir)
    inp = run_dir / "deck.inp"
    if not inp.is_file():
        inp = run_dir / "ccx" / "job.inp"
    if not inp.is_file():
        raise FileNotFoundError(f"no deck.inp or ccx/job.inp under {run_dir}")
    nodes, elements = parse_nodes_elements(inp)
    tip = parse_nset(inp, "far_face")
    heal = parse_nset(inp, "fixed_end")
    axis = 0 if long_axis == "x" else 1
    s0 = max(nodes[n][axis] for n in heal) if heal else min(v[axis] for v in nodes.values())
    plan = plot_planform(
        nodes,
        elements,
        run_dir / "post_planform.png",
        title=f"Undeformed mid-surface  ({inp.name})",
        tip_ids=tip,
        heal_ids=heal,
    )
    # only plot targets up to angles that make sense; include 90 always for saved 90 run
    side = plot_bend_side(
        nodes,
        run_dir / "post_bend_side.png",
        length_mm=arm_mm,
        s0=s0,
        long_axis=long_axis,
        angles_deg=angles_deg,
        title=f"Tip-drive circular poses  L={arm_mm:g} mm",
    )
    return {"planform": str(plan), "bend_side": str(side), "s0": s0, "n_nodes": len(nodes)}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Post tip U-bend: ELSE energy -> F(θ) CSV + SVG"
    )
    p.add_argument(
        "run_dir",
        type=Path,
        help="solve dir with ccx/job.dat (e.g. results/fin_ubend_90_saved)",
    )
    p.add_argument(
        "--arm-mm",
        type=float,
        required=True,
        help="undeformed free length / moment arm (mm)",
    )
    p.add_argument("--start-deg", type=float, default=1.0)
    p.add_argument("--step-deg", type=float, default=1.0)
    p.add_argument(
        "--end-deg",
        type=float,
        default=None,
        help="default: floor(tmax) from the .dat",
    )
    p.add_argument("--elset", default="blade")
    p.add_argument("--shapes", action="store_true",
                   help="also write planform + circular tip-drive side plots")
    p.add_argument("--long-axis", choices=("x", "y"), default="x",
                   help="span axis for shape plots (fin=x, strip=y)")
    args = p.parse_args(argv)

    meta = post_run(
        args.run_dir,
        arm_mm=args.arm_mm,
        start_deg=args.start_deg,
        step_deg=args.step_deg,
        end_deg=args.end_deg,
        elset=args.elset,
    )
    print(f"csv  {meta['csv']}")
    print(f"svg  {meta['svg']}")
    print(f"tmax {meta['tmax']:g}  angles={meta['n_angles']}")
    for k, v in meta["summary"].items():
        print(f"  {k} = {v:.6g}")
    if args.shapes:
        sh = shapes_run(
            args.run_dir, arm_mm=args.arm_mm, long_axis=args.long_axis
        )
        print(f"planform  {sh['planform']}")
        print(f"bend_side {sh['bend_side']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
