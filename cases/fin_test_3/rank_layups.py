#!/usr/bin/env python3
"""Rank a layup sweep's ``results.parquet`` against a target flex and kick point.

::

    dist = sqrt( ((flex_n - target_flex) / flex_tol)^2
                 + (kick_s_frac - target_frac)^2 )

``flex_tol`` (default 2 N) puts the two axes in comparable units -- a design one
``flex_tol`` off in flex is as far as one whose kick point is a full free-length
off, which never happens, so in practice flex dominates unless it is already
close. Rows with a kick-point ``warning`` (broad-plateau or bimodal curvature)
are dropped: "kick point" is not a single location for those blades. Rows that
errored or never reached 90 deg (``flex_n`` null) are dropped too.

No cost term yet. A later version ranks the *change* cost -- move a zone
boundary < change a ply count < swap a fabric < buy new stock.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

#: Centre of each kick bucket as a free-length fraction (heel <1/3 < mid < 2/3< tip).
BUCKET_FRAC = {"heel": 1.0 / 6.0, "mid": 0.5, "tip": 5.0 / 6.0}


def target_fraction(*, bucket: str | None, frac: float | None) -> float:
    if (bucket is None) == (frac is None):
        raise SystemExit("give exactly one of --target-kick or --target-kick-frac")
    if frac is not None:
        if not 0.0 <= frac <= 1.0:
            raise SystemExit("--target-kick-frac must be in [0, 1]")
        return float(frac)
    if bucket not in BUCKET_FRAC:
        raise SystemExit(f"--target-kick must be one of {sorted(BUCKET_FRAC)}")
    return BUCKET_FRAC[bucket]


def _warned(series: pd.Series) -> pd.Series:
    return series.notna() & (series.astype(str).str.strip() != "")


def drop_reasons(df: pd.DataFrame) -> dict[str, int]:
    """Count why rows are not rankable -- printed so a lost corner is visible."""
    errored = int((df["status"] != "ok").sum()) if "status" in df.columns else 0
    ok = df[df["status"] == "ok"] if "status" in df.columns else df
    no_ninety = 0
    warned = 0
    if "flex_n" in ok.columns:
        no_ninety = int(ok["flex_n"].isna().sum())
        ok = ok[ok["flex_n"].notna()]
    if "warning" in ok.columns:
        warned = int(_warned(ok["warning"]).sum())
    return {"errored": errored, "no_90_deg": no_ninety, "kick_warning": warned}


def rank(
    df: pd.DataFrame,
    *,
    target_flex: float,
    target_frac: float,
    flex_tol: float,
) -> pd.DataFrame:
    """Filtered, distance-scored, ascending. Pure; no I/O.

    Empty frame in -> empty frame out (an all-error sweep writes a parquet with
    no ``flex_n`` column); the caller turns that into a clear exit.
    """
    if flex_tol <= 0.0:
        raise SystemExit("--flex-tol must be > 0")
    keep = df.copy()
    if "status" in keep.columns:
        keep = keep[keep["status"] == "ok"]
    if "flex_n" not in keep.columns or "kick_s_frac" not in keep.columns:
        return keep.iloc[0:0]
    keep = keep[keep["flex_n"].notna()]
    if "warning" in keep.columns:
        keep = keep[~_warned(keep["warning"])]
    keep = keep.copy()
    keep["flex_dist"] = (keep["flex_n"] - target_flex) / flex_tol
    keep["kick_dist"] = keep["kick_s_frac"] - target_frac
    keep["dist"] = np.hypot(keep["flex_dist"], keep["kick_dist"])
    return keep.sort_values("dist", kind="stable").reset_index(drop=True)


def _scatter(
    ranked: pd.DataFrame,
    dropped: pd.DataFrame,
    *,
    target_flex: float,
    target_frac: float,
    top: int,
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    if not dropped.empty:
        ax.scatter(
            dropped["flex_n"], dropped["kick_s_frac"],
            facecolors="none", edgecolors="0.6", s=40,
            label="dropped (warned / no 90 deg)",
        )
    ax.scatter(
        ranked["flex_n"], ranked["kick_s_frac"], c=ranked["dist"],
        cmap="viridis", s=55, label="candidate",
    )
    ax.plot(
        [target_flex], [target_frac], "*", color="C3", ms=18, label="target"
    )
    for i, r in ranked.head(max(0, top)).iterrows():
        ax.annotate(
            f"{i + 1}. {str(r['fingerprint'])[:8]}",
            (r["flex_n"], r["kick_s_frac"]),
            textcoords="offset points", xytext=(6, 4), fontsize=7,
        )
    ax.set_xlabel("flex  W(90 deg)  (N)")
    ax.set_ylabel("kick point  s_kick / free length")
    ax.set_ylim(-0.02, 1.02)
    ax.axhspan(0, 1 / 3, color="0.93", zorder=0)
    ax.axhspan(2 / 3, 1, color="0.93", zorder=0)
    ax.legend(fontsize=8, loc="best")
    ax.set_title("layup sweep: distance to target")
    fig.tight_layout()
    fig.savefig(path, format="svg")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("parquet", type=Path, help="results.parquet from sweep_layups.py")
    p.add_argument("--target-flex", type=float, required=True, help="W(90 deg), N")
    p.add_argument(
        "--target-kick",
        choices=sorted(BUCKET_FRAC),
        default=None,
        help="aims at the bucket's centre fraction (heel 1/6, mid 1/2, tip "
        "5/6); ranks by distance to that, not by bucket membership",
    )
    p.add_argument("--target-kick-frac", type=float, default=None)
    p.add_argument("--flex-tol", type=float, default=2.0)
    p.add_argument("--top", type=int, default=10, help="rows to label on the plot")
    p.add_argument(
        "--out", type=Path, default=None, help="ranked CSV (default: next to parquet)"
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    target_frac = target_fraction(
        bucket=args.target_kick, frac=args.target_kick_frac
    )
    df = pd.read_parquet(args.parquet)
    ranked = rank(
        df,
        target_flex=args.target_flex,
        target_frac=target_frac,
        flex_tol=args.flex_tol,
    )
    reasons = drop_reasons(df)
    n_dropped = sum(reasons.values())
    if n_dropped:
        detail = ", ".join(f"{k}={v}" for k, v in reasons.items() if v)
        print(f"dropped {n_dropped} of {len(df)} rows: {detail}")
    if ranked.empty:
        raise SystemExit(
            "no rankable rows (all errored, warned, or short of 90 deg). "
            "A whole grid short of 90 deg usually means --uy-frac was too low "
            "in the sweep."
        )

    out = args.out or args.parquet.parent / "ranked.csv"
    ranked.to_csv(out, index=False)
    svg = out.with_suffix(".svg")
    kept_fps = set(ranked["fingerprint"]) if "fingerprint" in ranked else set()
    dropped = df[
        (df.get("status", "ok") == "ok")
        & ~df["fingerprint"].isin(kept_fps)
    ] if "fingerprint" in df.columns else df.iloc[0:0]
    _scatter(
        ranked, dropped,
        target_flex=args.target_flex, target_frac=target_frac,
        top=args.top, path=svg,
    )

    cols = [
        c
        for c in (
            "fingerprint", "flex_n", "kick_bucket", "kick_s_frac",
            "migration_delta_mm", "through_thickness_mm", "dist", "core",
        )
        if c in ranked.columns
    ]
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(ranked[cols].head(max(1, args.top)).to_string(index=False))
    print(f"\nwrote {out} and {svg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
