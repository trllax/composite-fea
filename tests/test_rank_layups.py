"""Ranker: filter warned / failed rows, score by distance, order ascending."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


rl = _module(ROOT / "cases/fin_test_3/rank_layups.py", "rank_layups_under_test")


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            # name        status  flex_n  kick_s_frac  warning
            ("near",      "ok",   12.2,   0.80,        None),
            ("far_flex",  "ok",   6.0,    0.79,        None),
            ("far_kick",  "ok",   12.0,   0.20,        None),
            ("warned",    "ok",   12.0,   0.80,        "bimodal curvature"),
            ("errored",   "error", math.nan, math.nan, None),
            ("no_ninety", "ok",   math.nan, 0.5,       None),
        ],
        columns=["fingerprint", "status", "flex_n", "kick_s_frac", "warning"],
    )


def test_target_fraction_maps_buckets_and_rejects_both_or_neither():
    assert rl.target_fraction(bucket="tip", frac=None) == pytest.approx(5 / 6)
    assert rl.target_fraction(bucket="heel", frac=None) == pytest.approx(1 / 6)
    assert rl.target_fraction(bucket=None, frac=0.42) == 0.42
    with pytest.raises(SystemExit):
        rl.target_fraction(bucket="tip", frac=0.5)
    with pytest.raises(SystemExit):
        rl.target_fraction(bucket=None, frac=None)


def test_rank_drops_warned_errored_and_short_rows():
    ranked = rl.rank(
        _frame(), target_flex=12.0, target_frac=0.8, flex_tol=2.0
    )
    kept = list(ranked["fingerprint"])
    assert kept == ["near", "far_kick", "far_flex"]
    assert "warned" not in kept
    assert "errored" not in kept
    assert "no_ninety" not in kept


def test_rank_distance_is_euclidean_in_scaled_axes():
    ranked = rl.rank(
        _frame(), target_flex=12.0, target_frac=0.8, flex_tol=2.0
    ).set_index("fingerprint")
    # near: flex_dist = 0.2/2 = 0.1, kick_dist = 0.0  -> 0.1
    assert ranked.loc["near", "dist"] == pytest.approx(0.1, abs=1e-9)
    # far_kick: flex_dist = 0, kick_dist = 0.2 - 0.8 = -0.6 -> 0.6
    assert ranked.loc["far_kick", "dist"] == pytest.approx(0.6, abs=1e-9)
    # far_flex: flex_dist = (6-12)/2 = -3, kick_dist = -0.01 -> ~3.0
    assert ranked.loc["far_flex", "dist"] == pytest.approx(
        math.hypot(3.0, 0.01), abs=1e-9
    )


def test_flex_tol_reweights_the_axes():
    tight = rl.rank(_frame(), target_flex=12.0, target_frac=0.8, flex_tol=2.0)
    # tight flex tolerance: the 6 N design is by far the worst
    assert list(tight["fingerprint"]) == ["near", "far_kick", "far_flex"]

    loose = rl.rank(_frame(), target_flex=12.0, target_frac=0.8, flex_tol=10.0)
    # loosen it and the 6 N flex miss (0.6 scaled) ties the kick miss (0.6),
    # so far_flex climbs level with far_kick instead of trailing far behind
    assert loose.iloc[0]["fingerprint"] == "near"
    assert set(loose["fingerprint"].iloc[1:]) == {"far_flex", "far_kick"}
    assert loose.set_index("fingerprint").loc["far_flex", "dist"] == pytest.approx(
        math.hypot(0.6, 0.01), abs=1e-9
    )


def test_drop_reasons_breaks_down_by_cause():
    reasons = rl.drop_reasons(_frame())
    assert reasons == {"errored": 1, "no_90_deg": 1, "kick_warning": 1}


def test_rank_on_an_all_error_frame_is_empty_not_a_keyerror():
    allerr = pd.DataFrame(
        [("a", "error"), ("b", "error")], columns=["fingerprint", "status"]
    )
    ranked = rl.rank(allerr, target_flex=10.0, target_frac=0.5, flex_tol=2.0)
    assert ranked.empty


def test_main_exits_cleanly_on_an_all_error_frame(tmp_path):
    parquet = tmp_path / "results.parquet"
    pd.DataFrame(
        [("a", "error"), ("b", "error")], columns=["fingerprint", "status"]
    ).to_parquet(parquet, index=False)
    with pytest.raises(SystemExit, match="no rankable rows"):
        rl.main([str(parquet), "--target-flex", "10", "--target-kick", "mid"])


def test_main_writes_ranked_csv_and_svg(tmp_path):
    parquet = tmp_path / "results.parquet"
    _frame().to_parquet(parquet, index=False)
    rc = rl.main(
        [str(parquet), "--target-flex", "12", "--target-kick", "tip"]
    )
    assert rc == 0
    out = tmp_path / "ranked.csv"
    assert out.is_file()
    assert (tmp_path / "ranked.svg").is_file()
    got = pd.read_csv(out)
    assert list(got["fingerprint"]) == ["near", "far_kick", "far_flex"]
    assert "dist" in got.columns
