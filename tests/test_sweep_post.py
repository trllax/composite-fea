"""sweep_post.py: ranking rules, schema drift, and degenerate sweeps.

Synthetic frames, following ``tests/test_post.py`` -- no ccx in the loop. The
two real sweeps on disk are exercised at the end behind a skipif, since
``results/`` is gitignored.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from compfea.sweep_post import (
    RATIO_COL,
    classify,
    deviation_note,
    label_axes,
    load_results,
    ok_rows,
    rank_designs,
    report_angles,
    sweep_post,
    varying_axes,
)


def _design_frame() -> pd.DataFrame:
    """Two designs differing in ply count; f_180 order is not the ratio order."""
    return pd.DataFrame(
        [
            {
                "stack": "[0/90]s", "fiber": "ud", "ply_mm": 0.1, "angles": "0/90",
                "zone_pairs": "1", "n_plies_root": 4, "n_plies_tip": 4,
                "thickness_root_mm": 0.4, "thickness_tip_mm": 0.4,
                "cache_key": "aaaa", "status": "ok",
                "f_90": 1.0, "f_180": 3.0, RATIO_COL: 3.0,
                "linearity_dev_f_90": 0.005, "linearity_dev_theta_f_90": 90.0,
                "linearity_dev_f_180": 0.02, "linearity_dev_theta_f_180": 175.0,
                "max_linearity_dev": 0.02,
            },
            {
                "stack": "[0/90]2s", "fiber": "ud", "ply_mm": 0.1, "angles": "0/90",
                "zone_pairs": "2", "n_plies_root": 8, "n_plies_tip": 8,
                "thickness_root_mm": 0.8, "thickness_tip_mm": 0.8,
                "cache_key": "bbbb", "status": "ok",
                "f_90": 2.0, "f_180": 5.0, RATIO_COL: 2.5,
                "linearity_dev_f_90": 0.004, "linearity_dev_theta_f_90": 90.0,
                "linearity_dev_f_180": 0.01, "linearity_dev_theta_f_180": 175.0,
                "max_linearity_dev": 0.01,
            },
        ]
    )


def test_ranking_never_uses_f_ratio_180_90():
    """Load-bearing. The ratio is ~2 by construction; ranking on it is noise.

    The rows are built so the ratio order (A then B) is the reverse of both
    objective orders (B then A), so a ranking that silently fell back to the
    ratio would be caught rather than coincidentally agreeing.
    """
    frame = _design_frame()
    assert list(frame.sort_values(RATIO_COL, ascending=False)["stack"]) == [
        "[0/90]s", "[0/90]2s",
    ]
    for by in ("f_90", "f_180"):
        assert list(rank_designs(frame, by=by)["stack"]) == ["[0/90]2s", "[0/90]s"]
    with pytest.raises(ValueError, match="diagnostic, not an objective"):
        rank_designs(frame, by=RATIO_COL)


def test_the_written_ranking_is_in_objective_order_end_to_end(tmp_path):
    """The unit guard covers rank_designs; this covers what actually lands.

    ``_design_frame`` is built so the ratio order is the reverse of both
    objective orders, so a ranking that fell back to ``f_ratio_180_90``
    anywhere between the parquet and the CSV shows up here as a reversal.
    """
    run = tmp_path / "run"
    run.mkdir()
    _design_frame().to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    ranked = pd.read_csv(summary["files"]["rank_csv"])
    assert list(ranked["stack"]) == ["[0/90]2s", "[0/90]s"]
    assert summary["ranked_on"] == "f_180"
    assert summary["best_f_180"]["design"] == "[0/90]2s"


def test_a_run_with_no_linearity_data_is_not_reported_as_clean(tmp_path):
    """Absent is not the same as measured-and-fine.

    ``NaN > warn`` is False, so a plain comparison stamps "not suspect" on a
    schema that never carried the deviation columns at all.
    """
    frame = pd.DataFrame([
        {"stack": "a", "cache_key": "k1", "status": "ok", "f_180": 3.0},
        {"stack": "b", "cache_key": "k2", "status": "ok", "f_180": 5.0},
    ])
    run = tmp_path / "run"
    run.mkdir()
    frame.to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    assert summary["best_f_180"]["suspect"] is None      # unknown, not False
    ranked = pd.read_csv(summary["files"]["rank_csv"])
    assert ranked["suspect_f_180"].isna().all()


def test_the_summary_json_is_strict_json_even_with_a_nan_force(tmp_path):
    frame = _design_frame()
    frame.loc[0, "f_90"] = float("nan")     # design solved to 180 but not 90
    run = tmp_path / "run"
    run.mkdir()
    frame.to_parquet(run / "results.parquet", index=False)
    sweep_post(run)
    text = (run / "sweep_post.json").read_text()
    assert "NaN" not in text
    assert json.loads(text)["n_ok"] == 2  # raises if not strict JSON


def test_the_label_uses_the_fewest_axes_that_tell_designs_apart():
    """A ply-count sweep moves six correlated columns; the label needs one."""
    frame = _design_frame()
    assert varying_axes(frame) == [
        "stack", "zone_pairs", "n_plies_root", "n_plies_tip",
        "thickness_root_mm", "thickness_tip_mm",
    ]
    assert label_axes(frame) == ["stack"]


def test_a_label_value_is_not_printed_with_float_noise(tmp_path):
    frame = pd.DataFrame([
        {"stack": "a", "thickness_root_mm": 1.2000000000000002,
         "cache_key": "k1", "status": "ok", "f_180": 3.0},
        {"stack": "a", "thickness_root_mm": 0.8,
         "cache_key": "k2", "status": "ok", "f_180": 5.0},
    ])
    run = tmp_path / "run"
    run.mkdir()
    frame.to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    assert summary["best_f_180"]["design"] == "thickness_root=0.8"


def test_report_angles_are_discovered_not_hardcoded():
    frame = pd.DataFrame(
        [{"f_45": 1.0, "f_22.5": 0.5, "f_120": 2.0, RATIO_COL: 2.0, "status": "ok"}]
    )
    assert report_angles(frame) == [22.5, 45.0, 120.0]
    # the ratio column ends in digits but is not an angle
    assert 180.0 not in report_angles(frame)


def test_error_rows_are_dropped_and_reported(tmp_path):
    frame = _design_frame()
    frame = pd.concat(
        [
            frame,
            pd.DataFrame([{
                "stack": "[45/-45]s", "cache_key": "cccc", "fiber": "ud",
                "ply_mm": 0.1, "angles": "45/-45", "zone_pairs": "1",
                "status": "error", "error": "NotConverged: stalled at inc 812",
            }]),
        ],
        ignore_index=True,
    )
    ok, bad = ok_rows(frame)
    assert len(ok) == 2 and len(bad) == 1

    run = tmp_path / "run"
    run.mkdir()
    frame.to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    assert summary["n_ok"] == 2 and summary["n_error"] == 1
    assert "NotConverged" in summary["errors"][0]["error"]
    # the failed design must not appear in the ranking
    ranked = pd.read_csv(summary["files"]["rank_csv"])
    assert "[45/-45]s" not in set(ranked["stack"])


def test_a_run_where_every_design_failed_exits_nonzero(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    pd.DataFrame([
        {"stack": "a", "status": "error", "error": "CcxFailed: exit 201"},
    ]).to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    assert summary["n_ok"] == 0
    assert "every design failed" in summary["note"]


def test_missing_linearity_columns_degrade_instead_of_keyerror(tmp_path):
    """The 18-column schema of results/ubend-demo-n1, which predates them."""
    frame = pd.DataFrame([
        {"n_pairs": 1, "n_plies": 4, "ply_mm": 0.1, "stack": "[0/90]s",
         "thickness_mm": 0.4, "cache_key": "old1", "status": "ok",
         "f_90": 1.76, "f_180": 3.45, RATIO_COL: 1.96},
        {"n_pairs": 2, "n_plies": 8, "ply_mm": 0.1, "stack": "[0/90]2s",
         "thickness_mm": 0.8, "cache_key": "old2", "status": "ok",
         "f_90": 3.5, "f_180": 6.9, RATIO_COL: 1.97},
    ])
    run = tmp_path / "run"
    run.mkdir()
    frame.to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    assert summary["mode"] == "design"
    assert "linearity_dev_f_90" in summary["missing_columns"]
    assert "max_linearity_dev" in summary["missing_columns"]
    assert Path(summary["files"]["rank_svg"]).is_file()


def test_missing_design_columns_fall_back_to_stack(tmp_path):
    frame = pd.DataFrame([
        {"stack": "[0/90]s", "cache_key": "k1", "status": "ok", "f_180": 3.0},
        {"stack": "[0/90]2s", "cache_key": "k2", "status": "ok", "f_180": 5.0},
    ])
    run = tmp_path / "run"
    run.mkdir()
    frame.to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    # stack alone tells the designs apart, so it is the whole label
    assert summary["best_f_180"]["design"] == "[0/90]2s"


def test_a_design_with_no_identifying_axes_falls_back_to_cache_key(tmp_path):
    """Nothing but the hash distinguishes these rows, so the label is the hash."""
    frame = pd.DataFrame([
        {"cache_key": "k1", "status": "ok", "f_180": 3.0},
        {"cache_key": "k2", "status": "ok", "f_180": 5.0},
    ])
    run = tmp_path / "run"
    run.mkdir()
    frame.to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    # no design axis varies -> read as a calibration, and labelled by hash
    assert summary["mode"] == "calibration"
    assert summary["best_f_180"]["design"] == "k2"


def test_a_run_that_varies_only_static_line_is_a_calibration_not_a_ranking(tmp_path):
    """The results/incr-calib shape: same design, five solver settings."""
    rows = []
    for i, (line, incs, secs) in enumerate([
        ("0.1, 1.0, 1.E-8, 0.5", 223, 25.2),
        ("0.05, 1.0, 1.E-8, 0.25", 288, 28.2),
        ("0.02, 1.0, 1.E-8, 0.1", 506, 42.3),
        ("0.01, 1.0, 1.E-8, 0.05", 865, 66.3),
        ("0.005, 1.0, 1.E-8, 0.02", 1908, 137.0),
    ]):
        rows.append({
            "stack": "[0/90]s", "fiber": "ud", "ply_mm": 0.1, "angles": "0/90",
            "zone_pairs": "1", "n_plies_root": 4, "n_plies_tip": 4,
            "thickness_root_mm": 0.4, "thickness_tip_mm": 0.4,
            "static_line": line, "increments": incs, "wall_time_s": secs,
            "cache_key": f"k{i}", "status": "ok",
            "f_90": 1.759634, "f_180": 3.451304, RATIO_COL: 1.961376,
        })
    run = tmp_path / "run"
    run.mkdir()
    pd.DataFrame(rows).to_parquet(run / "results.parquet", index=False)

    ok, _ = ok_rows(pd.DataFrame(rows))
    assert varying_axes(ok) == []
    assert classify(ok) == "calibration"

    summary = sweep_post(run)
    assert summary["mode"] == "calibration"
    assert "calibration" in summary["note"]
    assert not (run / "sweep_post_rank.svg").exists()
    table = pd.read_csv(summary["files"]["calibration_csv"])
    # identical forces across the grid: the deviation column is exactly zero
    assert table["dev_f_180"].abs().max() == 0.0
    assert table["speedup_vs_finest"].max() == pytest.approx(137.0 / 25.2)


def test_a_single_row_run_still_writes_a_summary(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _design_frame().head(1).to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    assert summary["mode"] == "single"
    assert not (run / "sweep_post_rank.svg").exists()
    assert json.loads((run / "sweep_post.json").read_text())["n_ok"] == 1


def test_a_high_linearity_dev_design_is_flagged_not_dropped(tmp_path):
    frame = _design_frame()
    frame.loc[0, "linearity_dev_f_180"] = 0.4
    run = tmp_path / "run"
    run.mkdir()
    frame.to_parquet(run / "results.parquet", index=False)
    summary = sweep_post(run)
    ranked = pd.read_csv(summary["files"]["rank_csv"])
    assert len(ranked) == 2                      # still ranked, not dropped
    flagged = ranked.loc[ranked["stack"] == "[0/90]s", "suspect_f_180"].iloc[0]
    assert bool(flagged) is True
    assert bool(ranked.loc[ranked["stack"] == "[0/90]2s",
                           "suspect_f_180"].iloc[0]) is False


def test_the_axis_label_names_the_angle_the_deviation_came_from():
    """On a 5-degree path the f_180 deviation is measured at 175, not 180."""
    frame = _design_frame()
    assert deviation_note(frame, 180.0) == "  (dev at 175°)"
    assert deviation_note(frame, 90.0) == ""


def test_a_run_dir_without_results_is_an_error(tmp_path):
    empty = tmp_path / "nothing"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        load_results(empty)


# --- against the real sweeps on disk ----------------------------------------

CAL = Path("results/incr-calib/results.parquet")
DEMO = Path("results/ubend-demo-n1/results.parquet")


@pytest.mark.skipif(not CAL.is_file(), reason=f"{CAL} not on disk")
def test_reads_the_real_incr_calib_parquet(tmp_path):
    """End to end on real output; the force matches cases/sweep_ubend/README."""
    summary = sweep_post("results/incr-calib", out_dir=tmp_path)
    assert summary["mode"] == "calibration"
    assert summary["n_ok"] == 5 and summary["n_error"] == 0
    assert summary["report_angles"] == [90.0, 180.0]
    table = pd.read_csv(summary["files"]["calibration_csv"])
    assert table["f_90"].iloc[0] == pytest.approx(1.759634, abs=5e-7)
    assert table["f_180"].iloc[0] == pytest.approx(3.451304, abs=5e-7)
    assert table["dev_f_90"].abs().max() == 0.0


@pytest.mark.skipif(not DEMO.is_file(), reason=f"{DEMO} not on disk")
def test_reads_the_real_old_schema_parquet(tmp_path):
    summary = sweep_post("results/ubend-demo-n1", out_dir=tmp_path)
    assert summary["mode"] == "single"
    assert "linearity_dev_f_90" in summary["missing_columns"]


# --- shapes from the cached .frd, when a sweep with DISP is on disk ---------

CAL_FRD = Path("results/incr-calib/cache/6418231dde201796/ccx/job.frd")
frd_only = pytest.mark.skipif(
    not CAL_FRD.is_file(), reason=f"{CAL_FRD} not on disk"
)


@frd_only
def test_designs_sharing_a_label_are_not_silently_merged(tmp_path):
    """incr-calib's five rows are one design, so every label is identical.

    Keying the shape comparison on the label alone kept only the last and wrote
    no comparison at all, without saying why.
    """
    from compfea.sweep_post import label_axes, ok_rows, rank_designs, write_shapes

    solved, _ = ok_rows(load_results("results/incr-calib"))
    ranked = rank_designs(solved, by="f_180")
    assert label_axes(solved) == []  # nothing distinguishes them but the hash
    files = write_shapes(
        Path("results/incr-calib"), ranked, tmp_path,
        n=3, angles_deg=(90.0, 180.0), axes=label_axes(solved), end_deg=180.0,
    )
    assert sum(1 for k in files if k.startswith("shape_") and k != "shape_compare") == 3
    assert "shape_compare" in files
    assert Path(files["shape_compare"]).is_file()


def _truncate_frd(src: Path, dest: Path, last_step: int) -> None:
    """Copy an frd, dropping every block from ``last_step + 1`` on."""
    out, drop = [], False
    for line in src.read_text().splitlines():
        if line.startswith("    1PSTEP"):
            parts = line.split()
            drop = len(parts) >= 4 and int(parts[3]) > last_step
        if not drop:
            out.append(line)
    dest.write_text("\n".join(out) + "\n")


@frd_only
def test_a_design_that_stops_short_is_left_out_of_the_comparison(tmp_path):
    """The discriminating case: one design has no DISP at the compared angle.

    Taking each design's own deepest pose instead of the requested one would
    draw this design's 90-degree shape on a plot titled 180 degrees -- the
    silent mislabel this repo exists to avoid. Built by truncating a real frd,
    because every design in the sweeps on disk reaches 180 and so cannot tell
    the two behaviours apart.
    """
    from compfea.sweep_post import label_axes, rank_designs, write_shapes

    run = tmp_path / "run"
    full_key, short_key = "fullfull0000000", "shortshort00000"
    for key, last in ((full_key, 36), (short_key, 18)):
        dest = run / "cache" / key / "ccx"
        dest.mkdir(parents=True)
        _truncate_frd(CAL_FRD, dest / "job.frd", last)

    ranked = rank_designs(
        pd.DataFrame([
            {"stack": "deep", "cache_key": full_key, "status": "ok", "f_180": 9.0},
            {"stack": "shallow", "cache_key": short_key, "status": "ok",
             "f_180": 4.0},
        ]),
        by="f_180",
    )
    axes = label_axes(ranked)
    files = write_shapes(
        run, ranked, tmp_path, n=2, angles_deg=(90.0, 180.0), axes=axes,
        end_deg=180.0,
    )

    # both designs get their own shape plot -- the short one at 90 only
    assert f"shape_{full_key}" in files and f"shape_{short_key}" in files
    # but only one design has 180, so there is no two-design comparison at 180
    assert "shape_compare" not in files
    assert not list(tmp_path.glob("sweep_post_shape_compare_90.svg"))


@frd_only
def test_the_comparison_is_named_for_the_requested_angle(tmp_path):
    from compfea.sweep_post import label_axes, ok_rows, rank_designs, write_shapes

    solved, _ = ok_rows(load_results("results/incr-calib"))
    ranked = rank_designs(solved, by="f_180")
    files = write_shapes(
        Path("results/incr-calib"), ranked, tmp_path,
        n=2, angles_deg=(90.0, 180.0), axes=label_axes(solved), end_deg=180.0,
    )
    assert Path(files["shape_compare"]).name == "sweep_post_shape_compare_180.svg"


@frd_only
def test_an_angle_with_no_disp_is_skipped_not_substituted():
    """There is no DISP below 90 degrees; asking for 45 must yield nothing."""
    from compfea.sweep_post import shape_curves

    curves = shape_curves(
        Path("results/incr-calib"), "6418231dde201796", (45.0, 90.0, 180.0),
        end_deg=180.0,
    )
    assert sorted(curves) == [90.0, 180.0]  # 45 dropped, not filled in
