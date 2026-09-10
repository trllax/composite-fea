"""Layup sweep: grid expansion, dispatch, caching, error rows.

The real per-design solve (``run_tipweight.py``, ~2 min) is stubbed here -- a
fake that writes a canned ``flex_kick.json``. What is under test is the sweep's
own logic: turning a designs JSON into ``TaperDesign``s, resolving + writing a
ply book per design, skipping one whose ``flex_kick.json`` is already on disk,
recording a failed design as a row rather than crashing, and assembling
``results.parquet``. The end-to-end multi-solve run on FIN_TEST_3 is a manual
step (see ``cases/fin_test_3/README.md``), not a suite test.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

from compfea.plybook_gen import TaperDesign, design_fingerprint
from compfea.shop_inventory import load_shop_inventory

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "materials" / "shop_inventory.csv"

TWILL = "hexcel_im7_twill_205"
UD = "hexcel_im2_uni_193"
BIAX = "hexcel_himax_biax_100"


def _module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


sl = _module(ROOT / "cases/fin_test_3/sweep_layups.py", "layup_sweep_under_test")


@pytest.fixture(scope="module")
def skus():
    return load_shop_inventory(INVENTORY)


BASE_DESIGN = {
    "through_zone": "FULL",
    "skins": [[TWILL, 0]],
    "core": [[UD, 0], [BIAX, 0]],
    "pads": {"HEAL": [[TWILL, 0]], "QUARTER": []},
    "pad_order": ["TIP", "MID", "QUARTER", "HEAL"],
}


def _canned_flex_kick(flex_n: float = 10.0, bucket: str = "mid", frac: float = 0.5):
    kp = {
        "s_kick_mm": frac * 640.0,
        "s_kick_frac": frac,
        "bucket": bucket,
        "s_peak_raw_mm": frac * 640.0,
        "s_rotmed_mm": 0.55 * 640.0,
        "s_rotmed_frac": 0.55,
        "plateau_frac": 0.2,
        "theta_deg": 90.0,
        "warning": None,
    }
    return {
        "flex_n": flex_n,
        "free_len_mm": 640.0,
        "headline": kp,
        "migration": {**kp, "theta_deg": 15.0},
        "migration_delta_mm": 1.0,
        "kick_bucket": bucket,
        "kick_s_mm": kp["s_kick_mm"],
        "kick_s_frac": frac,
    }


# --------------------------------------------------------------------------
# design list


def test_load_designs_expands_a_grid(tmp_path):
    spec = {
        "base": BASE_DESIGN,
        "axes": {
            "core": [[[UD, 0]], [[UD, 0], [BIAX, 0]], [[UD, 0], [UD, 0], [BIAX, 0]]],
            "pads.HEAL": [[], [[TWILL, 0]]],
        },
    }
    path = tmp_path / "grid.json"
    path.write_text(json.dumps(spec))
    designs = sl.load_designs(path)
    assert len(designs) == 6
    assert all(isinstance(d, TaperDesign) for d in designs)


def test_load_designs_list_and_single(tmp_path):
    as_list = tmp_path / "list.json"
    as_list.write_text(json.dumps([BASE_DESIGN, BASE_DESIGN]))
    assert len(sl.load_designs(as_list)) == 2

    as_single = tmp_path / "one.json"
    as_single.write_text(json.dumps(BASE_DESIGN))
    got = sl.load_designs(as_single)
    assert len(got) == 1 and isinstance(got[0], TaperDesign)


# --------------------------------------------------------------------------
# dispatch / collect / cache


def _three_designs():
    return [
        sl.design_from_dict({**BASE_DESIGN, "core": [[UD, 0]]}),
        sl.design_from_dict({**BASE_DESIGN, "core": [[UD, 0], [BIAX, 0]]}),
        sl.design_from_dict({**BASE_DESIGN, "core": [[UD, 0], [UD, 0], [BIAX, 0]]}),
    ]


def _fp_of(run_dir) -> str:
    """Cache dir name is '<fingerprint>.<solve_phash>'; recover the fingerprint."""
    return Path(run_dir).name.rsplit(".", 1)[0]


def _install_fake_solver(monkeypatch, *, calls: list, fail_fp: str | None = None):
    def fake(plybook_path, run_dir, **kw):
        fp = _fp_of(run_dir)
        calls.append(fp)
        if fail_fp is not None and fp == fail_fp:
            raise RuntimeError("solver said no")
        (Path(run_dir) / "flex_kick.json").write_text(
            json.dumps(_canned_flex_kick(flex_n=8.0 + len(calls)))
        )

    monkeypatch.setattr(sl, "_run_tipweight_solve", fake)


def _run(tmp_path, designs, **over):
    kw = dict(
        run_root=tmp_path / "run",
        jobs=2,
        materials=Path("materials/from_shop.csv"),
        inventory=INVENTORY,
        size_mm=16.0,
        uy_frac=0.6,
        bow_tip_mm=0.5,
        kick_bands=41,
        timeout_s=60.0,
    )
    kw.update(over)
    return sl.run_sweep(designs, **kw)


def test_run_sweep_dispatches_collects_and_writes_parquet(tmp_path, monkeypatch):
    calls: list[str] = []
    _install_fake_solver(monkeypatch, calls=calls)
    designs = _three_designs()

    frame = _run(tmp_path, designs)

    assert len(frame) == 3
    assert (frame["status"] == "ok").all()
    for col in (
        "fingerprint", "flex_n", "kick_bucket", "kick_s_frac", "kick_s_mm",
        "migration_delta_mm", "n_plies_total", "through_n_plies",
        "through_thickness_mm", "core",
    ):
        assert col in frame.columns, col
    assert len(calls) == 3  # one solve per design

    parquet = tmp_path / "run" / "results.parquet"
    assert parquet.is_file()
    status = json.loads((tmp_path / "run" / "status.json").read_text())
    assert status["state"] == "done" and status["ok"] == 3

    # n_plies_total is the whole book: skin*2 + core*2 (mirrored) + HEAL pad*2.
    by_core = frame.set_index("core")["n_plies_total"].to_dict()
    assert by_core[f"{UD}@0"] == 2 + 2 * 1 + 2
    assert by_core[f"{UD}@0/{BIAX}@0"] == 2 + 2 * 2 + 2
    # through_n_plies is the FULL stack only: skin*2 + core*2.
    by_core_n = frame.set_index("core")["through_n_plies"].to_dict()
    assert by_core_n[f"{UD}@0"] == 2 + 2 * 1
    # through_thickness is the FULL stack only: skin*2 + core*2, from inventory
    # (as-laid cured thickness = gsm/1000 * CURED_PLY_FACTOR).
    by_core_t = frame.set_index("core")["through_thickness_mm"].to_dict()
    assert by_core_t[f"{UD}@0"] == pytest.approx(2 * 0.2255 + 2 * 0.2123)


def test_cache_hit_skips_the_solve(tmp_path, monkeypatch):
    calls: list[str] = []
    _install_fake_solver(monkeypatch, calls=calls)
    designs = _three_designs()
    skus_local = load_shop_inventory(INVENTORY)

    # Pre-populate the cache dir for the middle design.
    cached_fp = design_fingerprint(designs[1], skus=skus_local)
    phash = sl.solve_phash(
        materials=Path("materials/from_shop.csv"), size_mm=16.0,
        uy_frac=0.6, bow_tip_mm=0.5, kick_bands=41,
    )
    cache_dir = tmp_path / "run" / "cache" / f"{cached_fp}.{phash}"
    cache_dir.mkdir(parents=True)
    (cache_dir / "flex_kick.json").write_text(
        json.dumps(_canned_flex_kick(flex_n=99.0, bucket="tip", frac=0.8))
    )

    frame = _run(tmp_path, designs)

    assert len(frame) == 3
    assert cached_fp not in calls  # never re-solved
    assert len(calls) == 2
    cached_row = frame[frame["fingerprint"] == cached_fp].iloc[0]
    assert bool(cached_row["cached"]) is True
    assert cached_row["flex_n"] == 99.0
    assert cached_row["kick_bucket"] == "tip"


def test_a_failing_solve_becomes_an_error_row(tmp_path, monkeypatch):
    designs = _three_designs()
    skus_local = load_shop_inventory(INVENTORY)
    bad_fp = design_fingerprint(designs[0], skus=skus_local)

    calls: list[str] = []
    _install_fake_solver(monkeypatch, calls=calls, fail_fp=bad_fp)

    frame = _run(tmp_path, designs)

    assert len(frame) == 3
    bad = frame[frame["fingerprint"] == bad_fp].iloc[0]
    assert bad["status"] == "error"
    assert "solver said no" in bad["error"]
    assert (frame["status"] == "ok").sum() == 2
    status = json.loads((tmp_path / "run" / "status.json").read_text())
    assert status["error"] == 1


def test_designs_that_resolve_alike_are_collapsed_to_one_solve(tmp_path, monkeypatch):
    calls: list[str] = []
    _install_fake_solver(monkeypatch, calls=calls)
    d = sl.design_from_dict({**BASE_DESIGN, "core": [[UD, 0]]})

    frame = _run(tmp_path, [d, d, d])

    assert len(calls) == 1  # one laminate -> one ccx run, no cache-dir race
    assert len(frame) == 1
    status = json.loads((tmp_path / "run" / "status.json").read_text())
    assert status["n_designs"] == 3 and status["n_unique"] == 1


def test_solve_phash_moves_with_kick_bands():
    common = dict(
        materials=Path("materials/from_shop.csv"),
        size_mm=16.0, uy_frac=0.6, bow_tip_mm=0.5,
    )
    assert sl.solve_phash(kick_bands=41, **common) != sl.solve_phash(
        kick_bands=81, **common
    )
    assert sl.solve_phash(kick_bands=41, **common) == sl.solve_phash(
        kick_bands=41, **common
    )


def test_solve_phash_moves_with_twist_probe():
    common = dict(
        materials=Path("materials/from_shop.csv"),
        size_mm=16.0, uy_frac=0.6, bow_tip_mm=0.5, kick_bands=41,
    )
    # a twist-probe run must not collide with a plain run's cache
    assert sl.solve_phash(**common) != sl.solve_phash(twist_probe=True, **common)
    # nor two twist runs with different couple forces
    assert sl.solve_phash(
        twist_probe=True, twist_cload_n=5.0, **common
    ) != sl.solve_phash(twist_probe=True, twist_cload_n=10.0, **common)
    assert sl.solve_phash(twist_probe=True, **common) == sl.solve_phash(
        twist_probe=True, **common
    )
