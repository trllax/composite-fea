"""compfea-build: the three-file front door, end to end without a solver."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from compfea.build import build, build_parser, layup_fingerprint, main

ROOT = Path(__file__).resolve().parents[1]
FIN2 = ROOT / "test_fin_2.step"
PLYBOOK = ROOT / "cases" / "fin_zoned" / "plybook.csv"
GENERIC = ROOT / "materials" / "generic.csv"

pytestmark = pytest.mark.skipif(
    not FIN2.is_file(), reason="test_fin_2.step not in repo root"
)

# 40 mm keeps these tests to a couple of seconds. cases/fin_zoned itself runs
# at 16; the mesh size does not change any laminate number reported here.
FAST_SIZE_MM = 40.0


def run(tmp_path: Path, **over) -> dict:
    kwargs = dict(
        step=FIN2,
        plybook=PLYBOOK,
        materials=GENERIC,
        long_axis="x",
        out=tmp_path,
        size_mm=FAST_SIZE_MM,
    )
    kwargs.update(over)
    return build(**kwargs)


def read(path: Path) -> list[dict]:
    return list(csv.DictReader(path.read_text().splitlines()))


def test_build_writes_every_report(tmp_path):
    manifest = run(tmp_path)
    for name in (
        "zone_report.csv", "stack_table.csv", "laminate_abd.csv",
        "build.json", "mesh.inp", "layup.inp",
    ):
        assert (tmp_path / name).is_file(), name
    assert json.loads((tmp_path / "build.json").read_text()) == manifest
    assert manifest["materials_used"] == ["cfrp", "cfrp_woven"]
    assert manifest["n_plies"] == 14


def test_the_deck_carries_two_materials(tmp_path):
    """The ply book names both cards; nothing in the repo did this before."""
    run(tmp_path)
    layup = (tmp_path / "layup.inp").read_text()
    assert layup.count("*MATERIAL, NAME=") == 2
    assert ", cfrp_woven, ori_p45" in layup
    assert ", cfrp, ori_p0" in layup


def elsets_from_inp(text: str) -> dict[str, set[int]]:
    """*ELSET blocks -> {name: element ids}. blade comes from *ELEMENT."""
    out: dict[str, set[int]] = {}
    current = None
    for line in text.splitlines():
        if line.startswith("*ELSET, ELSET="):
            current = line.split("=")[-1].strip()
            out.setdefault(current, set())
        elif line.startswith("*"):
            current = None
        elif current:
            out[current].update(int(v) for v in line.split(",") if v.strip())
    return out


def element_ids_from_inp(text: str) -> set[int]:
    ids: set[int] = set()
    reading = False
    for line in text.splitlines():
        if line.startswith("*ELEMENT"):
            reading = True
        elif line.startswith("*"):
            reading = False
        elif reading and line.strip():
            ids.add(int(line.split(",")[0]))
    return ids


def test_the_cov_sections_partition_the_part(tmp_path):
    """A partition, not a cover: an element in two sections is not a laminate,
    and an element in none goes into the deck with nothing on it."""
    manifest = run(tmp_path)
    mesh = (tmp_path / "mesh.inp").read_text()
    cov = {n: e for n, e in elsets_from_inp(mesh).items() if n.startswith("cov_")}
    assert len(cov) == manifest["n_stacks"] >= 2

    members = [e for ids in cov.values() for e in ids]
    assert len(members) == len(set(members)), "cov_* sets overlap"
    assert set(members) == element_ids_from_inp(mesh), "cov_* must cover the part"
    assert len(members) == manifest["elements"]

    sections = (tmp_path / "layup.inp").read_text()
    for name in cov:
        assert f"*SHELL SECTION, COMPOSITE, ELSET={name}\n" in sections


def test_blade_still_spans_the_part(tmp_path):
    """evaluate_design reads ELSE off 'blade'; it has to keep existing."""
    manifest = run(tmp_path)
    mesh = (tmp_path / "mesh.inp").read_text()
    assert "ELSET=blade" in mesh
    assert manifest["elements"] > 100


def test_zone_report_areas_match_the_stack_table_zones(tmp_path):
    run(tmp_path)
    zones = {r["zone"] for r in read(tmp_path / "zone_report.csv")}
    assert {"FULL", "HALF", "QUARTER", "TIP", "HEAL", "z_3_4ths", "blade"} == zones
    stacks = {r["zone"] for r in read(tmp_path / "stack_table.csv")}
    assert all(s.startswith("cov_") for s in stacks)


def test_every_zone_of_the_shipped_plybook_is_symmetric_and_balanced(tmp_path):
    """The case is meant to be a good laminate; say so as a test, not a comment."""
    run(tmp_path)
    rows = read(tmp_path / "laminate_abd.csv")
    assert len(rows) >= 5
    for row in rows:
        assert row["symmetric"] == "True", f"{row['zone']} {row['stack']}"
        assert row["balanced"] == "True", f"{row['zone']} {row['stack']}"


def test_stack_table_z_stations_are_contiguous_and_centred(tmp_path):
    run(tmp_path)
    rows = read(tmp_path / "stack_table.csv")
    by_zone: dict[str, list[dict]] = {}
    for row in rows:
        by_zone.setdefault(row["zone"], []).append(row)
    for zone, plies in by_zone.items():
        plies.sort(key=lambda r: int(r["ply"]))
        total = sum(float(r["thickness_mm"]) for r in plies)
        assert float(plies[0]["z_bot_mm"]) == pytest.approx(-0.5 * total), zone
        assert float(plies[-1]["z_top_mm"]) == pytest.approx(+0.5 * total), zone
        for lower, upper in zip(plies, plies[1:], strict=False):
            assert float(lower["z_top_mm"]) == pytest.approx(float(upper["z_bot_mm"]))


def test_the_fingerprint_tracks_the_resolved_stack_not_the_file(tmp_path):
    """A moved file must hit cache; an edited thickness must miss it."""
    rows = [
        {"zone": "cov_1", "ply": 1, "material": "cfrp",
         "angle_deg": 0.0, "thickness_mm": 0.15},
        {"zone": "cov_1", "ply": 2, "material": "cfrp",
         "angle_deg": 90.0, "thickness_mm": 0.15},
    ]
    base = layup_fingerprint(rows)
    assert layup_fingerprint(list(rows)) == base

    thicker = [dict(r) for r in rows]
    thicker[0]["thickness_mm"] = 0.16
    assert layup_fingerprint(thicker) != base

    rotated = [dict(r) for r in rows]
    rotated[0]["angle_deg"] = 45.0
    assert layup_fingerprint(rotated) != base

    swapped = [dict(r) for r in rows]
    swapped[0]["material"] = "cfrp_woven"
    assert layup_fingerprint(swapped) != base

    reordered = [dict(rows[1]), dict(rows[0])]
    assert layup_fingerprint(reordered) != base, "ply order must change the key"


def test_long_axis_has_no_default(tmp_path):
    """layup.py refuses to guess; the CLI must not guess on its behalf."""
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["--step", str(FIN2), "--plybook", str(PLYBOOK),
             "--materials", str(GENERIC), "--out", str(tmp_path)]
        )
    args = parser.parse_args(
        ["--step", str(FIN2), "--plybook", str(PLYBOOK),
         "--materials", str(GENERIC), "--out", str(tmp_path), "--long-axis", "x"]
    )
    assert args.long_axis == "x"


def test_main_runs_and_reports(tmp_path, capsys):
    assert main([
        "--step", str(FIN2), "--plybook", str(PLYBOOK),
        "--materials", str(GENERIC), "--long-axis", "x",
        "--out", str(tmp_path), "--size-mm", str(FAST_SIZE_MM),
    ]) == 0
    out = capsys.readouterr().out
    assert "zone" in out and "D11" in out and "cov_1" in out
    assert "wrote" in out


def test_a_zone_with_no_ply_of_its_own_is_a_warning_not_a_failure(tmp_path):
    """QUARTER on the equivalence book: nested, covered, and perfectly legal."""
    manifest = run(
        tmp_path, plybook=ROOT / "cases" / "step_fin" / "plybook_equivalent.csv"
    )
    assert any("carry no ply of their own" in w for w in manifest["warnings"])
    assert (tmp_path / "laminate_abd.csv").is_file()
