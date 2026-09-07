"""``compfea-build``: STEP + ply book + materials -> a deck and a laminate report.

The design front door. Three files you author, no flags describing the laminate::

    compfea-build --step test_fin_2.step \\
                  --plybook cases/fin_zoned/plybook.csv \\
                  --materials materials/generic.csv \\
                  --long-axis x --out results/fin_zoned

Writes ``deck.inp`` plus four reports, and solves nothing unless asked:

``zone_report.csv``
    Element count, area and extent per mesh zone. Compare against the CAD
    before spending a solve -- a zone you drew and a zone the deck carries are
    different objects.
``stack_table.csv``
    The resolved stack per zone, ply by ply, with z stations. Ply 1 is the -z
    ply.
``laminate_abd.csv``
    A, B, D, flexural moduli and the symmetric/balanced/bend-twist flags. These
    are the numbers to put beside ANSYS ACP's, and they need no solver.
``build.json``
    Inputs, hashes and the warnings anything raised.

``--long-axis`` has no default, for the reason ``layup.py`` gives: a 0-degree
ply runs along the part's long axis, this repo's cases disagree about which
global axis that is, and a default is silently wrong for half of them.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from compfea.abd import laminate_summary
from compfea.geometry import Mesh
from compfea.layup import LONG_AXES, coverages_from_mesh, mesh_elsets_for_stacks
from compfea.materials import load_materials
from compfea.plybook import (
    load_plybook,
    resolve,
    stack_table,
    unused_zones,
)
from compfea.step_mesh import mesh_step, zone_report


def write_csv(path: Path, rows: list[dict]) -> None:
    """Rows -> CSV. Column order is the first row's key order."""
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def layup_fingerprint(rows: list[dict]) -> str:
    """Hash of the **resolved** stacks, for caching.

    Deliberately not the ply book's path or mtime: a moved file must still hit
    cache, and an edited thickness must miss it. The material name alone is not
    enough either -- editing a modulus under a stable name once served stale
    rows, which is why sweep.py fingerprints all ten constants.
    """
    payload = "|".join(
        f"{r['zone']}:{r['ply']}:{r['material']}:"
        f"{float(r['angle_deg']):.6f}:{float(r['thickness_mm']):.9f}"
        for r in rows
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def build(
    *,
    step: Path,
    plybook: Path,
    materials: Path,
    long_axis: str,
    out: Path,
    size_mm: float = 16.0,
    clamp_coverage: str | None = "HEAL",
) -> dict:
    """Mesh, resolve, report. Returns the manifest that lands in build.json."""
    if long_axis not in LONG_AXES:
        raise ValueError(f"long_axis must be one of {LONG_AXES}, not {long_axis!r}")
    out.mkdir(parents=True, exist_ok=True)

    library = load_materials(materials)
    warnings = list(load_materials.warnings)
    rows = load_plybook(plybook)

    raw = mesh_step(step, size_mm=size_mm, clamp_coverage=clamp_coverage)
    zones = zone_report(raw)

    coverages = coverages_from_mesh(raw.elsets)
    layup, stacks = resolve(
        rows, library, coverages, long_axis=long_axis, elsets=raw.elsets
    )
    bare = unused_zones(rows, raw.elsets)
    if bare:
        warnings.append(
            f"mesh zone(s) {bare} carry no ply of their own. Legal if they are "
            "nested inside a zone whose plies already cover them, and exactly "
            "what a zone forgotten in the ply book looks like -- check the CAD."
        )

    mesh = Mesh(
        nodes=raw.nodes,
        elements=raw.elements,
        nsets=raw.nsets,
        elsets=mesh_elsets_for_stacks(stacks, all_elements=tuple(raw.elements)),
        heading=raw.heading,
    )

    table = stack_table(layup)
    abd_rows = laminate_summary(
        layup, {name: rec.kind for name, rec in library.items()}
    )
    write_csv(out / "zone_report.csv", zones)
    write_csv(out / "stack_table.csv", table)
    write_csv(out / "laminate_abd.csv", abd_rows)
    (out / "mesh.inp").write_text(mesh.to_inp() + "\n")
    (out / "layup.inp").write_text(layup.to_inp() + "\n")

    manifest = {
        "step": str(step),
        "plybook": str(plybook),
        "materials": str(materials),
        "long_axis": long_axis,
        "size_mm": size_mm,
        "elements": len(mesh.elements),
        "nodes": len(mesh.nodes),
        "mesh_zones": [r["zone"] for r in zones],
        "n_plies": len(rows),
        "n_stacks": len(layup.zones),
        "materials_used": [m.name for m in layup.materials],
        "layup_fingerprint": layup_fingerprint(table),
        "warnings": warnings,
    }
    (out / "build.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def _report(manifest: dict, out: Path) -> None:
    zones = list(csv.DictReader((out / "zone_report.csv").read_text().splitlines()))
    abd_rows = list(csv.DictReader((out / "laminate_abd.csv").read_text().splitlines()))
    print(
        f"{manifest['elements']} elements, {manifest['nodes']} nodes, "
        f"{manifest['n_plies']} plies -> {manifest['n_stacks']} stacks, "
        f"materials {manifest['materials_used']}"
    )
    print(f"\n{'zone':12s} {'elems':>6s} {'area_mm2':>12s} {'x_min':>9s} {'x_max':>9s}")
    for row in zones:
        print(
            f"{row['zone']:12s} {int(row['elements']):6d} "
            f"{float(row['area_mm2']):12.3f} {float(row['x_min_mm']):9.2f} "
            f"{float(row['x_max_mm']):9.2f}"
        )
    print(
        f"\n{'stack':12s} {'t_mm':>7s} {'D11':>12s} {'D16':>12s} "
        f"{'sym':>5s} {'bal':>5s}  layup"
    )
    for row in abd_rows:
        print(
            f"{row['zone']:12s} {float(row['thickness_mm']):7.3f} "
            f"{float(row['d11']):12.3f} {float(row['d16']):12.3f} "
            f"{row['symmetric']:>5s} {row['balanced']:>5s}  {row['stack']}"
        )
    for warning in manifest["warnings"]:
        print(f"\nwarning: {warning}")
    print(f"\nwrote {out}/  (fingerprint {manifest['layup_fingerprint']})")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="compfea-build",
        description="STEP + ply book + materials -> deck and laminate report",
    )
    parser.add_argument("--step", type=Path, required=True)
    parser.add_argument("--plybook", type=Path, required=True)
    parser.add_argument("--materials", type=Path, required=True)
    parser.add_argument(
        "--long-axis",
        required=True,
        choices=LONG_AXES,
        help="in-plane axis a 0-degree ply runs along; no default on purpose",
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--size-mm", type=float, default=16.0)
    parser.add_argument(
        "--clamp-coverage",
        default="HEAL",
        help="zone whose nodes form fixed_end; empty string for the min-x edge",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = build(
        step=args.step,
        plybook=args.plybook,
        materials=args.materials,
        long_axis=args.long_axis,
        out=args.out,
        size_mm=args.size_mm,
        clamp_coverage=args.clamp_coverage or None,
    )
    _report(manifest, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
