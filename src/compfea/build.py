"""``compfea-build``: STEP + ply book + materials -> a deck and a laminate report.

The design front door. Three files you author, no flags describing the laminate::

    compfea-build --step test_fin_2.step \\
                  --plybook cases/fin_zoned/plybook.csv \\
                  --materials materials/generic.csv \\
                  --long-axis x --out results/fin_zoned

Writes the deck's two pieces -- ``mesh.inp`` and ``layup.inp`` -- plus four
reports. It does not assemble a runnable deck and it never solves: boundary
conditions and steps belong to a load case, not to a laminate.

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
from compfea.materials import MaterialRecord, load_materials
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


def layup_fingerprint(
    rows: list[dict],
    materials: dict[str, MaterialRecord] | None = None,
    long_axis: str = "",
) -> str:
    """Hash of the **resolved laminate** -- stacks, cards, and the long axis.

    Covers everything that changes ``layup.to_inp()`` and nothing else. It is
    **not** a deck key: the mesh is not in it, so a cache keyed on this alone
    would serve a result computed on a different STEP or a different
    ``size_mm``. ``build.json`` records those separately; a cache wants both.

    Deliberately not the ply book's path or mtime: a moved file must still hit
    cache, and an edited thickness must miss it.

    The material *name* is not enough, and this function shipped believing it
    was while its own docstring said otherwise. Editing a modulus under a stable
    name once served stale rows in ``sweep.py``, which is why
    ``_material_fingerprint`` there interpolates all ten constants. Pass
    ``materials`` and it does the same here. ``long_axis`` is in the key too:
    it decides which global axis a 0-degree ply runs along, so the same ply book
    under ``x`` and under ``y`` are different laminates that would otherwise
    share a key.

    Both are optional only so the stack-shape behaviour can be tested on its
    own; ``build`` always passes them.
    """
    parts = [f"axis={long_axis}"]
    parts += [
        f"{r['zone']}:{r['ply']}:{r['material']}:"
        f"{float(r['angle_deg']):.6f}:{float(r['thickness_mm']):.9f}"
        for r in rows
    ]
    if materials is not None:
        for name in sorted({str(r["material"]) for r in rows}):
            ec = materials[name].constants
            parts.append(
                f"mat={name}:{ec.e1!r}:{ec.e2!r}:{ec.e3!r}:{ec.nu12!r}:"
                f"{ec.nu13!r}:{ec.nu23!r}:{ec.g12!r}:{ec.g13!r}:{ec.g23!r}:"
                f"{ec.density!r}"
            )
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


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

    raw = mesh_step(step, size_mm=size_mm, clamp_coverage=clamp_coverage, long_axis=long_axis)
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
    abd_rows = laminate_summary(layup)
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
        "layup_fingerprint": layup_fingerprint(table, library, long_axis),
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
