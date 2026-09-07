"""Ply book: an ACP-style CSV of plies in global stacking order -> a ``Layup``.

One row is one ply. The ``ply`` column is the **global** stacking index and
**ply 1 is the -z ply**, matching the first line of a ``*SHELL SECTION,
COMPOSITE`` card. ``zone`` names a coverage mask -- a mesh ELSET -- and an
element's stack is the subsequence of plies whose zone covers it::

    ply,zone,material,angle_deg,thickness_mm,kind
    1,FULL,cf_ud,0,0.125,UD
    2,FULL,cf_ud,45,0.125,UD
    3,SPAR,cf_ud,0,0.250,UD
    4,TIP,cf_woven,45,0.200,woven
    5,FULL,cf_ud,-45,0.125,UD

An element under ``FULL`` alone carries plies 1, 2, 5; one also under ``SPAR``
carries 1, 2, 3, 5. That is why the format is a ply book and not one stack per
zone: a per-zone table lets neighbouring zones disagree about which ply is which,
and a ply that stops and restarts across a boundary is not a laminate anybody
can lay up. Here a drop is a subsequence, so continuity holds by construction.

``kind`` is optional and is a **cross-check**, not an input: when present it must
match the library's kind for that material. It catches naming the UD card where
the woven one was meant, which is otherwise a silent 60% stiffness error.

Angles are degrees, thicknesses mm, per ``CLAUDE.md``. Nothing here has units to
convert -- the materials file is where that happens.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from compfea.layup import (
    Layup,
    Ply,
    ZoneLayup,
    canonical_angle,
    layup_from_coverage,
)
from compfea.materials import MaterialRecord
from compfea.step_mesh import _sanitize_elset

#: Zone value meaning "every element", i.e. Ply.coverage is None.
ALL_ZONES = ("all", "*", "")

REQUIRED_COLUMNS = ("ply", "zone", "material", "angle_deg", "thickness_mm")


class PlyBookError(ValueError):
    """A ply book that cannot be turned into a laminate anyone intended."""


@dataclass(frozen=True)
class PlyRow:
    """One row of the ply book, validated but not yet bound to a mesh."""

    ply: int
    zone: str | None
    material: str
    angle_deg: float
    thickness_mm: float
    kind: str | None = None


def load_plybook(path: str | Path) -> tuple[PlyRow, ...]:
    """Read a ply book CSV, in stacking order, ply 1 first (the -z ply)."""
    path = Path(path)
    body = [
        ln for ln in path.read_text().splitlines() if not ln.strip().startswith("#")
    ]
    reader = csv.DictReader(body)
    if reader.fieldnames is None:
        raise PlyBookError(f"{path}: no header row")
    header = [f.strip().lower() for f in reader.fieldnames]
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    if missing:
        raise PlyBookError(f"{path}: missing required columns {missing}")

    rows: list[PlyRow] = []
    for lineno, raw in enumerate(reader, start=2):
        row = {
            (k.strip().lower() if k else ""): (v or "").strip()
            for k, v in raw.items()
            if k is not None
        }
        if not any(row.get(c) for c in REQUIRED_COLUMNS):
            continue  # blank line
        where = f"{path}:{lineno}"
        try:
            index = int(row["ply"])
        except ValueError:
            raise PlyBookError(
                f"{where}: ply = {row['ply']!r} is not an integer"
            ) from None
        material = row["material"]
        if not material:
            raise PlyBookError(f"{where}: material is required")
        try:
            angle = canonical_angle(float(row["angle_deg"]))
        except ValueError as exc:
            raise PlyBookError(f"{where}: angle_deg -- {exc}") from None
        try:
            thickness = float(row["thickness_mm"])
        except ValueError:
            raise PlyBookError(
                f"{where}: thickness_mm = {row['thickness_mm']!r} is not a number"
            ) from None
        if not math.isfinite(thickness) or thickness <= 0.0:
            raise PlyBookError(
                f"{where}: thickness_mm must be finite and > 0, got {thickness!r}"
            )
        raw_zone = row["zone"]
        zone = None if raw_zone.lower() in ALL_ZONES else _sanitize_elset(raw_zone)
        kind = (row.get("kind") or "").strip().lower() or None
        rows.append(PlyRow(index, zone, material, angle, thickness, kind))

    if not rows:
        raise PlyBookError(f"{path}: no ply rows")

    seen = [r.ply for r in rows]
    expected = list(range(1, len(rows) + 1))
    if sorted(seen) != expected:
        duplicates = sorted({p for p in seen if seen.count(p) > 1})
        raise PlyBookError(
            f"{path}: ply numbers must be 1..{len(rows)} with no gaps or "
            f"duplicates, got {sorted(seen)}"
            + (f" (duplicated: {duplicates})" if duplicates else "")
            + ". A gap means a row was deleted and the stack is no longer the "
            "one that was drawn."
        )
    return tuple(sorted(rows, key=lambda r: r.ply))


def check_against_library(
    rows: tuple[PlyRow, ...], library: dict[str, MaterialRecord]
) -> None:
    """Every ply's material must exist, and any stated ``kind`` must agree."""
    unknown = sorted({r.material for r in rows} - set(library))
    if unknown:
        raise PlyBookError(
            f"ply material(s) {unknown} are not in the materials library "
            f"{sorted(library)}"
        )
    for row in rows:
        if row.kind is None:
            continue
        want = library[row.material].kind
        if row.kind != want:
            raise PlyBookError(
                f"ply {row.ply}: kind={row.kind!r} but material "
                f"{row.material!r} is {want!r} in the library. One of the two "
                "is wrong, and a UD card where a woven one was meant is a "
                "silent stiffness error, not a solver failure."
            )


def check_against_mesh(
    rows: tuple[PlyRow, ...],
    elsets: dict[str, tuple[int, ...]],
    *,
    whole: str = "blade",
) -> None:
    """Every zone a ply names must exist in the mesh.

    Only one direction is an error. A ply naming a zone the mesh does not have
    lands nowhere and is refused here.

    The other direction is **not** an error, and assuming it was is a mistake
    worth recording: on ``test_fin_2.step``, ``QUARTER`` carries no ply of its
    own, and it does not need one -- it is a nested region whose elements are
    covered by the ``FULL`` and ``HALF`` plies above them. Demanding that every
    named zone appear in the book would refuse a perfectly good laminate.

    What actually has to hold is that no element ends up with an empty stack,
    and that is checked where it can be seen: ``layup_from_coverage`` refuses an
    element matching no plies, and ``step_mesh`` refuses an element on a tile no
    shell covers. Use ``unused_zones`` to report the rest.
    """
    named = {name for name in elsets if name != whole}
    unknown = sorted({r.zone for r in rows if r.zone is not None} - named)
    if unknown:
        raise PlyBookError(
            f"ply book names zone(s) {unknown} that the mesh does not have. "
            f"Mesh zones are {sorted(named)} (names are sanitized for ccx, so "
            "'3_4ths' becomes 'z_3_4ths')."
        )


def unused_zones(
    rows: tuple[PlyRow, ...],
    elsets: dict[str, tuple[int, ...]],
    *,
    whole: str = "blade",
) -> list[str]:
    """Mesh zones no ply names -- informational, not an error.

    Legitimate for a nested region covered by outer plies. Worth printing,
    because it is also what a zone forgotten in the CSV looks like, and only
    the person who drew the CAD can tell the two apart.

    A whole-part ply (``zone=ALL``) does **not** suppress this. An earlier
    version returned early whenever one existed, on the reasoning that no zone
    could then be bare -- true, and beside the point: a zone forgotten in the
    ply book is still forgotten, it just is not empty. One ``ALL`` row silenced
    the warning for every zone in the file.
    """
    return sorted(
        {name for name in elsets if name != whole}
        - {r.zone for r in rows if r.zone is not None}
    )


def plies_from_rows(rows: tuple[PlyRow, ...]) -> tuple[Ply, ...]:
    """Ply book rows -> ``layup.Ply`` objects, stacking order preserved."""
    return tuple(
        Ply(r.thickness_mm, r.angle_deg, r.material, r.zone) for r in rows
    )


def resolve(
    rows: tuple[PlyRow, ...],
    library: dict[str, MaterialRecord],
    element_coverages: dict[int, frozenset[str]],
    *,
    long_axis: str,
    elsets: dict[str, tuple[int, ...]] | None = None,
) -> tuple[Layup, dict[str, tuple[int, ...]]]:
    """Ply book + materials + mesh coverage -> ``(Layup, exclusive elsets)``.

    Pass ``elsets`` (the mesh's, including ``blade``) to get the both-ways zone
    check; without it only the plies are validated.
    """
    check_against_library(rows, library)
    if elsets is not None:
        check_against_mesh(rows, elsets)
    used = sorted({r.material for r in rows})
    return layup_from_coverage(
        plies_from_rows(rows),
        element_coverages,
        long_axis=long_axis,
        materials=tuple(library[n].constants for n in used),
    )


def stack_table(layup: Layup) -> list[dict[str, object]]:
    """Resolved stacks, one row per ply per zone, with z stations.

    ``z`` is measured from the laminate mid-surface, which is where ccx puts it:
    the first ply line is the -z ply, so ``z_bot`` of ply 1 is ``-t/2``.
    """
    out: list[dict[str, object]] = []
    for zone in layup.zones:
        total = zone.thickness
        z = -0.5 * total
        for index, ply in enumerate(zone.plies, start=1):
            out.append(
                {
                    "zone": zone.elset,
                    "ply": index,
                    "material": ply.material,
                    "angle_deg": ply.angle_deg,
                    "thickness_mm": ply.thickness,
                    "z_bot_mm": z,
                    "z_top_mm": z + ply.thickness,
                }
            )
            z += ply.thickness
    return out


def stack_label(zone: ZoneLayup) -> str:
    """Compact ``[0/90/45]`` label for a zone's stack, bottom to top."""
    return "[" + "/".join(f"{p.angle_deg:g}" for p in zone.plies) + "]"
