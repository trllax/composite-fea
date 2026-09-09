"""Design vector -> ply book, for a tapered laminate over nested named zones.

A ``TaperDesign`` is the compact description of a tapered laminate:

- **skins**  -- ``n`` plies per face covering the through-going zone, the
  continuous outer surface;
- **core**   -- an ordered half-stack covering the through zone, mirrored about
  the mid-plane (adds uniform bending stiffness, ``EI ~ t^3``);
- **pads**   -- per inboard / tip zone, an ordered half-stack covering just that
  zone, mirrored (moves the stiffness taper, and so the location of peak
  curvature under load).

``design_rows`` turns one into the repo's existing ``plybook.PlyRow`` tuple, so
everything downstream (``plybook.plies_from_rows`` -> ``resolve`` ->
``layup_from_coverage``) is untouched. ``write_plybook`` emits a CSV that
``plybook.load_plybook`` reads back unchanged.

**Symmetry is by construction, not enforced.** The build mirrors ``core`` and
each pad half-stack about the mid-plane and keeps ply angles unchanged across
the mirror (``[0/45/-45]_s`` = ``0/45/-45/-45/45/0`` -- *not* sign-flipped, which
would be an antisymmetric stack), so ``B`` cancels to float noise. A negated
mirror is invisible to every force check in this repo -- see ``CLAUDE.md`` on
``*ORIENTATION`` and ``D16``. ``core_center`` (one ply at the mid-plane) keeps
``B`` at zero because its centroid is at ``z = 0``. What this model *cannot*
express is an odd ply count in a pad zone, or an arbitrary unsymmetric global
order -- real thin-laminate layups sometimes need that, and it is the job of a
future ``ExplicitDesign``. Nothing here rejects an unsymmetric result;
``abd.laminate_summary`` reports whatever ``B`` comes out.

Ply thickness is **taken from the shop inventory SKU**, never chosen: the
generator is handed the ``dict[str, ShopSku]`` from
``shop_inventory.load_shop_inventory`` and writes the ``kind`` column from the
SKU's own record, so ``plybook.check_against_library``'s cross-check passes.

Units: degrees and mm, per ``CLAUDE.md``. Nothing here converts units.
"""

from __future__ import annotations

import copy
import hashlib
import itertools
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from compfea.layup import canonical_angle
from compfea.plybook import ALL_ZONES, PlyBookError, PlyRow
from compfea.shop_inventory import ShopSku
from compfea.step_mesh import _sanitize_elset

REQUIRED_COLUMNS = ("ply", "zone", "material", "angle_deg", "thickness_mm", "kind")

# write_plybook uses these so the CSV carries the full precision that
# design_fingerprint hashes: angle to canonical_angle's 6 decimal places,
# thickness at full IEEE-double round-trip precision (``.17g`` -- an inventory
# thickness derived from oz/yd^2 is not a round decimal). A lossy ``:g`` here
# would let the solved deck disagree with the laminate the cache key names.
_ANGLE_FMT = ".6f"
_THICK_FMT = ".17g"


class PlyBookGenError(PlyBookError):
    """A design vector that cannot be turned into a ply book anyone intended."""


@dataclass(frozen=True)
class PlySpec:
    """One ply of a design vector: a stocked material and an angle in degrees.

    Thickness is not here -- it comes from the inventory SKU named by
    ``material``. The angle is passed through ``canonical_angle`` when the row
    is built.
    """

    material: str
    angle_deg: float

    def as_pair(self) -> list:
        return [self.material, self.angle_deg]

    @classmethod
    def from_pair(cls, pair: Sequence) -> PlySpec:
        try:
            material, angle_deg = pair
            return cls(str(material), float(angle_deg))
        except (TypeError, ValueError) as exc:
            raise PlyBookGenError(
                f"ply spec must be [material, angle_deg], got {pair!r}"
            ) from exc


#: One pad zone's half-stack: ``(zone_name, (PlySpec, ...))``.
PadHalfStack = tuple[str, tuple[PlySpec, ...]]


def _norm_specs(raw: object) -> tuple[PlySpec, ...]:
    if raw is None:
        return ()
    if isinstance(raw, PlySpec):
        return (raw,)
    return tuple(
        s if isinstance(s, PlySpec) else PlySpec.from_pair(s)
        for s in raw  # type: ignore[union-attr]
    )


def _check_named(zone: str, what: str) -> str:
    z = str(zone).strip()
    if not z or z.lower() in ALL_ZONES:
        raise PlyBookGenError(
            f"{what} must be a named zone, not {zone!r} -- this module cannot "
            "emit a whole-part ply; use a real ELSET name"
        )
    return z


@dataclass(frozen=True)
class TaperDesign:
    """A symmetric tapered laminate over an ordered list of nested zones.

    ``pads`` is normalised to sorted ``(zone, half_stack)`` items so the
    dataclass stays frozen, hashable and comparable regardless of construction
    path; ``pads_map`` gives the dict view. ``pad_order`` is the ``-z`` ->
    mid-plane order the pad groups are laid in, which *is* their
    through-thickness position and therefore their bending leverage. Every pad
    zone with a non-empty half-stack must appear in ``pad_order`` (checked
    here); names in ``pad_order`` with no pad are allowed and skipped, so a grid
    can list every zone unconditionally.
    """

    through_zone: str
    skins: tuple[PlySpec, ...] = ()
    core: tuple[PlySpec, ...] = ()
    core_center: PlySpec | None = None
    pads: tuple[PadHalfStack, ...] = ()
    pad_order: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        set_ = object.__setattr__
        set_(self, "through_zone", _check_named(self.through_zone, "through_zone"))
        set_(self, "skins", _norm_specs(self.skins))
        set_(self, "core", _norm_specs(self.core))
        if self.core_center is not None and not isinstance(self.core_center, PlySpec):
            set_(self, "core_center", PlySpec.from_pair(self.core_center))

        items = list(
            self.pads.items() if isinstance(self.pads, Mapping) else self.pads
        )
        pad_names = [_check_named(zone, "pad zone") for zone, _ in items]
        pad_dup = sorted({z for z in pad_names if pad_names.count(z) > 1})
        if pad_dup:
            raise PlyBookGenError(
                f"pads names zone(s) {pad_dup} more than once -- one half-stack "
                "per zone; merge them"
            )
        pads = tuple(
            sorted(
                (name, _norm_specs(half))
                for name, (_, half) in zip(pad_names, items, strict=True)
            )
        )
        set_(self, "pads", pads)
        set_(self, "pad_order", tuple(str(z).strip() for z in self.pad_order))

        dup = sorted({z for z in self.pad_order if self.pad_order.count(z) > 1})
        if dup:
            raise PlyBookGenError(f"pad_order has duplicate zone(s) {dup}")
        missing = sorted(
            {zone for zone, half in self.pads if half} - set(self.pad_order)
        )
        if missing:
            raise PlyBookGenError(
                f"pad zone(s) {missing} carry plies but are not in pad_order "
                f"{list(self.pad_order)} -- pad_order sets the through-thickness "
                "position of every pad group and cannot be guessed"
            )

    @property
    def pads_map(self) -> dict[str, tuple[PlySpec, ...]]:
        return {zone: half for zone, half in self.pads}


def design_from_dict(obj: Mapping) -> TaperDesign:
    """Parse a plain dict (JSON) into a validated ``TaperDesign``.

    Does not touch the inventory -- material names and angles are only checked
    against stock when ``design_rows`` runs. Structural validation
    (``through_zone`` named, pad zones in ``pad_order``) happens in
    ``TaperDesign.__post_init__``.
    """
    center = obj.get("core_center")
    return TaperDesign(
        through_zone=str(obj.get("through_zone") or ""),
        skins=_norm_specs(obj.get("skins")),
        core=_norm_specs(obj.get("core")),
        core_center=PlySpec.from_pair(center) if center else None,
        pads=tuple(
            (str(zone), _norm_specs(half))
            for zone, half in (obj.get("pads") or {}).items()
        ),
        pad_order=tuple(str(z) for z in (obj.get("pad_order") or ())),
    )


def design_to_dict(design: TaperDesign) -> dict:
    """``TaperDesign`` -> a plain JSON-round-trippable dict."""
    return {
        "through_zone": design.through_zone,
        "skins": [s.as_pair() for s in design.skins],
        "core": [s.as_pair() for s in design.core],
        "core_center": design.core_center.as_pair() if design.core_center else None,
        "pads": {zone: [s.as_pair() for s in half] for zone, half in design.pads},
        "pad_order": list(design.pad_order),
    }


def _row(idx: int, zone: str, spec: PlySpec, skus: Mapping[str, ShopSku]) -> PlyRow:
    try:
        sku = skus[spec.material]
    except KeyError:
        raise PlyBookGenError(
            f"ply {idx}: material {spec.material!r} is not stocked "
            f"(inventory has {sorted(skus)})"
        ) from None
    try:
        angle = canonical_angle(spec.angle_deg)
    except ValueError as exc:
        raise PlyBookGenError(f"ply {idx}: angle_deg -- {exc}") from None
    return PlyRow(
        ply=idx,
        zone=_sanitize_elset(zone),
        material=spec.material,
        angle_deg=angle,
        thickness_mm=sku.thickness_mm,
        kind=sku.record.kind,
    )


def design_rows(
    design: TaperDesign, *, skus: Mapping[str, ShopSku]
) -> tuple[PlyRow, ...]:
    """Design vector + shop inventory -> ``plybook.PlyRow`` tuple, ply 1 = -z.

    Stacking order, -z -> +z::

        skins -> pads (pad_order: outer -> inner) -> core half -> [core_center]
              -> core half reversed -> pads reversed -> skins reversed

    which is its own reverse, so the laminate is mirror-symmetric about the
    mid-plane for *every* element regardless of which pads cover it.

    Raises ``PlyBookGenError`` for an unstocked material, a non-finite angle, or
    a design with nothing covering ``through_zone`` (which would leave every
    element with an empty stack).
    """
    if not (design.skins or design.core or design.core_center):
        raise PlyBookGenError(
            f"design: nothing covers through_zone {design.through_zone!r} -- "
            "add at least one skin or core ply, or every element is bare"
        )

    pads_map = design.pads_map
    ordered_pads = [z for z in design.pad_order if pads_map.get(z)]

    rows: list[PlyRow] = []
    idx = 1

    def emit(zone: str, specs: Sequence[PlySpec]) -> None:
        nonlocal idx
        for spec in specs:
            rows.append(_row(idx, zone, spec, skus))
            idx += 1

    emit(design.through_zone, design.skins)
    for zone in ordered_pads:
        emit(zone, pads_map[zone])
    emit(design.through_zone, design.core)
    if design.core_center is not None:
        emit(design.through_zone, (design.core_center,))
    emit(design.through_zone, tuple(reversed(design.core)))
    for zone in reversed(ordered_pads):
        emit(zone, tuple(reversed(pads_map[zone])))
    emit(design.through_zone, tuple(reversed(design.skins)))

    if [r.ply for r in rows] != list(range(1, len(rows) + 1)):
        raise PlyBookGenError("internal: emitted ply numbering is not 1..N")
    return tuple(rows)


def design_fingerprint(
    design: TaperDesign, *, skus: Mapping[str, ShopSku]
) -> str:
    """16-hex hash of the **resolved laminate** -- the sweep cache key.

    Covers the ``PlyRow`` sequence (materials, angles, thicknesses, zones,
    order) *and* the ten elastic constants of every stocked lamina it uses, the
    same reason ``build.layup_fingerprint`` interpolates all of them: editing a
    modulus under a stable SKU name changes the laminate and must miss the
    cache. Not keyed on the dict form, so a reordered ``pads`` mapping that
    resolves identically still hits.
    """
    rows = design_rows(design, skus=skus)
    parts = [
        f"{r.ply}:{r.zone}:{r.material}:{r.angle_deg:{_ANGLE_FMT}}:"
        f"{r.thickness_mm:{_THICK_FMT}}:{r.kind}"
        for r in rows
    ]
    for name in sorted({r.material for r in rows}):
        ec = skus[name].record.constants
        parts.append(
            f"mat={name}:{ec.e1!r}:{ec.e2!r}:{ec.e3!r}:{ec.nu12!r}:{ec.nu13!r}:"
            f"{ec.nu23!r}:{ec.g12!r}:{ec.g13!r}:{ec.g23!r}:{ec.density!r}"
        )
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


def _comment_block(design: TaperDesign, note: str) -> list[str]:
    lines = ["# Generated by compfea.plybook_gen from a TaperDesign."]
    if note:
        lines += [f"# {ln}" for ln in note.splitlines()]
    lines.append(
        f"# through_zone={design.through_zone}  pad_order={list(design.pad_order)}"
    )
    lines.append(
        "# symmetric by construction (core + each pad mirrored about the "
        "mid-plane); B ~ 0"
    )
    return lines


def write_plybook(
    design: TaperDesign,
    path: str | Path,
    *,
    skus: Mapping[str, ShopSku],
    note: str = "",
) -> Path:
    """Write a ply book CSV that ``plybook.load_plybook`` reads back unchanged."""
    rows = design_rows(design, skus=skus)
    path = Path(path)
    out = [*_comment_block(design, note), ",".join(REQUIRED_COLUMNS)]
    for r in rows:
        out.append(
            f"{r.ply},{r.zone},{r.material},"
            f"{r.angle_deg:{_ANGLE_FMT}},{r.thickness_mm:{_THICK_FMT}},{r.kind}"
        )
    path.write_text("\n".join(out) + "\n")
    return path


def expand_grid(
    base: Mapping, axes: Mapping[str, Sequence]
) -> list[TaperDesign]:
    """Cartesian product of ``axes`` over a ``base`` design dict.

    An axis key is either a top-level field of the design dict (``core``,
    ``skins``, ``pad_order``, ``core_center``) or ``"pads.<ZONE>"`` to vary one
    zone's half-stack. Each combination is run through ``design_from_dict``, so
    a structurally broken point in the grid fails loudly here.
    """
    keys = list(axes)
    designs: list[TaperDesign] = []
    for combo in itertools.product(*(list(axes[k]) for k in keys)):
        obj = copy.deepcopy(dict(base))
        for key, value in zip(keys, combo, strict=True):
            if key.startswith("pads."):
                obj["pads"] = {
                    **(obj.get("pads") or {}),
                    key[len("pads.") :]: value,
                }
            else:
                obj[key] = value
        designs.append(design_from_dict(obj))
    return designs
