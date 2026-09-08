"""Shop inventory CSV -> lamina cards for ply books.

Stock is scarce and discrete: architecture, areal weight / thickness, and an
optional commercial id. Elastic constants are often missing on the roll. When
the nine lamina numbers are blank, this module fills them from the repo's
default UD card (``ANSYS_EPOXY_CARBON_UD_230``, Engineering Data
"Epoxy Carbon UD (230GPa) Prepreg") or ``ANSYS_EPOXY_CARBON_WOVEN_230_WET`` for ``0_90`` / ``biax``.
Closer to ANSYS ACP Engineering Data; still not a roll-specific datasheet.

Fibre-only labels (e.g. AS4C 33.5 Msi / 653 ksi) belong in ``fibre_e_msi`` /
``fibre_xt_ksi`` (or ``source``). They are **never** copied into ``e1``. A
fibre modulus is not a woven lamina card.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from compfea.layup import (
    ANSYS_EPOXY_CARBON_UD_230,
    ANSYS_EPOXY_CARBON_WOVEN_230_WET,
    EngineeringConstants,
    woven_from_ud,
)
from compfea.materials import (
    DENSITY_TO_TONNE_PER_MM3,
    STRESS_TO_MPA,
    MaterialError,
    MaterialRecord,
    check_admissible,
    parse_units_header,
)

ARCH_TO_KIND = {"ud": "ud", "0_90": "woven", "biax": "woven"}
ARCHITECTURES = tuple(ARCH_TO_KIND)
OZ_YD2_TO_GSM = 33.9057
GSM_TO_THICKNESS_MM = 1.0 / 1000.0  # 100 gsm -> 0.1 mm

_ELASTIC = ("e1", "e2", "e3", "nu12", "nu13", "nu23", "g12", "g13", "g23")
REQUIRED_COLUMNS = (
    "name",
    "architecture",
    "areal_weight_gsm",
    "areal_weight_oz_yd2",
    "thickness_mm",
)


class ShopInventoryError(MaterialError):
    """A shop inventory row that cannot become a trusted lamina card."""


@dataclass(frozen=True)
class ShopSku:
    """One stocked reinforcement, with a lamina card ready for the ply book."""

    name: str
    product: str
    architecture: str
    tow: str
    areal_weight_gsm: float
    thickness_mm: float
    width_in: float | None
    record: MaterialRecord
    used_default_lamina: bool
    fibre_e_msi: float | None = None
    fibre_xt_ksi: float | None = None


def default_lamina(architecture: str, name: str) -> EngineeringConstants:
    """ANSYS UD-230-based lamina card for a stock architecture.

    UD rows copy ``ANSYS_EPOXY_CARBON_UD_230``. ``0_90`` / ``biax`` use
    ``woven_from_ud`` of that card (membrane-equivalent weave stand-in).
    """
    if architecture not in ARCH_TO_KIND:
        raise ShopInventoryError(
            f"architecture must be one of {ARCHITECTURES}, got {architecture!r}"
        )
    if architecture == "ud":
        base = ANSYS_EPOXY_CARBON_UD_230
        return EngineeringConstants(
            e1=base.e1,
            e2=base.e2,
            e3=base.e3,
            nu12=base.nu12,
            nu13=base.nu13,
            nu23=base.nu23,
            g12=base.g12,
            g13=base.g13,
            g23=base.g23,
            density=base.density,
            name=name,
        )
    # 0_90 and biax: Engineering Data woven wet card (not woven_from_ud).
    # biax (+/-45) still uses this as a stand-in until a ±45 card is measured.
    base = ANSYS_EPOXY_CARBON_WOVEN_230_WET
    return EngineeringConstants(
        e1=base.e1,
        e2=base.e2,
        e3=base.e3,
        nu12=base.nu12,
        nu13=base.nu13,
        nu23=base.nu23,
        g12=base.g12,
        g13=base.g13,
        g23=base.g23,
        density=base.density,
        name=name,
    )


def _cell(row: dict[str, str], key: str) -> str:
    return (row.get(key) or "").strip()


def _optional_float(row: dict[str, str], key: str, where: str) -> float | None:
    raw = _cell(row, key)
    if not raw:
        return None
    try:
        value = float(raw)
    except ValueError as exc:
        raise ShopInventoryError(f"{where}: {key} = {raw!r} is not a number") from exc
    if not math.isfinite(value):
        raise ShopInventoryError(f"{where}: {key} must be finite")
    return value


def resolve_areal_weight_gsm(row: dict[str, str], where: str) -> float:
    """Prefer labeled gsm; else convert oz/yd2."""
    gsm = _optional_float(row, "areal_weight_gsm", where)
    oz = _optional_float(row, "areal_weight_oz_yd2", where)
    if gsm is not None and gsm <= 0:
        raise ShopInventoryError(f"{where}: areal_weight_gsm must be > 0")
    if oz is not None and oz <= 0:
        raise ShopInventoryError(f"{where}: areal_weight_oz_yd2 must be > 0")
    if gsm is not None and oz is not None:
        converted = oz * OZ_YD2_TO_GSM
        if abs(gsm - converted) / gsm > 0.05:
            raise ShopInventoryError(
                f"{where}: gsm={gsm:g} and oz/yd2={oz:g} disagree by more than "
                f"5% (oz implies {converted:g} gsm). Prefer the labeled gsm or "
                "fix the oz."
            )
        return gsm
    if gsm is not None:
        return gsm
    if oz is not None:
        return oz * OZ_YD2_TO_GSM
    raise ShopInventoryError(
        f"{where}: need areal_weight_gsm and/or areal_weight_oz_yd2"
    )


def resolve_thickness_mm(row: dict[str, str], gsm: float, where: str) -> float:
    """Explicit thickness, else gsm/1000."""
    t = _optional_float(row, "thickness_mm", where)
    derived = gsm * GSM_TO_THICKNESS_MM
    if t is None:
        return derived
    if t <= 0:
        raise ShopInventoryError(f"{where}: thickness_mm must be > 0")
    if abs(t - derived) / derived > 0.25:
        # Soft check: cured ply can differ; warn via load_shop_inventory.warnings
        pass
    return t


def _elastic_all_blank(row: dict[str, str]) -> bool:
    return all(not _cell(row, k) for k in _ELASTIC) and not _cell(row, "density")


def _elastic_from_row(
    row: dict[str, str],
    where: str,
    name: str,
    stress_scale: float,
    density_scale: float,
) -> EngineeringConstants:
    def req(key: str, scale: float) -> float:
        raw = _cell(row, key)
        if not raw:
            raise ShopInventoryError(
                f"{where}: {key} is empty. Fill all lamina elastic fields, or "
                "leave all of them blank to use the default lamina card."
            )
        try:
            return float(raw) * scale
        except ValueError as exc:
            raise ShopInventoryError(
                f"{where}: {key} = {raw!r} is not a number"
            ) from exc

    return EngineeringConstants(
        e1=req("e1", stress_scale),
        e2=req("e2", stress_scale),
        e3=req("e3", stress_scale),
        nu12=req("nu12", 1.0),
        nu13=req("nu13", 1.0),
        nu23=req("nu23", 1.0),
        g12=req("g12", stress_scale),
        g13=req("g13", stress_scale),
        g23=req("g23", stress_scale),
        density=req("density", density_scale),
        name=name,
    )


def load_shop_inventory(
    path: str | Path, *, on_warning: str = "collect"
) -> dict[str, ShopSku]:
    """Read shop stock into ``{name: ShopSku}`` with lamina cards attached.

    Blank lamina elastic fields -> ``default_lamina(architecture)``. Fibre
    columns never become ``e1``.
    """
    if on_warning not in ("collect", "raise", "ignore"):
        raise ValueError(f"on_warning must be collect/raise/ignore, got {on_warning!r}")
    path = Path(path)
    lines = path.read_text().splitlines()
    stress_scale, density_scale = parse_units_header(lines)
    body = [ln for ln in lines if not ln.strip().startswith("#")]
    reader = csv.DictReader(body)
    if reader.fieldnames is None:
        raise ShopInventoryError(f"{path}: no header row")
    header = [f.strip().lower() for f in reader.fieldnames if f]
    for col in ("name", "architecture"):
        if col not in header:
            raise ShopInventoryError(f"{path}: missing required column {col!r}")

    out: dict[str, ShopSku] = {}
    warnings: list[str] = []
    for lineno, raw_row in enumerate(reader, start=2):
        row = {
            (k.strip().lower() if k else ""): (v or "")
            for k, v in raw_row.items()
            if k is not None
        }
        name = _cell(row, "name")
        if not name:
            continue
        where = f"{path}:{lineno} ({name})"
        if name in out:
            raise ShopInventoryError(f"{where}: duplicate name")
        architecture = _cell(row, "architecture").lower()
        if architecture not in ARCH_TO_KIND:
            raise ShopInventoryError(
                f"{where}: architecture must be one of {ARCHITECTURES}, "
                f"got {architecture!r}"
            )
        kind = ARCH_TO_KIND[architecture]
        gsm = resolve_areal_weight_gsm(row, where)
        thickness_mm = resolve_thickness_mm(row, gsm, where)
        derived_t = gsm * GSM_TO_THICKNESS_MM
        if abs(thickness_mm - derived_t) / derived_t > 0.25:
            warnings.append(
                f"{where}: thickness_mm={thickness_mm:g} differs >25% from "
                f"gsm/1000={derived_t:g}; using the explicit thickness"
            )

        used_default = _elastic_all_blank(row)
        if used_default:
            constants = default_lamina(architecture, name)
            source = _cell(row, "source") or (
                f"default lamina ({architecture} -> {kind} from "
                f"{'ANSYS_EPOXY_CARBON_UD_230' if architecture == 'ud' else 'woven_from_ud'}); "
                "not a datasheet"
            )
            if "default lamina" not in source:
                source = f"{source}; default lamina card filled for blank E"
        else:
            constants = _elastic_from_row(
                row, where, name, stress_scale, density_scale
            )
            source = _cell(row, "source")
            if not source:
                raise ShopInventoryError(
                    f"{where}: source is required when lamina moduli are given"
                )

        fibre_e = _optional_float(row, "fibre_e_msi", where)
        fibre_xt = _optional_float(row, "fibre_xt_ksi", where)
        if fibre_e is not None and fibre_e <= 0:
            raise ShopInventoryError(f"{where}: fibre_e_msi must be > 0")
        if fibre_xt is not None and fibre_xt <= 0:
            raise ShopInventoryError(f"{where}: fibre_xt_ksi must be > 0")

        record = MaterialRecord(
            constants=constants, kind=kind, source=source
        )
        warnings.extend(check_admissible(constants))
        if record.kind == "woven":
            ratio = record.constants.e1 / record.constants.e2
            if not 0.8 <= ratio <= 1.25:
                warnings.append(
                    f"{where}: woven-like architecture but E1/E2 = {ratio:.3g}"
                )

        width = _optional_float(row, "width_in", where)
        sku = ShopSku(
            name=name,
            product=_cell(row, "product"),
            architecture=architecture,
            tow=_cell(row, "tow"),
            areal_weight_gsm=gsm,
            thickness_mm=thickness_mm,
            width_in=width,
            record=record,
            used_default_lamina=used_default,
            fibre_e_msi=fibre_e,
            fibre_xt_ksi=fibre_xt,
        )
        out[name] = sku

    if on_warning == "raise" and warnings:
        raise ShopInventoryError("; ".join(warnings))
    load_shop_inventory.warnings = warnings if on_warning == "collect" else []
    return out


load_shop_inventory.warnings = []  # type: ignore[attr-defined]


def materials_library(skus: dict[str, ShopSku]) -> dict[str, MaterialRecord]:
    """``{name: MaterialRecord}`` for ``plybook.resolve`` / ``compfea-build``."""
    return {name: sku.record for name, sku in skus.items()}


def write_materials_csv(skus: dict[str, ShopSku], path: str | Path) -> Path:
    """Emit a materials library file in repo units (MPa, tonne/mm^3)."""
    path = Path(path)
    with path.open("w", newline="") as fh:
        fh.write("# units: stress=MPa density=tonne/mm3\n")
        fh.write(
            "# Generated from shop inventory. Rows that used_default_lamina "
            "say so in source.\n"
        )
        writer = csv.writer(fh)
        writer.writerow(
            [
                "name",
                "kind",
                "e1",
                "e2",
                "e3",
                "nu12",
                "nu13",
                "nu23",
                "g12",
                "g13",
                "g23",
                "density",
                "xt",
                "xc",
                "yt",
                "yc",
                "s12",
                "source",
            ]
        )
        for sku in skus.values():
            c = sku.record.constants
            writer.writerow(
                [
                    c.name,
                    sku.record.kind,
                    c.e1,
                    c.e2,
                    c.e3,
                    c.nu12,
                    c.nu13,
                    c.nu23,
                    c.g12,
                    c.g13,
                    c.g23,
                    c.density,
                    "",
                    "",
                    "",
                    "",
                    "",
                    sku.record.source,
                ]
            )
    return path
