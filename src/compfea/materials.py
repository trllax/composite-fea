"""Engineering materials library: a CSV of lamina cards, checked on the way in.

One row is one lamina. The nine elastic constants and the density are what ccx
sees; the strength allowables ride along for post-processing and are never
written to a deck.

Units: **mm, N, tonne, MPa, s** once loaded, per ``CLAUDE.md``. Getting there is
the point of this module. ANSYS Engineering Data exports stiffness in Pa and
density in kg/m^3, this repo wants MPa and tonne/mm^3, and CalculiX has no unit
system to object with -- a Pa-valued modulus produces a converged, plausible,
wrong answer. So a library file **must** declare its units:

    # units: stress=Pa density=kg/m3
    name,kind,e1,e2,e3,nu12,nu13,nu23,g12,g13,g23,density,xt,xc,yt,yc,s12,source

A file with no ``# units:`` line is refused rather than guessed at.

``kind`` is ``ud`` or ``woven``. It is honest metadata: it gates validation and
it is echoed in reports. It **does not change the deck** -- a woven lamina is a
single orthotropic card either way, and nothing here expands one into two UD
plies. See ``layup.woven_from_ud`` for the derivation used when a woven card is
computed from a UD one instead of measured.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from compfea.layup import EngineeringConstants

#: Multiply a declared stress unit by this to reach MPa.
STRESS_TO_MPA: dict[str, float] = {
    "pa": 1e-6,
    "kpa": 1e-3,
    "mpa": 1.0,
    "gpa": 1e3,
    "psi": 6.894757293168361e-3,
    "ksi": 6.894757293168361,
}

#: Multiply a declared density unit by this to reach tonne/mm^3.
DENSITY_TO_TONNE_PER_MM3: dict[str, float] = {
    "kg/m3": 1e-12,
    "g/cm3": 1e-9,
    "tonne/mm3": 1.0,
    "t/mm3": 1.0,
}

#: Plausible modulus band in MPa. Basalt UD sits near 45 GPa at the bottom and
#: M60J-class carbon near 588 GPa at the top, so this is wide on purpose: it is
#: a unit check, not a materials opinion.
#:
#: Note what catches what. A file left in Pa lands at 1e11 and blows the ceiling
#: on every field. A file left in **GPa** is caught by the transverse and shear
#: moduli -- E2 at 9 and G23 at 3 are both under the floor -- and *not* by E1,
#: which lands at 121 and sits comfortably inside the band. The floor stays at
#: 100 rather than rising to catch E1 too, because a matrix-dominated G23 is
#: legitimately a few thousand MPa and a tighter floor would start refusing real
#: laminae. Every field is checked, so the GPa case is caught either way.
MODULUS_BAND_MPA = (100.0, 600_000.0)

#: Plausible density band in tonne/mm^3, i.e. 0.1 to 10 g/cm^3.
DENSITY_BAND = (1e-10, 1e-8)

KINDS = ("ud", "woven")

_ELASTIC = ("e1", "e2", "e3", "nu12", "nu13", "nu23", "g12", "g13", "g23")
_ALLOWABLES = ("xt", "xc", "yt", "yc", "s12")
REQUIRED_COLUMNS = ("name", "kind", *_ELASTIC, "density", "source")


class MaterialError(ValueError):
    """A library file that cannot be trusted to produce a correct deck."""


@dataclass(frozen=True)
class MaterialRecord:
    """One library row.

    ``constants`` is the part ccx consumes. Everything else is provenance and
    post-processing: no allowable reaches a deck, because ccx has nowhere to put
    one and inventing a failure criterion is not this module's job.
    """

    constants: EngineeringConstants
    kind: str
    source: str
    xt: float | None = None
    xc: float | None = None
    yt: float | None = None
    yc: float | None = None
    s12: float | None = None

    @property
    def name(self) -> str:
        return self.constants.name

    def allowables(self) -> dict[str, float | None]:
        return {k: getattr(self, k) for k in _ALLOWABLES}


def parse_units_header(lines: list[str]) -> tuple[float, float]:
    """``# units: stress=<u> density=<u>`` -> scales to (MPa, tonne/mm^3).

    Refuses a missing or partial declaration. This is the one place a factor of
    10^6 can enter the repo, so it does not get a default.
    """
    declared: dict[str, str] = {}
    for line in lines:
        stripped = line.strip()
        if not stripped.startswith("#"):
            break
        body = stripped.lstrip("#").strip()
        if not body.lower().startswith("units:"):
            continue
        for token in body[len("units:") :].replace(",", " ").split():
            if "=" not in token:
                raise MaterialError(
                    f"malformed units token {token!r}; expected key=value"
                )
            key, value = token.split("=", 1)
            declared[key.strip().lower()] = value.strip().lower()
    if not declared:
        raise MaterialError(
            "no '# units:' line. A materials file must declare its units --"
            " ANSYS exports Pa and kg/m3, this repo needs MPa and tonne/mm^3,"
            " and CalculiX will not object to either. Add e.g."
            " '# units: stress=Pa density=kg/m3' as the first line."
        )
    missing = {"stress", "density"} - set(declared)
    if missing:
        raise MaterialError(
            f"units line declares {sorted(declared)} but is missing "
            f"{sorted(missing)}; both are required"
        )
    try:
        stress = STRESS_TO_MPA[declared["stress"]]
    except KeyError:
        raise MaterialError(
            f"unknown stress unit {declared['stress']!r}; "
            f"known: {sorted(STRESS_TO_MPA)}"
        ) from None
    try:
        density = DENSITY_TO_TONNE_PER_MM3[declared["density"]]
    except KeyError:
        raise MaterialError(
            f"unknown density unit {declared['density']!r}; "
            f"known: {sorted(DENSITY_TO_TONNE_PER_MM3)}"
        ) from None
    return stress, density


def check_admissible(ec: EngineeringConstants) -> list[str]:
    """Refuse an impossible lamina; return warnings for a suspicious one.

    ccx accepts a thermodynamically inadmissible orthotropic card without
    complaint and returns a converged answer from it, so the check has to happen
    here or not at all.
    """
    for field in _ELASTIC[:3] + _ELASTIC[6:]:
        value = getattr(ec, field)
        if not (math.isfinite(value) and value > 0.0):
            raise MaterialError(
                f"{ec.name}: {field} must be finite and > 0, got {value!r}"
            )
        lo, hi = MODULUS_BAND_MPA
        if not lo <= value <= hi:
            raise MaterialError(
                f"{ec.name}: {field} = {value:g} MPa is outside the plausible "
                f"band [{lo:g}, {hi:g}] MPa. This is almost always a unit "
                "error -- check the '# units:' line against the file's numbers."
            )
    if not math.isfinite(ec.density) or not (
        DENSITY_BAND[0] <= ec.density <= DENSITY_BAND[1]
    ):
        raise MaterialError(
            f"{ec.name}: density = {ec.density:g} tonne/mm^3 is outside "
            f"[{DENSITY_BAND[0]:g}, {DENSITY_BAND[1]:g}] (0.1 to 10 g/cm^3)"
        )

    # Orthotropic admissibility. Each minor ratio follows from the major one.
    e = {"1": ec.e1, "2": ec.e2, "3": ec.e3}
    nu = {"12": ec.nu12, "13": ec.nu13, "23": ec.nu23}
    for pair, value in nu.items():
        if not math.isfinite(value):
            raise MaterialError(f"{ec.name}: nu{pair} must be finite, got {value!r}")
        i, j = pair
        limit = math.sqrt(e[i] / e[j])
        if abs(value) >= limit:
            raise MaterialError(
                f"{ec.name}: |nu{pair}| = {abs(value):g} must be < "
                f"sqrt(E{i}/E{j}) = {limit:g}; the compliance matrix is not "
                "positive definite and ccx will solve it anyway"
            )
    nu21 = ec.nu12 * ec.e2 / ec.e1
    nu31 = ec.nu13 * ec.e3 / ec.e1
    nu32 = ec.nu23 * ec.e3 / ec.e2
    delta = (
        1.0
        - ec.nu12 * nu21
        - ec.nu23 * nu32
        - ec.nu13 * nu31
        - 2.0 * nu21 * nu32 * ec.nu13
    )
    if delta <= 0.0:
        raise MaterialError(
            f"{ec.name}: 1 - nu12*nu21 - nu23*nu32 - nu13*nu31 - "
            f"2*nu21*nu32*nu13 = {delta:g} <= 0, so the 3-D orthotropic "
            "stiffness is not positive definite"
        )

    warnings: list[str] = []
    if ec.e1 / ec.e2 > 5.0 and ec.nu12 < 0.05:
        warnings.append(
            f"{ec.name}: nu12 = {ec.nu12:g} is very small for a lamina with "
            f"E1/E2 = {ec.e1 / ec.e2:.0f}. ANSYS lists the *major* ratio "
            "(nu_xy, typically 0.25-0.35) in the slot ccx wants; this looks "
            "like the minor ratio nu21."
        )
    return warnings


def _kind_warnings(record: MaterialRecord) -> list[str]:
    if record.kind != "woven":
        return []
    ratio = record.constants.e1 / record.constants.e2
    if not 0.8 <= ratio <= 1.25:
        return [
            f"{record.name}: kind=woven but E1/E2 = {ratio:.3g}. A balanced "
            "weave has near-equal warp and fill stiffness; this row looks like "
            "a UD lamina labelled woven."
        ]
    return []


def _cell(row: dict[str, str], key: str, name: str) -> str:
    if key not in row:
        raise MaterialError(f"{name}: missing column {key!r}")
    value = row[key]
    return "" if value is None else value.strip()


def _number(row: dict[str, str], key: str, name: str, scale: float) -> float:
    raw = _cell(row, key, name)
    if not raw:
        raise MaterialError(f"{name}: {key} is required and empty")
    try:
        return float(raw) * scale
    except ValueError:
        raise MaterialError(f"{name}: {key} = {raw!r} is not a number") from None


def _optional(row: dict[str, str], key: str, name: str, scale: float) -> float | None:
    raw = row.get(key)
    raw = "" if raw is None else raw.strip()
    if not raw:
        return None
    try:
        value = float(raw) * scale
    except ValueError:
        raise MaterialError(f"{name}: {key} = {raw!r} is not a number") from None
    if not math.isfinite(value) or value <= 0.0:
        raise MaterialError(
            f"{name}: {key} = {value:g} MPa must be a positive magnitude; "
            "compressive allowables are written as positive numbers here"
        )
    return value


def load_materials(
    path: str | Path, *, on_warning: str = "collect"
) -> dict[str, MaterialRecord]:
    """Read a materials CSV into ``{name: MaterialRecord}``, in file order.

    ``on_warning`` is ``"collect"`` (attach to ``load_materials.warnings``),
    ``"raise"``, or ``"ignore"``. Warnings are things that are legal and look
    wrong -- a minor Poisson ratio in the major slot, a UD card labelled woven.
    """
    if on_warning not in ("collect", "raise", "ignore"):
        raise ValueError(f"on_warning must be collect/raise/ignore, got {on_warning!r}")
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()
    stress_scale, density_scale = parse_units_header(lines)

    body = [ln for ln in lines if not ln.strip().startswith("#")]
    reader = csv.DictReader(body)
    if reader.fieldnames is None:
        raise MaterialError(f"{path}: no header row")
    header = [f.strip().lower() for f in reader.fieldnames]
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    if missing:
        raise MaterialError(f"{path}: missing required columns {missing}")

    out: dict[str, MaterialRecord] = {}
    warnings: list[str] = []
    for lineno, raw_row in enumerate(reader, start=2):
        row = {
            (k.strip().lower() if k else ""): (v or "")
            for k, v in raw_row.items()
            if k is not None
        }
        name = _cell(row, "name", f"{path}:{lineno}")
        if not name:
            continue  # blank line
        where = f"{path}:{lineno} ({name})"
        if name in out:
            raise MaterialError(
                f"{where}: duplicate material name; a ply line naming it would "
                "silently get whichever row ccx kept"
            )
        kind = _cell(row, "kind", where).lower()
        if kind not in KINDS:
            raise MaterialError(f"{where}: kind must be one of {KINDS}, got {kind!r}")
        source = _cell(row, "source", where)
        if not source:
            raise MaterialError(
                f"{where}: source is required. A modulus with no provenance is "
                "indistinguishable from a guess six months later."
            )
        constants = EngineeringConstants(
            e1=_number(row, "e1", where, stress_scale),
            e2=_number(row, "e2", where, stress_scale),
            e3=_number(row, "e3", where, stress_scale),
            nu12=_number(row, "nu12", where, 1.0),
            nu13=_number(row, "nu13", where, 1.0),
            nu23=_number(row, "nu23", where, 1.0),
            g12=_number(row, "g12", where, stress_scale),
            g13=_number(row, "g13", where, stress_scale),
            g23=_number(row, "g23", where, stress_scale),
            density=_number(row, "density", where, density_scale),
            name=name,
        )
        record = MaterialRecord(
            constants=constants,
            kind=kind,
            source=source,
            **{k: _optional(row, k, where, stress_scale) for k in _ALLOWABLES},
        )
        warnings.extend(check_admissible(constants))
        warnings.extend(_kind_warnings(record))
        out[name] = record

    if not out:
        raise MaterialError(f"{path}: no material rows")
    if warnings and on_warning == "raise":
        raise MaterialError("; ".join(warnings))
    load_materials.warnings = [] if on_warning == "ignore" else warnings
    return out


load_materials.warnings = []


def constants_for(
    library: dict[str, MaterialRecord], names: list[str] | tuple[str, ...]
) -> tuple[EngineeringConstants, ...]:
    """The cards for ``names``, deduplicated and sorted -- deck-ready.

    Sorted so a deck is byte-identical across runs; only referenced materials
    are emitted so an unused library row cannot drift into a cache key.
    """
    unknown = sorted(set(names) - set(library))
    if unknown:
        raise MaterialError(
            f"material(s) {unknown} not in the library {sorted(library)}"
        )
    return tuple(library[n].constants for n in sorted(set(names)))
