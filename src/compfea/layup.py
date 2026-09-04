"""Layup design vectors -> CalculiX materials, orientations, COMPOSITE sections.

Ply angles become *ORIENTATION names (ccx field 4), never bare degrees.
Units: mm, N, tonne, MPa, s (density tonne/mm^3).

Two things here are load-bearing and both fail silently if they are wrong:

- **Names must be one-to-one with angles.** Every distinct ply angle gets its
  own *ORIENTATION card, and the ply lines refer to those cards by name. If two
  distinct angles ever collapse onto one name, the deck carries two cards with
  the same NAME, ccx keeps one, and both plies get a fibre direction that was
  never asked for. Name and direction cosines are therefore both derived from a
  single canonicalised angle, and `Layup.to_inp` re-checks the mapping before
  emitting.
- **The angle sign convention is a rotation, not a mirror.** Positive angle is
  counter-clockwise about +z (the shell normal) from the 0-degree axis, viewed
  down +z, which is the standard composites convention. The mirrored form is
  exactly angle negation -- it emits every +theta ply as -theta. A stack of only
  0 and 90 degree plies cannot detect that at all, since +90 and -90 are the same
  fibre direction; every other stack shows it only in the *sign* of bend-twist
  coupling, never in the magnitude of a reaction. So no force check anywhere in
  this repo will catch it, and the direction cosines are asserted directly in
  tests/test_layup_deck.py instead.

There is no default long axis. `CLAUDE.md` states +x, every case deck in this
repo so far uses +y, and a default that is right for one and silently wrong for
the other is exactly the kind of plausible wrong number this repo exists to
avoid. Callers say which, and each case README records it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Ply angles are rounded to this many decimal places before anything is derived
# from them. It sets both the resolution of the *ORIENTATION name and the
# resolution at which two plies count as sharing an angle -- one constant, so
# the two can never disagree.
_ANGLE_DECIMALS = 6

LONG_AXES = ("x", "y")


def canonical_angle(angle_deg: float) -> float:
    """Ply angle at the one resolution names and direction cosines both use.

    A sweep hands this module computed design vectors, so a non-finite angle is
    reachable. It is rejected here rather than emitted: nan would otherwise reach
    the deck as an *ORIENTATION named ori_pnan whose direction cosines are all
    nan, and ccx does not object to that.
    """
    angle = float(angle_deg)
    if not math.isfinite(angle):
        raise ValueError(f"ply angle must be finite, got {angle_deg!r}")
    angle = round(angle, _ANGLE_DECIMALS)
    return 0.0 if angle == 0.0 else angle  # collapse -0.0


@dataclass(frozen=True)
class EngineeringConstants:
    """Orthotropic lamina card (*ELASTIC, TYPE=ENGINEERING CONSTANTS)."""

    e1: float
    e2: float
    e3: float
    nu12: float
    nu13: float
    nu23: float
    g12: float
    g13: float
    g23: float
    density: float = 1.6e-9
    name: str = "cfrp"

    def to_inp(self) -> str:
        return "\n".join(
            [
                f"*MATERIAL, NAME={self.name}",
                "*ELASTIC, TYPE=ENGINEERING CONSTANTS",
                (
                    f"{self.e1}, {self.e2}, {self.e3}, "
                    f"{self.nu12}, {self.nu13}, {self.nu23}, "
                    f"{self.g12}, {self.g13},"
                ),
                f"{self.g23},",
                "*DENSITY",
                f"{self.density}",
            ]
        )


# Placeholder UD CFRP used in strip cases (MPa, tonne/mm^3).
PLACEHOLDER_CFRP = EngineeringConstants(
    e1=135000.0,
    e2=9000.0,
    e3=9000.0,
    nu12=0.30,
    nu13=0.30,
    nu23=0.45,
    g12=4500.0,
    g13=4500.0,
    g23=3000.0,
)


@dataclass(frozen=True)
class Ply:
    """One lamina. ``angle_deg`` is fiber angle about the shell normal."""

    thickness: float
    angle_deg: float
    material: str = "cfrp"


@dataclass(frozen=True)
class ZoneLayup:
    """Stack for one ELSET. ``plies`` are bottom -> top (first ply at -z)."""

    elset: str
    plies: tuple[Ply, ...]

    def __post_init__(self) -> None:
        if not self.plies:
            raise ValueError(f"zone {self.elset!r} has no plies")
        if not self.elset:
            raise ValueError("elset name is required")

    @property
    def thickness(self) -> float:
        return sum(p.thickness for p in self.plies)


def orientation_name(angle_deg: float) -> str:
    """*ORIENTATION name for a fiber angle, one-to-one with the angle.

    ``ori_p0``, ``ori_m45``, ``ori_p22p5``, ``ori_p30p02``. Sign is the ``p``/``m``
    prefix and the decimal point is ``p``, so the name is a valid ccx identifier
    and two angles that differ at all -- to `_ANGLE_DECIMALS` places -- differ in
    the name.
    """
    angle = canonical_angle(angle_deg)
    digits = f"{abs(angle):.{_ANGLE_DECIMALS}f}".rstrip("0").rstrip(".")
    return f"ori_{'m' if angle < 0 else 'p'}{digits.replace('.', 'p')}"


def _fmt(value: float) -> str:
    """Direction cosine, with sin/cos round-off snapped to a clean zero.

    Keeps the card readable at the axis-aligned angles that make up most decks:
    without it a 90-degree ply prints -0.0000000000 for cos(pi/2). The 10 places
    are chosen to resolve `_ANGLE_DECIMALS`: printing fewer would let two angles
    that earn distinct *ORIENTATION names carry identical direction cosines.
    """
    return f"{0.0 if abs(value) < 1e-12 else value:.10f}"


def orientation_card(angle_deg: float, *, long_axis: str) -> str:
    """Rectangular orientation: 0 deg along the long axis, +z the shell normal.

    Positive ``angle_deg`` rotates the fibre counter-clockwise about +z. The card
    carries two in-plane vectors, ``a`` (the fibre direction) and ``b`` (``a``
    turned +90 deg about +z), so ``a x b = +z`` for every angle and both axis
    choices.

    Note the 90-degree card here is ``-1, 0, 0, 0, -1, 0`` where the hand-written
    ``cases/cantilever_ansys/cantilever_89deg.inp`` writes ``1, 0, 0, 0, 1, 0``.
    Those are 180 degrees apart about +z, which is the same fibre direction for an
    orthotropic ply, and both are correct. Do not "fix" one to match the other by
    dropping a sign -- that is how the mirror gets reintroduced.
    """
    if long_axis not in LONG_AXES:
        raise ValueError(f"long_axis must be one of {LONG_AXES}, not {long_axis!r}")
    angle = canonical_angle(angle_deg)
    phi = math.radians(angle)
    s, c = math.sin(phi), math.cos(phi)
    if long_axis == "x":
        # 0 deg -> +x, +90 deg -> +y.
        a1, a2 = c, s
        b1, b2 = -s, c
    else:
        # 0 deg -> +y, +90 deg -> -x. Rotating the same way about +z as the x
        # case; writing (sin, cos) here instead would mirror it.
        a1, a2 = -s, c
        b1, b2 = -c, -s
    cosines = ", ".join(_fmt(v) for v in (a1, a2, 0.0, b1, b2, 0.0))
    return "\n".join(
        [
            f"*ORIENTATION, NAME={orientation_name(angle)}, SYSTEM=RECTANGULAR",
            cosines,
        ]
    )


@dataclass(frozen=True)
class Layup:
    """Materials + zone stacks ready to emit into a deck.

    ``long_axis`` is the in-plane axis a 0-degree ply runs along, ``"x"`` or
    ``"y"``. It has no default on purpose -- see the module docstring.
    """

    materials: tuple[EngineeringConstants, ...]
    zones: tuple[ZoneLayup, ...]
    long_axis: str

    def __post_init__(self) -> None:
        if not self.materials:
            raise ValueError("at least one material is required")
        if not self.zones:
            raise ValueError("at least one zone layup is required")
        if self.long_axis not in LONG_AXES:
            raise ValueError(
                f"long_axis must be one of {LONG_AXES}, not {self.long_axis!r}"
            )
        known = {m.name for m in self.materials}
        for zone in self.zones:
            for ply in zone.plies:
                if ply.material not in known:
                    raise ValueError(
                        f"ply material {ply.material!r} not in "
                        f"materials {sorted(known)}"
                    )

    @classmethod
    def uniform(
        cls,
        plies: list[Ply] | tuple[Ply, ...],
        *,
        long_axis: str,
        elset: str = "blade",
        material: EngineeringConstants = PLACEHOLDER_CFRP,
    ) -> Layup:
        """Single-zone layup (one ELSET, one material library entry)."""
        fixed = tuple(Ply(p.thickness, p.angle_deg, material.name) for p in plies)
        return cls(
            materials=(material,),
            zones=(ZoneLayup(elset, fixed),),
            long_axis=long_axis,
        )

    def angles(self) -> tuple[float, ...]:
        """Distinct canonical ply angles, in first-seen order."""
        seen: dict[float, None] = {}
        for zone in self.zones:
            for ply in zone.plies:
                seen.setdefault(canonical_angle(ply.angle_deg), None)
        return tuple(seen)

    def to_inp(self) -> str:
        angles = self.angles()
        names = [orientation_name(a) for a in angles]
        if len(set(names)) != len(angles):
            # Should be unreachable while names are derived from the canonical
            # angle, and checked anyway: the failure it guards is a deck that
            # solves and reports a ply angle nobody asked for.
            raise ValueError(
                f"orientation names collide across distinct angles: "
                f"{sorted(zip(angles, names, strict=True))}"
            )
        blocks: list[str] = [m.to_inp() for m in self.materials]
        blocks += [orientation_card(a, long_axis=self.long_axis) for a in angles]
        for zone in self.zones:
            lines = [f"*SHELL SECTION, COMPOSITE, ELSET={zone.elset}"]
            # Bottom (-z) ply first: ccx reads the first ply line as the -z ply.
            for ply in zone.plies:
                ori = orientation_name(ply.angle_deg)
                lines.append(f"{ply.thickness}, , {ply.material}, {ori}")
            blocks.append("\n".join(lines))
        return "\n".join(blocks)


def cross_ply_symmetric(
    ply_thickness: float,
    *,
    long_axis: str,
    elset: str = "blade",
    material: EngineeringConstants = PLACEHOLDER_CFRP,
) -> Layup:
    """[0/90]s with four plies of equal thickness, bottom -> top."""
    t = ply_thickness
    return Layup.uniform(
        [
            Ply(t, 0.0),
            Ply(t, 90.0),
            Ply(t, 90.0),
            Ply(t, 0.0),
        ],
        long_axis=long_axis,
        elset=elset,
        material=material,
    )
