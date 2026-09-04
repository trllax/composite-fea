"""Layup design vectors -> CalculiX materials, orientations, COMPOSITE sections.

Ply angles become *ORIENTATION names (ccx field 4), never bare degrees.
Units: mm, N, tonne, MPa, s (density tonne/mm^3).
"""

from __future__ import annotations

import math
from dataclasses import dataclass


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
    """Stable *ORIENTATION name for a fiber angle in degrees."""
    # Avoid "-0" and long floats; keep one decimal when needed.
    if abs(angle_deg - round(angle_deg)) < 1e-9:
        tag = f"{int(round(angle_deg))}"
    else:
        tag = f"{angle_deg:.1f}".replace("-", "m").replace(".", "p")
        return f"ori_a{tag}"
    if angle_deg < 0:
        return f"ori_m{abs(int(round(angle_deg)))}"
    return f"ori_p{tag}"


def orientation_card(
    angle_deg: float,
    *,
    zero_along: str = "y",
) -> str:
    """Rectangular orientation: 0° along +x or +y in the shell plane, +z normal.

    ``zero_along='y'`` matches the strip/U-bend convention (long axis +y).
    ``zero_along='x'`` matches CLAUDE.md's "+x along the part long axis".
    """
    name = orientation_name(angle_deg)
    phi = math.radians(angle_deg)
    s, c = math.sin(phi), math.cos(phi)
    if zero_along == "y":
        # phi=0 -> (0,1,0); phi=90 -> (1,0,0)
        a1, a2, a3 = s, c, 0.0
        b1, b2, b3 = -c, s, 0.0
    elif zero_along == "x":
        # phi=0 -> (1,0,0); phi=90 -> (0,1,0)
        a1, a2, a3 = c, s, 0.0
        b1, b2, b3 = -s, c, 0.0
    else:
        raise ValueError("zero_along must be 'x' or 'y'")
    return "\n".join(
        [
            f"*ORIENTATION, NAME={name}, SYSTEM=RECTANGULAR",
            f"{a1:.8f}, {a2:.8f}, {a3:.8f}, {b1:.8f}, {b2:.8f}, {b3:.8f}",
        ]
    )


@dataclass(frozen=True)
class Layup:
    """Materials + zone stacks ready to emit into a deck."""

    materials: tuple[EngineeringConstants, ...]
    zones: tuple[ZoneLayup, ...]
    zero_along: str = "y"

    def __post_init__(self) -> None:
        if not self.materials:
            raise ValueError("at least one material is required")
        if not self.zones:
            raise ValueError("at least one zone layup is required")
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
        elset: str = "blade",
        material: EngineeringConstants = PLACEHOLDER_CFRP,
        zero_along: str = "y",
    ) -> Layup:
        """Single-zone layup (one ELSET, one material library entry)."""
        fixed = tuple(
            Ply(p.thickness, p.angle_deg, material.name) if isinstance(p, Ply) else p
            for p in plies
        )
        # Ensure material name on plies
        fixed = tuple(
            Ply(p.thickness, p.angle_deg, material.name) for p in fixed
        )
        return cls(
            materials=(material,),
            zones=(ZoneLayup(elset, fixed),),
            zero_along=zero_along,
        )

    def angles(self) -> tuple[float, ...]:
        seen: dict[float, None] = {}
        for zone in self.zones:
            for ply in zone.plies:
                key = round(ply.angle_deg, 9)
                seen.setdefault(key, None)
        return tuple(seen.keys())

    def to_inp(self) -> str:
        blocks: list[str] = [m.to_inp() for m in self.materials]
        for angle in self.angles():
            blocks.append(orientation_card(angle, zero_along=self.zero_along))
        for zone in self.zones:
            lines = [f"*SHELL SECTION, COMPOSITE, ELSET={zone.elset}"]
            for ply in zone.plies:
                ori = orientation_name(ply.angle_deg)
                lines.append(f"{ply.thickness}, , {ply.material}, {ori}")
            blocks.append("\n".join(lines))
        return "\n".join(blocks)


def cross_ply_symmetric(
    ply_thickness: float,
    *,
    elset: str = "blade",
    material: EngineeringConstants = PLACEHOLDER_CFRP,
    zero_along: str = "y",
) -> Layup:
    """[0/90]s with four plies of equal thickness."""
    t = ply_thickness
    return Layup.uniform(
        [
            Ply(t, 0.0),
            Ply(t, 90.0),
            Ply(t, 90.0),
            Ply(t, 0.0),
        ],
        elset=elset,
        material=material,
        zero_along=zero_along,
    )
