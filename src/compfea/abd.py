"""Classical lamination theory on a resolved stack: A, B, D per zone.

This is the cheap half of a cross-check against ANSYS. Every number here is one
ACP also reports for a laminate, and none of it needs a solver or a mesh -- so
you can establish that ccx and ANSYS were handed the *same laminate* before
anyone argues about a force.

Units, and they are a new place to go wrong::

    Q, Qbar   MPa            (N/mm^2)
    A         N/mm           membrane
    B         N              membrane-bending coupling
    D         N*mm           bending

with ply z stations in mm measured from the laminate mid-surface, ply 1 (the
first ``*SHELL SECTION`` line) at ``-t/2``.

Sign convention follows ``layup.orientation_card``: positive angle rotates the
fibre counter-clockwise about +z. That matters for exactly one thing here, and
it is the reason this module is worth having -- see ``bend_twist_sign``.
"""

from __future__ import annotations

import math

import numpy as np

from compfea.layup import (
    EngineeringConstants,
    Layup,
    ZoneLayup,
    _plane_stress_q,
    canonical_angle,
)

#: Below this a coupling term is float noise, not physics. B for a symmetric
#: stack cancels in exact arithmetic; in floating point it lands near 1e-12 of
#: the A scale.
COUPLING_REL_TOL = 1e-9


def qbar(ec: EngineeringConstants, angle_deg: float) -> np.ndarray:
    """Reduced stiffness of one lamina rotated to the laminate axes, MPa.

    3x3 in Voigt order (11, 22, 12). The on-axis Q comes from
    ``layup._plane_stress_q`` rather than being re-derived, so a change to the
    lamina reduction reaches both the deck and this report.
    """
    q11, q12, q22, q66 = _plane_stress_q(ec)
    phi = math.radians(canonical_angle(angle_deg))
    c, s = math.cos(phi), math.sin(phi)
    c2, s2 = c * c, s * s
    c4, s4 = c2 * c2, s2 * s2
    cs = c * s
    return np.array(
        [
            [
                q11 * c4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * s4,
                (q11 + q22 - 4.0 * q66) * s2 * c2 + q12 * (s4 + c4),
                (q11 - q12 - 2.0 * q66) * c2 * cs
                + (q12 - q22 + 2.0 * q66) * s2 * cs,
            ],
            [
                (q11 + q22 - 4.0 * q66) * s2 * c2 + q12 * (s4 + c4),
                q11 * s4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * c4,
                (q11 - q12 - 2.0 * q66) * s2 * cs
                + (q12 - q22 + 2.0 * q66) * c2 * cs,
            ],
            [
                (q11 - q12 - 2.0 * q66) * c2 * cs
                + (q12 - q22 + 2.0 * q66) * s2 * cs,
                (q11 - q12 - 2.0 * q66) * s2 * cs
                + (q12 - q22 + 2.0 * q66) * c2 * cs,
                (q11 + q22 - 2.0 * q12 - 2.0 * q66) * s2 * c2
                + q66 * (s4 + c4),
            ],
        ]
    )


def ply_z_stations(zone: ZoneLayup) -> list[tuple[float, float]]:
    """``(z_bot, z_top)`` per ply in mm, mid-surface origin, ply 1 at -t/2."""
    z = -0.5 * zone.thickness
    out: list[tuple[float, float]] = []
    for ply in zone.plies:
        out.append((z, z + ply.thickness))
        z += ply.thickness
    return out


def abd(
    zone: ZoneLayup, materials: dict[str, EngineeringConstants]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(A, B, D)`` for one zone. N/mm, N, N*mm."""
    a = np.zeros((3, 3))
    b = np.zeros((3, 3))
    d = np.zeros((3, 3))
    for ply, (z0, z1) in zip(zone.plies, ply_z_stations(zone), strict=True):
        try:
            ec = materials[ply.material]
        except KeyError:
            raise KeyError(
                f"zone {zone.elset!r} uses material {ply.material!r}, which is "
                f"not among {sorted(materials)}"
            ) from None
        q = qbar(ec, ply.angle_deg)
        a += q * (z1 - z0)
        b += q * (z1**2 - z0**2) / 2.0
        d += q * (z1**3 - z0**3) / 3.0
    return a, b, d


def flexural_moduli(
    d: np.ndarray, thickness: float
) -> tuple[float, float]:
    """``(Ef_x, Ef_y)`` in MPa from the bending stiffness of a laminate.

    The standard reduction ``12 / (t^3 * d_inv[i, i])``, which is the modulus of
    an equivalent homogeneous plate in bending. Not the same as ``A``'s
    membrane modulus for anything but a homogeneous stack.
    """
    inv = np.linalg.inv(d)
    factor = 12.0 / thickness**3
    return factor / inv[0, 0], factor / inv[1, 1]


def bend_twist_sign(d: np.ndarray) -> int:
    """``sign(D16)`` -- the one check in this repo that sees the angle convention.

    ``layup.py`` states that writing the mirrored form of ``*ORIENTATION``
    "flips the sign of bend-twist coupling", that a cross-ply or balanced
    +/-45 stack cannot detect it at all, and that no force check anywhere here
    will catch it. ``D16`` will: it is zero for those blind stacks and changes
    sign with the convention for any unbalanced one. This is the first quantity
    in the repo that can fail on a mirrored deck.
    """
    scale = max(abs(d[0, 0]), abs(d[1, 1]), 1e-300)
    if abs(d[0, 2]) <= COUPLING_REL_TOL * scale:
        return 0
    return 1 if d[0, 2] > 0 else -1


def is_symmetric(b: np.ndarray, a: np.ndarray) -> bool:
    """``B == 0`` to within float noise, scaled by the membrane stiffness.

    ``B`` has units of N and ``A`` of N/mm, so the comparison needs a length.
    Using ``max|A|`` times a millimetre is crude but it is the right order and
    it does not need the thickness threaded through.
    """
    return bool(np.max(np.abs(b)) <= COUPLING_REL_TOL * max(np.max(np.abs(a)), 1e-300))


def is_balanced(zone: ZoneLayup) -> bool:
    """Every ``+theta`` has a matching ``-theta`` of equal thickness.

    0 and 90 are their own opposites and never unbalance a stack.
    """
    net: dict[float, float] = {}
    for ply in zone.plies:
        angle = canonical_angle(ply.angle_deg)
        if angle in (0.0, 90.0, -90.0):
            continue
        net[abs(angle)] = net.get(abs(angle), 0.0) + math.copysign(
            ply.thickness, angle
        )
    return all(abs(v) <= 1e-12 for v in net.values())


def laminate_summary(layup: Layup) -> list[dict[str, object]]:
    """One row per zone: thickness, A/B/D terms, flags -- the ACP comparison.

    Flat keys rather than nested matrices so this drops straight into a CSV and
    a dataframe. ``a11`` .. ``d66`` are the six independent terms of each 3x3.
    """
    materials = {m.name: m for m in layup.materials}
    rows: list[dict[str, object]] = []
    for zone in layup.zones:
        a, b, d = abd(zone, materials)
        thickness = zone.thickness
        ef_x, ef_y = flexural_moduli(d, thickness)
        row: dict[str, object] = {
            "zone": zone.elset,
            "n_plies": len(zone.plies),
            "thickness_mm": thickness,
            "stack": "[" + "/".join(f"{p.angle_deg:g}" for p in zone.plies) + "]",
            "symmetric": is_symmetric(b, a),
            "balanced": is_balanced(zone),
            "bend_twist_sign": bend_twist_sign(d),
            "ef_x_mpa": ef_x,
            "ef_y_mpa": ef_y,
        }
        for name, matrix in (("a", a), ("b", b), ("d", d)):
            for label, (i, j) in (
                ("11", (0, 0)), ("12", (0, 1)), ("16", (0, 2)),
                ("22", (1, 1)), ("26", (1, 2)), ("66", (2, 2)),
            ):
                row[f"{name}{label}"] = float(matrix[i, j])
        rows.append(row)
    return rows


def beam_ei_free_edge(d: np.ndarray, width_mm: float) -> float:
    """``b (D11 - D12^2 / D22)`` -- a strip with free long edges, N*mm^2.

    The reduction ``CLAUDE.md`` singles out. ``D11 * b`` is the wide-plate form
    and is the wrong one for a strip; the two differ by 0.247% on the
    smoke_cantilever laminate, and taking the wrong one is what once got
    explained away as a solver property.
    """
    return width_mm * (d[0, 0] - d[0, 1] ** 2 / d[1, 1])
