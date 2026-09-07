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
    """``(Ef_1, Ef_2)`` in MPa from the bending stiffness of a laminate.

    Directions 1 and 2 are the **laminate** axes -- 1 along the 0-degree ply --
    not global x and y. See ``laminate_summary``.

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


def is_balanced(a: np.ndarray) -> bool:
    """``A16 == A26 == 0`` to within float noise -- ACP's definition of balance.

    Read off the membrane stiffness rather than by pairing plies, because
    pairing gets it wrong in both directions and this report exists to be laid
    beside ACP's:

    - **False balanced.** Pairing on angle and thickness ignores *which
      lamina*, so ``[+30 cfrp / -30 cfrp_woven]`` of equal thickness cancels in
      the bookkeeping and not in ``A``: measured ``A16`` is 21% of ``A11``.
    - **False unbalanced.** Fibre direction has period 180 degrees, so ``135``
      and ``-45`` are the same ply. ``qbar`` knows that; a pairing rule that
      compares raw angles does not, and calls ``[45/135]`` unbalanced when its
      ``A16`` is exactly zero.

    It also handles a weave for free, with no material metadata. ``Q11 == Q22``
    reduces ``Qbar16`` to ``(Q11 - Q12 - 2 Q66) * cs (c^2 - s^2)``, which
    vanishes at 0, 45 and 90 degrees -- so a woven ply laid at 45 needs no
    partner, while one at 22.5 genuinely is unbalanced, its tows lying at 22.5
    and 112.5. An earlier version took a ``kinds`` map and called *every* woven
    ply self-balancing, which is right only at those three angles.
    """
    scale = max(abs(a[0, 0]), abs(a[1, 1]), 1e-300)
    return bool(
        abs(a[0, 2]) <= COUPLING_REL_TOL * scale
        and abs(a[1, 2]) <= COUPLING_REL_TOL * scale
    )


def laminate_summary(layup: Layup) -> list[dict[str, object]]:
    """One row per zone: thickness, A/B/D terms, flags -- the ACP comparison.

    Flat keys rather than nested matrices so this drops straight into a CSV and
    a dataframe. ``a11`` .. ``d66`` are the six independent terms of each 3x3.

    Axis naming: ``1`` is the laminate's 0-degree direction, which is the
    global axis ``layup.long_axis`` names. The columns are ``ef_1``/``ef_2``,
    not ``ef_x``/``ef_y``, because for ``long_axis="y"`` -- half this repo's
    cases -- direction 1 *is* global y, and a column called ``ef_x`` would be
    read as the global-x modulus. ``long_axis`` rides along in every row so the
    mapping travels with the CSV.
    """
    materials = {m.name: m for m in layup.materials}
    rows: list[dict[str, object]] = []
    for zone in layup.zones:
        a, b, d = abd(zone, materials)
        thickness = zone.thickness
        ef_1, ef_2 = flexural_moduli(d, thickness)
        row: dict[str, object] = {
            "zone": zone.elset,
            "n_plies": len(zone.plies),
            "thickness_mm": thickness,
            "stack": "[" + "/".join(f"{p.angle_deg:g}" for p in zone.plies) + "]",
            "symmetric": is_symmetric(b, a),
            "balanced": is_balanced(a),
            "bend_twist_sign": bend_twist_sign(d),
            "long_axis": layup.long_axis,
            "ef_1_mpa": ef_1,
            "ef_2_mpa": ef_2,
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
