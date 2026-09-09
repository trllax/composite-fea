"""Tip twist stiffness: a comparative torsional metric on the buckled blade.

The tip-weight bench (``cases/fin_test_3``) drives a blade past its Euler load
and reports flex and kick point. Both are bending quantities; neither sees how
the loaded blade resists a *twist*, which is what biax / +/-45 / transverse plies
(``D66``, ``A66``) buy and what separates two layups that match on flex and kick.

This module builds the node sets and reads the result for a small twist
perturbation appended to that solve as two extra ``*STATIC`` steps: a force
couple ``+F`` / ``-F`` on the leading / trailing chord extremes of the tip edge
(``deck.twist_couple_body``) twists the tip section about its own tangent (global
z at a 90 degree tip tangent). A couple rather than a prescribed rotation because
an unsymmetric layup has *rolled* the tip edge by the time it is buckled over, so
its nodes sit at unknown non-uniform axial positions -- a prescribed total
displacement there cannot be formed without first solving, but a force is
incremental by construction.

``k_twist = couple / twist_angle`` in N.mm/rad. The twist angle is the ``U`` swing
between the ``+F`` and ``-F`` steps divided by the *deformed* chord width (so the
swing removes the roll's contribution to the angle and the deformed width removes
its effect on the lever arm). It is **comparative only**: no closed form,
absolute value not meaningful, larger is stiffer. The couple is a moment about
global z, which is the twist axis only near a 90 degree tip tangent -- away from
there it also drives lateral bending, so ``k_twist`` is comparable only across
layups that reach a similar ``theta_at_probe``, which is recorded with every
result.
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from compfea.geometry import Mesh
from compfea.kick import axis_index, read_sets_disp
from compfea.run import parse_dat_totals

LE_NAME = "twist_le"
TE_NAME = "twist_te"

_ARM_MISMATCH_TOL = 0.20


def tip_chord_nsets(
    mesh: Mesh,
    *,
    long_axis: str = "y",
    rib_nset: str = "far_face",
    chord_frac_min: float = 0.55,
    k: int = 4,
    exclude: Iterable[int] | None = None,
    le_name: str = LE_NAME,
    te_name: str = TE_NAME,
) -> dict[str, tuple[int, ...]]:
    """Split the tip edge into leading- and trailing-edge node sets of equal size.

    The repo's first chordwise node grouping -- every other node set here
    collapses to mid-chord. From ``mesh.nsets[rib_nset]`` (the tip edge) minus
    ``exclude``, keep the chord extremes (``|x - x_mid| >= chord_frac_min *
    max|x - x_mid|``), split by ``sign(x - x_mid)``, take up to ``k`` per side
    nearest the extreme, then trim both to the same count so a ``+F`` / ``-F``
    couple applied one force per node has exactly zero net force.

    ``exclude`` drops node ids -- the ``clip`` patch that carries the held axial
    drive, and ``fixed_end``. Raises if either edge has fewer than two nodes or
    the two edges' mean arms differ by more than 20% (a lopsided pick whose
    couple would not be centred on the chord line).
    """
    chord_axis = 1 - axis_index(long_axis)
    drop = set(exclude or ())

    rib = [n for n in mesh.nsets[rib_nset] if n not in drop]
    if len(rib) < 4:
        raise ValueError(
            f"tip_chord_nsets: {rib_nset!r} has {len(rib)} usable nodes after "
            "exclude (need >= 4); refine the chord mesh"
        )

    chord = {n: mesh.nodes[n][chord_axis] for n in rib}
    x_mid = 0.5 * (min(chord.values()) + max(chord.values()))
    half = max(abs(c - x_mid) for c in chord.values())
    if half <= 0.0:
        raise ValueError("tip_chord_nsets: tip edge has no chord extent")

    le = sorted(
        (n for n in rib if chord[n] - x_mid >= chord_frac_min * half),
        key=lambda n: -(chord[n] - x_mid),
    )[:k]
    te = sorted(
        (n for n in rib if x_mid - chord[n] >= chord_frac_min * half),
        key=lambda n: -(x_mid - chord[n]),
    )[:k]
    n_side = min(len(le), len(te))
    if n_side < 2:
        raise ValueError(
            f"tip_chord_nsets: chord split gave le={len(le)} te={len(te)} nodes "
            f"(need >= 2 each); lower chord_frac_min or refine the chord mesh"
        )
    le, te = le[:n_side], te[:n_side]

    arm_le = sum(chord[n] - x_mid for n in le) / n_side
    arm_te = sum(x_mid - chord[n] for n in te) / n_side
    if abs(arm_le - arm_te) / (0.5 * (arm_le + arm_te)) > _ARM_MISMATCH_TOL:
        raise ValueError(
            f"tip_chord_nsets: le/te mean arms {arm_le:.2f}/{arm_te:.2f} mm differ "
            f"by more than {_ARM_MISMATCH_TOL * 100:.0f}%; the chord pick is lopsided"
        )

    return {le_name: tuple(sorted(le)), te_name: tuple(sorted(te))}


def section_arms(
    mesh: Mesh, le_nset: str, te_nset: str, *, long_axis: str = "y"
) -> tuple[float, float, float]:
    """``(x_mid, arm_le, arm_te)`` from undeformed chord coords; arms are positive.

    ``x_mid`` is the midpoint of the chord extent spanned by the two edge sets --
    the line the couple twists the section about. ``arm_le`` / ``arm_te`` are the
    mean chord offsets of the leading / trailing edge nodes from it.
    """
    chord_axis = 1 - axis_index(long_axis)
    le = [mesh.nodes[n][chord_axis] for n in mesh.nsets[le_nset]]
    te = [mesh.nodes[n][chord_axis] for n in mesh.nsets[te_nset]]
    x_mid = 0.5 * (min(le + te) + max(le + te))
    arm_le = sum(c - x_mid for c in le) / len(le)
    arm_te = sum(x_mid - c for c in te) / len(te)
    return x_mid, arm_le, arm_te


def couple_nmm(
    mesh: Mesh, le_nset: str, te_nset: str, *, force: float, long_axis: str = "y"
) -> float:
    """The moment about the chord line of ``+force`` on every ``le_nset`` node
    and ``-force`` on every ``te_nset`` node, on the *undeformed* mesh:
    ``|force| * sum |x_n - x_mid|``.

    ``twist_stiffness`` recomputes this on the deformed chord so a rolled tip
    section is not credited a longer lever arm than it has; this stays as the
    undeformed reference.
    """
    chord_axis = 1 - axis_index(long_axis)
    x_mid, _, _ = section_arms(mesh, le_nset, te_nset, long_axis=long_axis)
    arms = sum(
        abs(mesh.nodes[n][chord_axis] - x_mid)
        for n in (*mesh.nsets[le_nset], *mesh.nsets[te_nset])
    )
    return abs(force) * arms


@dataclass(frozen=True)
class TwistResult:
    """The comparative tip twist stiffness and its self-checks.

    ``k_twist_nmm_per_rad`` is ``None`` when the swing was below resolution or
    non-finite -- a ranker must treat that as "no twist number", not as an
    infinitely stiff design.
    """

    k_twist_nmm_per_rad: float | None
    force_n: float
    couple_nmm: float
    chord_span_def_mm: float
    phi_deg: float
    phi_plus_deg: float
    phi_minus_deg: float
    theta_at_probe_deg: float
    le_te_asymmetry: float
    rf_le_n: float
    rf_te_n: float
    heel_rf_swing_n: float
    warning: str | None


def _disp_at(disp: dict, nset: str, t: float) -> tuple[float, float, float]:
    per_time = disp[nset]
    hit = [v for tt, v in per_time.items() if abs(tt - t) <= 1e-6]
    if len(hit) != 1:
        raise ValueError(
            f"twist read: {len(hit)} U blocks for {nset!r} at time {t} "
            "(need exactly 1); did the drive step print this set, and did the "
            "twist step reach its end time?"
        )
    return hit[0]


def _rf_at(totals, nset: str, comp: str, t: float) -> float:
    hit = totals[
        (totals["nset"] == nset.lower()) & ((totals["time"] - t).abs() <= 1e-6)
    ]
    if len(hit) != 1:
        raise ValueError(
            f"twist read: {len(hit)} RF rows for {nset!r} at time {t} (need 1)"
        )
    return float(hit[comp].iloc[0])


def twist_stiffness(
    dat_path: str | Path,
    mesh: Mesh,
    *,
    le_nset: str = LE_NAME,
    te_nset: str = TE_NAME,
    read_nset: str = "fixed_end",
    drive_end_time: float,
    plus_time: float,
    minus_time: float,
    force_n: float,
    theta_at_probe_deg: float,
    n_side: int | None = None,
    long_axis: str = "y",
    asym_tol: float = 0.35,
    theta_lo_deg: float = 45.0,
    theta_hi_deg: float = 135.0,
    phi_lo_deg: float = 0.03,
    phi_hi_deg: float = 3.0,
) -> TwistResult:
    """Read ``k_twist`` from a solved ``+F`` / ``-F`` twist-couple swing.

    ``drive_end_time`` / ``plus_time`` / ``minus_time`` are the total times at
    the end of the drive step and the two couple steps. ``force_n`` is the
    per-node force magnitude of the ``+F`` step; ``n_side`` the node count per
    edge (defaults to ``len(mesh.nsets[le_nset])``). ``theta_at_probe_deg`` is
    the tip tangent at ``drive_end_time`` (from the caller's canonical
    clip->band fit) -- stored and range-checked, not recomputed.

    The tip-edge ``U`` at ``drive_end_time`` is the baseline the two couple
    states are measured against (so the drive step must print ``le_nset`` /
    ``te_nset``). Everything is done on the *deformed* chord, whose width between
    the two edge groups is ``s = (x_le + ux_le) - (x_te + ux_te)`` at the
    baseline::

        phi(t)  = ((uy_le(t) - uy_le(0)) - (uy_te(t) - uy_te(0))) / s   rad about z
        phi     = 1/2 * (phi(+F) - phi(-F))                             half-swing
        couple  = |force_n| * n_side * s
        k_twist = couple / |phi|

    ``k_twist`` about global z -- which is the twist axis only when the tip
    tangent is ~90 degrees; away from there the couple also drives lateral
    bending, so the number is comparable only across layups that reach a
    *similar* ``theta_at_probe_deg``. Returned ``None`` if the swing is below
    ``phi_lo_deg`` or non-finite. The blade never reaches here with no positive
    torsional stiffness -- that step does not converge.
    """
    axis = axis_index(long_axis)
    chord_axis = 1 - axis
    comp = ("fx", "fy", "fz")[axis]
    if n_side is None:
        n_side = len(mesh.nsets[le_nset])

    totals = parse_dat_totals(dat_path)
    disp = read_sets_disp(dat_path, [le_nset, te_nset])

    x_mid, arm_le_u, arm_te_u = section_arms(
        mesh, le_nset, te_nset, long_axis=long_axis
    )
    x_le_u = x_mid + arm_le_u
    x_te_u = x_mid - arm_te_u

    base_le = _disp_at(disp, le_nset, drive_end_time)
    base_te = _disp_at(disp, te_nset, drive_end_time)
    # deformed chord width between the edge groups at the baseline
    span = (x_le_u + base_le[chord_axis]) - (x_te_u + base_te[chord_axis])
    if span <= 0.0:
        raise ValueError(
            f"twist read: deformed chord span {span:.3f} mm <= 0; the tip "
            "section has rolled past edge-on, the couple is not a twist here"
        )
    couple = abs(force_n) * n_side * span

    def phi_at(t: float) -> tuple[float, float, float]:
        d_le = _disp_at(disp, le_nset, t)[axis] - base_le[axis]
        d_te = _disp_at(disp, te_nset, t)[axis] - base_te[axis]
        return (d_le - d_te) / span, d_le, d_te

    phi_plus, dle_p, dte_p = phi_at(plus_time)
    phi_minus, _, _ = phi_at(minus_time)
    phi_half = 0.5 * (phi_plus - phi_minus)
    phi_deg = math.degrees(abs(phi_half))

    k_twist: float | None
    if math.isfinite(phi_half) and phi_deg >= phi_lo_deg:
        k_twist = couple / abs(phi_half)
    else:
        k_twist = None

    le_te_asym = (
        abs(dle_p + dte_p) / abs(dle_p - dte_p)
        if abs(dle_p - dte_p) > 1e-12
        else 999.0
    )
    # ccx reports the applied *CLOAD in RF for a free loaded node, so le/te RF is
    # not a free-vs-constrained discriminator -- a constrained twist edge instead
    # shows up as no twist swing, which the phi_lo_deg check below catches. These
    # are kept as reported diagnostics only.
    rf_le = _rf_at(totals, le_nset, comp, plus_time)
    rf_te = _rf_at(totals, te_nset, comp, plus_time)
    heel_swing = _rf_at(totals, read_nset, comp, plus_time) - _rf_at(
        totals, read_nset, comp, minus_time
    )
    applied_side = abs(force_n) * n_side

    warning = None
    if not (theta_lo_deg <= theta_at_probe_deg <= theta_hi_deg):
        warning = (
            f"probe tip tangent {theta_at_probe_deg:.1f} deg is outside "
            f"{theta_lo_deg:g}..{theta_hi_deg:g}; the z-couple is more lateral "
            "bending than twist here -- k_twist is not usable"
        )
    elif k_twist is None:
        warning = (
            f"twist angle {phi_deg:.3f} deg is below {phi_lo_deg} deg or "
            "non-finite -- raise --twist-cload-n; k_twist is None"
        )
    elif phi_deg > phi_hi_deg:
        warning = (
            f"twist angle {phi_deg:.3f} deg exceeds {phi_hi_deg} deg -- lower "
            "--twist-cload-n to keep the perturbation linear"
        )
    elif le_te_asym > asym_tol:
        warning = (
            f"le/te edge-motion asymmetry {le_te_asym:.2f} > {asym_tol}; strong "
            "bend-twist coupling or roll -- the comparison is soft"
        )
    elif abs(heel_swing) > 0.25 * applied_side:
        warning = (
            f"heel RF swing {heel_swing:.1f} N is large next to the applied "
            f"{applied_side:.1f} N/side; the perturbation is not a clean twist"
        )

    return TwistResult(
        k_twist_nmm_per_rad=k_twist,
        force_n=force_n,
        couple_nmm=couple,
        chord_span_def_mm=span,
        phi_deg=phi_deg,
        phi_plus_deg=math.degrees(phi_plus),
        phi_minus_deg=math.degrees(phi_minus),
        theta_at_probe_deg=theta_at_probe_deg,
        le_te_asymmetry=le_te_asym,
        rf_le_n=rf_le,
        rf_te_n=rf_te,
        heel_rf_swing_n=heel_swing,
        warning=warning,
    )


def write_twist_json(res: TwistResult, path: str | Path) -> None:
    """Serialise a :class:`TwistResult` to ``path`` (sorted keys, 2-space indent).

    ``allow_nan=False``: a non-finite slipped into the JSON would be written as
    bare ``Infinity`` / ``NaN``, which strict readers reject.
    """
    Path(path).write_text(
        json.dumps(asdict(res), indent=2, sort_keys=True, allow_nan=False)
    )
