"""Closed-form large-deflection cantilever under an axial (compressive) tip load.

The reference the fin tip-weight gate is checked against. Like
``cases/smoke_cantilever/elastica.py`` it lives in the case, not in
``src/compfea/``: a reference that shares a module with the code under test can
be wrong together with it. Only ``numpy`` is used and nothing here knows what
CalculiX is.

## The apparatus this models

A blade is clamped at one end and held so its axis lies along a hung weight's
line of action -- vertical, tip up. The dead weight ``W`` pulls straight down,
i.e. *along* the undeformed axis, in compression. The blade is a slender column
past its buckling load, so it lies over; the operator nudges the tip to pick the
direction. "Tip at 90 degrees" means the tip *tangent* has turned a quarter turn
from vertical to horizontal.

This is not the transverse-tip-load elastica of ``smoke_cantilever`` (curvature
concentrated at the root, load perpendicular to the undeformed axis) and it is
not a circular-arc fold (uniform curvature). It is the clamped-free buckling
column, and its small-angle limit is the Euler load.

## The elastica

Inextensible Euler-Bernoulli beam, arc length ``s`` from the clamp (``s = 0``) to
the free tip (``s = L``). ``phi(s)`` is the tangent angle from the undeformed
axis (the load line), so ``dX/ds = sin(phi)``, ``dY/ds = cos(phi)`` with ``Y``
along the load.

The only external force on the segment ``[s, L]`` is the tip load ``(0, -W)`` at
the tip point, so the internal moment is ``M(s) = -W (X_tip - X(s))``.
Differentiating ``EI phi' = -M_z``:

    EI phi'' = -W sin(phi)        =>   phi'' + (W/EI) sin(phi) = 0

the pendulum equation, as a buckling column must give. Clamped at the root,
``phi(0) = 0``; free at the tip, ``phi'(L) = 0``. First integral with
``phi'(L) = 0`` and tip angle ``phi_L = phi(L)``:

    (phi')^2 = (2W/EI) (cos(phi) - cos(phi_L))

    L     = sqrt(EI / 2W) * I,   I = int_0^{phi_L} dphi / sqrt(cos phi - cos phi_L)
    Y_tip = sqrt(EI / 2W) * J,   J = int_0^{phi_L} cos(phi) dphi / sqrt(...)

so ``W = EI * I(phi_L)^2 / (2 L^2)``, and the axial shortening
``delta = L - Y_tip`` has ``delta / L = 1 - J / I``, a function of ``phi_L``
alone. As ``phi_L -> 0``, ``I -> pi / sqrt(2)`` and ``W -> pi^2 EI / (4 L^2)``,
the clamped-free Euler load: the small-angle limit is the bifurcation, exactly as
it should be. ``I`` and ``J`` both rise monotonically with ``phi_L`` up to
``pi``, so ``W`` does too -- there is no limit point below 180 degrees of tip
rotation and a dead weight can hold any tip angle up to there.

The integrand is singular at the upper limit
(``cos phi - cos phi_L -> sin(phi_L) (phi_L - phi)``); the substitution
``phi = phi_L - u^2`` removes it, leaving a smooth integrand Gauss-Legendre
handles in a few dozen nodes. It degenerates as ``phi_L -> 180 deg`` where
``sin(phi_L) -> 0``; the design point is 90 deg, well clear.

Units are the caller's, used consistently: N and mm here, so ``EI`` is N.mm^2 and
``W`` is N.
"""

from __future__ import annotations

import numpy as np

# The phi = phi_L - u^2 substitution needs sin(phi_L) bounded away from zero.
# 90 deg (tip tangent horizontal) is the design point; allow headroom past it for
# the fin run, but refuse the 180 deg neighbourhood where the quadrature
# degenerates.
_MAX_PHI_L = np.radians(150.0)
# Below this, cos(phi) - cos(phi_L) cancels to round-off and the integrals return
# nan/inf. The bifurcation itself (phi_L -> 0) is the Euler load; use euler_load()
# there, not this. delta/L = 1e-6 already gives phi_L ~ 0.6 deg, well clear.
_MIN_PHI_L = np.radians(0.05)


def _integrals(phi_l: float, nodes: int = 200) -> tuple[float, float]:
    """``(I, J)`` for tip angle ``phi_l`` in radians.

    ``I = int_0^{phi_l} dphi / sqrt(cos phi - cos phi_l)``
    ``J = int_0^{phi_l} cos(phi) dphi / sqrt(cos phi - cos phi_l)``

    The ``phi = phi_l - u^2`` substitution folds the ``dphi`` Jacobian and the
    singular denominator into a bounded ``2u / sqrt(cos(phi_l - u^2) - cos phi_l)``
    that tends to ``2 / sqrt(sin phi_l)`` at ``u -> 0``.
    """
    if not _MIN_PHI_L <= phi_l <= _MAX_PHI_L:
        raise ValueError(
            f"phi_l must be in [{np.degrees(_MIN_PHI_L):g}, "
            f"{np.degrees(_MAX_PHI_L):g}] deg; use euler_load() near phi_l = 0"
        )
    x, w = np.polynomial.legendre.leggauss(nodes)
    u_max = np.sqrt(phi_l)
    u = 0.5 * u_max * (x + 1.0)
    weights = 0.5 * u_max * w
    phi = phi_l - u**2
    integrand = 2.0 * u / np.sqrt(np.cos(phi) - np.cos(phi_l))
    return (
        float(np.sum(weights * integrand)),
        float(np.sum(weights * integrand * np.cos(phi))),
    )


def shortening_ratio(phi_l: float) -> float:
    """Axial shortening ``delta / L`` for tip angle ``phi_l`` (radians).

    ``delta`` is how far the tip moves *along* the load line -- the quantity the
    FE model prescribes. Monotonic in ``phi_l``, zero at the bifurcation.
    """
    i, j = _integrals(phi_l)
    return 1.0 - j / i


def tip_angle(delta_over_l: float) -> float:
    """Tip angle (radians) whose axial shortening is ``delta_over_l``. Bisected."""
    lo_ratio = shortening_ratio(_MIN_PHI_L)
    hi = shortening_ratio(_MAX_PHI_L)
    if not lo_ratio <= delta_over_l <= hi:
        raise ValueError(
            f"delta/L out of range [{lo_ratio:g}, {hi:g}]: {delta_over_l}"
        )
    lo, hi_phi = _MIN_PHI_L, _MAX_PHI_L
    for _ in range(100):
        mid = 0.5 * (lo + hi_phi)
        if shortening_ratio(mid) < delta_over_l:
            lo = mid
        else:
            hi_phi = mid
    return 0.5 * (lo + hi_phi)


def axial_load_at_angle(phi_l: float, *, ei: float, length: float) -> float:
    """Compressive tip load ``W`` that bends the tip to ``phi_l`` (radians).

    ``W = EI * I(phi_l)^2 / (2 L^2)``.
    """
    i, _ = _integrals(phi_l)
    return i**2 * ei / (2.0 * length**2)


def axial_load(delta_over_l: float, *, ei: float, length: float) -> float:
    """Compressive tip load ``W`` that shortens the column by ``delta_over_l``.

    As ``delta/L -> 0`` this approaches the clamped-free Euler load
    ``pi^2 EI / (4 L^2)`` and rises monotonically above it.
    """
    return axial_load_at_angle(tip_angle(delta_over_l), ei=ei, length=length)


def euler_load(*, ei: float, length: float) -> float:
    """Clamped-free buckling load ``pi^2 EI / (4 L^2)`` -- the ``delta/L -> 0``
    limit of :func:`axial_load`."""
    return float(np.pi**2 * ei / (4.0 * length**2))


# --------------------------------------------------------------------------
# Curvature along the span -- the reference for the kick-point metric.
#
# From the first integral ``(phi')^2 = (2W/EI)(cos phi - cos phi_L)``:
#
#     kappa(s) = phi'(s) = sqrt(2W/EI) * sqrt(cos phi - cos phi_L)
#     s(phi)   = sqrt(EI/2W) * int_0^phi dpsi / sqrt(cos psi - cos phi_L)
#
# so ``kappa`` is largest at the clamp (``s = 0``, ``phi = 0``) and zero at the
# free tip (``s = L``, ``phi = phi_L``). A uniform column therefore "kicks" at
# the heel; a stiffness taper that thins the tip moves the peak outboard. Same
# ``phi = phi_L - v^2`` substitution as :func:`_integrals` for the s-quadrature.
#
# ``k = sqrt(2W/EI) = I(phi_L)/L`` -- the ``EI`` inside ``W`` cancels, so the
# *shape* of ``kappa(s)`` at a given tip angle does not depend on ``EI`` (only
# ``W`` does). ``ei`` is still taken for signature parity with
# :func:`axial_load_at_angle`; it does not change the result.
# --------------------------------------------------------------------------


def _wavenumber(phi_l: float, *, ei: float, length: float) -> float:
    """``k = sqrt(2W/EI) = I(phi_l) / L`` in 1/mm, so ``kappa = k sqrt(cos..)``.

    ``ei`` is accepted and ignored -- it cancels (see the section note above).
    """
    del ei
    i, _ = _integrals(phi_l)
    return i / length


def curvature_root(phi_l: float, *, ei: float, length: float) -> float:
    """Curvature at the clamp, ``kappa(0) = k sqrt(1 - cos phi_l)`` in 1/mm.

    The maximum of :func:`curvature_profile`; the simplest closed-form target
    for a kick-point check on a uniform column.
    """
    k = _wavenumber(phi_l, ei=ei, length=length)
    return float(k * np.sqrt(1.0 - np.cos(phi_l)))


def curvature_profile(
    phi_l: float, *, ei: float, length: float, n: int = 200
) -> tuple[np.ndarray, np.ndarray]:
    """``(s_mm[n], kappa_1pmm[n])`` from clamp (``s = 0``) to tip (``s = L``).

    ``kappa`` is monotone decreasing: maximal at ``s = 0``, zero at ``s = L``.
    """
    if not _MIN_PHI_L <= phi_l <= _MAX_PHI_L:
        raise ValueError(
            f"phi_l must be in [{np.degrees(_MIN_PHI_L):g}, "
            f"{np.degrees(_MAX_PHI_L):g}] deg; use euler_load() near phi_l = 0"
        )
    if n < 2:
        raise ValueError(f"n must be >= 2, got {n}")
    k = _wavenumber(phi_l, ei=ei, length=length)

    # Fine grid in v = sqrt(phi_l - phi); phi = phi_l - v^2 runs the tip (v = 0)
    # to the clamp (v = v_max). The s integrand 2v / sqrt(cos(phi_l - v^2) -
    # cos phi_l) is bounded, tending to 2 / sqrt(sin phi_l) as v -> 0. The grid
    # always over-resolves the requested n, so refining n really refines.
    v_max = np.sqrt(phi_l)
    v = np.linspace(0.0, v_max, max(4001, 2 * n + 1))
    diff = np.cos(phi_l - v**2) - np.cos(phi_l)
    with np.errstate(invalid="ignore", divide="ignore"):
        g = 2.0 * v / np.sqrt(diff)
    g[0] = 2.0 / np.sqrt(np.sin(phi_l))
    cum = np.concatenate(
        [[0.0], np.cumsum(np.diff(v) * 0.5 * (g[:-1] + g[1:]))]
    )
    # cum(v) = int_phi^{phi_l} dpsi/sqrt(..) = I - F(phi), so s(phi) = F/k =
    # (I - cum)/k. v = 0 -> s = L (tip); v = v_max -> s = 0 (clamp).
    s_of_v = (cum[-1] - cum) / k
    kappa_of_v = k * np.sqrt(np.clip(diff, 0.0, None))

    # Reverse to clamp -> tip (s ascending), then decimate to n points.
    s = s_of_v[::-1]
    kappa = kappa_of_v[::-1]
    idx = np.linspace(0, s.size - 1, n).round().astype(int)
    return s[idx], kappa[idx]


def rotation_median_s(phi_l: float, *, length: float) -> float:
    """Arc length from the clamp to where the tangent angle reaches ``phi_l / 2``.

    The station that splits the tip rotation in half -- the robustness
    cross-check for the kick point. ``s = L * F(phi_l/2) / I``; the integrand is
    singular only at ``phi_l`` so the half-angle integral is taken directly.
    """
    i, _ = _integrals(phi_l)
    x, w = np.polynomial.legendre.leggauss(200)
    half = 0.5 * phi_l
    psi = 0.5 * half * (x + 1.0)
    wts = 0.5 * half * w
    f_half = float(np.sum(wts / np.sqrt(np.cos(psi) - np.cos(phi_l))))
    return length * f_half / i
