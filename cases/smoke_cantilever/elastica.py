"""Closed-form large-deflection cantilever under a transverse tip load.

The reference the NLGEOM gate is checked against. It lives in the case rather
than in `src/compfea/` on purpose: a reference that imports the code under test,
or that shares a module with it, can agree with it and still be wrong together.
Only `numpy` is used here, and nothing in this file knows what CalculiX is.

Inextensible Euler-Bernoulli beam, clamped at s = 0, transverse point load P at
the free end whose direction stays fixed as the beam bends. With phi(s) the
tangent angle from the undeformed axis and phi0 the tip value:

    EI phi'' = -P cos(phi)                       moment balance
    (phi')^2 = (2P/EI) (sin(phi0) - sin(phi))    integrated, phi'(L) = 0

    L     = sqrt(EI/2P) * D,   D = int_0^phi0 dphi / sqrt(sin phi0 - sin phi)
    delta = sqrt(EI/2P) * N,   N = int_0^phi0 sin(phi) dphi / sqrt(...)

so delta/L = N/D depends on phi0 alone, and P = D^2 EI / (2 L^2). The integrand
is singular at phi0; the substitution phi = phi0 - u^2 removes it exactly
(sin phi0 - sin phi -> u^2 cos phi0 as u -> 0), leaving something smooth that
Gauss-Legendre nails in a few dozen nodes. It degenerates as phi0 -> 90 deg,
where cos phi0 -> 0; this case drives 26 deg.

Units are whatever the caller uses consistently: N and mm here, so EI is N.mm^2.
"""

from __future__ import annotations

import numpy as np

# Above this tip angle the phi = phi0 - u^2 substitution loses its footing, and
# a vertical tip load cannot reach 90 deg anyway (delta/L ~ 0.92 gets 89.9 deg).
_MAX_PHI0 = np.radians(85.0)


def _integrals(phi0: float, nodes: int = 200) -> tuple[float, float]:
    """(D, N) for tip angle ``phi0`` in radians."""
    if not 0.0 < phi0 <= _MAX_PHI0:
        raise ValueError(f"phi0 must be in (0, {np.degrees(_MAX_PHI0):g}) deg")
    x, w = np.polynomial.legendre.leggauss(nodes)
    u_max = np.sqrt(phi0)
    u = 0.5 * u_max * (x + 1.0)
    weights = 0.5 * u_max * w
    phi = phi0 - u**2
    integrand = 2.0 * u / np.sqrt(np.sin(phi0) - np.sin(phi))
    return (
        float(np.sum(weights * integrand)),
        float(np.sum(weights * integrand * np.sin(phi))),
    )


def deflection_ratio(phi0: float) -> float:
    """delta/L for a given tip angle in radians."""
    d, n = _integrals(phi0)
    return n / d


def tip_angle(delta_over_l: float) -> float:
    """Tip angle in radians that produces ``delta/L``. Monotonic, so bisect."""
    if not 0.0 < delta_over_l <= deflection_ratio(_MAX_PHI0):
        raise ValueError(f"delta/L out of range: {delta_over_l}")
    lo, hi = 1e-9, _MAX_PHI0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if deflection_ratio(mid) < delta_over_l:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def tip_load(delta_over_l: float, *, ei: float, length: float) -> float:
    """Transverse tip load P that deflects the tip to ``delta/L``.

    Reduces to the small-deflection ``P = 3 EI delta / L^3`` as delta/L -> 0 and
    stiffens above it: at delta/L = 0.3 it is 10% higher.
    """
    d, _ = _integrals(tip_angle(delta_over_l))
    return d**2 * ei / (2.0 * length**2)
