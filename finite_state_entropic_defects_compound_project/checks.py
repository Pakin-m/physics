"""Minimal numerical checks for the finite-state entropic-defect draft.

No external dependencies.  The script verifies the algebraic threshold equations
used in the manuscript and the Z_p periodicity of the Aharonov--Bohm phase.
All quantities are in consistent natural/geometrized units chosen by the caller.
"""

from __future__ import annotations

import cmath
import math


def nfw_activation_x(rho_s: float, rho_c: float, tol: float = 1e-13) -> float:
    """Solve x(1+x)^2 = rho_s/rho_c for the unique positive x."""
    if rho_s <= 0 or rho_c <= 0:
        raise ValueError("rho_s and rho_c must be positive")
    target = rho_s / rho_c

    def f(x: float) -> float:
        return x * (1.0 + x) ** 2 - target

    lo, hi = 0.0, max(1.0, target ** (1.0 / 3.0) + 1.0)
    while f(hi) < 0:
        hi *= 2.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol * max(1.0, mid):
            break
    return 0.5 * (lo + hi)


def schwarzschild_activation_radius(M: float, m0: float, MG: float, alphaG: float) -> float:
    """r_G from 48 alpha_G M^2/(M_G^2 r^6) = m0^2."""
    if min(M, m0, MG, alphaG) <= 0:
        raise ValueError("arguments must be positive")
    return (48.0 * alphaG * M * M / (m0 * m0 * MG * MG)) ** (1.0 / 6.0)


def gauss_bonnet_schwarzschild(M: float, r: float) -> float:
    return 48.0 * M * M / r**6


def outside_horizon(M: float, m0: float, MG: float, alphaG: float) -> bool:
    direct = schwarzschild_activation_radius(M, m0, MG, alphaG) > 2.0 * M
    inequality = (m0 * m0 * MG * MG * M**4) < (0.75 * alphaG)
    assert direct == inequality
    return direct


def hidden_flux(n: int, p: int, gE: float) -> float:
    if p < 2 or gE <= 0:
        raise ValueError("p >= 2 and gE > 0 required")
    return 2.0 * math.pi * n / (p * gE)


def ab_phase(n: int, p: int) -> complex:
    if p < 2:
        raise ValueError("p >= 2 required")
    return cmath.exp(2j * math.pi * n / p)


def horizon_cell_area(p: int, ell_p: float = 1.0) -> float:
    """Consistency matching a_* = 4 ell_P^2 ln p; not a new derivation."""
    if p < 2 or ell_p <= 0:
        raise ValueError("p >= 2 and ell_p > 0 required")
    return 4.0 * ell_p * ell_p * math.log(p)


def run_checks() -> None:
    # NFW threshold identity.
    for ratio in (1e-6, 0.1, 1.0, 10.0, 1e6):
        x = nfw_activation_x(rho_s=ratio, rho_c=1.0)
        residual = x * (1.0 + x) ** 2 - ratio
        assert abs(residual) <= 1e-10 * max(1.0, ratio)

    # Schwarzschild threshold identity and M^(1/3) scaling.
    M, m0, MG, alphaG = 3.0, 0.2, 5.0, 0.7
    r = schwarzschild_activation_radius(M, m0, MG, alphaG)
    lhs = alphaG * gauss_bonnet_schwarzschild(M, r) / (MG * MG)
    assert math.isclose(lhs, m0 * m0, rel_tol=1e-12, abs_tol=1e-12)
    r8 = schwarzschild_activation_radius(8.0 * M, m0, MG, alphaG)
    assert math.isclose(r8 / r, 2.0, rel_tol=1e-12)
    outside_horizon(M, m0, MG, alphaG)

    # Z_p selection rule: n and n+p have identical discrete holonomy.
    for p in (2, 3, 5, 7, 11):
        for n in range(-2 * p, 2 * p + 1):
            assert abs(ab_phase(n + p, p) - ab_phase(n, p)) < 1e-12

    # Flux is reduced by p for the minimal n=1 defect.
    assert math.isclose(hidden_flux(1, 5, 2.0), math.pi / 5.0, rel_tol=1e-15)

    # Area matching is positive and monotone in p.
    vals = [horizon_cell_area(p) for p in (2, 3, 5, 7)]
    assert vals == sorted(vals) and vals[0] > 0

    print("All analytic consistency checks passed.")
    print("Example NFW activation x for rho_s/rho_c=10:", nfw_activation_x(10.0, 1.0))
    print("Example Schwarzschild activation radius:", r)
    print("Example Z_5 phase for n=1:", ab_phase(1, 5))


if __name__ == "__main__":
    run_checks()
