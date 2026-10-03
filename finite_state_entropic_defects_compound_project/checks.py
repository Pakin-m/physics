"""Consistency checks for the finite-state entropic-defect draft.

Dependency-free.  The script checks the analytic thresholds, Kerr invariants and
slow-spin deformation, the Z_p holonomy, and a Rayleigh-Ritz sufficient test for
one Schwarzschild tachyonic parameter point.

All quantities are in internally consistent natural/geometrized units.
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
    """Local eligibility radius from 48 alpha_G M^2/(M_G^2 r^6) = m0^2."""
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


def kerr_gauss_bonnet(M: float, a: float, r: float, theta: float) -> float:
    """Gauss-Bonnet density on vacuum Kerr (equal to Kretschmann scalar)."""
    c = math.cos(theta)
    rho2 = r * r + a * a * c * c
    numerator = r**6 - 15*a*a*r**4*c*c + 15*a**4*r*r*c**4 - a**6*c**6
    return 48.0 * M * M * numerator / rho2**6


def kerr_pontryagin(M: float, a: float, r: float, theta: float) -> float:
    """Kerr Pontryagin density; overall sign depends on orientation convention."""
    c = math.cos(theta)
    rho2 = r * r + a * a * c * c
    return (-96.0 * a * M*M * r * c *
            (r*r - 3*a*a*c*c) * (3*r*r - a*a*c*c) / rho2**6)


def slow_spin_eta(alphaO: float, m0: float, MG: float, alphaG: float, MO: float) -> float:
    if min(alphaO, m0, MG, alphaG, MO) <= 0:
        raise ValueError("positive parameters required")
    return alphaO * m0*m0 * MG**4 / (alphaG*alphaG * MO**6)


def slow_spin_weather_radius(
    M: float, a: float, theta: float,
    m0: float, MG: float, alphaG: float,
    alphaO: float, MO: float,
) -> float:
    """O(a^2) local Kerr eligibility surface from the manuscript."""
    r0 = schwarzschild_activation_radius(M, m0, MG, alphaG)
    eta = slow_spin_eta(alphaO, m0, MG, alphaG, MO)
    return r0 * (1.0 + (a*a * math.cos(theta)**2 / r0**2) * (-3.5 + 6.0*eta))



def exact_kerr_local_mass2(M: float, a: float, r: float, theta: float,
                           m0: float, MG: float, alphaG: float,
                           alphaO: float, MO: float) -> float:
    return (m0*m0
            - alphaG*kerr_gauss_bonnet(M, a, r, theta)/(MG*MG)
            - alphaO*kerr_pontryagin(M, a, r, theta)**2/(MO**6))


def exact_kerr_weather_root(M: float, a: float, theta: float,
                            m0: float, MG: float, alphaG: float,
                            alphaO: float, MO: float) -> float:
    """Root near the Schwarzschild local radius, for a sufficiently small spin."""
    r0 = schwarzschild_activation_radius(M, m0, MG, alphaG)
    lo, hi = 0.5*r0, 1.5*r0
    flo = exact_kerr_local_mass2(M, a, lo, theta, m0, MG, alphaG, alphaO, MO)
    fhi = exact_kerr_local_mass2(M, a, hi, theta, m0, MG, alphaG, alphaO, MO)
    if flo*fhi > 0:
        raise ValueError("root not bracketed")
    for _ in range(120):
        mid = 0.5*(lo+hi)
        fm = exact_kerr_local_mass2(M, a, mid, theta, m0, MG, alphaG, alphaO, MO)
        if flo*fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5*(lo+hi)

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


def simpson_integral(fn, a: float, b: float, n: int = 200_000) -> float:
    """Composite Simpson rule, n forced even."""
    if n % 2:
        n += 1
    h = (b - a) / n
    total = fn(a) + fn(b)
    for i in range(1, n):
        total += (4.0 if i % 2 else 2.0) * fn(a + i*h)
    return total * h / 3.0


def schwarzschild_rayleigh(mu: float, gamma: float, ell: int = 0,
                           q: float = 2.0, xmax: float = 30.0) -> float:
    """Rayleigh quotient for u=(x-1) exp[-q(x-1)].

    x=r/(2M), f=1-1/x, B=l(l+1)/x^2+1/x^3+mu^2-gamma/x^6.
    A negative quotient is a sufficient variational proof of a negative mode
    within the numerical quadrature accuracy.
    """
    if mu < 0 or gamma < 0 or ell < 0 or q <= 0 or xmax <= 1:
        raise ValueError("invalid parameters")

    eps = 1e-8
    lo = 1.0 + eps

    def u(x: float) -> float:
        z = x - 1.0
        return z * math.exp(-q*z)

    def ux(x: float) -> float:
        z = x - 1.0
        return math.exp(-q*z) * (1.0 - q*z)

    def numerator_integrand(x: float) -> float:
        f = 1.0 - 1.0/x
        B = ell*(ell+1.0)/x**2 + 1.0/x**3 + mu*mu - gamma/x**6
        return f*ux(x)**2 + B*u(x)**2

    def denominator_integrand(x: float) -> float:
        f = 1.0 - 1.0/x
        return u(x)**2 / f

    num = simpson_integral(numerator_integrand, lo, xmax)
    den = simpson_integral(denominator_integrand, lo, xmax)
    return num / den


def run_checks() -> None:
    # NFW local threshold identity.
    for ratio in (1e-6, 0.1, 1.0, 10.0, 1e6):
        x = nfw_activation_x(rho_s=ratio, rho_c=1.0)
        residual = x * (1.0 + x) ** 2 - ratio
        assert abs(residual) <= 1e-10 * max(1.0, ratio)

    # Schwarzschild local threshold identity and M^(1/3) scaling.
    M, m0, MG, alphaG = 3.0, 0.2, 5.0, 0.7
    r = schwarzschild_activation_radius(M, m0, MG, alphaG)
    lhs = alphaG * gauss_bonnet_schwarzschild(M, r) / (MG * MG)
    assert math.isclose(lhs, m0 * m0, rel_tol=1e-12, abs_tol=1e-12)
    r8 = schwarzschild_activation_radius(8.0 * M, m0, MG, alphaG)
    assert math.isclose(r8 / r, 2.0, rel_tol=1e-12)
    outside_horizon(M, m0, MG, alphaG)

    # Kerr invariants reduce correctly: P=0 for Schwarzschild and on equator;
    # G becomes 48 M^2/r^6 at a=0.
    th = 0.73
    assert math.isclose(kerr_gauss_bonnet(M, 0.0, r, th),
                        gauss_bonnet_schwarzschild(M, r), rel_tol=1e-14)
    assert abs(kerr_pontryagin(M, 0.0, r, th)) < 1e-15
    assert abs(kerr_pontryagin(M, 0.7, r, math.pi/2.0)) < 1e-12

    # Numerical small-spin coefficients against exact Kerr invariants.
    Mx, rx, ax, tx = 1.3, 7.0, 1e-4, 0.41
    c2 = math.cos(tx)**2
    G0 = 48.0*Mx*Mx/rx**6
    G2_expected = -1008.0*Mx*Mx*c2/rx**8
    G2_numeric = (kerr_gauss_bonnet(Mx, ax, rx, tx)-G0)/(ax*ax)
    assert math.isclose(G2_numeric, G2_expected, rel_tol=2e-6)
    P2_expected = 82944.0*Mx**4*c2/rx**14
    P2_numeric = kerr_pontryagin(Mx, ax, rx, tx)**2/(ax*ax)
    assert math.isclose(P2_numeric, P2_expected, rel_tol=2e-6)


    # Slow-spin weather-surface formula agrees with the exact local root at O(a^2).
    pars = dict(M=1.0, a=1e-2, theta=0.4, m0=0.1, MG=1.0, alphaG=1.0, alphaO=0.2, MO=1.0)
    rex = exact_kerr_weather_root(**pars)
    rap = slow_spin_weather_radius(**pars)
    assert abs(rex-rap) / rex < 2e-8, (rex, rap)

    # Z_p selection rule: n and n+p have identical discrete holonomy.
    for p in (2, 3, 5, 7, 11):
        for n in range(-2 * p, 2 * p + 1):
            assert abs(ab_phase(n + p, p) - ab_phase(n, p)) < 1e-12
    assert math.isclose(hidden_flux(1, 5, 2.0), math.pi / 5.0, rel_tol=1e-15)

    # Area matching remains diagnostic only.
    vals = [horizon_cell_area(p) for p in (2, 3, 5, 7)]
    assert vals == sorted(vals) and vals[0] > 0

    # True global instability: one explicit Rayleigh-Ritz sufficient point.
    rq = schwarzschild_rayleigh(mu=0.2, gamma=12.0, ell=0, q=2.0)
    assert rq < -0.06, rq

    print("All analytic/numerical consistency checks passed.")
    print("Example NFW eligibility x for rho_s/rho_c=10:", nfw_activation_x(10.0, 1.0))
    print("Example Schwarzschild local eligibility radius:", r)
    print("Rayleigh quotient (mu=0.2, gamma=12, l=0):", rq)
    print("Slow-spin polar-shape sign flips at eta_Omega = 7/12 =", 7.0/12.0)
    print("Example Z_5 phase for n=1:", ab_phase(1, 5))


if __name__ == "__main__":
    run_checks()
