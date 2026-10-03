# Finite-State Entropic Defects — compact compound project

## Current research kernel

The project uses the "universe weather" idea only as a source of invariant questions. It does **not** identify gravity with literal fluid flow or identify an orbit with a spacetime vortex.

The retained EFT has one compact mechanism:

- a hidden charge-`p` Higgs field leaves a residual local `Z_p` gauge phase when it condenses;
- the local Higgs mass is shifted by Gauss-Bonnet curvature, a parity-even rotational-curvature term `Pontryagin^2`, and a dark-sector stress-tensor trace;
- the finite arithmetic is physical only through the exact holonomy `exp(2π i n/p)`;
- bulk "entropic pressure" is calculated from the hidden-sector stress tensor;
- a negative local mass is treated only as **phase eligibility**, not as proof of condensation.

Nearby ingredients are already known: curvature/spin scalarization, quadratic Pontryagin scalar couplings, symmetron/environment-dependent defects, dark-matter vortices, and `Z_p` gauge theory. The novelty hypothesis is only their specific finite-state gauge synthesis plus its derived predictions.

## Files

- `main.tex` — main manuscript and hard novelty/consistency audit.
- `references.bib` — compact peer-reviewed literature basis.
- `checks.py` — dependency-free algebra, Kerr-invariant, slow-spin, and variational-instability checks.
- `README.md` — project map and kill criteria.

## Strongest derived results

### 1. Exact finite-state observable

```text
Phi_E(n) = 2π n/(p g_E)
exp(i g_E Phi_E) = exp(2π i n/p)
```

### 2. Covariant environmental mass

```text
m_E^2 = m_0^2
        - alpha_G G/M_G^2
        - alpha_Omega P^2/M_Omega^6
        - alpha_D D/M_D^2
```

`G` is the Gauss-Bonnet invariant, `P` the gravitational Pontryagin density, and `D=-T^(D)mu_mu` the dark-sector trace.

### 3. Schwarzschild: local threshold vs true instability

The local curvature-eligibility radius is

```text
r_0 = [48 alpha_G M^2/(m_0^2 M_G^2)]^(1/6).
```

But the actual onset is the lowest eigenvalue of

```text
[-d^2/dr_*^2 + V_l(r)] u_l = omega^2 u_l,
```

not merely `m_E^2<0`. `checks.py` evaluates a Rayleigh-Ritz trial function and finds an explicit parameter point with a negative quotient, establishing a genuine tachyonic mode for that test point.

### 4. Kerr "weather" deformation

For slow spin,

```text
r_w(theta)/r_0
 = 1 + [a^2 cos^2(theta)/r_0^2] [-7/2 + 6 eta_Omega] + O(a^4),

eta_Omega = alpha_Omega m_0^2 M_G^4/(alpha_G^2 M_Omega^6).
```

The leading polar deformation changes sign at

```text
eta_Omega = 7/12.
```

Notably, this shape criterion is independent of black-hole mass at this order.

### 5. Dark-halo local eligibility

For an NFW profile,

```text
x_c(1+x_c)^2 = rho_s/rho_c,
rho_c = m_0^2 M_D^2/alpha_D.
```

This predicts only a locally favored phase region; the global field equation must still be solved.

## Reproduce checks

```bash
python3 checks.py
```

## Build manuscript

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Current publication status

Not journal-ready yet. The minimum upgrade path is:

1. compute the Schwarzschild critical spectral curve and the Kerr 2D instability boundary;
2. solve at least one nonlinear defect configuration and its stress/pressure profile;
3. obtain one observable that depends explicitly on `p`.

A short arXiv research note becomes more defensible once the spectral boundary is solved and independently checked.

## Kill criteria

Downgrade or abandon the entropic-weather interpretation if any of these occurs:

- an exact literature collision is found;
- no global unstable mode exists inside the EFT-controlled region;
- the `Pontryagin^2` operator is only relevant beyond its cutoff;
- viable defect/cosmological bounds remove the observable parameter space;
- `p` never enters an accessible observable;
- no independent microscopic entropy interpretation can be derived.

If only the last item fails, the equations may still define a hidden discrete-gauge defect EFT, but the word "entropic" should be removed.
