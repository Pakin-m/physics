# Finite-State Entropic Defects — compact research package

## Core claim under test

The project no longer treats gravity as literal atmospheric flow.  Its narrow candidate novelty is a hidden, environmentally activated gauge phase:

- a charge-`p` hidden Higgs field breaks `U(1) -> Z_p` when a covariant effective mass becomes negative;
- activation can be driven jointly by strong curvature and by the dark-sector stress-tensor trace;
- the residual finite arithmetic is physical because a unit hidden charge encircling a defect acquires the exact phase `exp(2π i n/p)`;
- bulk “entropic pressure” is computed from the hidden-sector stress tensor, not postulated.

The manuscript explicitly **does not** claim that these ingredients individually are new.  Density-dependent gauge strings, environment-dependent defects, curvature-triggered scalar phases, discrete gauge theories, and quantized dark-matter vortices all have prior literature.  The novelty hypothesis is the specific high-density/curvature activation plus residual `Z_p` topological observable in one covariant EFT.

## Files

- `main.tex` — main manuscript; includes derivations, collision audit, falsification conditions, and next calculations.
- `references.bib` — compact bibliography, prioritizing peer-reviewed sources.
- `checks.py` — dependency-free checks of the NFW threshold, Schwarzschild activation radius, horizon condition, and `Z_p` periodicity.

## Reproduce the algebra checks

```bash
python3 checks.py
```

## Build the manuscript

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Current strongest results

1. **Exact finite-arithmetic observable**
   \[
   \Phi_E(n)=\frac{2\pi n}{p g_E},\qquad
   e^{i g_E\Phi_E}=e^{2\pi i n/p}.
   \]

2. **Covariant phase threshold**
   \[
   m_E^2=m_0^2-\frac{\alpha_G}{M_G^2}\mathcal G
              -\frac{\alpha_D}{M_D^2}\mathcal D<0.
   \]

3. **Schwarzschild activation radius**
   \[
   r_G=\left(\frac{48\alpha_G M^2}{m_0^2M_G^2}\right)^{1/6}.
   \]

4. **NFW halo activation condition**
   \[
   x_c(1+x_c)^2=\rho_s/\rho_c,
   \quad \rho_c=m_0^2M_D^2/\alpha_D.
   \]

## Kill criteria

The project should be downgraded or abandoned if: an exact literature collision is found; no controlled EFT parameter window survives; defect/cosmological bounds eliminate the observable regime; the `Z_p` phase cannot affect any measurable dark-sector or gravitational observable; or no independent entropy derivation exists.  If the last point fails, the same equations may still describe a hidden discrete-gauge defect model, but the word “entropic” should be removed.
