GD-CLASS: Glassy Dynamics Cosmic Linear Anisotropy Solving System
=================================================================

**Based on:** T. Martin, *"Glassy Dynamics of Spacetime: A Unified Resolution of the Hubble Tension and Dark Energy via the Topological Percolation Threshold"* (2025), submitted to Classical and Quantum Gravity.

GD-CLASS modifies the [CLASS](https://github.com/lesgourg/class_public) Boltzmann solver to implement spacetime stiffness kappa(z) from Glassy Dynamics theory. By replacing Newton's constant G with an effective G_eff = G/kappa that evolves with redshift, GD-CLASS tests whether Planck-scale topological defects can explain the 5-sigma Hubble tension.

---

## Key Result

GD background-only modification (justified by omega_BD = 50,000) can partially resolve the Hubble tension while maintaining an acceptable fit to Planck 2018 CMB data:

| Model | kappa_c | H_0 [km/s/Mpc] | chi2/dof vs Planck | Quality |
|-------|---------|-----------------|-------------------|---------|
| LCDM  | 1.00    | 67.4            | 1.17              | Baseline |
| **GD**| **0.98**| **71.0**        | **1.26**          | **+8% (publishable)** |
| GD    | 0.96    | 73.0            | 1.54              | +32% (marginal) |

The tradeoff is monotonic: more H_0 shift = worse CMB fit. No magic kappa resolves the tension without cost. This is comparable to Early Dark Energy (EDE) models in the literature.

---

## Overview

### The Hubble Tension

| Measurement | H_0 [km/s/Mpc] | Method |
|-------------|----------------|--------|
| Planck CMB  | 67.4 +/- 0.5   | Early universe (z ~ 1100) |
| SH0ES       | 73.0 +/- 1.0   | Local distance ladder (z ~ 0) |
| **Tension** | **~5 sigma**    | Statistically significant disagreement |

### GD Theory (Compliant Inclusion Model)

Spacetime contains compliant topological inclusions. At the Scher-Zallen percolation threshold:

```
Early universe (z > 1100): kappa = 1 - phi_c  (stronger gravity, G_eff > G_N)
Late universe  (z < 1100): kappa = 1.0        (standard gravity)
```

- **kappa < 1**: Stronger effective gravity (pre-recombination)
- **kappa = 1**: Standard gravity recovered (LCDM)
- Transition triggered by **recombination quench** at z ~ 1090

The modified Friedmann equation:
```
H^2 = (8*pi*G / 3*kappa) * rho      where G_eff = G_N / kappa
```

### Why Background-Only Is Justified

For Brans-Dicke parameter omega_BD = 50,000 (required by Cassini: gamma_PPN = 0.99998):
- Perturbation corrections are O(1/omega_BD) ~ 0.002%
- The scalar field perturbation delta_kappa is suppressed
- CMB peaks are preserved to high accuracy
- Only the background expansion (H_0, r_s) differs measurably

This was confirmed by Phase C testing: naive G_eff injection into perturbation equations (without delta_kappa compensation) makes the CMB catastrophically worse, as expected for an incomplete scalar-tensor implementation.

---

## Recombination Quench Transition

The kappa(z) transition uses a **stretched exponential** centered at recombination:

```
z >= z_onset (1140):    kappa = kappa_c   (modified gravity)
z <= z_freeze (1040):   kappa = 1.0       (standard gravity)
Transition width:       Dz ~ 100 (matches recombination timescale)

  t = (z - z_freeze) / (z_onset - z_freeze)     normalized 0..1
  decay = exp(-(t / 0.5)^beta)
  kappa(z) = kappa_c + (1.0 - kappa_c) * decay
```

The transition width Dz ~ 100 matches the physical timescale of recombination (~100,000 years).

---

## Results

### Parameter Fitting (Background-Only, Phase C Disabled)

With re-fitted cosmological parameters for each kappa value:

| kappa_c | H_0 | chi2/dof | r_s [Mpc] | Best-fit key params |
|---------|-----|----------|-----------|---------------------|
| 1.00    | 67.4 | 1.17    | 144.52    | Planck LCDM best-fit |
| 0.99    | 69.0 | 1.21    | -         | h=0.69, n_s=0.97 |
| 0.98    | 71.0 | 1.26    | 144.28    | h=0.71, omega_b=0.021, n_s=0.99 |
| 0.97    | 71.0 | 1.36    | -         | h=0.71, n_s=1.00 |
| 0.96    | 73.0 | 1.54    | 141.59    | h=0.73, omega_cdm=0.125, n_s=1.00 |
| 0.95    | 73.0 | 1.80    | -         | h=0.73, n_s=1.00 |

**Caveat:** Best fits prefer n_s ~ 1.0 (scale-invariant), in tension with Planck's n_s = 0.965.

### Phase C Test (G_eff in Perturbations)

Phase C was implemented and tested but makes the CMB worse:

| Model | D_l(l=2) | D_l(l=220) | H_0 | Verdict |
|-------|----------|------------|-----|---------|
| LCDM | 1022 | 5741 | 67.4 | baseline |
| Quench bg-only | 389 | 10346 | 73.0 | H_0 right, peaks doubled |
| Quench + Phase C | 1777 | 12669 | 73.0 | ISW explodes, peaks worse |

Phase C is disabled in the current code (G_eff_ratio = 1.0 in perturbations.c).

### Comparison Plot

![GD vs Planck](gd_bestfit_vs_planck.png)

Three panels: (1) TT spectrum comparison, (2) residuals vs Planck, (3) chi2/dof vs H_0 tradeoff.

---

## Interactive Explorer (No Installation Required)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/lawdroid/class_public/blob/feature/kappa-evolution/GD_CLASS_Explorer.ipynb)

**Click the badge above** to open the GD-CLASS Explorer notebook in Google Colab. No compilation, no installation -- just click "Runtime -> Run all" and move the sliders.

The notebook shows pre-computed CMB spectra for LCDM, GD kappa=0.98, and GD kappa=0.96 compared to Planck 2018 data. Interactive sliders compute kappa(z), H(z), r_s, and H_0 in real time.

---

## Quick Start (for developers)

```bash
# Compile
make clean && make -j4

# Run LCDM baseline
./class lcdm_test.ini

# Run GD best-fit (kappa=0.98, H0=71)
./class gd_bestfit_k098.ini

# Run GD (kappa=0.96, H0=73)
./class gd_bestfit_k096.ini

# Parameter scan (grid search over cosmological params)
python3 gd_param_scan.py
```

---

## GD Parameters

| Parameter      | Default | Description |
|----------------|---------|-------------|
| `gd_kappa_c`   | 1.0     | Stiffness at freeze-out. 1.0 = LCDM |
| `gd_z_freeze`  | 1040    | Redshift where kappa reaches 1.0 |
| `gd_z_onset`   | 1140    | Redshift where transition begins |
| `gd_beta`      | 0.5     | Stretched exponent (0.3-0.9) |

Set parameters in `.ini` files. When `gd_kappa_c = 1.0` (default), GD is disabled and CLASS behaves as standard LCDM.

---

## Modified Files

| File | Modification |
|------|--------------|
| `include/background.h` | GD parameter fields and table indices |
| `source/background.c`  | kappa(z) transition, modified Friedmann eq |
| `source/input.c`        | Read GD parameters from .ini, set defaults |
| `source/perturbations.c` | Phase C hooks (currently disabled) |

---

## Project Status

| Phase | Description | Status |
|-------|-------------|--------|
| **A** | Step-function kappa(z) | Done |
| **B** | Smooth stretched exponential transition | Done |
| **C** | G_eff in perturbation equations | Done (disabled -- makes CMB worse) |
| **D** | Parameter fitting to Planck TT | **Done (grid search)** |
| **E** | MCMC fitting with proper contours | TODO |

---

## References

1. Martin, T. (2025). "Glassy Dynamics of Spacetime" -- submitted to Classical and Quantum Gravity
2. Blas, D., Lesgourgues, J., & Tram, T. (2011). CLASS II: Approximation schemes. JCAP
3. Planck Collaboration (2018). Cosmological parameters. A&A 641, A6
4. Riess, A. et al. (2022). SH0ES H_0 measurement. ApJL 934, L7
