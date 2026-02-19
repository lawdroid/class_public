GD-CLASS: Glassy Dynamics Cosmic Linear Anisotropy Solving System
=================================================================

**Based on:** T. Martin, *"Glassy Dynamics of Spacetime: A Unified Resolution of the Hubble Tension and Dark Energy via the Topological Percolation Threshold"* (2025), submitted to Classical and Quantum Gravity.

GD-CLASS modifies the [CLASS](https://github.com/lesgourg/class_public) Boltzmann solver to implement spacetime stiffness kappa(z) from Glassy Dynamics theory. By replacing Newton's constant G with an effective G_eff = G/kappa that evolves with redshift, GD-CLASS tests whether Planck-scale topological defects can explain why early-universe and late-universe measurements of the Hubble constant disagree by 5 sigma.

---

## Overview

### The Hubble Tension

| Measurement | H_0 [km/s/Mpc] | Method |
|-------------|----------------|--------|
| Planck CMB  | 67.4 +/- 0.5   | Early universe (z ~ 1100) |
| SH0ES       | 73.0 +/- 1.0   | Local distance ladder (z ~ 0) |
| **Tension** | **~5 sigma**    | Statistically significant disagreement |

### GD Theory

Spacetime has a **stiffness parameter kappa(z)** arising from topological defects at the Planck scale:

```
Standard:  H^2 = (8*pi*G / 3) * rho
GD Theory: H^2 = (8*pi*G / 3*kappa) * rho      where G_eff = G_N / kappa
```

- **kappa > 1**: Weaker effective gravity (early universe)
- **kappa = 1**: Standard gravity recovered (LCDM)
- **kappa_c = 1.176**: Scher-Zallen percolation threshold (phi_c = 0.15)

---

## Stretched Exponential Transition

The kappa(z) transition uses a **stretched exponential** form with three regimes:

```
z >= z_onset (10^6):     kappa = 1        (standard gravity)
z <= z_freeze (1100):    kappa = kappa_c   (fully modified)
z_freeze < z < z_onset:  stretched exponential interpolation

  t = (z - z_freeze) / (z_onset - z_freeze)     normalized 0..1
  decay = exp(-(t / 0.5)^beta)
  kappa(z) = 1 + (kappa_c - 1) * decay
```

### Parameters: z_i, z_f, and beta

- **z_onset (z_i ~ 10^6):** The redshift where kappa begins departing from 1. Setting z_i = 10^6 clears Big Bang Nucleosynthesis (BBN at z ~ 10^9) by three orders of magnitude, so nucleosynthesis is completely standard. Higher z_i risks conflicting with BBN; lower z_i (e.g., 10^4) concentrates the transition near recombination, leaving detectable features in the CMB damping tail.

- **z_freeze (z_f = 1100):** The redshift where kappa freezes at kappa_c. With z_f at recombination, the acoustic peak structure is preserved, but the sound horizon r_s is computed with the modified expansion history throughout the transition, which directly shifts the CMB-inferred H_0. If z_f were lower (e.g., z ~ 100), kappa would still be evolving after recombination, modifying the ISW effect and structure growth — constrained by observations.

- **beta (0.3-0.9):** Controls the transition shape. Small beta (0.3-0.5) means most of the change in kappa happens early, with a long gradual tail — the CMB sees an almost-standard cosmology with a slightly shifted G. Large beta (0.7-0.9) concentrates the change near the middle of the interval, leaving sharper imprints on the damping tail. In practice, beta is currently unconstrained: all values from 0.3 to 0.9 give nearly identical results (see Beta Parameter Scan below).

### Why Stretched Exponential?

In condensed matter, glasses relax not with a single timescale (simple exponential) but with a broad distribution of timescales (stretched exponential). Spacetime's topological defects freeze the same way: as the universe cools through recombination, defect clusters freeze out over a range of redshifts rather than at a single instant.

A tanh or step-function transition creates an abrupt "gravitational jolt" — kappa changes over a narrow redshift range, producing a spike in the Integrated Sachs-Wolfe (ISW) effect at low multipoles. The stretched exponential spreads the transition over z = 10^6 down to z = 1100, with beta controlling how gradual the spread is. This is physically motivated: real glass transitions are never instantaneous.

### Cosmological Consequences

- **Sound horizon increases by 7.8%** -- The frozen acoustic ruler (r_s) is larger because weaker effective gravity (G_eff = G/kappa < G) lets sound waves travel farther before recombination at z = 1100.
- **H_0 = 72.6 km/s/Mpc** -- When the larger ruler is re-matched to Planck's observed angular scale (theta_* ~ 1 degree), the inferred Hubble constant shifts from 67.4 to 72.6, resolving the tension with local SH0ES measurements.
- **ISW excess at l = 2-30 (~10x LCDM)** -- Tom's original step-function test (kappa = 1.08) showed a +20% low-l boost. He predicted the stretched exponential would cure this "gravitational jolt." The smooth transition does help slightly (~2% variation across beta), but at the full kappa_c = 1.176 the ISW excess is +921% — the dominant driver is the magnitude of kappa_c, not the transition shape. This confirms that background-only modifications are insufficient. Phase C (perturbation-level G_eff) is essential to suppress this excess.

---

## Results

### kappa = 1.176 vs LCDM (Planck 2018 baseline)

| Quantity           | LCDM    | GD beta=0.5 | Change  |
|--------------------|---------|-------------|---------|
| Sound horizon r_s  | 144.5 Mpc | 155.8 Mpc | +7.8%  |
| Implied H_0        | 67.4    | 72.6        | +7.7%  |
| Age of universe    | 13.80 Gyr | 14.59 Gyr | +5.7%  |
| Quadrupole (l=2)   | 1022 uK^2 | 10441 uK^2 | +921% |
| First peak (l=220) | 5741 uK^2 | 4530 uK^2  | -21%  |

**Key findings:**
- H_0 = 72.6 km/s/Mpc -- on target for resolving Hubble tension
- ISW spike at low-l is ~10x LCDM -- Phase C (perturbation G_eff) needed to fix
- Beta parameter barely affects ISW: all beta values give similar excess

### Beta Parameter Scan

| beta | r_s [Mpc] | Implied H_0 | D_l(l=2) uK^2 |
|------|-----------|-------------|----------------|
| 0.3  | 154.5     | 72.0        | 10619          |
| 0.5  | 155.8     | 72.6        | 10441          |
| 0.7  | 156.2     | 72.8        | 10411          |
| 0.9  | 156.5     | 72.9        | 10405          |

---

## Interactive Explorer (No Installation Required)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/lawdroid/class_public/blob/feature/kappa-evolution/GD_CLASS_Explorer.ipynb)

**Click the badge above** to open the GD-CLASS Explorer notebook in Google Colab. No compilation, no installation — just click "Runtime → Run all" and move the sliders to explore how GD parameters affect the expansion history, sound horizon, and CMB power spectrum.

The notebook computes kappa(z), H(z), r_s, and H_0 in real time for any parameter combination, and shows pre-computed CMB spectra compared to Planck 2018 data.

---

## Quick Start (for developers)

```bash
# Compile
make clean && make -j4

# Run LCDM baseline
./class lcdm_test.ini

# Run GD with stretched exponential transition (kappa = 1.176)
./class gd_test_KR.ini

# Compare with Planck 2018 data
python3 compare_planck.py
```

The comparison script produces `compare_planck.png` with three panels:
1. Full TT power spectrum (LCDM vs GD vs Planck data)
2. Residuals relative to LCDM
3. Low-l ISW region zoom

---

## GD Parameters

| Parameter      | Default | Description |
|----------------|---------|-------------|
| `gd_kappa_c`   | 1.0     | Stiffness at freeze-out. 1.0 = LCDM, 1.176 = Scher-Zallen |
| `gd_z_freeze`  | 0.8     | Redshift where kappa freezes (use 1100 for recombination) |
| `gd_z_onset`   | 10^6    | Redshift where transition begins (clears BBN) |
| `gd_beta`      | 0.5     | Stretched exponent (0.3-0.9) |
| `gd_omega_BD`  | 50000   | Brans-Dicke parameter (reserved for Phase C) |
| `gd_phi_c`     | 0.15    | Critical defect density (reserved for Phase C) |

Set parameters in `.ini` files. When `gd_kappa_c = 1.0` (default), GD is disabled and CLASS behaves as standard LCDM.

---

## Modified Files

| File | Modification |
|------|--------------|
| `include/background.h` | GD parameter fields (7 fields) and table indices |
| `source/background.c`  | kappa(z) stretched exponential transition, modified Friedmann eq, verbose output |
| `source/input.c`        | Read 6 GD parameters from .ini, set defaults |

Total: ~140 lines of GD code across 3 files. All other CLASS files are unmodified.

---

## Project Status

| Phase | Description | Status |
|-------|-------------|--------|
| **A** | Step-function kappa(z) | Done |
| **B** | Smooth stretched exponential kappa(z) transition | **Done** |
| **C** | Consistent perturbations (G_eff in growth equations) | TODO |
| **D** | Screening mechanisms (Vainshtein/Chameleon) | TODO |
| **E** | MCMC fitting to Planck + BAO + SN data | TODO |

### Phase C Targets (perturbations.c)
- Replace G_N with G_eff = G_N/kappa in Poisson equation
- Add scalar field perturbation delta_kappa
- Implement anisotropic stress from kappa gradient

---

## Documentation

- [README_GD.md](README_GD.md) -- Full implementation details
- [TESTING_RESULTS.md](TESTING_RESULTS.md) -- Simulation data
- [NEXT_STEPS.md](NEXT_STEPS.md) -- Development roadmap
- [GD_IMPLEMENTATION_STATUS.md](GD_IMPLEMENTATION_STATUS.md) -- Milestone tracking

---

## References

1. Martin, T. (2025). "Glassy Dynamics of Spacetime" -- submitted to Classical and Quantum Gravity
2. Blas, D., Lesgourgues, J., & Tram, T. (2011). CLASS II: Approximation schemes. JCAP
3. Planck Collaboration (2018). Cosmological parameters. A&A 641, A6
4. Riess, A. et al. (2022). SH0ES H_0 measurement. ApJL 934, L7

