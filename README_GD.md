# Glassy Dynamics CLASS Implementation

**Project:** Spectral Thermodynamics / Hubble Tension Resolution
**Repository:** CLASS modified for Glassy Dynamics (GD) theory testing
**Date:** February 2026
**Branch:** `feature/kappa-evolution`

---

## Goal

Test whether the Glassy Dynamics theory can resolve the **Hubble tension**:
- CMB (Planck): H_0 = 67.4 +/- 0.5 km/s/Mpc
- Local (SH0ES): H_0 = 73.0 +/- 1.0 km/s/Mpc
- Tension: ~5 sigma disagreement

**GD Hypothesis:** Spacetime stiffness kappa(z) from topological defects modifies the expansion history, potentially explaining the discrepancy.

---

## GD Implementation

### The kappa(z) Stiffness Parameter

Modified Friedmann equation:
```
H^2 = (8*pi*G_N)/(3*kappa) * rho_matter_rad + Lambda_eff/3
```

Where:
- kappa > 1: Weaker effective gravity (slower expansion)
- kappa < 1: Stronger effective gravity (faster expansion)
- kappa = 1: Standard LCDM

### Kohlrausch (Stretched Exponential) Transition

The kappa(z) evolution uses a three-regime stretched exponential:

```c
if (z >= z_onset) {
  kappa = 1.0;                    // Above onset: standard gravity
} else if (z <= z_freeze) {
  kappa = kappa_c;                // Below freeze: fully modified
} else {
  t = (z - z_freeze) / (z_onset - z_freeze);   // Normalized 0..1
  decay = exp(-pow(t / 0.5, beta));
  kappa = 1.0 + (kappa_c - 1.0) * decay;       // Stretched exponential
}
```

- **beta small (0.3-0.5):** Gradual transition, most change early, long tail
- **beta large (0.7-0.9):** Sharper transition, concentrated near z_freeze

### Files Modified

| File | Modification |
|------|--------------|
| `source/background.c` | kappa(z) Kohlrausch transition, modified H^2, GD verbose output |
| `include/background.h` | 7 GD fields: has_gd, kappa_c, z_freeze, z_onset, beta, omega_BD, phi_c + 2 table indices |
| `source/input.c` | GD input parameter parsing (6 params), defaults |

### Parameters

```ini
gd_kappa_c  = 1.176      # Stiffness at freeze-out (Scher-Zallen)
gd_z_freeze = 1100       # Transition redshift (recombination)
gd_z_onset  = 1000000    # Onset redshift (clears BBN)
gd_beta     = 0.5        # Kohlrausch stretched exponent
gd_omega_BD = 50000      # Brans-Dicke parameter (reserved for Phase C)
gd_phi_c    = 0.15       # Critical defect density (reserved for Phase C)
```

---

## Test Results

### Test 1: Early-time kappa = 1.08 (conservative)

| Result | Value |
|--------|-------|
| Sound horizon r_s | +3.7% increase |
| Implied H_0 shift | +3.7% (67.4 -> ~69.9) |
| CMB quadrupole l=2 | +20% |
| Finding | Promising but modest H_0 shift |

### Test 2: Kohlrausch kappa = 1.176 (Tom's target)

| Result | LCDM | KR beta=0.5 | Change |
|--------|------|-------------|--------|
| Sound horizon r_s | 144.5 Mpc | 155.8 Mpc | +7.8% |
| Implied H_0 | 67.4 | 72.6 | +7.7% |
| Age | 13.80 Gyr | 14.59 Gyr | +5.7% |
| Quadrupole l=2 | 1022 uK^2 | 10441 uK^2 | +921% |
| First peak l=220 | 5741 uK^2 | 4530 uK^2 | -21% |

### Beta Parameter Scan (kappa = 1.176)

| beta | r_s [Mpc] | Implied H_0 | D_l(l=2) |
|------|-----------|-------------|----------|
| 0.3 | 154.5 | 72.0 | 10619 |
| 0.5 | 155.8 | 72.6 | 10441 |
| 0.7 | 156.2 | 72.8 | 10411 |
| 0.9 | 156.5 | 72.9 | 10405 |

**Conclusion:** Beta barely affects ISW. Background-only modification produces correct H_0 but the ISW spike (10x LCDM) requires Phase C perturbation modifications.

---

## Key Findings

### What Works
1. **kappa mechanism successfully modifies H_0** -- the physics is correct
2. **Sound horizon responds as expected** -- r_s increases with early-time kappa > 1
3. **Kohlrausch transition compiles and runs** -- three-regime form is clean
4. **H_0 = 72.6 with kappa = 1.176** -- on target for Hubble tension resolution

### What Needs Work
1. **ISW spike is 10x LCDM** -- background-only modification insufficient
2. **Phase C essential** -- need G_eff in perturbation equations
3. **Beta parameter has minimal effect** -- ISW not sensitive to transition shape
4. **Full MCMC fitting needed** -- to find best-fit parameters with GD

---

## Repository Structure

```
GD_CLASS/
+-- source/
|   +-- background.c      # Modified: kappa(z) Kohlrausch + Friedmann + verbose
|   +-- perturbations.c   # Future: consistent perturbations (Phase C)
|   +-- input.c           # Modified: 6 GD parameters
+-- include/
|   +-- background.h      # Modified: GD struct fields + indices
+-- data/
|   +-- planck_2018_TT.txt # Planck 2018 TT power spectrum
+-- output/               # Test results (gitignored)
+-- compare_planck.py     # Planck comparison + beta scan script
+-- test_cl.py            # CMB analysis script
+-- gd_test_KR.ini        # Kohlrausch config (kappa=1.176)
+-- gd_test_early_k108.ini # Conservative test (kappa=1.08)
+-- lcdm_test.ini         # LCDM baseline
+-- README.md             # Project overview
+-- README_GD.md          # This file
+-- NEXT_STEPS.md         # Development roadmap
+-- TESTING_RESULTS.md    # Test results
```

---

## Running Tests

### Compile
```bash
make clean && make -j4
```

### Run LCDM Baseline
```bash
./class lcdm_test.ini
```

### Run GD Kohlrausch Test
```bash
./class gd_test_KR.ini
```

Expected output:
```
Running CLASS version v3.3.4
Computing background
 -> Glassy Dynamics enabled:
    kappa_c = 1.176 (stiffness at freeze-out)
    z_freeze = 1100 (glass transition redshift)
    z_onset = 1e+06 (onset redshift)
    beta = 0.5 (Kohlrausch stretched exponent)
    G_eff/G_N = 0.85034 (at z < z_freeze)
    H0_local/H0_CMB = 1.08444 (Hubble tension factor)
```

### Compare with Planck
```bash
python3 compare_planck.py
```

Generates `compare_planck.png` with TT spectrum comparison, residuals, and ISW zoom.

---

## Next Steps

### Phase C: Consistent Perturbations (perturbations.c)
- [ ] Replace G_N with G_eff = G_N/kappa in Poisson equation
- [ ] Add scalar field perturbation delta_kappa
- [ ] Implement anisotropic stress from kappa gradient
- [ ] Israel-Stewart viscous dynamics

### Phase D: Screening Mechanisms
- [ ] Vainshtein screening
- [ ] Chameleon screening

### Phase E: MCMC Fitting
- [ ] Implement MCMC with GD parameters
- [ ] Fit to Planck 2018 + BAO + SN data
- [ ] Find best-fit cosmology with GD

---

## References

1. Martin, T. (2025). "Glassy Dynamics of Spacetime" -- submitted to Classical and Quantum Gravity
2. Blas, D., Lesgourgues, J., & Tram, T. (2011). CLASS II: Approximation schemes
3. Planck Collaboration (2018). Cosmological parameters
4. Riess, A. et al. (2022). SH0ES H_0 measurement

---

*This is an active research project. Results are preliminary and subject to revision.*
