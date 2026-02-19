# Glassy Dynamics CLASS Implementation Status

**Date:** February 11, 2026
**Phase:** B (Step Function κ(z)) - COMPLETED
**Repository:** lawdroid/class_public (local: week7/GD_CLASS/)

---

## Implementation Summary

### Files Modified

| File | Modification | Status |
|------|--------------|--------|
| `include/background.h` | Added GD struct members (kappa_c, z_freeze, omega_BD, phi_c) | ✅ Done |
| `source/background.c` | Implemented step-function κ(z), modified Friedmann H² | ✅ Done |
| `source/input.c` | Added GD parameter reading from .ini files | ✅ Done |
| `gd_test.ini` | Test configuration with Planck 2018 baseline | ✅ Done |

### GD Parameters

```ini
gd_kappa_c = 1.176    # κ_c = 1/(1-φ_c) = 1/(1-0.15)
gd_z_freeze = 0.8     # Glass transition redshift
gd_omega_BD = 50000   # Brans-Dicke parameter (Cassini constraint)
gd_phi_c = 0.15       # Scher-Zallen percolation threshold
```

---

## Test Results

### Background Comparison (z=0)

| Quantity | ΛCDM | GD+ (κ=1.176) | Ratio | Notes |
|----------|------|---------------|-------|-------|
| H(z=0) [1/Mpc] | 2.247e-04 | 2.193e-04 | 0.976 | ✓ Expected: sqrt(Ω_m/κ + Ω_Λ) |
| Age [Gyr] | 13.797 | 14.083 | 1.021 | Older universe with GD |

**Physical interpretation:**
- H_GD/H_ΛCDM = sqrt((Ω_m/κ + Ω_Λ)/(Ω_m + Ω_Λ)) = sqrt(0.315/1.176 + 0.685) = sqrt(0.953) = 0.976 ✓

### CMB Power Spectrum Comparison

#### Low-ℓ (Quadrupole Region)

| ℓ | ΛCDM C_ℓ | GD+ C_ℓ | Change |
|---|----------|---------|--------|
| 2 | 1.376e-10 | 1.565e-10 | **+13.7%** |
| 3 | 1.303e-10 | 1.441e-10 | +10.6% |
| 4 | 1.234e-10 | 1.338e-10 | +8.5% |
| 5 | 1.182e-10 | 1.264e-10 | +7.0% |

**Issue:** Low-ℓ power INCREASES instead of decreasing. This is opposite to the desired quadrupole suppression.

**Reason:** Step function at z=0.8 creates abrupt ISW contribution that enhances low-ℓ power.

#### High-ℓ (Acoustic Peaks)

| ℓ | ΛCDM C_ℓ | GD+ C_ℓ | Change |
|---|----------|---------|--------|
| 200 | 7.553e-10 | 7.523e-10 | -0.4% |
| 500 | 3.300e-10 | 3.251e-10 | -1.5% |
| 1000 | 1.383e-10 | 1.416e-10 | +2.4% |
| 1500 | 9.507e-11 | 1.007e-10 | **+6.0%** |

**Issue:** High-ℓ deviations up to 6% exceed the 1% target for CMB compatibility.

---

## Physical Interpretation

### What the Step Function Does

```
κ(z) = { 1.0      for z > 0.8   (early universe)
       { 1.176    for z < 0.8   (late universe)
```

At z = 0.8 (about 7 Gyr ago):
- Gravity suddenly weakens: G_eff = G_N/1.176 ≈ 0.85 G_N
- Expansion rate drops: H² ∝ G_eff × ρ
- Gravitational potentials start decaying faster
- Enhanced late-ISW effect → MORE low-ℓ power (not less!)

### Why Quadrupole is NOT Suppressed

The step function creates a **discontinuity** in H(z) that:
1. Creates a sudden change in potential evolution
2. Generates additional ISW contribution at z ≈ 0.8
3. This ISW power adds to the low-ℓ spectrum

**Solution needed:** Smooth Israel-Stewart evolution that transitions gradually.

---

## Next Steps

### Phase C: Smooth Israel-Stewart Evolution

Replace step function with differential equation:
```
τ_κ κ̈ + (1 + 3H τ_κ) κ̇ + Γ(κ - κ_c) = source
```

This requires:
1. Adding κ, κ̇ as dynamical variables in background
2. Integrating alongside other background quantities
3. Choosing appropriate τ_κ and Γ parameters

### Phase D: Perturbations Modifications

For full scalar-tensor consistency:
1. Modify Poisson equation: ∇²Φ = 4πG_eff ρ δ
2. Add scalar field perturbation δκ
3. Include anisotropic stress from κ field

### Parameter Studies

Test different configurations:
- z_freeze = 1100 (recombination) vs z_freeze = 0.8
- Transition width Δz
- κ evolution direction (κ < 1 early → κ > 1 late, or vice versa)

---

## Commands to Reproduce

```bash
cd /Volumes/Backup/Upwork/spectralThermodynamics/week7/GD_CLASS

# Build
make clean && make

# Run ΛCDM baseline
./class lcdm_test.ini

# Run GD+ test
./class gd_test.ini

# Compare outputs
paste output/lcdm_test_00_cl.dat output/gd_test_00_cl.dat | \
  awk 'NR>11 {printf "%d %.2f%%\n", $1, ($9/$2-1)*100}'
```

---

## Key Findings

1. **Implementation works:** κ(z) correctly modifies H(z) as expected from physics
2. **Step function too simple:** Creates artifacts that enhance rather than suppress low-ℓ
3. **Need smooth evolution:** Israel-Stewart formalism is essential for CMB physics
4. **Perturbations matter:** Background-only modification insufficient for accurate C_ℓ

---

*Status report generated: February 11, 2026*
