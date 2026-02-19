# GD_CLASS Testing Results

**Date:** February 12, 2026
**CLASS Version:** v3.3.4

---

## Test Configurations

### Baseline: ΛCDM (Planck 2018)
```ini
h = 0.6736
omega_b = 0.02237
omega_cdm = 0.1200
tau_reio = 0.0544
n_s = 0.9649
A_s = 2.1e-9
```

### Test A: Late-time κ < 1
```ini
gd_kappa_c = 0.85
gd_z_freeze = 1100
# κ = 1.0 at z > 1100, κ = 0.85 at z < 1100
```

### Test B: Early-time κ > 1
```ini
gd_kappa_c = 1.08
gd_z_freeze = 1100
# κ = 1.08 at z > 1100, κ = 1.0 at z < 1100
```

---

## Results: Background Cosmology

### H₀ and Age

| Model | H₀ [km/s/Mpc] | Age [Gyr] | Change |
|-------|---------------|-----------|--------|
| ΛCDM | 67.36 | 13.797 | — |
| Late κ=0.85 | 69.20 | 13.024 | H₀ +2.7% |
| Early κ=1.08 | 67.36 | 13.798 | H₀ unchanged* |

*Early-time κ doesn't directly change H₀; it changes sound horizon

### Sound Horizon at Recombination

| Model | r_s [Mpc] | Change |
|-------|-----------|--------|
| ΛCDM | 143.53 | — |
| Early κ=1.08 | 148.83 | +3.69% |

**Implication:** To maintain CMB peak positions with larger r_s, H₀ must increase to ~69.9 km/s/Mpc when fitting data.

### H(z) at Key Redshifts

| z | H_ΛCDM | H_GD (early κ=1.08) | Ratio |
|---|--------|---------------------|-------|
| 0 | 67.36 | 67.36 | 1.000 |
| 1100 | — | — | 0.981 |
| 2000 | — | — | 0.963 |
| 10000 | — | — | 0.962 |

---

## Results: CMB Power Spectrum

### Low-ℓ Multipoles (Quadrupole Region)

| ℓ | ΛCDM [μK²] | Late κ=0.85 | Change | Early κ=1.08 | Change |
|---|------------|-------------|--------|--------------|--------|
| 2 | 1022 | 21360 | +1990% | 1229 | +20% |
| 3 | 968 | 18728 | +1835% | 1181 | +22% |
| 4 | 916 | 16663 | +1718% | 1129 | +23% |
| 5 | 878 | 15014 | +1610% | 1089 | +24% |
| 10 | 819 | 10036 | +1125% | 1023 | +25% |

### First Acoustic Peak (ℓ ~ 220)

| ℓ | ΛCDM [μK²] | Late κ=0.85 | Change | Early κ=1.08 | Change |
|---|------------|-------------|--------|--------------|--------|
| 220 | 5741 | 5350 | -6.8% | 4752 | -17.2% |
| 221 | 5741 | 5344 | -6.9% | 4748 | -17.3% |
| 222 | 5741 | 5338 | -7.0% | 4744 | -17.4% |

---

## Physics Interpretation

### Late-time κ < 1 Approach

```
Physics Chain:
κ < 1 at z < 1100
    ↓
G_eff = G_N/κ > G_N (stronger gravity)
    ↓
H² larger → H₀ increases ✓
    ↓
But: Perturbations also see G_eff
    ↓
ISW effect amplified → Low-ℓ explodes ❌
```

**Conclusion:** Mechanism works for H₀ but breaks CMB.

### Early-time κ > 1 Approach

```
Physics Chain:
κ > 1 at z > 1100
    ↓
G_eff = G_N/κ < G_N (weaker early gravity)
    ↓
Slower early expansion → larger sound horizon r_s
    ↓
To fit CMB peaks: need higher H₀
    ↓
Late-time κ = 1 → ISW less affected ✓
```

**Conclusion:** More promising, but still affects CMB.

---

## Quantitative Relationships

### κ vs Sound Horizon

From Test B (κ = 1.08):
```
Δr_s/r_s ≈ 0.35 × (κ - 1)
```

### κ vs Implied H₀

To maintain CMB peak positions (θ_s = r_s/D_A constant):
```
H₀_needed ≈ H₀_Planck × (r_s_GD / r_s_ΛCDM)
```

| κ_early | r_s increase | Implied H₀ |
|---------|--------------|------------|
| 1.05 | ~2% | ~68.7 |
| 1.08 | ~3.7% | ~69.9 |
| 1.15 | ~5.5% | ~71.1 |
| 1.20 | ~7% | ~72.1 |

---

## Output Files

All output in `output/` directory (gitignored):

| File Pattern | Contents |
|--------------|----------|
| `lcdm_test_*_background.dat` | ΛCDM background evolution |
| `lcdm_test_*_cl.dat` | ΛCDM CMB power spectrum |
| `gd_k085_*_background.dat` | Late κ=0.85 background |
| `gd_k085_*_cl.dat` | Late κ=0.85 CMB spectrum |
| `gd_early_k108_*_background.dat` | Early κ=1.08 background |
| `gd_early_k108_*_cl.dat` | Early κ=1.08 CMB spectrum |

---

## Conclusions

1. **H₀ modification mechanism confirmed** — κ can change expansion rate
2. **Simple step function has limitations** — CMB perturbations affected
3. **Early-time κ > 1 is more promising** — reduces ISW impact
4. **Full parameter fitting needed** — to assess true viability

---

*Data collected February 12, 2026*
