# Milestone 01: Late-time κ < 1 Test Results

**Date:** February 12, 2026
**Branch:** feature/kappa-evolution

---

## What We Tested

Step-function κ(z) with late-time modification:
- z > 1100: κ = 1.0 (standard)
- z < 1100: κ = 0.85 (modified)

Physics: G_eff = G_N/κ = 1.176 G_N at late times

---

## Results

| Parameter | ΛCDM | GD (κ=0.85) | Change |
|-----------|------|-------------|--------|
| H₀ [km/s/Mpc] | 67.36 | 69.20 | +2.7% |
| Age [Gyr] | 13.80 | 13.02 | -5.6% |
| Quadrupole ℓ=2 [μK²] | 1022 | 21360 | +2000% |
| First peak ℓ=220 [μK²] | 5741 | 5350 | -7% |

---

## Conclusion

**H₀ mechanism works** — increased from 67.4 toward 73.

**CMB badly broken** — low-ℓ exploded due to enhanced ISW effect.

**Root cause:** G_eff modification affects both:
1. Background expansion (intended) ✓
2. Perturbation growth / ISW (unintended) ❌

---

## Next Step

Try **early-time κ > 1** approach:
- z > 1100: κ > 1 (weaker early gravity)
- z < 1100: κ = 1 (standard late-time)

Rationale: CMB forms with modified physics but late-time ISW protected.

---

## Files

- `source/background.c` — κ implementation (fixed for κ < 1)
- `gd_test_k085.ini` — test configuration (in .gitignore)
- `test_cl.py` — CMB power spectrum analysis script
- `output/gd_k085_*` — raw output data (in .gitignore)

---

*Milestone recorded for future reference*
