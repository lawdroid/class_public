# GD_CLASS: Next Steps

**Date:** February 12, 2026

---

## Current Status

| Component | Status |
|-----------|--------|
| κ(z) in background | ✅ Implemented |
| Step function transition | ✅ Tested |
| H₀ mechanism | ✅ Confirmed working |
| CMB consistency | ⚠️ Needs work |
| Perturbation consistency | 🔲 Not yet implemented |

---

## Recommended Path Forward

### Phase 1: Parameter Exploration (1 week)

**Goal:** Map the full κ parameter space

Tasks:
- [ ] Test κ = 1.10, 1.12, 1.15, 1.18, 1.20
- [ ] Record r_s, H₀(implied), CMB changes for each
- [ ] Find κ needed for H₀ = 73 km/s/Mpc
- [ ] Document scaling relationships

**Deliverable:** Parameter table showing κ vs observables

---

### Phase 2: CMB Parameter Fitting (2 weeks)

**Goal:** Find if CMB can be restored with adjusted parameters

Tasks:
- [ ] Vary (ω_b, ω_cdm, n_s, A_s) along with κ
- [ ] Check if CMB peaks can be restored
- [ ] Calculate χ² against Planck data
- [ ] Identify parameter degeneracies

**Approach:**
```
For each κ value:
1. Run CLASS with baseline parameters
2. Compare CMB to Planck
3. Adjust parameters to minimize residuals
4. Record best-fit and χ²
```

**Deliverable:** Best-fit parameters for GD model, χ² comparison

---

### Phase 3: Consistent Perturbations (2-3 weeks)

**Goal:** Properly implement κ in perturbation equations

Current issue:
```
Background: H² = (8πG/3κ) × ρ    ✓ implemented
Perturbations: δ̈ + 2Hδ̇ = 4πGρδ  ← still uses G, not G/κ
```

Tasks:
- [ ] Modify `source/perturbations.c` for consistent G_eff
- [ ] Check energy-momentum conservation
- [ ] Verify causality preserved
- [ ] Test with Israel-Stewart formalism

**Key files:**
- `source/perturbations.c` — growth equations
- `source/thermodynamics.c` — recombination

---

### Phase 4: Screening Mechanism (2-3 weeks)

**Goal:** Implement scale-dependent κ

Physics:
```
Large scales (Hubble): κ ≠ 1 → modified expansion
Small scales (CMB):    κ → 1 → standard perturbations
```

Options:
1. **Vainshtein-type:** κ depends on curvature
2. **Chameleon-type:** κ depends on density
3. **k-dependent:** κ(k) varies with wavenumber

---

### Phase 5: Observational Comparison (ongoing)

**Goal:** Test against real data

Datasets:
- [ ] Planck 2018 CMB
- [ ] BAO measurements
- [ ] SH0ES H₀
- [ ] DES/KiDS weak lensing (σ₈)
- [ ] Pantheon+ supernovae

---

## Technical Improvements

### Code Quality
- [ ] Add unit tests for κ functions
- [ ] Validate against known limits (κ→1 = ΛCDM)
- [ ] Document all modified functions

### Python Interface
- [ ] Expose GD parameters to `classy`
- [ ] Create analysis notebooks
- [ ] Automate comparison plots

### Performance
- [ ] Profile code for bottlenecks
- [ ] Optimize if needed for MCMC runs

---

## Timeline Estimate

| Phase | Duration | Cumulative |
|-------|----------|------------|
| Parameter exploration | 1 week | Week 1 |
| CMB fitting | 2 weeks | Week 3 |
| Perturbations | 2-3 weeks | Week 6 |
| Screening | 2-3 weeks | Week 9 |
| Observational tests | Ongoing | Week 10+ |

---

## Decision Points

### After Phase 1:
> "Can any κ value give H₀ = 73 with reasonable CMB changes?"
> - If YES → proceed to Phase 2
> - If NO → reconsider theory

### After Phase 2:
> "Can parameter adjustment restore CMB fit?"
> - If YES → viable with refitting
> - If NO → need Phase 3/4 modifications

### After Phase 3:
> "Does consistent perturbation treatment help?"
> - If YES → continue development
> - If NO → screening mechanism essential

---

## Resources Needed

| Resource | Purpose |
|----------|---------|
| Computing | MCMC parameter scans |
| Planck data | CMB likelihood |
| Literature review | Screening mechanisms |
| Collaboration | Cross-check results |

---

## Success Criteria

The GD model is viable if it can:

1. ✅ Shift H₀ from 67 → 73 km/s/Mpc
2. 🔲 Preserve CMB peaks within Planck errors
3. 🔲 Keep low-ℓ within cosmic variance (~10%)
4. 🔲 Satisfy BBN constraints
5. 🔲 Pass solar system tests (Cassini)
6. 🔲 Provide testable predictions

---

---

*Roadmap prepared February 12, 2026*
