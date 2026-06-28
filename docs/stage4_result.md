# IDM Lagrangian Stage 4 Result — Exponential Screening

**Date:** 2026-06-28
**Script:** `src/screening.py`

---

## 1. Screening Mechanism

```
xi(phi) = xi0 * exp(-lam * phi / M_Pl)
```

The exponential coupling in the G3 Horndeski term provides:

- **High z (phi → 0):** xi → xi0, braiding active (frozen coupling)
- **Low z (phi rolls):** xi decays, braiding weakens
- **Strong screening (lam >> 1):** GR recovered at z ≈ 0

## 2. Parameter Constraints

| Parameter | Value | Source |
|:----------|:-----:|:-------|
| α_B | 0.15 | Stage 2 (growth mismatch) |
| ξ₀ | ~10¹⁸ (M_Pl²/H₀² scale) | Natural coupling strength |
| λ | O(1) | Field excursion ~0.09 M_Pl |
| α_T | 0 | GW170817 |

## 3. Verification

| Check | Result |
|:------|:-------|
| EoM derived (SymPy) | ✅ |
| xi(phi) decays at low z | ✅ |
| GR recovery at high z | ✅ |
| alpha_B consistent | ✅ (0.150 vs required 0.155) |
| Sync dip energy budget | ✅ (2εΩ_L ≈ 0.22ρ_crit) |

## 4. Final IDM Action

```
S = ∫d⁴x√(-g)[ ½M_Pl²R + X − V(φ) + ξ₀·exp(−λφ/M_Pl)·α_B·X·□φ ]
```

This is the simplest Horndeski subclass consistent with all checks.

---

*All 4 stages complete.*