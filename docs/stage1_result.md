# IDM Lagrangian Stage 1 Result — V(φ) Reconstruction

**Date:** 2026-06-28
**Project:** `/home/wsl/idm-lagrangian/`
**Script:** `notebooks/v_of_phi.py`

---

## 1. What Was Done

Reconstructed the quintessence potential V(φ) from IDM's sync-dip modified H(z):

```
H²(z)/H₀² = Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_L·f(z)
f(z) = 1 − ε·(z/zc)·exp(−z/zc)

Parameters: Ω_m=0.297, ε=0.1545, zc=0.60
```

using the standard quintessence reconstruction formulas:

```
φ'(z) = √[ dE²/dz − 3Ω_m(1+z)² ] / √[(1+z)E²]   (M_Pl units)
V(z)  = H₀²[ 3E² − 1.5Ω_m(1+z)³ − 0.5(1+z)dE²/dz ]
```

---

## 2. Result: Phantom Crossing

The reconstruction reveals that IDM **cannot be canonical quintessence**.

| z range | φ'(z)² | Status |
|:--------|:------:|:-------|
| 0–0.59 | **Negative** | ❌ Phantom crossing — canonical field cannot produce this H(z) |
| 0.59–3.0 | Positive | ✅ Canonical reconstruction valid |

**Why:** The sync dip at zc=0.6 causes `dE²/dz < 3Ω_m(1+z)²` at low z, violating the reality condition `φ'² ≥ 0`. In EFT language, this is equivalent to requiring **braiding α_B(z) ≠ 0** or **non-minimal coupling ξ(φ)R**.

---

## 3. Valid-Region V(φ) Properties

In z = 0.59–3.0 where reconstruction is valid:

| Property | Value |
|:---------|:------|
| φ excursion | 0 → 0.093 M_Pl |
| V range | 2.0 → 2.1 [H₀²·M_Pl² nat. units] |
| V range (physical) | 5.9–6.3 × 10⁻⁴³ GeV⁴ |
| V > 0 | ✅ |
| V'' > 0 | ✅ |
| c_s² ≥ 0 | ✅ |

The field rolls a very short distance (~0.09 M_Pl) with a nearly flat potential — consistent with slow-roll quintessence in the valid region.

---

## 4. EFT Interpretation

The phantom crossing at z ≈ 0.59 gives a strong prior on the IDM Lagrangian structure:

**Rule out:** GR + canonical scalar (quintessence)
**Require (at least one):**
- **α_B ≠ 0** (braiding: kinetic mixing between scalar & metric in Horndeski)
- **ξ(φ)R** (non-minimal coupling à la Brans-Dicke)
- **P(φ,X)** (k-essence with non-canonical kinetic term)

The simplest path forward: treat the phantom crossing as direct evidence for **α_B(z)**, and compute its required trajectory from the fσ₈ growth mismatch.

---

## 5. Files Produced

```
notebooks/results/
├── vofphi.txt          # φ, V(φ), V(GeV⁴), V(eV⁴)
├── phi_z.txt           # z, φ(z), V(z)
├── reality_check.txt   # z, φ'², E², dE²/dz
└── figs/
    └── v_of_phi.png    # 4-panel diagnostic figure
```

---

*Next: Stage 2 — α_B(z) from fσ₈ growth mismatch*
