# IDM Lagrangian Stage 2 Result — α_B(z) from Growth Mismatch

**Date:** 2026-06-28
**Script:** `src/eft_mapping.py`

---

## 1. Method

Compared the IDM growth rate f(z) (with full sync dip H(z)) against ΛCDM
growth rate (same Ω_m, Ω_L). The ratio f_IDM/f_ΛCDM encodes the required
braiding α_B(z) via the Horndeski quasi-static expression:

```
G_eff/G_N = [1 + α_B/(1+α_B)] / [1 + α_K/(6(1+α_B)²)]
```

Assumptions: α_T=0 (GW170817), α_M=0, α_K=1.

---

## 2. Results

| z | f_IDM/f_ΛCDM | G_eff/G_N | α_B |
|:-:|:------------:|:---------:|:---:|
| 0.0 | 0.9965 | 0.993 | +0.14 |
| 0.6 | 1.0000 | 1.000 | +0.15 |
| 1.0 | 1.0061 | 1.012 | +0.16 |

α_B ≈ +0.14–0.16 across all z — **small but non-zero braiding**.

---

## 3. Combined Stages 1+2 Classification

| Criterion | Result |
|:----------|:-------|
| Can IDM be GR + canonical scalar? | ❌ No (phantom crossing at z<0.59) |
| α_T = 0? | ✅ Required (GW170817) |
| α_B ≠ 0? | ✅ Yes, +0.14–0.16 |
| α_K = O(1)? | ✅ Required (not canonical) |
| α_M = 0? | ⚠️ Assumed — relaxable |

**IDM candidate subclass:** Horndeski with α_B ≈ +0.15, α_T = 0,
finite α_K (O(1)).

---

## 4. Next Step

Stage 3: SymPy EoM — verify a Horndeski action with α_B(z) ≈ +0.15
produces the IDM field equations.