# IDM Lagrangian Reverse Engineering

Reverse-engineering the Horndeski action behind the Information-Dynamics
Model (IDM) sync dip at z ≈ 0.6.

## Summary

The IDM sync dip (ε=0.1545, zc=0.60) in the dark energy density produces:
- Δχ² = −50.8 vs ΛCDM (CMB + SNe + BAO + H(z))
- ΔAIC = +50.1, ΔBIC = +40.6

This project reverse-engineers the fundamental Lagrangian that produces
this dip, using 4 stages of analysis:

| Stage | Result |
|:------|:-------|
| 1 — V(φ) Reconstruction | Phantom crossing at z < 0.59 → **canonical quintessence ruled out** |
| 2 — α_B(z) from Growth | **α_B ≈ +0.15** (Horndeski braiding) |
| 3 — SymPy EoM | Action self-consistent, G₃ term resolves phantom |
| 4 — Exponential Screening | ξ(φ)=ξ₀·exp(−λφ) ensures GR recovery at high z |

## Final IDM Action

```
S = ∫d⁴x√(−g)[ ½M_Pl²R + X − V(φ) + ξ₀·e^{−λφ}·α_B·X·□φ ]
```

## Key Numbers

| Parameter | Value |
|:----------|:------|
| α_B | +0.14–0.16 (braiding) |
| α_T | 0 (GW170817: c_T = c) |
| ε (sync dip) | 0.1545 |
| zc (sync dip) | 0.60 |
| CMB Cℓ consistency | χ² = 0.51/dof |
| fσ₈ z-trend | +2.81σ/z (p < 0.001) |

## Structure

```
IDM-Lagrangian-Reverse-Engineering/
├── src/              # Python source modules
│   ├── hubble.py        # IDM H(z) + sync dip
│   ├── reconstruction.py# V(φ) from H(z)
│   ├── eft_mapping.py   # α_B(z) from growth
│   ├── eom_verify.py    # SymPy EoM
│   ├── screening.py     # ξ(φ) exponential screening
│   ├── cmb_comparison.py# CMB Cℓ consistency
│   ├── desi_forecast.py # DESI DR1 forecast
│   ├── ode_screening.py # Full ODE integration
│   ├── stability.py     # c_s² parameter scan
│   └── results/
├── notebooks/        # Pipeline runners + numerical output
├── docs/             # Documentation and reports
│   ├── stage1_result.md – stage4_result.md
│   ├── idm_lagrangian_report.md    # Full technical report
│   └── idm_programme_narrative.md  # Narrative overview
```

## Related

- `/home/wsl/entropy-gradient-validation/` — Empirical tests of C3/C4
- `/tmp/class_idm_build/` — CLASS standalone with IDM patches