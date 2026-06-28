# IDM-Lagrangian-Reverse-Engineering

**GitHub:** github.com/LuciferNg/IDM-Lagrangian-Reverse-Engineering

The sync dip cross-correlation + V(φ) reconstruction follows the
curator.md packaging rules. The Cℓ data is tracked in results/; CLASS
build stays in /tmp/class_idm_build/ (not tracked).

## Key files
- `src/hubble.py` — IDM background with sync dip
- `src/reconstruction.py` — V(φ) from H(z)
- `src/eft_mapping.py` — α_B(z) from growth
- `src/eom_verify.py` — SymPy EoM
- `docs/idm_lagrangian_report.md` — Full technical report