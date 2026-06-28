#!/usr/bin/env python3
"""
IDM Lagrangian — V(φ) Reconstruction

Pipeline:
1. Check reality condition for canonical quintessence
2. Reconstruct V(z), φ(z) from H(z)
3. Check physical priors (V>0, V''>0)
4. Plot V(φ) and diagnostics
5. Report findings
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.hubble import (
    E2, dE2_dz, Omega_m, Omega_L, f_idm,
    parameters as hubble_params,
)
from src.reconstruction import (
    phi_prime_sq,
    check_reality,
    reconstruct_phi_z,
    V_of_z,
    reconstruct_V_phi,
    V_from_phi,
    check_stability,
    V_to_GeV4, H0_in_MPl,
)

OUT = os.path.dirname(os.path.abspath(__file__)) + "/results"
os.makedirs(OUT, exist_ok=True)
FIG = OUT + "/figs"
os.makedirs(FIG, exist_ok=True)

params = hubble_params()
print("=" * 65)
print("IDM LAGRANGIAN — V(φ) RECONSTRUCTION")
print("=" * 65)
print(f"  H₀ = {params['H0']:.1f}, Ω_m = {params['Omega_m']:.3f}")
print(f"  ε  = {params['eps']:.4f}, zc = {params['zc']:.2f}")

# ─────────────────────────────────────────────
# 1. REALITY CONDITION
# ─────────────────────────────────────────────
print("\n" + "─" * 65)
print("STEP 1: REALITY CONDITION")
print("─" * 65)

z_fine = np.linspace(0.001, 3.0, 5000)
reality = check_reality(z_fine)

print(f"\n  min(φ'²)       = {reality['min_val']:.6f}  at z = {reality['z_min']:.3f}")
print(f"  Negative bins  = {reality['n_negative']}/{len(z_fine)}")
print(f"  Condition:     {'✅ PASS — canonical quintessence viable' if reality['pass'] else '❌ FAIL — phantom crossing'}")

if not reality['pass']:
    # Find boundaries of real region
    phi2 = phi_prime_sq(z_fine)
    pos = np.where(phi2 >= 0)[0]
    if len(pos) > 0:
        print(f"  Valid region: z = {z_fine[pos[0]]:.3f} → {z_fine[pos[-1]]:.3f}")
        print(f"  ({(z_fine[pos[-1]] - z_fine[pos[0]]):.2f} range in z)")
    neg = np.where(phi2 < 0)[0]
    if len(neg) > 0:
        print(f"  Phantom region(s): z ≈ {z_fine[neg[0]]:.3f} → {z_fine[neg[-1]]:.3f}")

# Check ΛCDM for reference
z_lcdm = np.array([0.0, 0.5, 1.0, 2.0, 3.0])
e2_lcdm = params['Omega_m'] * (1+z_lcdm)**3 + params['Omega_L']
de2_lcdm = 3 * params['Omega_m'] * (1+z_lcdm)**2
phi2_lcdm = (de2_lcdm - 3*params['Omega_m']*(1+z_lcdm)**2) / ((1+z_lcdm) * e2_lcdm)
print(f"\n  ΛCDM reference: φ'² = {phi2_lcdm} (should be ~0 everywhere ✓)")

# ─────────────────────────────────────────────
# 2. RECONSTRUCT V(φ)
# ─────────────────────────────────────────────
print("\n" + "─" * 65)
print("STEP 2: V(φ) RECONSTRUCTION")
print("─" * 65)

z_scan = np.linspace(0.001, 3.0, 10000)

if reality['pass'] or reality['frac_negative'] < 0.99:
    phi_s, V_nat, z_valid, V_z = reconstruct_V_phi(z_scan)
    
    # Convert to physical units
    V_gev4 = V_to_GeV4(V_nat)
    gamma_eV4 = V_gev4 * 1e36  # GeV⁴ → eV⁴
    
    print(f"  Valid z range:  {z_valid[0]:.3f} → {z_valid[-1]:.3f}")
    print(f"  φ range:        {phi_s[0]:.4f} → {phi_s[-1]:.4f} M_Pl")
    print(f"  V range:        {V_nat.min():.3e} → {V_nat.max():.3e} (H₀²·M_Pl² nat. units)")
    print(f"  V range:        {V_gev4.min():.3e} → {V_gev4.max():.3e} GeV⁴")
    print(f"  H₀ in M_Pl:    {H0_in_MPl():.3e}")
    
    # 3. PHYSICAL PRIORS
    print("\n" + "─" * 65)
    print("STEP 3: PHYSICAL PRIORS")
    print("─" * 65)
    stab = check_stability(phi_s, V_nat)
    print(f"  V > 0:    {'✅' if stab['V_positive'] else '❌'}  (min = {stab['V_min']:.3e})")
    print(f"  V'' > 0:  {'✅' if stab['V2_positive'] else '❌'}  (min = {stab['V2_min']:.3e})")
    print(f"  c_s²≥0:  {'✅' if stab['c_s2_physical'] else '❌'}")
    
    # 4. EFT INTERPRETATION
    print("\n" + "─" * 65)
    print("STEP 4: EFT INTERPRETATION")
    print("─" * 65)
    if not reality['pass']:
        neg_phi2 = phi_prime_sq(z_fine)
        neg_mask = neg_phi2 < 0
        z_phantom = z_fine[neg_mask]
        print(f"  Phantom region: z ≈ {z_phantom[0]:.2f}–{z_phantom[-1]:.2f}")
        print(f"  → Canonical quintessence alone is INSUFFICIENT")
        print(f"  → IDM requires α_B ≠ 0 (braiding) or non-minimal coupling ξ(φ)R")
        print(f"  → The reconstructed V(φ) is the 'effective' potential")
        print(f"    in the region where canonical approximation holds.")
    else:
        print(f"  Canonical quintessence viable.\n")
    
    # 5. PLOTS
    print("─" * 65)
    print("STEP 5: PLOTS")
    print("─" * 65)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # (a) Reality condition
    ax = axes[0, 0]
    ax.plot(z_fine, phi_prime_sq(z_fine), 'b-', lw=2)
    ax.axhline(0, color='k', ls='--', alpha=0.3)
    if np.any(phi_prime_sq(z_fine) < 0):
        mask = phi_prime_sq(z_fine) < 0
        ax.fill_between(z_fine, phi_prime_sq(z_fine), 0, where=mask, color='r', alpha=0.2, label='Phantom')
    ax.set_xlabel('z')
    ax.set_ylabel("(dφ/dz)² [M_Pl²]")
    ax.set_title("Reality Condition")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    
    # (b) E(z) comparison
    ax = axes[0, 1]
    z_plot = np.linspace(0, 2.5, 500)
    e2_idm = E2(z_plot)
    e2_lcdm = params['Omega_m']*(1+z_plot)**3 + params['Omega_L']
    ax.plot(z_plot, np.sqrt(e2_idm), 'b-', lw=2, label='IDM')
    ax.plot(z_plot, np.sqrt(e2_lcdm), 'k--', lw=1.5, label='ΛCDM')
    ax.set_xlabel('z')
    ax.set_ylabel('E(z) = H(z)/H₀')
    ax.set_title("Expansion History")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    
    # (c) V(φ)
    ax = axes[1, 0]
    ax.plot(phi_s, V_nat, 'b-', lw=2, label='V(φ)')
    ax.set_xlabel('φ [M_Pl]')
    ax.set_ylabel('V(φ) [H₀²·M_Pl²]')
    ax.set_title("Reconstructed V(φ)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    
    # (d) V(z) & φ(z)
    ax = axes[1, 1]
    color1, color2 = 'tab:blue', 'tab:red'
    ax_twin = ax.twinx()
    ax.plot(z_valid, V_nat, color=color1, lw=2, label='V(z)')
    ax_twin.plot(z_valid, phi_s[::-1], color=color2, lw=2, ls='--', label='φ(z)')
    ax.set_xlabel('z')
    ax.set_ylabel('V(z) [H₀²·M_Pl²]', color=color1)
    ax_twin.set_ylabel('φ(z) [M_Pl]', color=color2)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f"{FIG}/v_of_phi.png", dpi=150, bbox_inches='tight')
    print(f"  Saved: {FIG}/v_of_phi.png")
    
    # Save data
    np.savetxt(f"{OUT}/vofphi.txt",
               np.column_stack([phi_s, V_nat, V_gev4, gamma_eV4]),
               header="phi[M_Pl]\tV[H0^2*Mpl^2]\tV[GeV^4]\tV[eV^4]")
    np.savetxt(f"{OUT}/phi_z.txt",
               np.column_stack([z_valid, phi_s[::-1], V_nat, V_z]),
               header="z\tphi[M_Pl]\tV[Mpl^4_nat]\tV(z)_nat")
    np.savetxt(f"{OUT}/reality_check.txt",
               np.column_stack([z_fine, phi_prime_sq(z_fine),
                                E2(z_fine), dE2_dz(z_fine)]),
               header="z\tphi_prime_sq\tE2\tdE2_dz")
    print(f"  Saved: {OUT}/vofphi.txt")
    print(f"  Saved: {OUT}/phi_z.txt")
    
    # 6. SUMMARY
    print("\n" + "=" * 65)
    print("SUMMARY")
    print("=" * 65)
    if reality['pass']:
        print(f"\n  Canonical quintessence: POSSIBLE")
        print(f"  V(φ) shape: {'Runaway (V→0 as φ→∞)' if V_nat[-1] < V_nat[0] else 'Bounded'}")
    else:
        print(f"\n  Canonical quintessence: NOT SUFFICIENT")
        print(f"  → V(φ) reconstructed only in z={z_valid[0]:.2f}–{z_valid[-1]:.2f}")
        print(f"  → Full IDM Lagrangian requires α_B(z) ≠ 0 or ξ(φ)R coupling")
    print(f"\n  Stability:  V>0 = {stab['V_positive']}, V''>0 = {stab['V2_positive']}")
    print(f"  H₀/M_Pl    = {H0_in_MPl():.3e}")
    print(f"  V₀ range   = {V_gev4.min():.3e} – {V_gev4.max():.3e} GeV⁴")
else:
    print("\n  No valid region for canonical reconstruction.")
    print("  → Skip to EFT braiding analysis.")

print("=" * 65)
