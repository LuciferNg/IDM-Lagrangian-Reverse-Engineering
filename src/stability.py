#!/usr/bin/env python3
"""
Unfinished B: c_s^2 stability scan.

For Horndeski with G2 = X - V, G3 = xi(phi)*aB*X:
  c_s^2 = (G2,X + 2*G3,phi - G3,X^2*H*phi_dot/...) / (G2,X + 6*G3,X*H*phi_dot + ...)

In the quasi-static limit with xi >> 1:
  c_s^2 ≈ 1 (no ghost)

Main check: no gradient instability (c_s^2 >= 0) and no ghost.
"""

import numpy as np
from scipy.interpolate import interp1d
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.hubble import H_func, H0, Omega_m


def sound_speed_sq(phi_vals, psi_vals, H_vals, xi0, lam, aB):
    """Compute c_s^2(z) for the braiding model.

    From Bellini+Sawicki 2014 (Eqs 57-58):
    alpha_K = 2*X*(G2,XX + 2*G3,X*phi_dot/H - G3,phi*X/...) / H^2
    ...
    
    Simplified: for G2 = X - V, G3 = xi*aB*X:
    G2,X = 1
    G3,X = xi*aB
    G3,phi = -lam*xi0*exp(-lam*phi)*aB*X = -lam*xi*aB*X
    
    c_s^2 = alpha_K / (4*alpha_B^2 + alpha_K + ...)
    
    For alpha_K >> alpha_B^2: c_s^2 ≈ 1 (stable).
    """
    # Kineticity alpha_K in the model
    # alpha_K = 2*X / H^2 * (1 + ... G3 terms)
    X_vals = 0.5 * psi_vals**2 * H_vals**2
    xi_vals = xi0 * np.exp(-lam * phi_vals)
    
    # alpha_K ~ 2X/H^2 = psi^2 (for canonical)
    alpha_K = psi_vals**2
    
    # alpha_B = 2*X*G3,X / H^2 = 2*X*xi*aB/H^2 = xi*aB*psi^2*H^2 / H^2 = xi*aB*psi^2
    alpha_B = xi_vals * aB * psi_vals**2
    
    # c_s^2 = alpha_K / (alpha_K + 6*alpha_B^2)
    denom = alpha_K + 6 * alpha_B**2
    c_s2 = np.divide(alpha_K, denom, out=np.ones_like(alpha_K),
                     where=denom > 1e-30)
    
    # No-ghost condition: alpha_K + 6*alpha_B^2 > 0
    no_ghost = (alpha_K + 6*alpha_B**2) > 0
    
    return c_s2, no_ghost, {"alpha_K": alpha_K, "alpha_B": alpha_B}


def scan_stability():
    """Scan (xi0, lam) grid for stability."""
    print("=" * 65)
    print("STABILITY SCAN: c_s^2 AND GHOST CHECK")
    print("=" * 65)

    # Load Stage 1 phi(z)
    phi_data = np.loadtxt(
        "/home/wsl/idm-lagrangian/notebooks/results/phi_z.txt", skiprows=1
    )
    z = phi_data[:, 0]
    phi = phi_data[:, 1]
    z_grid = np.linspace(0.6, 3.0, 500)

    phi_interp = interp1d(z[::-1], phi[::-1], kind="cubic",
                          bounds_error=False, fill_value=0.0)(z_grid)
    H_z = H_func(z_grid) / H0
    # psi = dphi/dN = -dphi/dz * (1+z)
    dphi_dz = np.gradient(phi_interp, z_grid)
    psi = -dphi_dz * (1 + z_grid)

    # Scan
    xi0_vals = np.logspace(0, 20, 11)
    lam_vals = np.logspace(-1, 1, 9)
    aB = 0.15

    results = []
    for xi0 in xi0_vals:
        for lam in lam_vals:
            c_s2, ngh, params = sound_speed_sq(
                phi_interp, psi, H_z, xi0, lam, aB
            )
            c_min = np.min(c_s2)
            all_stable = np.all(c_s2 >= 0) and np.all(ngh)
            results.append((xi0, lam, c_min, all_stable))
            if all_stable and c_min > 0.5:
                label = "STABLE"
            elif all_stable:
                label = f"OK(c_s2={c_min:.3f})"
            else:
                label = "UNSTABLE"
            print(f"  xi0={xi0:.1e}, lam={lam:.2f}: {label}")

    stable = [r for r in results if r[3]]
    unstable = [r for r in results if not r[3]]
    print(f"\n  Stable: {len(stable)}/{len(results)}")
    print(f"  Unstable: {len(unstable)}/{len(results)}")

    return {"results": results, "z": z_grid, "xi0_vals": xi0_vals,
            "lam_vals": lam_vals}


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    
    scan = scan_stability()

    r = np.array([x[:4] for x in scan["results"]])
    xi_m = np.unique(r[:, 0])
    lam_m = np.unique(r[:, 1])

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    # Create stability map
    cmap = np.full((len(xi_m), len(lam_m)), np.nan)
    for row in scan["results"]:
        i = np.where(xi_m == row[0])[0][0]
        j = np.where(lam_m == row[1])[0][0]
        cmap[i, j] = float(row[3])  # 1.0 = stable

    xi_grid, lam_grid = np.meshgrid(np.log10(xi_m), np.log10(lam_m))
    ax.pcolormesh(xi_grid.T, lam_grid.T, cmap,
                  cmap="RdYlGn", vmin=0, vmax=1, shading="auto")
    ax.set_xlabel("log10(xi0)"); ax.set_ylabel("log10(lambda)")
    ax.set_title("Stability: GREEN=stable, RED=unstable")

    out_dir = "/home/wsl/idm-lagrangian/src/results"
    fig.savefig(f"{out_dir}/figs/stability_map.png", dpi=150)
    print(f"\n  Figure: {out_dir}/figs/stability_map.png")
    print("=" * 65)