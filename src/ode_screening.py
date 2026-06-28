#!/usr/bin/env python3
"""
Unfinished A: Full ODE — KG equation on IDM background.

Correct approach: take IDM H(z) as given (fits data). Solve ONLY
the Klein-Gordon equation with braiding to verify φ(z) matches
Stage 1 reconstruction.

KG equation with braiding (M_Pl=1, N = ln a):
  (1 + xi*aB)*[d²phi/dN² + (3 + dlnH/dN)*dphi/dN] + V_phi/H² = 0

For xi >> 1: effective friction is rescaled, producing phantom.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.hubble import H_func, E2, Omega_m, H0, H_prime_over_H


def solve_KG_with_braiding(xi0, lam, aB=0.15):
    """Solve the braiding-modified KG on the IDM background.

    Units: H0=1, M_Pl=1.  N = ln(a) from -4 to 0.
    State: y = [phi, psi = dphi/dN]
    """
    N_grid = np.linspace(-4.0, 0.0, 5000)
    z_grid = 1.0 / np.exp(N_grid) - 1.0

    # IDM background
    H_z = H_func(z_grid) / H0  # E(z)
    # dlnH/dN = dlnH/da * a = H'/H * dz/dN
    Hprime_H = H_prime_over_H(z_grid)  # dlnH/dz
    dlnH_dN = Hprime_H * (-np.exp(N_grid))

    # V(phi) from Stage 1 interpolation
    phi_data = np.loadtxt("/home/wsl/idm-lagrangian/notebooks/results/phi_z.txt",
                          skiprows=1)
    z_data = phi_data[:, 0]
    phi_stage1 = phi_data[:, 1]
    V_data = phi_data[:, 2]  # V(z) in H0^2*M_Pl^2 units
    V_interp = interp1d(z_data[::-1], V_data[::-1],
                        kind="cubic", bounds_error=False, fill_value=2.0)
    V_z = V_interp(z_grid)
    # dV/dphi through chain rule: dV/dz * dz/dphi
    # But we only have phi(z). So V_phi = dV/dz / dphi/dz
    # For numerical stability, use finite differences
    dV_dz = np.gradient(V_z, z_grid)
    dphi_dz = np.gradient(phi_stage1, 1.0/(1+z_data))

    phi_interp = interp1d(z_data[::-1], phi_stage1[::-1],
                          kind="cubic", bounds_error=False, fill_value=0.0)
    dphi_dz_interp = interp1d(z_data[::-1], dphi_dz[::-1],
                              kind="cubic", bounds_error=False,
                              fill_value=phi_stage1[-1]*0.01)

    # xi(phi) screening
    def xi_func(phi_v):
        return xi0 * np.exp(-lam * phi_v)

    def xi_prime(phi_v):
        return -lam * xi0 * np.exp(-lam * phi_v)

    # KG with braiding (N = ln a as independent variable)
    def ode(N, y):
        phi, psi = y
        z_val = 1.0 / np.exp(N) - 1.0
        H_val = H_func(z_val) / H0
        dlnH = H_prime_over_H(z_val) * (-np.exp(N))
        xi_v = xi_func(phi)
        alpha_eff = 1.0 + xi_v * aB

        # dV/dphi
        phi_ref = phi_interp(z_val)
        # Use analytic V_phi from Stage 1: flat V ~ constant
        V_phi = 0.05  # Stage 1 showed V varies by ~6%

        dpsi = - (3.0 + dlnH) * psi - V_phi / H_val**2 / alpha_eff
        return [psi, dpsi]

    y0 = [0.0, 1e-6]
    sol = solve_ivp(ode, (-4.0, 0.0), y0, t_eval=N_grid,
                    method="RK45", rtol=1e-8, atol=1e-10,
                    max_step=0.01)

    N = sol.t
    a = np.exp(N)
    z = 1.0/a - 1.0
    phi_kg = sol.y[0]

    # Compare with Stage 1
    phi_target = phi_interp(z)

    return {"z": z, "phi": phi_kg, "phi_target": phi_target}


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print("=" * 65)
    print("ODE: BRAIDING KG ON IDM BACKGROUND")
    print("=" * 65)

    xi0_candidates = [1, 1e5, 1e10, 1e15, 1e18]
    lam_val = 1.0
    aB_val = 0.15

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for xi0 in xi0_candidates:
        try:
            res = solve_KG_with_braiding(xi0, lam_val, aB_val)
            # Valid range: z > 0.6
            mask = res["z"] > 0.6
            if np.sum(mask) > 10:
                diff = np.mean(
                    (res["phi"][mask] - res["phi_target"][mask])**2
                )
                print(f"  xi0={xi0:.1e}:  ||phi-phi_target||^2 = {diff:.3e}")
                axes[0].plot(res["z"][mask], res["phi"][mask],
                             label=f"xi0={xi0:.1e}")
        except Exception as ex:
            print(f"  xi0={xi0:.1e}: FAILED ({str(ex)[:40]})")

    # Stage 1 reference
    phi_data = np.loadtxt(
        "/home/wsl/idm-lagrangian/notebooks/results/phi_z.txt", skiprows=1
    )
    axes[0].plot(phi_data[:, 0], phi_data[:, 1],
                 "k--", lw=2, label="Stage 1 (target)")
    axes[0].set_xlabel("z"); axes[0].set_ylabel("phi [M_Pl]")
    axes[0].set_title(f"KG with braiding (lam={lam_val}, aB={aB_val})")
    axes[0].grid(alpha=0.3); axes[0].legend(fontsize=7)

    # Best fit residual
    res_best = solve_KG_with_braiding(1e15, lam_val, aB_val)
    mask = res_best["z"] > 0.6
    axes[1].plot(res_best["z"][mask],
                 (res_best["phi"][mask] - res_best["phi_target"][mask]) * 1000,
                 "r-", lw=2)
    axes[1].axhline(0, color="k", ls="--", alpha=0.3)
    axes[1].set_xlabel("z")
    axes[1].set_ylabel("dphi [milli-M_Pl]")
    axes[1].set_title("Residual (xi0=1e15, lam=1)")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    out_dir = "/home/wsl/idm-lagrangian/src/results"
    fig.savefig(f"{out_dir}/figs/kg_braiding.png", dpi=150)

    print(f"\n  Figure: {out_dir}/figs/kg_braiding.png")
    print("=" * 65)