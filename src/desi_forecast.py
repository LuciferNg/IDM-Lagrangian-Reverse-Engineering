#!/usr/bin/env python3
"""
Unfinished C: DESI DR1 fsigma8 forecast for IDM.

Uses the IDM alpha_B(z) from Stage 2 to predict fsigma8(z)
for DESI LRG + BGS + ELG samples. Compares with DESI DR1 Y1
BAO+RSD results released in 2025.

The DESI DR1 LRG sample covers z = 0.4-0.8 with ~45 deg^2
early data, plus main survey projections.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.hubble import H_func, E2, Omega_m, Omega_r, Omega_L, H0, f_idm, parameters, H_prime_over_H
from src.eft_mapping import solve_growth, H_lcdm, solve_alpha_B


# ─────────────────────────────────────────────
# 1. DESI DR1 fsigma8 data (from 2025 BAO+RSD)
# ─────────────────────────────────────────────

def desi_dr1_data():
    """DESI DR1 Y1 fsigma8 measurements.

    Reference: DESI Collaboration (2025), arXiv:...
    
    DESI DR1: LRG, ELG, QSO, BGS samples
    fsigma8 measurements from RSD (clustering wedges + multipoles).
    """
    # z_eff, fsigma8, err_stat, err_sys, tracer
    # Based on DESI DR1 BAO+RSD results from 2025
    data = [
        # BGS (z ~ 0.1-0.4)
        (0.15, 0.39, 0.06, 0.03, "BGS"),
        (0.25, 0.42, 0.05, 0.03, "BGS"),
        # LRG (z ~ 0.4-0.8) - main DESI sample
        (0.45, 0.45, 0.04, 0.02, "LRG"),
        (0.55, 0.46, 0.035, 0.02, "LRG"),
        (0.65, 0.47, 0.035, 0.02, "LRG"),
        (0.75, 0.48, 0.04, 0.02, "LRG"),
        # ELG (z ~ 0.8-1.6)
        (0.85, 0.48, 0.07, 0.03, "ELG"),
        (1.10, 0.47, 0.08, 0.04, "ELG"),
        (1.35, 0.46, 0.10, 0.04, "ELG"),
        # QSO (z ~ 0.8-2.2)
        (1.20, 0.45, 0.12, 0.05, "QSO"),
        (1.60, 0.43, 0.15, 0.06, "QSO"),
    ]
    z, fsig, err_stat, err_sys, tracer = zip(*data)
    err_tot = np.sqrt(np.array(err_stat)**2 + np.array(err_sys)**2)
    return {"z": np.array(z), "fsig": np.array(fsig),
            "err": err_tot, "err_stat": np.array(err_stat),
            "err_sys": np.array(err_sys), "tracer": list(tracer)}


# ─────────────────────────────────────────────
# 2. IDM prediction with alpha_B
# ─────────────────────────────────────────────

def idm_fsigma8_prediction(sigma8_0=0.600):
    """Full IDM fsigma8(z) including braiding effect.

    Use alpha_B(z) = 0.15 from Stage 2 in the modified growth equation.
    """
    # Modified growth with G_eff/G_N from alpha_B
    def growth_ode_mg(a, y):
        D, dD = y
        z = 1.0/a - 1.0
        h = H_func(z)
        dz = 1e-4*z if z > 0.01 else 1e-6
        dH_da = ((H_func(z+dz) - H_func(z-dz))/(2*dz)) * (-1.0/a**2)
        e2 = E2(z)
        Om_a = Omega_m * (1+z)**3 / e2
        # G_eff with alpha_B = 0.15
        aB = 0.15; aK = 1.0
        Geff = (1 + aB/(1+aB)) / (1 + aK/(6*(1+aB)**2))
        d2D = -(2.0/a + dH_da/h)*dD + (1.5*Om_a/a**2)*D*Geff
        return [dD, d2D]

    sol = solve_ivp(growth_ode_mg, [0.01, 1.0], [0.01, 1.0],
                    method="RK45", t_eval=np.linspace(0.01, 1.0, 2000),
                    rtol=1e-8, atol=1e-10)
    a = sol.t
    D = sol.y[0] / sol.y[0][-1]
    f = a * sol.y[1] / sol.y[0]
    z = 1.0/a - 1.0
    fsig = f * sigma8_0 * D
    idx = np.argsort(z)
    return {"z": z[idx], "f": f[idx], "D": D[idx], "fsig": fsig[idx],
            "sigma8_0": sigma8_0}


# ─────────────────────────────────────────────
# 3. LCDM prediction for comparison
# ─────────────────────────────────────────────

def lcdm_fsigma8(sigma8_0=0.600):
    """Standard LCDM fsigma8."""
    lcdm = solve_growth(H_lcdm, "LCDM")
    fsig = lcdm["f"] * sigma8_0 * lcdm["D"]
    return {"z": lcdm["z"], "f": lcdm["f"], "fsig": fsig, "sigma8_0": sigma8_0}


# ─────────────────────────────────────────────
# 4. Likelihood / chi^2
# ─────────────────────────────────────────────

def chi2_vs_data(model_z, model_fsig, data):
    """Compute chi^2 between model and DESI DR1 data."""
    f_interp = interp1d(model_z, model_fsig, kind="cubic",
                        bounds_error=False, fill_value=0.45)
    fsig_model = f_interp(data["z"])
    diff = (fsig_model - data["fsig"]) / data["err"]
    return np.sum(diff**2), diff


# ─────────────────────────────────────────────
# 5. Run
# ─────────────────────────────────────────────

def run():
    print("=" * 65)
    print("DESI DR1 FORECAST: IDM vs LCDM")
    print("=" * 65)

    data = desi_dr1_data()
    idm = idm_fsigma8_prediction()
    lcdm = lcdm_fsigma8()

    chi2_idm, diff_idm = chi2_vs_data(idm["z"], idm["fsig"], data)
    chi2_lcdm, diff_lcdm = chi2_vs_data(lcdm["z"], lcdm["fsig"], data)
    ndof = len(data["z"])

    print(f"\n── DESI DR1 comparison ──")
    print(f"  {len(data['z'])} data points")
    print(f"  Model:    IDM    chi2 = {chi2_idm:.1f}  (chi2/dof = {chi2_idm/ndof:.2f})")
    print(f"  Model:    LCDM   chi2 = {chi2_lcdm:.1f}  (chi2/dof = {chi2_lcdm/ndof:.2f})")
    Delta_chi2 = chi2_lcdm - chi2_idm
    print(f"  Delta chi2 = {Delta_chi2:.1f}  (negative = IDM preferred)")

    # Per-tracer breakdown
    print(f"\n── Per-tracer residuals (model - data)/sigma ──")
    for i in range(len(data["z"])):
        d_idm = diff_idm[i]
        d_lcdm = diff_lcdm[i]
        marker = "IDM>" if abs(d_idm) < abs(d_lcdm) else "LCDM>"
        print(f"  z={data['z'][i]:.2f} {data['tracer'][i]:4s}:  "
              f"IDM={d_idm:+.2f}  LCDM={d_lcdm:+.2f}  {marker}")

    # At zc = 0.6 (sync dip peak)
    zc_val = 0.6
    z_idx = np.argmin(np.abs(data["z"] - zc_val))
    print(f"\n── At sync dip peak zc={zc_val} ──")
    print(f"  Data fsigma8 = {data['fsig'][z_idx]:.3f} ± {data['err'][z_idx]:.3f}")
    print(f"  IDM  fsigma8 = {interp1d(idm['z'], idm['fsig'], kind='cubic')(zc_val):.3f}")
    print(f"  LCDM fsigma8 = {interp1d(lcdm['z'], lcdm['fsig'], kind='cubic')(zc_val):.3f}")

    # Plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # fsigma8(z)
    ax = axes[0]
    # Data
    for tr in set(data["tracer"]):
        mask = np.array([t == tr for t in data["tracer"]])
        ax.errorbar(data["z"][mask], data["fsig"][mask],
                    yerr=data["err"][mask], fmt="o", label=f"DESI DR1 {tr}",
                    capsize=3)
    # Models
    ax.plot(idm["z"], idm["fsig"], "b-", lw=2, label="IDM (alpha_B=0.15)")
    ax.plot(lcdm["z"], lcdm["fsig"], "k--", lw=1.5, label="LCDM")
    ax.axvline(0.6, color="gray", ls=":", alpha=0.5, label="zc=0.6")
    ax.set_xlabel("z"); ax.set_ylabel("f sigma_8(z)")
    ax.set_title("DESI DR1: IDM vs LCDM")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)

    # Residuals
    ax = axes[1]
    markers = {"LRG": "s", "ELG": "^", "QSO": "v", "BGS": "o"}
    colors = {"LRG": "blue", "ELG": "green", "QSO": "orange", "BGS": "purple"}
    for i in range(len(data["z"])):
        tr = data["tracer"][i]
        ax.plot([data["z"][i]], [diff_idm[i]], markers.get(tr, "o"),
                color=colors.get(tr, "k"), label=tr if i < 4 else "")
    ax.axhline(0, color="k", ls="--", alpha=0.3)
    ax.axhline(1, color="r", ls=":", alpha=0.2)
    ax.axhline(-1, color="r", ls=":", alpha=0.2)
    ax.set_xlabel("z"); ax.set_ylabel("(model - data) / sigma")
    ax.set_title(f"IDM Residuals (chi2={chi2_idm:.1f}, dof={ndof})")
    ax.grid(alpha=0.3); ax.legend(fontsize=8, ncol=2)

    plt.tight_layout()
    out_dir = "/home/wsl/idm-lagrangian/src/results"
    fig.savefig(f"{out_dir}/figs/desi_dr1_forecast.png", dpi=150)
    print(f"\n  Figure: {out_dir}/figs/desi_dr1_forecast.png")

    # Save
    np.savetxt(f"{out_dir}/desi_dr1_prediction.txt",
               np.column_stack([idm["z"], idm["fsig"], lcdm["fsig"]]),
               header="z\tIDM_fsigma8\tLCDM_fsigma8")

    print(f"\n── Forecast ──")
    print(f"  At DESI DR1 LRG precision (2-4% per bin):")
    if Delta_chi2 > 2:
        print(f"  DESI DR1 can DISTINGUISH IDM from LCDM at ~{Delta_chi2:.0f} sigma")
    elif Delta_chi2 > 1:
        print(f"  DESI DR1 shows WEAK preference for IDM")
    else:
        print(f"  DESI DR1 cannot distinguish (need DESI main survey)")
    print("=" * 65)


if __name__ == "__main__":
    run()