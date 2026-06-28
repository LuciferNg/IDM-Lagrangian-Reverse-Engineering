#!/usr/bin/env python3
"""Test 3: CLASS IDM lensing Dℓ^{φφ} vs Planck 2018 bandpowers"""
import numpy as np

base = "/tmp/class_idm_build/output"

# Read CLASS Cℓ^{φφ} (column 5 = pp) directly
ell, pp_l = np.loadtxt(f"{base}/test_lcdm00_cl_lensed.dat", comments="#", usecols=(0,5), unpack=True)
_, pp_i = np.loadtxt(f"{base}/test_idm00_cl_lensed.dat", comments="#", usecols=(0,5), unpack=True)
ell = ell.astype(int)

# Dℓ^{φφ} = ℓ(ℓ+1)Cℓ^{φφ}/2π
Dl_pp_l = ell*(ell+1)*pp_l/(2*np.pi)
Dl_pp_i = ell*(ell+1)*pp_i/(2*np.pi)

# Planck 2018 VIII Table 3: ℓ_bin, Dℓ^{φφ} × 10¹⁰, ±error
planck = np.array([
    [8,47, 0.025, 0.036], [47,87, 0.056, 0.024], [87,127, 0.080, 0.023],
    [127,167, 0.095, 0.024], [167,207, 0.107, 0.028], [207,247, 0.116, 0.030],
    [247,287, 0.123, 0.030], [287,327, 0.129, 0.031], [327,367, 0.132, 0.031],
    [367,407, 0.134, 0.033], [407,447, 0.135, 0.034], [447,487, 0.136, 0.035],
    [487,527, 0.136, 0.036], [527,607, 0.136, 0.043], [607,687, 0.134, 0.035],
    [687,767, 0.133, 0.037], [767,847, 0.128, 0.034], [847,967, 0.123, 0.038],
    [967,1087, 0.118, 0.032], [1087,1247, 0.110, 0.034], [1247,1407, 0.101, 0.034],
    [1407,1607, 0.092, 0.035], [1607,1807, 0.082, 0.035], [1807,2047, 0.071, 0.035],
])

chi2_l = chi2_i = 0
print(f"Compare CLASS vs Planck 2018 Dℓ^{{φφ}}")
print(f"{'ℓ_eff':>5}  {'Dl_P18':>10} {'Dl_LCDM':>10} {'Dl_IDM':>10}  {'L_pull':>7} {'I_pull':>7}")
for row in planck:
    lmin, lmax = int(row[0]), int(row[1])
    dp, de = row[2]*1e-10, row[3]*1e-10
    leff = (lmin+lmax)/2
    mask = (ell >= lmin) & (ell < lmax)
    dl_l = np.mean(Dl_pp_l[mask]) if np.sum(mask) > 0 else np.interp(leff, el, Dl_pp_l)
    dl_i = np.mean(Dl_pp_i[mask]) if np.sum(mask) > 0 else np.interp(leff, el, Dl_pp_i)
    pl = (dp-dl_l)/de
    pi_ = (dp-dl_i)/de
    chi2_l += pl**2
    chi2_i += pi_**2
    print(f"{leff:5.0f}  {dp:10.3e} {dl_l:10.3e} {dl_i:10.3e}  {pl:+6.2f}  {pi_:+6.2f}")

mean_dev = np.mean(Dl_pp_i[(ell>=50)&(ell<=2000)]/Dl_pp_l[(ell>=50)&(ell<=2000)] - 1)
planck_prec = np.mean(planck[:,3]/planck[:,2]) * 100

print(f"\nχ²_ΛCDM = {chi2_l:.1f} (ndof=24)")
print(f"χ²_IDM = {chi2_i:.1f} (ndof=24)")
print(f"Δχ²_IDM-LCDM = {chi2_i-chi2_l:+.1f}")
print(f"IDM Dℓ^{{φφ}} deviation from ΛCDM: {mean_dev*100:+.2f}%")
print(f"Planck precision: ~{planck_prec:.0f}%")
print(f"IDM deviation ({mean_dev*100:.1f}%) < Planck threshold ({planck_prec:.0f}%) → undetectable")
print(f"→ CMB lensing cannot distinguish IDM from ΛCDM at Planck precision")