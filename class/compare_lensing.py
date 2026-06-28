#!/usr/bin/env python3
"""Compare CLASS lensing: LCDM vs IDM"""
import numpy as np

base = "/tmp/class_idm_build/output"

def read_cl(path):
    data = np.loadtxt(path, comments="#")
    ell = data[:,0].astype(int)
    cl_pp = data[:,5]  # Cℓ^{φφ}
    cl_tp = data[:,6]  # Cℓ^{Tφ}
    return ell, cl_pp, cl_tp

ell_l, pp_l, tp_l = read_cl(f"{base}/test_lcdm00_cl_lensed.dat")
ell_i, pp_i, tp_i = read_cl(f"{base}/test_idm00_cl_lensed.dat")

# Cℓ^{κκ} = [ℓ(ℓ+1)/2]² × Cℓ^{φφ}
kk_l = (ell_l*(ell_l+1)/2)**2 * pp_l
kk_i = (ell_i*(ell_i+1)/2)**2 * pp_i
Dl_l = ell_l*(ell_l+1)*kk_l/(2*np.pi)
Dl_i = ell_i*(ell_i+1)*kk_i/(2*np.pi)

print(f"{'ℓ':>6}  {'Dℓ_LCDM':>12}  {'Dℓ_IDM':>12}  {'I/L':>8}  {'dev%':>8}")
for L in [40, 100, 200, 400, 600, 800, 1000, 1500, 2000]:
    idx = np.searchsorted(ell_l, L)
    r = Dl_i[idx]/Dl_l[idx]
    print(f"{ell_l[idx]:6d}  {Dl_l[idx]:12.4e}  {Dl_i[idx]:12.4e}  {r:7.5f}  {(r-1)*100:+7.2f}%")

mask = (ell_l >= 50) & (ell_l <= 2000)
mean = np.mean(Dl_i[mask]/Dl_l[mask])
print(f"\nMean IDM/LCDM ℓ=50-2000: {mean:.5f} ({(mean-1)*100:+.2f}%)")

# Save for Test 3
np.savetxt(f"{base}/cl_kk_lcdm.txt", 
           np.column_stack([ell_l, kk_l, Dl_l]),
           fmt="%d %.6e %.6e", header="ell Ckk, Dkk")
np.savetxt(f"{base}/cl_kk_idm.txt",
           np.column_stack([ell_i, kk_i, Dl_i]),
           fmt="%d %.6e %.6e", header="ell Ckk, Dkk")
print(f"\nFiles saved:")
print(f"  {base}/cl_kk_lcdm.txt")
print(f"  {base}/cl_kk_idm.txt")