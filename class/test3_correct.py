"""Test 3 with correct Planck bandpowers"""
import numpy as np

out = "/tmp/class_idm_build/output"
ell, pp_l = np.loadtxt(f"{out}/test_lcdm00_cl_lensed.dat", comments="#", usecols=(0,5), unpack=True)
_, pp_i = np.loadtxt(f"{out}/test_idm00_cl_lensed.dat", comments="#", usecols=(0,5), unpack=True)

# Planck 2018 VIII Table 3: lmin, lmax, Dl_x_1e10, stat, sys, total
planck = np.array([
    [2,8,-0.050,0.021,0.016,0.026],[8,47,1.414,0.130,0.039,0.136],
    [47,84,4.621,0.328,0.075,0.337],[84,125,7.108,0.387,0.084,0.396],
    [125,163,8.913,0.465,0.089,0.474],[163,203,10.063,0.534,0.094,0.542],
    [203,244,10.824,0.619,0.110,0.629],[244,285,11.418,0.698,0.118,0.708],
    [285,327,11.885,0.793,0.129,0.804],[327,370,12.087,0.882,0.142,0.893],
    [370,414,12.101,0.965,0.155,0.977],[414,458,12.102,1.033,0.162,1.046],
    [458,504,11.910,1.056,0.164,1.069],[504,593,11.565,1.161,0.131,1.168],
    [593,679,11.131,1.075,0.108,1.080],[679,767,10.543,1.025,0.105,1.031],
    [767,858,9.906,0.981,0.100,0.986],[858,1000,8.681,0.880,0.093,0.885],
    [1000,1144,7.363,0.778,0.099,0.785],[1144,1290,6.160,0.695,0.097,0.702],
    [1290,1461,5.022,0.608,0.091,0.615],[1461,1635,4.147,0.533,0.086,0.540],
    [1635,1830,3.340,0.462,0.079,0.469],[1830,2048,2.657,0.397,0.073,0.404],
])

chi2_l = chi2_i = 0.0
print("{:>8}  {:>11}  {:>11}  {:>11}  {:>8}  {:>8}".format("l_bin","Dl_P18","Dl_LCDM","Dl_IDM","P_l","P_i"))
for row in planck:
    lmin, lmax = int(row[0]), int(row[1])
    dp, de = row[2]*1e-10, row[5]*1e-10
    mask = (ell >= lmin) & (ell < lmax)
    dl_l = np.mean(pp_l[mask])
    dl_i = np.mean(pp_i[mask])
    pl = (dp-dl_l)/de
    pi_ = (dp-dl_i)/de
    chi2_l += pl**2
    chi2_i += pi_**2
    print("{:3d}-{:<3d}  {:10.3f}  {:10.3f}  {:10.3f}  {:7.2f}  {:7.2f}".format(
        lmin, lmax, dp/1e-10, dl_l/1e-10, dl_i/1e-10, pl, pi_))

ndof = len(planck)
print(f"\nchi2_LCDM = {chi2_l:.1f} (ndof={ndof})")
print(f"chi2_IDM  = {chi2_i:.1f} (ndof={ndof})")
print(f"Delta_chi2 = {chi2_i-chi2_l:+.1f}")
ratio = np.mean(pp_i[(ell>=50)&(ell<=2000)]/pp_l[(ell>=50)&(ell<=2000)])
print(f"IDM/LCDM Dl ratio l=50-2000: {ratio:.4f} ({(ratio-1)*100:+.2f}%)")
prec = np.mean(planck[:,5]/planck[:,2])*100
print(f"Planck precision: {prec:.0f}%")
if abs(ratio-1)*100 < prec:
    print("IDM suppression UNDETECTABLE at Planck precision")
