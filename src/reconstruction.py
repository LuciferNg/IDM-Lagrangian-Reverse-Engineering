"""V(φ) reconstruction from IDM background H(z).

Corrected derivation, M_Pl = 1 units (reduced Planck mass).

From Friedmann equations:
    H² = (ρ_m + ρ_φ)/3
    2Ȟ + 3H² = -P_φ

For canonical quintessence:
    ρ_φ = ½φ̇² + V
    P_φ = ½φ̇² - V
    ⇒ φ̇² = ρ_φ + P_φ
    ⇒ V = (ρ_φ - P_φ)/2

Expressing in redshift (d/dt = -(1+z)H d/dz, Ȟ = -(1+z)HH'):

    φ'² = [dE²/dz - 3Ω_m(1+z)²] / [(1+z) E²]
    V(z) = H₀²/2 · [ (dE²/dz)/(1+z) + (6 - 1/E²·dE²/dz)·E² - 3Ω_m(1+z)³ ]

    where E²(z) = H²(z)/H₀²

Reality condition: dE²/dz - 3Ω_m(1+z)² ≥ 0
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
from src.hubble import (
    E2, dE2_dz, H_func, H_prime_over_H,
    H0, Omega_m, parameters as hubble_params,
)
from typing import Tuple


# ─────────────────────────────────────────────
# Corrected φ'(z)
# ─────────────────────────────────────────────

def phi_prime_sq(z: np.ndarray) -> np.ndarray:
    """(dφ/dz)² in M_Pl² units.

    φ'² = [dE²/dz - 3Ω_m(1+z)²] / [(1+z) E²]

    For ΛCDM (dE²/dz = 3Ω_m(1+z)²): φ'² = 0 ✓
    """
    e2 = E2(z)
    de2 = dE2_dz(z)
    numerator = de2 - 3.0 * Omega_m * (1.0 + z) ** 2
    return numerator / ((1.0 + z) * e2)


def check_reality(z_grid: np.ndarray) -> dict:
    """Check φ'² ≥ 0 for canonical quintessence viability."""
    phi2 = phi_prime_sq(z_grid)
    n_neg = int(np.sum(phi2 < 0))
    min_val = float(np.min(phi2))
    z_min = float(z_grid[np.argmin(phi2)])
    return {
        "pass": n_neg == 0,
        "min_val": min_val,
        "z_min": z_min,
        "n_negative": n_neg,
        "frac_negative": n_neg / len(z_grid),
    }


# ─────────────────────────────────────────────
# Corrected V(z)
# ─────────────────────────────────────────────

def V_of_z(z: np.ndarray) -> np.ndarray:
    """Potential V(z) in M_Pl⁴ (natural) units.

    Derivation in M_Pl=1:
    V = (ρ_φ - P_φ)/2
    ρ_φ = 3H² - ρ_m = 3H₀²E² - 3H₀²Ω_m(1+z)³
    P_φ = -(2Ȟ + 3H²)
    Ȟ = dH/dt = -(1+z)HH' = -(1+z)H₀√E · H₀·½dE²/dz/√E
      = -½(1+z)H₀²·dE²/dz/E²·E² ... 
      
    Let me simplify differently.
    
    ρ_φ = 3H₀²E² - 3H₀²Ω_m(1+z)³
    P_φ = -2Ȟ - 3H²
    
    Ȟ = -(1+z)H·H' = -(1+z)(H₀√E)(H₀·½dE²/dz/√E) = -½(1+z)H₀²·dE²/dz
    
    So P_φ = (1+z)H₀²·dE²/dz - 3H₀²E²
    
    V = (ρ_φ - P_φ)/2
      = [3H₀²E² - 3H₀²Ω_m(1+z)³ - (1+z)H₀²·dE²/dz + 3H₀²E²]/2
      = H₀²[6E² - 3Ω_m(1+z)³ - (1+z)·dE²/dz]/2
      = H₀²[3E² - 1.5Ω_m(1+z)³ - 0.5(1+z)·dE²/dz]

    In M_Pl=1 units, H₀ is a dimensionless number.
    Let's keep numerical values and convert to GeV⁴ at output.
    """
    e2 = E2(z)
    de2 = dE2_dz(z)
    return 3.0 * e2 - 1.5 * Omega_m * (1.0 + z) ** 3 - 0.5 * (1.0 + z) * de2


# ─────────────────────────────────────────────
# φ(z) integration
# ─────────────────────────────────────────────

def reconstruct_phi_z(z_grid: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute φ(z) by integrating φ'(z) from high z to low z.

    Parameters
    ----------
    z_grid : ndarray
        Redshift grid (any order, will sort descending).

    Returns
    -------
    z : ndarray
        Sorted descending redshift (trimmed to real region).
    phi : ndarray
        Field value φ(z) in M_Pl.
    """
    z = np.sort(z_grid)[::-1]
    phi2 = phi_prime_sq(z)

    # Trim to real region
    if np.any(phi2 < 0):
        # Find first negative from high z
        neg = np.where(phi2 < 0)[0]
        if len(neg) > 0:
            first_neg = neg[0]
            if first_neg > 0:
                z = z[:first_neg]
                phi2 = phi2[:first_neg]
            else:
                raise ValueError("φ'² < 0 at highest z — cannot integrate canonical field")

    phi_prime = np.sqrt(phi2)
    phi = cumulative_trapezoid(phi_prime, z, initial=0.0)
    # Normalize: φ → 0 at high z (field rolls from 0)
    phi -= phi[-1]
    return z, phi


# ─────────────────────────────────────────────
# V(φ) from z interpolation
# ─────────────────────────────────────────────

def reconstruct_V_phi(
    z_grid: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Full V(φ) reconstruction pipeline.

    Returns
    -------
    phi : ndarray
        φ(z) in M_Pl.
    V_nat : ndarray
        V(φ) in H₀²·M_Pl² units (natural dimensionless).
    z_valid : ndarray
        Redshift where reconstruction is valid.
    V_z : ndarray
        V(z) in same units.
    """
    z_valid, phi = reconstruct_phi_z(z_grid)
    V_z = V_of_z(z_valid)

    # Sort by φ
    phi_s, V_s = V_from_phi(phi, V_z)

    return phi_s, V_s, z_valid, V_z


def V_from_phi(phi: np.ndarray, V: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Sort V(φ) by monotonic φ."""
    idx = np.argsort(phi)
    return phi[idx], V[idx]


# ─────────────────────────────────────────────
# Unit conversion factors
# ─────────────────────────────────────────────

def H0_in_MPl() -> float:
    """H₀ in reduced Planck mass units.

    69 km/s/Mpc → 2.133e-33 GeV → 8.76e-61 M_Pl
    """
    H0_GeV = 69.0 * 3.241e-42  # km/s/Mpc → GeV
    M_Pl_GeV = 2.435e18
    return H0_GeV / M_Pl_GeV


def V_to_GeV4(V_nat: np.ndarray) -> np.ndarray:
    """Convert V from H₀²·M_Pl² units to GeV⁴.

    V_phys = V_nat · H₀² · M_Pl²
    = V_nat · (H0_in_MPl · M_Pl)² · M_Pl²
    = V_nat · H0_in_MPl² · M_Pl⁴
    """
    h0_mpl = H0_in_MPl()
    M_Pl_GeV = 2.435e18
    return V_nat * h0_mpl ** 2 * M_Pl_GeV ** 4


# ─────────────────────────────────────────────
# Physical priors
# ─────────────────────────────────────────────

def check_stability(phi: np.ndarray, V: np.ndarray) -> dict:
    """Check V(φ) physical priors."""
    dV = np.gradient(V, phi)
    d2V = np.gradient(dV, phi)

    return {
        "V_min": float(np.min(V)),
        "V_positive": bool(np.all(V > 0)),
        "V2_min": float(np.min(d2V)),
        "V2_positive": bool(np.all(d2V > -1e-10)),  # Allow tiny numerical noise
        "c_s2_physical": bool(np.all(V > 0) and np.all(d2V > -1e-10)),
    }