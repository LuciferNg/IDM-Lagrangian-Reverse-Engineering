"""IDM background H(z) implementation.

Provides the sync-dip modified expansion history following the
Information-Dynamics Model v4.0 parameters.

IDM sync dip form:
    f(z) = 1 - ε·(z/zc)·exp(-z/zc)

    H²(z)/H₀² = Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_L·f(z)
"""

import numpy as np
from typing import Tuple

# IDM v4.0 MAP parameters
H0 = 69.0          # km/s/Mpc
Omega_m = 0.297
Omega_r = 4.15e-5 / (H0 / 100) ** 2
Omega_L = 1.0 - Omega_m - Omega_r
eps = 0.1545       # sync dip amplitude (MAP)
zc = 0.60          # sync dip peak redshift
sigma = 0.3        # dip width (not used in the x*exp(-x) form, kept for reference)

# Planck units
M_Pl = 2.435e18    # GeV, reduced Planck mass
M_Pl_GeV = M_Pl
M_Pl_eV = M_Pl_GeV * 1e9


# ─────────────────────────────────────────────
# IDM H(z)
# ─────────────────────────────────────────────

def f_idm(z: np.ndarray) -> np.ndarray:
    """Sync dip function: f(z) = 1 - ε·(z/zc)·exp(-z/zc).

    Parameters
    ----------
    z : array_like
        Redshift.

    Returns
    -------
    array_like
        f(z), the sync dip modulation factor.
    """
    x = z / zc
    return 1.0 - eps * x * np.exp(-x)


def E2(z: np.ndarray) -> np.ndarray:
    """Normalized Hubble rate squared: E²(z) = H²(z)/H₀².

    Parameters
    ----------
    z : array_like
        Redshift.

    Returns
    -------
    array_like
        E²(z) = Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_L·f(z).
    """
    return Omega_m * (1 + z) ** 3 + Omega_r * (1 + z) ** 4 + Omega_L * f_idm(z)


def H_func(z: np.ndarray) -> np.ndarray:
    """Hubble rate H(z) in km/s/Mpc.

    Parameters
    ----------
    z : array_like
        Redshift.

    Returns
    -------
    array_like
        H(z).
    """
    return H0 * np.sqrt(E2(z))


# ─────────────────────────────────────────────
# Derivatives (analytic)
# ─────────────────────────────────────────────

def df_idm_dz(z: np.ndarray) -> np.ndarray:
    """Derivative of sync dip function df/dz.

    f(z) = 1 - ε·(z/zc)·exp(-z/zc)
    df/dz = -ε/zc · exp(-z/zc) + ε·(z/zc)·(1/zc)·exp(-z/zc)
          = ε/zc · exp(-z/zc) · (z/zc - 1)
    """
    x = z / zc
    return (eps / zc) * np.exp(-x) * (x - 1.0)


def dE2_dz(z: np.ndarray) -> np.ndarray:
    """Derivative of E²(z) w.r.t. z."""
    return (
        3.0 * Omega_m * (1 + z) ** 2
        + 4.0 * Omega_r * (1 + z) ** 3
        + Omega_L * df_idm_dz(z)
    )


def dH_dz(z: np.ndarray) -> np.ndarray:
    """Derivative of H(z) w.r.t. z."""
    return H0 * 0.5 * dE2_dz(z) / np.sqrt(E2(z))


def H_prime_over_H(z: np.ndarray) -> np.ndarray:
    """Logarithmic derivative: H'(z) / H(z)."""
    return 0.5 * dE2_dz(z) / E2(z)


# ─────────────────────────────────────────────
# Cosmological functions
# ─────────────────────────────────────────────

def Omega_m_z(z: np.ndarray) -> np.ndarray:
    """Matter density parameter at redshift z."""
    return Omega_m * (1 + z) ** 3 / E2(z)


def Omega_L_z(z: np.ndarray) -> np.ndarray:
    """Dark energy density parameter at redshift z (including sync dip)."""
    return Omega_L * f_idm(z) / E2(z)


# ─────────────────────────────────────────────
# Convenience
# ─────────────────────────────────────────────

def parameters() -> dict:
    """Return current IDM parameters as a dict."""
    return {
        "H0": H0,
        "Omega_m": Omega_m,
        "Omega_r": Omega_r,
        "Omega_L": Omega_L,
        "eps": eps,
        "zc": zc,
        "sigma": sigma,
    }
