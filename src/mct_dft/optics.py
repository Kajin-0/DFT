"""Optical property conversions: dielectric function -> n, k, alpha, R, A.

All formulas in SI unless noted; wavelengths accepted in metres, with *_um
wrappers for micrometre convenience.
"""
from __future__ import annotations

import cmath
import math

import numpy as np


def dielectric_to_nk(eps1: float, eps2: float) -> tuple[float, float]:
    """Complex refractive index N = n + i kappa from eps = eps1 + i eps2.

    Uses a proper complex square root and fixes the sign so that kappa >= 0
    (passive medium convention).
    """
    n_complex = cmath.sqrt(complex(eps1, eps2))
    n, kappa = n_complex.real, n_complex.imag
    if kappa < 0:
        n, kappa = -n, -kappa
    return n, kappa


def nk_to_dielectric(n: float, kappa: float) -> tuple[float, float]:
    """eps1 = n^2 - k^2 ; eps2 = 2 n k  (consistency cross-check helper)."""
    return n * n - kappa * kappa, 2.0 * n * kappa


def absorption_coefficient_per_m(kappa: float, wavelength_m: float) -> float:
    """alpha = 4*pi*kappa / lambda   [1/m]."""
    if wavelength_m <= 0:
        raise ValueError("wavelength must be positive")
    return 4.0 * math.pi * kappa / wavelength_m


def normal_reflectance(n: float, kappa: float, n_ambient: float = 1.0) -> float:
    """Normal-incidence reflectance R = |(N - n0)/(N + n0)|^2."""
    z = complex(n, kappa)
    r = (z - n_ambient) / (z + n_ambient)
    return abs(r) ** 2


def absorptance_zero_refl(alpha_per_m: float, thickness_m: float) -> float:
    """Single-pass ideal absorber: A = 1 - exp(-alpha * d)."""
    return 1.0 - math.exp(-alpha_per_m * thickness_m)


def external_quantum_efficiency(
    n: float, kappa: float, wavelength_m: float, thickness_m: float
) -> float:
    """Ideal external QE proxy: eta = (1 - R) * [1 - exp(-alpha d)]."""
    alpha = absorption_coefficient_per_m(kappa, wavelength_m)
    return (1.0 - normal_reflectance(n, kappa)) * absorptance_zero_refl(
        alpha, thickness_m
    )


# numpy-vectorised variants ----------------------------------------------------
def dielectric_to_nk_vec(eps1, eps2):
    e = np.sqrt(np.asarray(eps1) + 1j * np.asarray(eps2))
    return e.real, np.abs(e.imag)


def absorption_crosscheck_um(wavelength_um: float, kappa: float) -> float:
    """alpha in 1/cm from lambda in um, handy tabulated form."""
    return absorption_coefficient_per_m(kappa, wavelength_um * 1e-6) / 100.0
