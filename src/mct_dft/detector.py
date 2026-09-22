"""HgCdTe empirical bandgap relation and ideal detector figures of merit.

The empirical relation implemented here is the Hansen / Schmit / Casselman
polynomial (as specified for this project):

    Eg(x, T) = -0.302 + 1.93 x - 0.81 x^2 + 0.832 x^3
               + 5.35e-4 (1 - 2x) T          [eV], T in kelvin

This is an EMPIRICAL BENCHMARK, not a DFT result. It must never be presented
as first-principles output.
"""
from __future__ import annotations

import math

from .constants import HC_EV_UM, C_LIGHT, ELEMENTARY_CHARGE, PLANCK


def hansen_eg(x: float, temperature: float = 0.0) -> float:
    """Empirical Hg(1-x)Cd(x)Te band gap in eV.

    Parameters
    ----------
    x : Cd mole fraction (0 <= x <= 1)
    temperature : lattice temperature in kelvin

    Returns
    -------
    Gap in eV. May be <= 0 (semimetallic / inverted regime); callers must
    treat non-positive gaps as "no normal semiconductor gap".
    """
    eg = (
        -0.302
        + 1.93 * x
        - 0.81 * x**2
        + 0.832 * x**3
        + 5.35e-4 * (1.0 - 2.0 * x) * temperature
    )
    return eg


def cutoff_wavelength_um(eg_ev: float) -> float | None:
    """Cutoff wavelength in micrometres for a *positive* gap.

    Returns None when eg_ev <= 0: a cutoff wavelength is undefined for a
    zero/negative-gap (semimetallic or band-inverted) system.
    """
    if eg_ev <= 0.0:
        return None
    return HC_EV_UM / eg_ev


def ideal_responsivity_a_per_w(eta: float, wavelength_m: float) -> float:
    """Ideal photovoltaic responsivity, R = eta * q * lambda / (h c)  [A/W].

    eta is the external quantum efficiency (0..1), wavelength_m in metres.
    Assumes perfect collection of photogenerated carriers and unity gain.
    """
    return eta * ELEMENTARY_CHARGE * wavelength_m / (PLANCK * C_LIGHT)


def ideal_responsivity_from_um(eta: float, wavelength_um: float) -> float:
    """Same as :func:`ideal_responsivity_a_per_w` with lambda in micrometres.

    Equals 0.80655 * eta * lambda[um]  (A/W).
    """
    return ideal_responsivity_a_per_w(eta, wavelength_um * 1e-6)


def empirical_series(x: float, temperatures=(0.0, 77.0, 300.0)) -> list[dict]:
    """Convenience table of the empirical benchmark at several temperatures."""
    rows = []
    for t in temperatures:
        eg = hansen_eg(x, t)
        rows.append(
            {
                "x": x,
                "T_K": t,
                "Eg_eV": eg,
                "lambda_c_um": cutoff_wavelength_um(eg),
                "note": "EMPIRICAL INPUT (Hansen et al.), not a DFT result",
            }
        )
    return rows


if __name__ == "__main__":  # quick manual check
    for x in (0.199, 0.200):
        for row in empirical_series(x):
            lam = row["lambda_c_um"]
            print(
                f"x={row['x']:.3f} T={row['T_K']:6.1f}  "
                f"Eg={row['Eg_eV']:.8f} eV  lambda_c="
                f"{lam if lam is None else f'{lam:.5f} um'}"
            )
