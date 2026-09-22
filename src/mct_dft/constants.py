"""Physical constants and unit conversions.

All constants come from scipy.constants (CODATA). Every numerical quantity in
the rest of the package is in SI or in the explicitly documented unit of the
function; conversions happen only here.
"""
from __future__ import annotations

from scipy.constants import (
    c as C_LIGHT,          # m/s, speed of light in vacuum
    e as ELEMENTARY_CHARGE,  # C, elementary charge (= J per eV)
    electron_mass as M_E,   # kg
    h as PLANCK,            # J s
    hbar as HBAR,           # J s
)

# Convenient derived constants ------------------------------------------------
EV_TO_J = ELEMENTARY_CHARGE          # 1 eV in joules
RY_TO_EV = 13.605693122994           # 1 Ry in eV (CODATA Rydberg * hc / e)
BOHR_TO_ANG = 0.529177210903         # 1 Bohr radius in Angstrom
ANG_TO_M = 1e-10                     # 1 Angstrom in metres
HC_EV_UM = 1.239841984               # h*c in eV*um (used for lambda_c)

__all__ = [
    "C_LIGHT",
    "ELEMENTARY_CHARGE",
    "M_E",
    "PLANCK",
    "HBAR",
    "EV_TO_J",
    "RY_TO_EV",
    "BOHR_TO_ANG",
    "ANG_TO_M",
    "HC_EV_UM",
    "ev_to_joule",
    "joule_to_ev",
    "inv_ang_to_inv_m",
    "inv_m_to_inv_ang",
]


def ev_to_joule(e_ev: float) -> float:
    """Convert electron-volts to joules."""
    return e_ev * EV_TO_J


def joule_to_ev(e_j: float) -> float:
    """Convert joules to electron-volts."""
    return e_j / EV_TO_J


def inv_ang_to_inv_m(k_inv_ang: float) -> float:
    """Convert inverse Angstrom to inverse metres."""
    return k_inv_ang / ANG_TO_M


def inv_m_to_inv_ang(k_inv_m: float) -> float:
    """Convert inverse metres to inverse Angstrom."""
    return k_inv_m * ANG_TO_M
