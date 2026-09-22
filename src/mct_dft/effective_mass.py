"""Effective-mass extraction from E(k) band data.

Near a band extremum:  E(k) = E0 + hbar^2 k^2 / (2 m*)
=> d^2E/dk^2 = hbar^2 / m*  =>  m* = hbar^2 / (d^2E/dk^2)  (SI throughout).

Public API works with k in 1/Angstrom and E in eV (plotting units); the
conversion to SI happens in exactly one place (:func:`_si_curvature`).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .constants import HBAR, M_E, EV_TO_J, ANG_TO_M


@dataclass
class ParabolicFit:
    e0_ev: float               # extremum energy
    curvature_ev_ang2: float   # d2E/dk2 in eV * Angstrom^2
    curvature_j_m2: float      # d2E/dk2 in J * m^2
    mass_kg: float             # |m*| in kg (sign carried by `sign`)
    mass_me: float             # |m*| in electron masses
    sign: int                  # +1 electron-like (min), -1 hole-like (max)
    r_squared: float           # quality of quadratic fit
    k_center_inv_ang: float
    window_inv_ang: float      # half-width of fitting window around k_center
    n_points: int
    parabolic_ok: bool         # heuristic validity flag
    notes: list = field(default_factory=list)


def _si_curvature(curvature_ev_ang2: float) -> float:
    """eV*Angstrom^2 -> J*m^2."""
    return curvature_ev_ang2 * EV_TO_J * ANG_TO_M**2


def mass_from_curvature(curvature_ev_ang2: float) -> tuple[float, float]:
    """Return (|m*| in kg, |m*| in m_e) from d2E/dk2 in eV*Angstrom^2."""
    curv_si = abs(_si_curvature(curvature_ev_ang2))
    if curv_si == 0.0:
        return math.inf, math.inf
    m_kg = HBAR**2 / curv_si
    return m_kg, m_kg / M_E


def fit_parabola(
    k_inv_ang: np.ndarray,
    e_ev: np.ndarray,
    k_center: float | None = None,
    window: float | None = None,
    r2_warn: float = 0.9,
) -> ParabolicFit:
    """Least-squares quadratic fit E = a k^2 + b k + c over a k-window.

    If k_center/window are None, the extremum is found from the data and a
    window covering all passed points is used.
    """
    k = np.asarray(k_inv_ang, dtype=float)
    e = np.asarray(e_ev, dtype=float)
    if k.size < 5:
        raise ValueError("need at least 5 (k, E) points for a stable fit")
    order = np.argsort(k)
    k, e = k[order], e[order]

    if k_center is None:
        k_center = float(k[np.argmin(np.abs(e - np.median(e)))])
        # better: use index of local extremum of the raw data
        i_ext = int(np.argmin(e)) if e.std() else 0
        k_center = float(k[i_ext])
    if window is None:
        window = float(0.5 * (k.max() - k.min()))

    mask = np.abs(k - k_center) <= window
    if mask.sum() < 5:
        raise ValueError(
            f"fewer than 5 points inside window +/-{window} around {k_center}"
        )
    kk, ee = k[mask], e[mask]

    a, b, c = np.polyfit(kk, ee, 2)  # E = a k^2 + b k + c
    curvature = 2.0 * a              # eV / Angstrom^2
    e0 = c - b * b / (4.0 * a) if a != 0 else float("nan")

    # R^2 of the quadratic model
    e_pred = a * kk**2 + b * kk + c
    ss_res = float(np.sum((ee - e_pred) ** 2))
    ss_tot = float(np.sum((ee - ee.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    m_kg, m_me = mass_from_curvature(curvature)
    sign = 1 if a > 0 else -1  # minimum vs maximum
    notes = []
    if r2 < r2_warn:
        notes.append(
            f"low R^2={r2:.4f}: band is not parabolic over this window; "
            "mass is indicative only"
        )
    return ParabolicFit(
        e0_ev=float(e0),
        curvature_ev_ang2=float(curvature),
        curvature_j_m2=float(_si_curvature(curvature)),
        mass_kg=float(m_kg),
        mass_me=float(m_me),
        sign=sign,
        r_squared=r2,
        k_center_inv_ang=k_center,
        window_inv_ang=float(window),
        n_points=int(mask.sum()),
        parabolic_ok=r2 >= r2_warn,
        notes=notes,
    )


def window_scan(
    k_inv_ang,
    e_ev,
    windows=(0.01, 0.02, 0.03),
    k_center: float | None = None,
) -> list[ParabolicFit]:
    """Fit over several window half-widths to expose curvature sensitivity."""
    fits = []
    for w in windows:
        try:
            fits.append(
                fit_parabola(k_inv_ang, e_ev, k_center=k_center, window=float(w))
            )
        except ValueError:
            continue
    return fits
