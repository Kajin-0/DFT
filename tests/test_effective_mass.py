"""Effective-mass fitting: synthetic parabolas must round-trip m*."""
import numpy as np
import pytest

from mct_dft.constants import HBAR, M_E, EV_TO_J, ANG_TO_M
from mct_dft.effective_mass import fit_parabola, window_scan, mass_from_curvature


def _synthetic_parabola(m_over_me, e0=0.0, nk=101, kmax=0.05):
    """Exact free-electron-like parabola in (1/Angstrom, eV) units."""
    k_ang = np.linspace(-kmax, kmax, nk)
    k_m = k_ang / ANG_TO_M                        # 1/m
    e_j = HBAR**2 * k_m**2 / (2.0 * m_over_me * M_E)
    e_ev = e0 + e_j / EV_TO_J
    return k_ang, e_ev


@pytest.mark.parametrize("m_over_me", [0.05, 0.1, 0.5, 1.0])
def test_parabola_recovers_mass(m_over_me):
    k, e = _synthetic_parabola(m_over_me)
    fit = fit_parabola(k, e, k_center=0.0, window=0.04)
    assert fit.mass_me == pytest.approx(m_over_me, rel=2e-3)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-12)
    assert fit.sign == +1 and fit.parabolic_ok


def test_hole_sign():
    k, e = _synthetic_parabola(0.2)
    fit = fit_parabola(k, -e, k_center=0.0, window=0.04)  # band maximum
    assert fit.sign == -1
    assert fit.mass_me == pytest.approx(0.2, rel=2e-3)


def test_window_scan_consistency():
    k, e = _synthetic_parabola(0.15, nk=201)
    fits = window_scan(k, e, windows=(0.01, 0.02, 0.03), k_center=0.0)
    assert len(fits) == 3
    for f in fits:
        assert f.mass_me == pytest.approx(0.15, rel=5e-3)
        assert f.n_points >= 5


def test_nonparabolic_flagged():
    k = np.linspace(-0.05, 0.05, 101)
    e = 50.0 * k**4  # quartic: quadratic fit must be flagged bad
    fit = fit_parabola(k, e, k_center=0.0, window=0.04)
    assert not fit.parabolic_ok or fit.r_squared < 0.999


def test_mass_units_si():
    # curvature for m* = m_e: d2E/dk2 = hbar^2/m_e in SI
    curv_ev_ang2 = (HBAR**2 / M_E) / (EV_TO_J * ANG_TO_M**2)
    m_kg, m_me = mass_from_curvature(curv_ev_ang2)
    assert m_kg == pytest.approx(M_E, rel=1e-9)
    assert m_me == pytest.approx(1.0, rel=1e-9)
