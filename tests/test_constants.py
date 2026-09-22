"""Unit conversion checks — explicit units are a hard project requirement."""
import pytest

from mct_dft import constants as C


def test_ev_joule_roundtrip():
    e_ev = 1.602176634e-19 / C.EV_TO_J
    assert C.ev_to_joule(1.0) == pytest.approx(1.602176634e-19, rel=1e-12)
    assert C.joule_to_ev(C.ev_to_joule(3.7)) == pytest.approx(3.7, rel=1e-15)


def test_ry_ev():
    assert C.RY_TO_EV == pytest.approx(13.605693122, rel=1e-6)


def test_inverse_length_units():
    k = 0.02  # 1/Angstrom
    assert C.inv_ang_to_inv_m(k) == pytest.approx(2e8)   # 1/m
    assert C.inv_m_to_inv_ang(C.inv_ang_to_inv_m(k)) == pytest.approx(k)


def test_hc_ev_um():
    assert C.HC_EV_UM == pytest.approx(1.239841984, rel=1e-9)
