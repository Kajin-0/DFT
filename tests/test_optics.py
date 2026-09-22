"""Optical conversions round-trip and physical sanity checks."""
import math

import numpy as np
import pytest

from mct_dft import optics as O


def test_dielectric_nk_roundtrip():
    for n, k in [(3.5, 0.0), (4.0, 0.5), (2.5, 2.0), (1.5, 0.01)]:
        e1, e2 = O.nk_to_dielectric(n, k)
        n2, k2 = O.dielectric_to_nk(e1, e2)
        assert n2 == pytest.approx(n, rel=1e-12)
        assert k2 == pytest.approx(k, rel=1e-12)


def test_absorption_positive():
    a = O.absorption_coefficient_per_m(0.3, 5e-6)
    assert a > 0
    # 4 pi k / lambda
    assert a == pytest.approx(4 * math.pi * 0.3 / 5e-6)


def test_reflectance_known_limits():
    # lossless n: R = ((n-1)/(n+1))^2 ; n=4 -> 0.36
    assert O.normal_reflectance(4.0, 0.0) == pytest.approx(0.36)
    # impedance matched / vacuum: R -> 0
    assert O.normal_reflectance(1.0, 0.0) == pytest.approx(0.0)


def test_absorptance_bounds():
    assert O.absorptance_zero_refl(0.0, 1e-5) == 0.0
    assert O.absorptance_zero_refl(1e7, 1e-5) == pytest.approx(1.0)


def test_external_qe_bounds():
    eta = O.external_quantum_efficiency(3.5, 1.0, 10e-6, 10e-6)
    assert 0.0 <= eta <= 1.0


def test_vec_api():
    # build eps from known (n,k) so the vectorised inversion is exact
    e1, e2 = O.nk_to_dielectric(3.5, 0.4)
    n, k = O.dielectric_to_nk_vec(np.array([e1, 10.0]), np.array([e2, 0.0]))
    assert n.shape == (2,)
    assert n[0] == pytest.approx(3.5) and k[0] == pytest.approx(0.4)
    assert n[1] == pytest.approx(np.sqrt(10.0)) and k[1] == pytest.approx(0.0)
