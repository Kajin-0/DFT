"""Empirical HgCdTe benchmark: Hansen-style Eg(x,T) and cutoff wavelength."""
import math

import pytest

from mct_dft import detector


@pytest.mark.parametrize(
    "x,T,eg_expected,lam_expected",
    [
        (0.199, 0.0, 0.05654985, 21.9248),
        (0.199, 77.0, 0.08134924, 15.2410),
        (0.199, 300.0, 0.15317085, 8.09450),
        (0.200, 0.0, 0.058256, 21.28265),
        (0.200, 77.0, 0.082973, 14.94272),
        (0.200, 300.0, 0.154556, 8.02196),
    ],
)
def test_hansen_benchmark_values(x, T, eg_expected, lam_expected):
    eg = detector.hansen_eg(x, T)
    assert eg == pytest.approx(eg_expected, rel=1e-4, abs=2e-8)
    lam = detector.cutoff_wavelength_um(eg)
    assert lam == pytest.approx(lam_expected, rel=3e-4)


def test_cutoff_rejects_nonpositive_gap():
    assert detector.cutoff_wavelength_um(0.0) is None
    assert detector.cutoff_wavelength_um(-0.25) is None


def test_hansen_signs():
    # x=0 (HgTe): empirical gap is negative (semimetallic) at low T
    assert detector.hansen_eg(0.0, 0.0) < 0.0
    # positive temperature coefficient for x < 0.5
    # sign flips above x = 0.5: dEg/dT = 5.35e-4 * (1-2x)
    assert detector.hansen_eg(0.2, 300.0) > detector.hansen_eg(0.2, 0.0)
    assert detector.hansen_eg(0.7, 300.0) < detector.hansen_eg(0.7, 0.0)


def test_ideal_responsivity_shortcut():
    # R [A/W] = 0.80655 * eta * lambda[um]
    lam_um = 10.0
    eta = 0.6
    exact = detector.ideal_responsivity_a_per_w(eta, lam_um * 1e-6)
    assert exact == pytest.approx(0.80655 * eta * lam_um, rel=1e-4)
    for eta in (0.0, 1.0):
        r = detector.ideal_responsivity_from_um(eta, 5.0)
        assert math.isfinite(r)
