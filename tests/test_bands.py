"""band_edges must reflect computed ordering honestly (no fake gaps)."""
import numpy as np
import pytest

from mct_dft.bands import band_edges, cutoff_um_from_gap


def _grid(n=21):
    # 1D toy: k along one axis, two parabolas
    k = np.linspace(-0.5, 0.5, n)
    kpts = np.column_stack([k, np.zeros(n), np.zeros(n)])
    return kpts, k


def test_direct_gap_detected():
    kpts, k = _grid()
    ev = np.column_stack([-k**2 * 20, 0.7 + k**2 * 40])  # vb, cb
    ed = band_edges(kpts, ev, reference_ev=0.3)
    assert ed.gap_ev == pytest.approx(0.7)
    assert ed.direct and "direct" in ed.status
    assert cutoff_um_from_gap(ed.gap_ev) == pytest.approx(1.239841984 / 0.7)
    assert np.allclose(ed.k_vbm, 0.0)


def test_inverted_overlap_reported_not_forced():
    kpts, k = _grid()
    # cb dips below vb at centre: indirect overlap -> negative 'gap'
    ev = np.column_stack([-k**2 * 5, -0.1 + k**2 * 10])
    # band-count definition: n_valence=1 -> VBM from band 0, CBM from band 1
    ed = band_edges(kpts, ev, reference_ev=0.05, n_valence=1)
    assert ed.gap_ev == pytest.approx(-0.1)
    assert "inverted" in ed.status or "overlap" in ed.status
    assert cutoff_um_from_gap(ed.gap_ev) is None


def test_bandcount_direct_gap():
    kpts, k = _grid()
    ev = np.column_stack([-k**2 * 20, 0.7 + k**2 * 40])
    ed = band_edges(kpts, ev, n_valence=1)
    assert ed.gap_ev == pytest.approx(0.7)
    assert ed.direct


def test_metallic_flag():
    kpts, k = _grid()
    ev = np.column_stack([np.full(k.shape, -0.1)])  # all below reference
    ed = band_edges(kpts, ev, reference_ev=0.0)
    assert ed.metallic
