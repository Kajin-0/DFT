"""pw.x output parser checks on stored text fixtures (no QE required)."""
import numpy as np
import pytest

from mct_dft.parser import parse_pw_output, parse_bands_output

SCF_OK = """
     Program PWSCF v.7.5 starts on 22Sep2026
     calculation = 'scf'
     Lattice parameter (alat)  = 12.2457  a.u.
     unit-cell volume          =    251.1230 (a.u.)^3
     number of atoms           = 2
     number of Kohn-Sham states=   20
     number of k points=    29
     With spin-orbit
     noncollinear calculation
     total cpu time spent up to now is        1.2 secs

     total energy              =     -63.12345678 Ry
     Harris-Foulkes estimate   =     -63.12345679 Ry
     estimated scf accuracy    <          2.0E-11 Ry

!    total energy              =     -63.12345678 Ry
     Harris-Foulkes estimate   =     -63.12345679 Ry
     estimated scf accuracy    <          1.0E-13 Ry

     convergence has been achieved in   9 iterations

     highest occupied level (ev):     8.2341

     JOB DONE.
"""

SCF_FAILED = SCF_OK.replace(
    "convergence has been achieved in   9 iterations",
    "convergence NOT achieved after 200 iterations: stopping",
).replace("JOB DONE.", "JOB ABORTED.")

BANDS_TXT = """
     Program PWSCF v.7.5 starts

          k =(   0.000000   0.000000   0.000000  ) ,  P = 0.0

     bands (ev):

    -5.3535  -5.3534   0.0000   0.0000   0.0001   0.0001   1.6020   1.6021

          k =(   0.100000   0.000000   0.000000  ) ,  P = 0.0

     bands (ev):

    -5.3000  -5.2999  -0.1000  -0.0999   0.0500   0.0501   1.5000   1.5001

     JOB DONE.
"""


def test_scf_success_parsed():
    r = parse_pw_output(SCF_OK)
    assert r.ok and r.converged and r.job_done
    assert r.total_energy_ry == pytest.approx(-63.12345678)
    assert r.total_energy_ev == pytest.approx(-63.12345678 * 13.605693122994)
    assert r.nscf_iterations == 9
    assert r.nat == 2 and r.n_bands == 20 and r.n_kpts == 29
    assert r.noncollinear and r.spin_orbit
    assert r.errors == []
    assert r.highest_occupied_ev == pytest.approx(8.2341)


def test_failed_scf_rejected():
    r = parse_pw_output(SCF_FAILED)
    assert not r.ok
    assert not r.converged
    assert any("NOT achieved" in e for e in r.errors)


def test_empty_output_rejected():
    r = parse_pw_output("")
    assert not r.ok and r.errors


NSCF_TXT = """
          k = 0.0442 0.0442 0.0442 (  1959 PWs)   bands (ev):

    -6.0318  -3.1379  -3.1379  -3.1298  -2.8952  -2.8952   4.6557   5.0930
     5.0930   6.2969   9.7106   9.8006   9.8006

     occupation numbers
     1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000
     1.0000   0.0000   0.0000   0.0000   0.0000

          k = 0.1326 0.1326-0.0442 (  1963 PWs)   bands (ev):

    -5.9568  -3.1263  -3.1190  -3.0703  -2.8988  -2.8925   3.5678   4.7155
     5.0070   7.0574   9.6357   9.9555  10.3033

     occupation numbers
     1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000
     1.0000   0.0000   0.0000   0.0000   0.0000

     the Fermi energy is     5.4414 ev

     JOB DONE.
"""


def test_bands_parsing():
    kpts, eigs = parse_bands_output(BANDS_TXT)
    assert kpts.shape == (2, 3)
    assert eigs.shape == (2, 8)
    assert np.allclose(kpts[1], [0.1, 0.0, 0.0])
    assert eigs[0, 0] == pytest.approx(-5.3535)
    assert eigs[1, -1] == pytest.approx(1.5001)


def test_nscf_bands_parsing():
    """NSCF listing: 'k = x y z ( N PWs) bands (ev):' on one line, possibly
    abutting floats, followed by an occupation-numbers block that must not
    be swallowed."""
    kpts, eigs = parse_bands_output(NSCF_TXT)
    assert kpts.shape == (2, 3)
    assert eigs.shape == (2, 13)
    assert np.allclose(kpts[1], [0.1326, 0.1326, -0.0442])
    assert eigs[0, 0] == pytest.approx(-6.0318)
    # occupation numbers were excluded
    assert eigs.max() > 10.3 - 1e-9
    assert not (np.isclose(eigs, 1.0, atol=1e-6)).any()
