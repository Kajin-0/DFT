"""SOC pseudopotential validation on synthetic UPF headers."""
import textwrap

import pytest

from mct_dft.pseudopotentials import parse_upf_header, validate_soc_set


FR_HEADER = textwrap.dedent("""\
    <UPF version="2.0.1">
    <PP_INFO>
    test
    </PP_INFO>
    <PP_HEADER element="Hg" pseudo_type="US" relativistic="full"
     is_ultrasoft="T" is_paw="F" is_coulomb="F" has_so="T"
     has_wfc="F" has_gipaw="F" core_correction="T" functional="PBE"
     z_valence="2.000000000000000E+001" total_psenergy="-1.0E+002"
     wfc_cutoff="5.0E+001" rho_cutoff="4.0E+002" l_max="3"
     l_local="-1" mesh_size="99" number_of_wfc="7" number_of_proj="16"/>
    </UPF>
""")


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_parse_fully_relativistic(tmp_path):
    p = _write(tmp_path, "Hg.UPF", FR_HEADER)
    m = parse_upf_header(p)
    assert m.element == "Hg"
    assert m.relativistic == "full"
    assert m.has_so is True
    assert m.is_ultrasoft is True
    assert m.pseudo_type == "US"
    assert m.z_valence == pytest.approx(20.0)
    assert m.suggested_ecutwfc_ry == pytest.approx(50.0)
    assert m.suggested_ecutrho_ry == pytest.approx(400.0)
    assert m.soc_capable()


def test_soc_validation_accepts_fr_set(tmp_path):
    p = _write(tmp_path, "Hg.UPF", FR_HEADER)
    problems = validate_soc_set([parse_upf_header(p)], xc_expected="PBE")
    assert problems == []


SR_HEADER = FR_HEADER.replace('relativistic="full"', 'relativistic="scalar"') \
                     .replace('has_so="T"', 'has_so="F"')


def test_soc_validation_rejects_scalar(tmp_path):
    p = _write(tmp_path, "Cd.UPF", SR_HEADER)
    problems = validate_soc_set([parse_upf_header(p)])
    assert any("relativistic" in s for s in problems)
    assert any("has_so" in s for s in problems)


def test_soc_validation_rejects_wrong_functional(tmp_path):
    p = _write(tmp_path, "Cd.UPF", FR_HEADER.replace('functional="PBE"',
                                                     'functional="LDA"'))
    problems = validate_soc_set([parse_upf_header(p)], xc_expected="PBE")
    assert any("functional" in s for s in problems)
