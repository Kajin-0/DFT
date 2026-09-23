"""QE integration smoke test — marked `qe`, excluded from default pytest.

Runs a *tiny* real pw.x SCF (CdTe primitive, coarse parameters, 1 rank) and
checks the parser verdict chain end-to-end. Requires .local/qe-env.
Skipped automatically if QE isn't installed.
"""
import os
import shutil
from pathlib import Path

import pytest

pytestmark = pytest.mark.qe

ROOT = Path(__file__).resolve().parents[1]
QE = ROOT / ".local/qe-env/bin/pw.x"


@pytest.mark.skipif(not QE.exists(), reason="project-local QE not installed")
def test_tiny_scf_roundtrip():
    # run INSIDE the workspace (.tmp/) — containment guard bans /tmp
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from mct_dft.structures import zincblende_primitive
    from mct_dft.qe_inputs import pw_input, write_text
    from mct_dft.parser import parse_pw_output
    from mct_dft.runner import run_command

    tmp_path = ROOT / ".tmp" / "qe_integration"
    if tmp_path.exists():
        shutil.rmtree(tmp_path)
    tmp_path.mkdir(parents=True)

    man = __import__("json").loads(
        (ROOT / "pseudopotentials/manifest.json").read_text())
    pseudo = {e["element"]: Path(e["local_path"]).name
              for e in man["potentials"] if e["role"] == "nosoc"}

    atoms = zincblende_primitive("Cd", "Te", 6.482)
    txt = pw_input(atoms, pseudopotentials=pseudo, ecutwfc_ry=25.0,
                   ecutrho_ry=150.0, calculation="scf", kgrid=(2, 2, 2),
                   pseudo_dir=str(ROOT / "pseudopotentials"),
                   outdir=str(tmp_path / "out"), conv_thr=1e-8)
    write_text(tmp_path / "pw.in", txt, root=ROOT)
    r = run_command([str(QE), "-in", "pw.in"], tmp_path, tmp_path / "pw.out")
    assert r["returncode"] == 0
    res = parse_pw_output((tmp_path / "pw.out").read_text(errors="replace"))
    assert res.ok and res.job_done and res.total_energy_ry is not None
    assert res.nat == 2
