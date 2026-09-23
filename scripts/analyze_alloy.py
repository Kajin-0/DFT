#!/usr/bin/env python3
"""Hg0.8Cd0.2Te SQS analysis: band-count edges from the SCF/nscf eigenvalue
listings of a run dir, DOS figure, comparison vs the empirical T=0 benchmark.

Band edges use n_valence = nelec (SOC) or nelec/2 (collinear) derived from
the run's own pw output. A negative or zero band-count gap is reported
honestly (no fabricated cutoff wavelength).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mct_dft.bands import band_edges, n_valence_from_electrons, cutoff_um_from_gap  # noqa
from mct_dft.parser import parse_bands_output  # noqa: E402
from mct_dft.detector import hansen_eg, cutoff_wavelength_um  # noqa: E402


def analyze_run(run_dir: Path) -> dict:
    res = json.loads((run_dir / "result.json").read_text())
    soc = bool(res["soc"])
    out_text = (run_dir / "pw.out").read_text(errors="replace")
    import re
    m = re.search(r"number of electrons\s*=\s*([\d.]+)", out_text)
    nelec = float(m.group(1)) if m else None
    nv = n_valence_from_electrons(nelec, soc) if nelec else None
    kpts, eigs = parse_bands_output(out_text)
    edges = None
    if eigs.size and nv and eigs.shape[1] > nv:
        e_fermi = res["pw"]["fermi_ev"]
        edges = band_edges(kpts, eigs, reference_ev=e_fermi, n_valence=nv)
    return {
        "run": run_dir.name, "soc": soc, "n_electrons": nelec,
        "n_valence": nv,
        "ok": res["pw"]["ok"], "E_tot_Ry": res["pw"]["total_energy_ry"],
        "E_F_eV": res["pw"]["fermi_ev"],
        "band_edges": (None if edges is None else {
            "vbm_eV": edges.vbm_ev, "cbm_eV": edges.cbm_ev,
            "gap_eV": edges.gap_ev, "status": edges.status,
            "cutoff_um": cutoff_um_from_gap(edges.gap_ev)}),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dirs", nargs="+", type=Path)
    ap.add_argument("--out-prefix", default="hgcdte_x020")
    args = ap.parse_args()

    rows = [analyze_run(ROOT / d) for d in args.run_dirs]
    print(json.dumps(rows, indent=2))

    # store
    out = ROOT / "results/tables" / f"{args.out_prefix}_scf_edges.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "empirical_T0_gap_eV": hansen_eg(0.2, 0.0),
        "empirical_77K_gap_eV": hansen_eg(0.2, 77.0),
        "empirical_300K_gap_eV": hansen_eg(0.2, 300.0),
        "empirical_cutoff_T0_um": cutoff_wavelength_um(hansen_eg(0.2, 0.0)),
        "note": ("Empirical rows are benchmarks (Hansen). The T=0K number "
                 "is the static-lattice comparison partner; PBE error on "
                 "narrow-gap HgCdTe can exceed the gap itself."),
        "runs": rows,
    }
    out.write_text(json.dumps(payload, indent=2))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
