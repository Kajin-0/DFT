#!/usr/bin/env python3
"""Aggregate production results into results/tables/material_summary.csv.

Every row is assembled ONLY from committed run artefacts (result.json,
bands.npz). Missing pieces stay empty — never filled with placeholders.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mct_dft.bands import band_edges, cutoff_um_from_gap, \
    n_valence_from_electrons  # noqa: E402


def _nval(scf_dir: Path, soc: bool) -> int | None:
    import re
    out = scf_dir / "scf_placeholder"
    p = scf_dir / "pw.out"
    if not p.exists():
        return None
    m = re.search(r"number of electrons\s*=\s*([\d.]+)",
                  p.read_text(errors="replace"))
    if not m:
        return None
    return n_valence_from_electrons(float(m.group(1)), soc)


def load_json(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


def gamma_splitting(run_dir: Path, n_valence: int | None, soc: bool) -> dict:
    """Gamma-point spin-orbit splitting: mean(Gamma8 quartet) - mean(Gamma7
    doublet) within the valence manifold. Only meaningful with SOC on; the
    no-SOC crystal-field splitting is a different quantity and is NOT
    reported as delta_so.
    """
    if not soc or n_valence is None:
        return {"delta_so_ev": None,
                "note_dso": "needs SOC bands + known valence band count"}
    z = np.load(run_dir / "bands.npz")
    kpts, eigs = z["kpts"], z["eigs_ev"]
    ig = int(np.argmin(np.linalg.norm(kpts, axis=1)))
    valence = np.sort(eigs[ig][:n_valence])  # first n_valence spinor bands
    if len(valence) < 6:
        return {"delta_so_ev": None}
    dso = float(np.mean(valence[-4:]) - np.mean(valence[-6:-4]))
    return {"delta_so_ev": dso,
            "k_gamma_index": ig,
            "vbm_gamma_eV": float(valence[-1]),
            "splitoff_gamma_eV": float(np.mean(valence[-6:-4]))}


def analyze_system(system: str) -> list[dict]:
    base = ROOT / "calculations" / system
    rows = []
    for soc in (False, True):
        tag = "soc" if soc else "nosoc"
        row = {
            "system": system,
            "soc": soc,
            "composition_x": {"cdte": 1.0, "hgte": 0.0}.get(system),
            "structure": "zincblende",
            "xc": "PBE",
        }
        scf = load_json(base / f"production_scf_{tag}" / "result.json")
        relax = load_json(base / f"relax_{tag}" / "result.json")
        bandsd = base / f"production_bands_{tag}"
        dos = load_json(base / f"production_dos_{tag}" / "result.json")

        if relax and relax["pw"]["ok"]:
            cell = np.asarray(relax["pw"]["final_cell_ang"], float)
            row["lattice_parameter_A"] = float(np.linalg.norm(cell[0])
                                               * np.sqrt(2.0))
            row["atom_count"] = relax["structure_check"]["n_atoms"]
            row["scf_converged"] = relax["pw"]["converged"]
        if scf and scf["pw"]["ok"]:
            row["total_energy_Ry"] = scf["pw"]["total_energy_ry"]
            row["ecutwfc_Ry"] = scf["ecutwfc_ry"]
            row["ecutrho_Ry"] = scf["ecutrho_ry"]
            row["kgrid"] = "x".join(map(str, scf["kgrid"]))
            row["scf_converged"] = scf["pw"]["converged"]
            row["input_sha256"] = json.loads(
                (base / f"production_scf_{tag}" / "provenance.json"
                 ).read_text())["input_sha256"].get("pw.in")
        if (bandsd / "bands.npz").exists():
            zb = np.load(bandsd / "bands.npz")
            bres = load_json(bandsd / "result.json")
            ef = (bres["pw"]["highest_occupied_ev"]
                  or bres.get("charge_reference_occupancy_ev"))
            nv = _nval(base / f"production_scf_{tag}", soc)
            ed = band_edges(zb["kpts"], zb["eigs_ev"], reference_ev=ef,
                            n_valence=nv)
            row.update({
                "vbm_eV": ed.vbm_ev, "cbm_eV": ed.cbm_ev,
                "gap_eV": ed.gap_ev,
                "gap_type": ("direct" if ed.direct else "indirect"),
                "gap_status": ed.status,
                "metallic": ed.metallic,
                "k_vbm": np.array2string(ed.k_vbm, precision=4,
                                         separator=","),
                "k_cbm": np.array2string(ed.k_cbm, precision=4,
                                         separator=","),
                "cutoff_um": cutoff_um_from_gap(ed.gap_ev),
                "method": "PBE" + ("+SOC" if soc else ""),
            })
            row.update(gamma_splitting(bandsd, nv, soc))
        if dos and dos["pw"]["ok"]:
            row["dos_gap_ev"] = dos["pw"].get("gap_ev")
        rows.append(row)
    return rows


FIELDS = ["system", "composition_x", "structure", "method", "xc", "soc",
          "atom_count", "lattice_parameter_A", "ecutwfc_Ry", "ecutrho_Ry",
          "kgrid", "scf_converged", "total_energy_Ry", "vbm_eV", "cbm_eV",
          "gap_eV", "gap_type", "gap_status", "metallic", "k_vbm", "k_cbm",
          "cutoff_um", "delta_so_ev", "dos_gap_ev", "input_sha256",
          "electron_mass_me", "hole_mass_me", "notes"]


def main() -> int:
    rows = []
    for system in ("cdte", "hgte"):
        rows += analyze_system(system)
    out = ROOT / "results/tables/material_summary.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(out)
    for r in rows:
        print(json.dumps(r, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
