#!/usr/bin/env python3
"""Publication-quality figures from committed/available run outputs.

Subcommands
-----------
bands    : band structure from a run dir containing bands.npz + seekpath.json
dos      : total DOS from dos output .dat files in a run dir
socc     : SOC vs no-SOC comparison overlay for a system
dosedge  : DOS near the gap (band-edge region)
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

from mct_dft.bands import band_edges, n_valence_from_electrons  # noqa: E402
from mct_dft.parser import parse_pw_output  # noqa: E402


def load_run(run_dir: Path) -> dict:
    return json.loads((run_dir / "result.json").read_text())


def n_valence_for(run_dir: Path) -> int | None:
    """valence-band count from the SCF run this band/dos run drew charge from"""
    res = load_run(run_dir)
    cf = res.get("charge_from")
    if not cf:
        # convention fallback for runs pre-dating charge_from recording
        sysname = res.get("system")
        tag = "soc" if res.get("soc") else "nosoc"
        guess = ROOT / "calculations" / str(sysname) / f"production_scf_{tag}"
        cf = str(guess.relative_to(ROOT)) if guess.exists() else None
    if not cf:
        return None
    scf_out = (ROOT / cf / "pw.out")
    if not scf_out.exists():
        return None
    # cheap targeted search instead of full parse
    import re
    t = scf_out.read_text(errors="replace")
    m = re.search(r"number of electrons\s*=\s*([\d.]+)", t)
    if not m:
        return None
    return n_valence_from_electrons(float(m.group(1)),
                                    bool(res.get("soc")))


def label_axis(ax, seekpath: dict, kpts: np.ndarray) -> np.ndarray:
    """x-axis = cumulative |dk|; vertical lines + labels at path vertices."""
    sp_kpts = np.asarray(seekpath["explicit_kpoints_rel"])
    labels = {tuple(np.round(k, 6)): lab for k, lab in
              zip(sp_kpts, seekpath["explicit_kpoints_labels"])}
    # cumulative distance of our fine k-mesh
    d = np.linalg.norm(np.diff(kpts, axis=0), axis=1)
    x = np.concatenate([[0.0], np.cumsum(d)])
    ticks, names = [], []
    for i, k in enumerate(kpts):
        lab = labels.get(tuple(np.round(k, 6)))
        if lab is not None and (not names or names[-1] != lab):
            ticks.append(x[i]); names.append(lab)
    ax.set_xticks(ticks)
    ax.set_xticklabels([n.replace("GAMMA", r"$\Gamma$") for n in names])
    for t in ticks:
        ax.axvline(t, color="0.6", lw=0.5, zorder=0)
    ax.set_xlim(x[0], x[-1])
    return x


def fig_bands(run_dir: Path, out: Path, title: str) -> dict:
    z = np.load(run_dir / "bands.npz")
    kpts, eigs = z["kpts"], z["eigs_ev"]  # eV eigenvalues, absolute
    res = load_run(run_dir)
    ef = res["pw"].get("fermi_ev") or res.get("charge_reference_occupancy_ev")
    nv = n_valence_for(run_dir)
    edges = band_edges(kpts, eigs, reference_ev=ef, n_valence=nv)
    ref_ev = edges.vbm_ev if not edges.metallic else (ef or edges.vbm_ev)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    sp = json.loads((run_dir / "seekpath.json").read_text())
    x = label_axis(ax, sp, kpts)
    for ib in range(eigs.shape[1]):
        ax.plot(x, eigs[:, ib] - ref_ev, color="#0b5394", lw=0.8)
    ax.axhline(0.0, color="k", lw=0.6, ls="--")
    ax.set_ylabel("E - VBM (eV)")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return {
        "figure": str(out),
        "vbm_eV": edges.vbm_ev, "cbm_eV": edges.cbm_ev,
        "gap_eV": edges.gap_ev, "direct": edges.direct,
        "status": edges.status,
        "n_valence": nv,
        "k_vbm": edges.k_vbm.tolist(), "k_cbm": edges.k_cbm.tolist(),
        "reference": "VBM" if not edges.metallic else "E_F",
    }


def read_dos_dat(path: Path):
    arr = np.loadtxt(path, comments="#")
    return arr[:, 0], arr[:, 1]  # E (eV), DOS (states/eV)


def fig_dos(run_dir: Path, out: Path, title: str, emin=-8, emax=6,
            highlight_edge=True) -> dict:
    # find the dos file QE wrote (name from our dos.in: prefix dependent)
    cands = sorted(run_dir.glob("*.dat")) + sorted((run_dir / "out").glob(
        "*.dos")) if (run_dir / "out").exists() else sorted(run_dir.glob("*.dat"))
    dos_file = None
    for c in cands:
        if "dos" in c.name and "pdos" not in c.name:
            dos_file = c
            break
    if dos_file is None:
        raise FileNotFoundError(f"no DOS file found in {run_dir}")
    e, d = read_dos_dat(dos_file)
    res = load_run(run_dir)
    fermi = res["pw"]["fermi_ev"] or res["pw"]["highest_occupied_ev"] or 0.0
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    m = (e - fermi > emin) & (e - fermi < emax)
    ax.plot(d[m], e[m] - fermi, color="#38761d", lw=1.0)
    ax.fill_betweenx(e[m] - fermi, 0, d[m], color="#38761d", alpha=0.25)
    ax.axhline(0, color="k", lw=0.6, ls="--")
    ax.set_xlabel("DOS (states / eV / cell)")
    ax.set_ylabel("E - E_ref (eV)")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return {"figure": str(out), "dos_file": str(dos_file),
            "reference_eV": fermi}


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("bands")
    p.add_argument("run_dir", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--title", default="")

    p = sub.add_parser("dos")
    p.add_argument("run_dir", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--title", default="")
    p.add_argument("--emin", type=float, default=-8.0)
    p.add_argument("--emax", type=float, default=6.0)

    args = ap.parse_args()
    if args.cmd == "bands":
        print(json.dumps(fig_bands(args.run_dir, args.out, args.title), indent=2))
    else:
        print(json.dumps(fig_dos(args.run_dir, args.out, args.title,
                                 args.emin, args.emax), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
