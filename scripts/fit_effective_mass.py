#!/usr/bin/env python3
"""Directional effective-mass extraction near a band extremum.

Runs a dedicated dense NSCF along short segments around Gamma in several
crystallographic directions ([100],[110],[111]) and fits parabolas at a few
window half-widths. Writes a JSON + a diagnostic plot.
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

from mct_dft.effective_mass import window_scan  # noqa: E402
from mct_dft.parser import parse_bands_output   # noqa: E402

DIRECTIONS = {
    "[100]": [1.0, 0.0, 0.0],
    "[110]": [1.0, 1.0, 0.0],
    "[111]": [1.0, 1.0, 1.0],
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("system", choices=["cdte", "hgte"])
    ap.add_argument("--band", choices=["cb", "vb"], default="cb",
                    help="conduction band (minimum at G) or valence band")
    ap.add_argument("--maxk", type=float, default=0.05,
                    help="segment length in crystal reciprocal units")
    ap.add_argument("--npts", type=int, default=41)
    ap.add_argument("--ecut", type=float, default=120.0)
    ap.add_argument("--ecutrho", type=float, default=720.0)
    ap.add_argument("--lattice", type=float, default=None,
                    help="lattice constant (default: config lattice_start_ang)")
    ap.add_argument("--charge-from", required=True,
                    help="repo-relative dir of the SCF run with charge density")
    ap.add_argument("--np", type=int, default=6)
    ap.add_argument("--soc", action="store_true")
    args = ap.parse_args()

    import subprocess
    which = "soc" if args.soc else "nosoc"
    tag = f"mass_{which}_{args.band}_e{args.ecut:g}"
    run_dir = ROOT / "calculations" / args.system / tag
    run_dir.mkdir(parents=True, exist_ok=True)

    n = args.npts
    # reciprocal lattice for |k| in 1/Angstrom
    import yaml as _yaml
    from mct_dft.structures import zincblende_primitive
    cfg = _yaml.safe_load(
        (ROOT / "config/systems" / f"{args.system}.yaml").read_text())
    sp = cfg["species"]
    lat = args.lattice or cfg["lattice_start_ang"]
    atoms = zincblende_primitive(sp[0], sp[1], lat)
    recip = atoms.cell.reciprocal()  # 1/Angstrom (crystal coords -> 1/A)

    # single pw.x bands run with all segments concatenated
    cat = np.vstack([np.hstack([np.outer(
        np.linspace(-args.maxk, args.maxk, n),
        np.asarray(DIRECTIONS[d], float)
        / np.linalg.norm(np.asarray(DIRECTIONS[d], float))),
        np.ones((n, 1))]) for d in DIRECTIONS])

    from mct_dft.qe_inputs import pw_input, write_text
    import json as _json
    ppath = ROOT / "pseudopotentials" / "manifest.json"
    manifest = _json.loads(ppath.read_text())
    role = "soc" if args.soc else "nosoc"
    pseudo = {e["element"]: Path(e["local_path"]).name
              for e in manifest["potentials"] if e["role"] == role}

    text = pw_input(atoms, pseudopotentials=pseudo, ecutwfc_ry=args.ecut,
                    ecutrho_ry=args.ecutrho, calculation="bands",
                    kpoints_explicit=cat, soc=args.soc,
                    pseudo_dir=str(ROOT / "pseudopotentials"),
                    outdir=str(run_dir / "out"))
    write_text(run_dir / "pw_mass.in", text)

    # reuse charge density
    import shutil
    src = (ROOT / args.charge_from / "out")
    if not src.exists():
        print(f"charge density missing: {src}")
        return 2
    if (run_dir / "out").exists():
        shutil.rmtree(run_dir / "out")
    shutil.copytree(src, run_dir / "out")

    from mct_dft.runner import run_command
    import os
    env = dict(os.environ, OMP_NUM_THREADS="1", TMPDIR=str(ROOT / ".tmp"))
    r = run_command([str(ROOT / ".local/qe-env/bin/mpirun"), "-np",
                     str(args.np), str(ROOT / ".local/qe-env/bin/pw.x"),
                     "-in", "pw_mass.in"], run_dir, run_dir / "pw_mass.out",
                    env=env)
    out_text = (run_dir / "pw_mass.out").read_text(errors="replace")
    kpts, eigs = parse_bands_output(out_text)
    if eigs.size == 0:
        print("no eigenvalues parsed; job failed?")
        return 1

    # reference: occupied/unoccupied split from the median heuristic
    ref = float(np.median(eigs))
    results = {}
    nk = n
    fig, axes = plt.subplots(1, len(DIRECTIONS), figsize=(4 * len(DIRECTIONS), 4),
                             sharey=False)
    for i, (label, vec) in enumerate(DIRECTIONS.items()):
        seg_k = kpts[i * nk:(i + 1) * nk]
        seg_e = eigs[i * nk:(i + 1) * nk]
        # Cartesian |k| relative to Gamma, signed along the direction
        kcart = seg_k @ recip
        s = np.sign(np.linspace(-1, 1, nk))
        kd = np.linalg.norm(kcart, axis=1) * s   # signed 1/Angstrom

        e_gamma = seg_e[nk // 2]
        if args.band == "cb":
            above = np.where(e_gamma > ref)[0]
            bi = above[np.argmin(e_gamma[above])]
        else:
            below = np.where(e_gamma < ref)[0]
            bi = below[np.argmax(e_gamma[below])]
        band_e = seg_e[:, bi]

        fits = window_scan(kd, band_e, windows=(0.01, 0.02, 0.03))
        results[label] = [
            {"window_inv_ang": f.window_inv_ang, "mass_me": f.mass_me,
             "sign": f.sign, "r2": f.r_squared, "ok": f.parabolic_ok,
             "n_points": f.n_points} for f in fits
        ]
        ax = axes[i]
        ax.plot(kd, band_e - band_e[nk // 2], "o", ms=3, color="#0b5394")
        for f in fits:
            m = np.abs(kd - f.k_center_inv_ang) <= f.window_inv_ang
            x = kd[m] - f.k_center_inv_ang
            ax.plot(kd[m], f.e0_ev + 0.5 * f.curvature_ev_ang2 * x**2
                    - band_e[nk // 2], ls="--", lw=0.8,
                    label=f"w={f.window_inv_ang}: m*={f.mass_me:.3f}")
        ax.set_title(label)
        ax.set_xlabel("k (1/Angstrom, signed)")
        ax.legend(fontsize=7)
    axes[0].set_ylabel(f"E - E(extremum) (eV) [{args.band}]")
    fig.suptitle(f"{args.system} {'+SOC' if args.soc else 'no-SOC'} "
                 f"effective mass fit ({args.band})")
    fig.tight_layout()
    fig_path = ROOT / "results/figures" / f"{args.system}_mass_{args.band}_{which}.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)

    out = {"system": args.system, "band": args.band, "soc": args.soc,
           "run_dir": str(run_dir), "figure": str(fig_path),
           "fits": results}
    jf = run_dir / "mass_fit.json"
    jf.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    main()
