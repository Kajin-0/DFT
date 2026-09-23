#!/usr/bin/env python3
"""Optical + ideal detector chain from epsilon.x outputs.

epsilon.x writes epsr/epsi files (real/imag dielectric vs energy, eV), per
Cartesian component. We aggregate the diagonal (zincblende is cubic:
isotropic diagnostic), derive n,k, alpha, R, ideal absorptance, ideal PV
responsivity, with SI-unit conversions via mct_dft.optics/constants.
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
from mct_dft.optics import (dielectric_to_nk, absorption_coefficient_per_m,  # noqa
                            normal_reflectance, external_quantum_efficiency)
from mct_dft.detector import ideal_responsivity_a_per_w  # noqa: E402
from mct_dft.constants import HC_EV_UM  # noqa: E402


def load_epsilon(run_dir: Path):
    """Returns arrays (energy_eV, eps1_avg, eps2_avg) over diagonal comps."""
    epsr, epsi = None, None
    for cand in sorted(run_dir.glob("epsr*.dat")) + sorted(
            run_dir.glob("*epsr*")):
        epsr = np.loadtxt(cand)
        break
    for cand in sorted(run_dir.glob("epsi*.dat")) + sorted(
            run_dir.glob("*epsi*")):
        epsi = np.loadtxt(cand)
        break
    if epsr is None or epsi is None:
        raise FileNotFoundError(f"epsr/epsi dat files not found in {run_dir}")
    # columns: E, xx, yy, zz (sometimes extra)
    e = epsr[:, 0]
    eps1 = epsr[:, 1:4].mean(axis=1) if epsr.shape[1] >= 4 else epsr[:, 1]
    eps2 = epsi[:, 1:4].mean(axis=1) if epsi.shape[1] >= 4 else epsi[:, 1]
    return e, eps1, eps2


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path,
                    help="run dir of an eps stage (e.g. calculations/cdte/optics_*)")
    ap.add_argument("--thickness-um", type=float, default=10.0)
    ap.add_argument("--tag", default="cdte")
    args = ap.parse_args()

    e, e1, e2 = load_epsilon(args.run_dir)
    lam_m = HC_EV_UM / e * 1e-6
    n, k = dielectric_to_nk_vec(e1, e2)
    alpha = absorption_coefficient_per_m(k, lam_m)      # 1/m
    with np.errstate(invalid="ignore"):
        R = (np.abs((n + 1j * k - 1) / (n + 1j * k + 1))) ** 2
    d = args.thickness_um * 1e-6
    eta = np.array([external_quantum_efficiency(nn, kk, lm, d)
                    for nn, kk, lm in zip(n, k, lam_m)])
    resp = ideal_responsivity_a_per_w(eta, lam_m)

    out_csv = ROOT / "results/tables" / f"{args.tag}_optics.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out_csv, np.column_stack(
        [e, lam_m * 1e6, e1, e2, n, k, alpha, R, eta, resp]),
        header=("E_eV lambda_um eps1 eps2 n k alpha_1perm R "
                "ideal_QE ideal_R_A_per_W"), comments="# ")
    print(f"wrote {out_csv}")

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    ax = axes[0, 0]
    ax.plot(e, n, label="n")
    ax.plot(e, k, label="k")
    ax.set_xlabel("E (eV)"); ax.set_ylabel("index")
    ax.set_title("complex refractive index"); ax.legend()

    ax = axes[0, 1]
    ax.plot(lam_m * 1e6, alpha / 100.0)
    ax.set_xlabel("wavelength (um)"); ax.set_ylabel("alpha (1/cm)")
    ax.set_title("absorption coefficient")
    ax.set_xlim(0, min(15, (lam_m * 1e6)[-1]))

    ax = axes[1, 0]
    ax.plot(lam_m * 1e6, R, label="R")
    ax.plot(lam_m * 1e6, eta, label="ideal QE (1-R)(1-e^{-ad})")
    ax.set_xlabel("wavelength (um)")
    ax.set_title("reflectance / ideal QE")
    ax.legend()
    ax.set_xlim(0, min(15, (lam_m * 1e6)[-1]))

    ax = axes[1, 1]
    ax.plot(lam_m * 1e6, resp)
    ax.set_xlabel("wavelength (um)")
    ax.set_ylabel("R (A/W)")
    ax.set_title("IDEAL PV responsivity (perfect collection)")
    ax.set_xlim(0, min(15, (lam_m * 1e6)[-1]))
    fig.suptitle(f"{args.tag} optical chain — d = {args.thickness_um} um")
    fig.tight_layout()
    out_png = ROOT / "results/figures" / f"{args.tag}_optics.png"
    fig.savefig(out_png, dpi=200)
    print(f"wrote {out_png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())


def dielectric_to_nk_vec(e1, e2):
    v = np.sqrt(np.asarray(e1) + 1j * np.asarray(e2))
    return v.real, np.abs(v.imag)
