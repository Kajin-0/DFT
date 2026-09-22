#!/usr/bin/env python3
"""Generate Hg0.8Cd0.2Te SQS supercells with icet (cation sublattice only).

Cd substitutes Hg on the FCC cation sublattice; the Te anion sublattice
stays pure. Target Cd fraction on the cation sublattice: exactly 0.20.
Cluster-space quality (deviation from the ideal random-alloy cluster
vector) is evaluated and stored together with the random seed, so the
structure is fully reproducible.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def primitive_binary_template(a_ang: float):
    """ZnS-type primitive cell with mixed sublattice declaration for icet."""
    from ase import Atoms
    pos = [(0.0, 0.0, 0.0), (0.25, 0.25, 0.25)]
    cell = 0.5 * a_ang * np.array([[0.0, 1, 1], [1, 0, 1], [1, 1, 0]])
    # placeholder species; chemical_symbols given to ClusterSpace control
    # which sites mix: site 0 (cation) {Hg,Cd}, site 1 (anion) {Te}
    return Atoms("HgTe", cell=cell, scaled_positions=pos, pbc=True)


def random_reference_vector(cs, n_atoms, cation_fraction, n_samples, rng):
    """Average cluster vector of randomly-occupied supercells (= target)."""
    from ase import Atoms
    # big cubic supercell of the FCC cation lattice is subsumed by icet's
    # internal machinery; here we sample random occupations of the primitive
    # repetitions directly using icet's enumeration-free approach:
    import itertools
    from ase.build import make_supercell
    # 5x5x5 of the primitive cell gives 125 cations -> plenty for statistics
    P = np.eye(3, dtype=int) * 5
    prim = build_template_atoms(cs)
    base = make_supercell(prim, P)
    idx_cat = [i for i, s in enumerate(base.get_chemical_symbols())
               if s in ("Hg", "Cd")]
    vecs = []
    for _ in range(n_samples):
        occ = np.array(["Hg"] * len(idx_cat))
        occ[rng.random(len(idx_cat)) < cation_fraction] = "Cd"
        cand = base.copy()
        syms = cand.get_chemical_symbols()
        for j, i in enumerate(idx_cat):
            syms[i] = occ[j]
        cand.set_chemical_symbols(syms)
        vecs.append(cs.get_cluster_vector(cand))
    return np.mean(vecs, axis=0)


def build_template_atoms(cs):
    return cs.primitive_structure


def sqs_metrics(cs, sqs, rng, x_cd=0.2, n_ref=300):
    ref = random_reference_vector(cs, None, x_cd, n_ref, rng)
    got = cs.get_cluster_vector(sqs)
    d = np.asarray(got) - np.asarray(ref)
    return {
        "rms_deviation": float(np.sqrt(np.mean(d**2))),
        "max_abs_deviation": float(np.max(np.abs(d))),
        "n_clusters": int(len(d)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-start", type=float, default=None,
                    help="zincblende lattice parameter; default = Vegard of "
                         "configs' starting constants (0.8*6.453+0.2*6.482)")
    ap.add_argument("--total-atoms", type=int, default=40,
                    help="SQS supercell size (must contain N cations divisible by 5)")
    ap.add_argument("--seed", type=int, default=20260922)
    ap.add_argument("--pair-cutoff", type=float, default=8.0)
    ap.add_argument("--triplet-cutoff", type=float, default=5.0)
    ap.add_argument("--n-steps", type=int, default=60000)
    ap.add_argument("--tag", default="sqs40")
    args = ap.parse_args()

    if args.total_atoms % 10 != 0:
        raise SystemExit("total atoms must be a multiple of 10 for x=0.20 exact")

    from icet import ClusterSpace
    from icet.tools.structure_generation import generate_sqs

    a_start = args.a_start or (0.8 * 6.453 + 0.2 * 6.482)
    prim = primitive_binary_template(a_start)
    cs = ClusterSpace(prim, cutoffs=[args.pair_cutoff, args.triplet_cutoff],
                      chemical_symbols=[["Hg", "Cd"], ["Te"]])

    rng = np.random.default_rng(args.seed)
    seed = int(rng.integers(0, 2**31 - 1))
    # icet's max_size is measured in primitive-cell units for this binary
    # template (verified empirically: max_size=40 -> 80-atom supercell).
    sqs = generate_sqs(cluster_space=cs,
                       max_size=args.total_atoms // 2,
                       target_concentrations={"Hg": 0.80, "Cd": 0.20},
                       include_smaller_cells=False,
                       n_steps=args.n_steps,
                       random_seed=seed,
                       optimality_weight=1.0)

    n_hg = sum(1 for s in sqs if s.symbol == "Hg")
    n_cd = sum(1 for s in sqs if s.symbol == "Cd")
    n_te = sum(1 for s in sqs if s.symbol == "Te")
    assert n_te == args.total_atoms // 2, "anion sublattice must be pure Te"
    x = n_cd / (n_hg + n_cd)

    metrics = sqs_metrics(cs, sqs, np.random.default_rng(args.seed + 1))

    outdir = ROOT / "calculations" / "hgcdte" / args.tag
    outdir.mkdir(parents=True, exist_ok=True)
    from ase.io import write as ase_write
    ase_write(outdir / "sqs.cif", sqs)
    ase_write(outdir / "sqs.extxyz", sqs)
    record = {
        "x_cd_cation_sublattice": x,
        "counts": {"Hg": n_hg, "Cd": n_cd, "Te": n_te},
        "total_atoms": len(sqs),
        "lattice_param_used_ang": a_start,
        "cluster_space_cutoffs": [args.pair_cutoff, args.triplet_cutoff],
        "random_seed_external": args.seed,
        "random_seed_icet": seed,
        "n_monte_carlo_steps": args.n_steps,
        "sqs_quality": metrics,
        "cell_angstrom": sqs.cell.array.tolist(),
        "note": ("Cd substitutes Hg on the FCC cation sublattice only; Te "
                 "anion sublattice is untouched. Cluster-vector deviations "
                 "vs. the ideal random alloy are given above."),
    }
    (outdir / "sqs_record.json").write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
