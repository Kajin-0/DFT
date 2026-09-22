#!/usr/bin/env python3
"""Run one QE stage for one system with full provenance.

Examples
--------
  # CdTe cutoff-convergence SCF
  python scripts/run_stage.py cdte scf --ecutwfc 60 --ecutrho 480 --kgrid 8 8 8
  # CdTe SOC band structure along the automatic SeekPath path
  python scripts/run_stage.py cdte bands --soc --ecutwfc 60 --ecutrho 480
  # HgTe DOS on a dense grid
  python scripts/run_stage.py hgte dos --soc --ecutwfc 60 --ecutrho 480 --kgrid 16 16 16

Every materialised artefact is contained inside the repository and each run
directory receives: input, raw stdout/stderr, parsed result.json, and a
provenance.json (git commit, QE/python versions, PP sha256, MPI ranks...).
Existing results are NOT overwritten unless --force is given.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mct_dft.parser import parse_pw_output, parse_bands_output  # noqa: E402
from mct_dft.provenance import (package_versions, git_commit,      # noqa: E402
                                write_provenance, sha256_file)
from mct_dft.qe_inputs import (pw_input, bands_x_input, dos_x_input,  # noqa: E402
                               projwfc_x_input, write_text)
from mct_dft.runner import ensure_within_workspace, run_command  # noqa: E402
from mct_dft.structures import (zincblende_primitive, validate_structure,  # noqa
                                LATTICE_START)

QE_BIN = ROOT / ".local/qe-env/bin"
PSEUDO_DIR = ROOT / "pseudopotentials"


def default_np() -> int:
    cores = os.cpu_count() or 1
    return max(1, int(cores * 0.75))


def load_manifest(soc: bool) -> dict[str, str]:
    """element -> UPF filename, for the requested relativistic treatment."""
    role = "soc" if soc else "nosoc"
    m = json.loads((PSEUDO_DIR / "manifest.json").read_text())
    return {e["element"]: Path(e["local_path"]).name
            for e in m["potentials"] if e["role"] == role}


def kpath_seekpath(atoms, npoints_per_segment: int = 40):
    """High-symmetry k path from SeeK-path (works for the primitive zb cell)."""
    from seekpath import get_explicit_k_path
    cell = (atoms.cell.array, atoms.get_scaled_positions(),
            atoms.get_atomic_numbers())
    res = get_explicit_k_path(cell, reference_distance=0.1)
    return res  # dict with explicit_kpoints_rel, path segments, labels


def segment_kpoints(pts: list[list[float]], nseg: int = 40):
    """Explicit K_POINTS list (crystal) walking segments, unit weights."""
    rows = []
    for a, b in zip(pts[:-1], pts[1:]):
        a, b = np.asarray(a, float), np.asarray(b, float)
        for i in range(nseg):
            t = i / nseg
            rows.append([*(a + t * (b - a)), 1.0])
    rows.append([*np.asarray(pts[-1], float), 1.0])
    return np.asarray(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", choices=["cdte", "hgte", "hgcdte"])
    ap.add_argument("stage", choices=["scf", "relax", "vc-relax", "nscf",
                                      "bands", "dos", "pdos"])
    ap.add_argument("--ecutwfc", type=float, default=None)
    ap.add_argument("--ecutrho", type=float, default=None)
    ap.add_argument("--kgrid", type=int, nargs=3, default=[8, 8, 8])
    ap.add_argument("--lattice", type=float, default=None,
                    help="lattice parameter in Angstrom; default = config start")
    ap.add_argument("--soc", action="store_true")
    ap.add_argument("--np", type=int, default=None, help="MPI ranks (default 3/4 cores)")
    ap.add_argument("--nbnd", type=int, default=None)
    ap.add_argument("--tag", default=None, help="override run-directory name")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--charge-from", default=None,
                    help="repo-relative run dir whose out/ (SCF charge density) "
                         "is copied in before a bands/nscf/dos/pdos stage")
    ap.add_argument("--input-sha", default=None,
                    help="reuse structure/pseudos exactly like a previous run")
    ap.add_argument("--struct", default=None,
                    help="repo-relative ASE-readable structure file (e.g. an "
                         "SQS extxyz). Overrides the zincblende primitive "
                         "builder; required for --system hgcdte")
    draw = ap.add_argument_group("bands")
    draw.add_argument("--path-seg-points", type=int, default=40)
    args = ap.parse_args()

    cfg = yaml.safe_load(
        (ROOT / "config/systems" / f"{args.system}.yaml").read_text())
    species: list[str] = cfg.get("species", [])
    a0 = args.lattice or cfg.get("lattice_start_ang") or (
        LATTICE_START.get(cfg["name"]) if cfg.get("name") in LATTICE_START else None)
    if args.struct is None and (not species or a0 is None):
        print("ERROR: system requires --struct (no primitive builder defined)")
        return 2
    ecut = args.ecutwfc or 60.0
    ecutrho = args.ecutrho or 4 * ecut  # placeholder; override explicitly

    tag = args.tag
    if tag is None:
        bits = [args.stage, f"e{ecut:g}", f"k{'x'.join(map(str, args.kgrid))}"]
        if args.soc:
            bits.append("soc")
        if args.lattice:
            bits.append(f"a{args.lattice:.4f}")
        tag = "_".join(bits)
    run_dir = ensure_within_workspace(ROOT / "calculations" / args.system / tag)
    result_json = run_dir / "result.json"
    if result_json.exists() and not args.force:
        try:
            prev = json.loads(result_json.read_text())
            if prev.get("pw", {}).get("ok"):
                print(f"[skip] {run_dir} already has a successful result.json"
                      " (use --force)")
                return 0
            print(f"[retry] {run_dir} has a failed result.json; re-running")
        except Exception:
            pass
    run_dir.mkdir(parents=True, exist_ok=True)
    outdir = run_dir / "out"
    if outdir.exists():
        shutil.rmtree(outdir)
    if args.charge_from:
        src = ensure_within_workspace(ROOT / args.charge_from) / "out"
        if not src.is_dir():
            print(f"ERROR: --charge-from dir has no out/: {src}")
            return 2
        shutil.copytree(src, outdir)
        print(f"[charge] copied out/ from {src}")

    if args.struct:
        from ase.io import read as ase_read
        atoms = ase_read(str(ensure_within_workspace(ROOT / args.struct)))
        atoms.set_pbc(True)
        struct_check = validate_structure(atoms)
        struct_check["source_structure"] = args.struct
    else:
        atoms = zincblende_primitive(species[0], species[1], a0)
        struct_check = validate_structure(atoms)
    pseudo = load_manifest(args.soc)

    np_ranks = args.np or default_np()
    np_ranks = max(1, min(np_ranks, os.cpu_count() or 1))

    env = dict(os.environ)
    env.update({
        "OMP_NUM_THREADS": "1",
        "TMPDIR": str(ROOT / ".tmp"),
        "ESPRESSO_PSEUDO": str(PSEUDO_DIR),
    })

    runs: list[dict] = []

    def record(argv, cwd, stdout, record_extra=None):
        r = run_command(argv, cwd, stdout, env=env)
        if record_extra:
            r.update(record_extra)
        runs.append(r)
        return r

    common = dict(
        pseudopotentials=pseudo, ecutwfc_ry=ecut, ecutrho_ry=ecutrho,
        soc=args.soc, pseudo_dir=str(PSEUDO_DIR), outdir=str(outdir),
    )

    if args.stage in ("scf", "relax", "vc-relax"):
        text = pw_input(
            atoms, kgrid=tuple(args.kgrid), nbnd=args.nbnd,
            relax=args.stage == "relax", vc_relax=args.stage == "vc-relax",
            **common)
        write_text(run_dir / "pw.in", text)
        record([str(QE_BIN / "mpirun"), "-np", str(np_ranks), str(QE_BIN / "pw.x"),
                "-in", "pw.in"], run_dir, run_dir / "pw.out")

    elif args.stage == "nscf":
        text = pw_input(atoms, calculation="nscf", occupations="tetrahedra",
                        kgrid=tuple(args.kgrid), nbnd=args.nbnd, **common)
        write_text(run_dir / "pw.in", text)
        record([str(QE_BIN / "mpirun"), "-np", str(np_ranks), str(QE_BIN / "pw.x"),
                "-in", "pw.in"], run_dir, run_dir / "pw.out")

    elif args.stage == "bands":
        sp = kpath_seekpath(atoms)
        # write dense walk of the path
        seg_points = segment_kpoints(sp["explicit_kpoints_rel"],
                                     args.path_seg_points)
        text = pw_input(atoms, calculation="bands", kpoints_explicit=seg_points,
                        nbnd=args.nbnd, **common)
        write_text(run_dir / "pw_bands.in", text)
        record([str(QE_BIN / "mpirun"), "-np", str(np_ranks), str(QE_BIN / "pw.x"),
                "-in", "pw_bands.in"], run_dir, run_dir / "pw_bands.out")
        write_text(run_dir / "bands.in", bands_x_input())
        record([str(QE_BIN / "bands.x"), "-in", "bands.in"], run_dir,
               run_dir / "bands.out")
        with open(run_dir / "seekpath.json", "w") as fh:
            json.dump(sp, fh, indent=2, default=lambda o: np.asarray(o).tolist())
    elif args.stage == "dos":
        text = pw_input(atoms, calculation="nscf", occupations="tetrahedra",
                        kgrid=tuple(args.kgrid), nbnd=args.nbnd, **common)
        write_text(run_dir / "pw.in", text)
        record([str(QE_BIN / "mpirun"), "-np", str(np_ranks), str(QE_BIN / "pw.x"),
                "-in", "pw.in"], run_dir, run_dir / "pw.out")
        write_text(run_dir / "dos.in", dos_x_input(de=0.005))
        record([str(QE_BIN / "dos.x"), "-in", "dos.in"], run_dir,
               run_dir / "dos.out")
    elif args.stage == "pdos":
        text = pw_input(atoms, calculation="nscf", occupations="tetrahedra",
                        kgrid=tuple(args.kgrid), nbnd=args.nbnd, **common)
        write_text(run_dir / "pw.in", text)
        record([str(QE_BIN / "mpirun"), "-np", str(np_ranks), str(QE_BIN / "pw.x"),
                "-in", "pw.in"], run_dir, run_dir / "pw.out")
        write_text(run_dir / "projwfc.in", projwfc_x_input())
        record([str(QE_BIN / "mpirun"), "-np", str(np_ranks),
                str(QE_BIN / "projwfc.x"), "-in", "projwfc.in"], run_dir,
               run_dir / "projwfc.out")

    # ---------------- parse + provenance ----------------
    pw_out_path = run_dir / ("pw_bands.out" if args.stage == "bands" else "pw.out")
    text_out = pw_out_path.read_text(errors="replace") if pw_out_path.exists() else ""
    parsed = parse_pw_output(text_out)
    n_fail = sum(1 for r in runs if r["returncode"] != 0)

    result = {
        "system": args.system, "stage": args.stage, "soc": args.soc,
        "ecutwfc_ry": ecut, "ecutrho_ry": ecutrho, "kgrid": args.kgrid,
        "lattice_ang_in": a0,
        "structure_check": struct_check,
        "mpi_ranks": np_ranks,
        "runs": runs,
        "pw": {
            "ok": parsed.ok, "converged": parsed.converged,
            "job_done": parsed.job_done,
            "total_energy_ry": parsed.total_energy_ry,
            "total_energy_ev": parsed.total_energy_ev,
            "fermi_ev": parsed.fermi_ev,
            "highest_occupied_ev": parsed.highest_occupied_ev,
            "lowest_unoccupied_ev": parsed.lowest_unoccupied_ev,
            "n_iterations": parsed.nscf_iterations,
            "n_bands": parsed.n_bands, "n_kpts": parsed.n_kpts,
            "gap_ev": None,
            "noncollinear": parsed.noncollinear, "spin_orbit": parsed.spin_orbit,
            "calculation": parsed.calculation,
            "final_cell_ang": parsed.final_cell_ang,
            "final_positions_crystal": parsed.final_positions_crystal,
            "errors": parsed.errors,
        },
        "returncode_failures": n_fail,
    }
    if args.stage == "bands" and parsed.job_done:
        kpts, eigs = parse_bands_output(text_out)
        if eigs.size:
            np.savez(run_dir / "bands.npz", kpts=kpts, eigs_ev=eigs)
    # occupancy-based band edges from the eigenvalue listing when available
    if parsed.fermi_ev is not None:
        _k, _e = parse_bands_output(text_out)
        if _e.size:
            occ = _e[_e < parsed.fermi_ev]
            uno = _e[_e >= parsed.fermi_ev]
            if occ.size:
                result["pw"]["highest_occupied_ev"] = float(occ.max())
            if uno.size:
                result["pw"]["lowest_unoccupied_ev"] = float(uno.min())
            if occ.size and uno.size:
                result["pw"]["gap_ev"] = float(uno.min() - occ.max())

    result_json.write_text(json.dumps(result, indent=2, default=str))

    # input hash
    in_files = [p for p in run_dir.glob("*.in")]
    hashes = {p.name: sha256_file(str(p)) for p in in_files}
    prov = {
        "git_commit": git_commit(ROOT),
        "python_packages": package_versions(),
        "qe_bin_dir": str(QE_BIN),
        "mpi_ranks": np_ranks,
        "input_sha256": hashes,
        "pseudopotentials": json.loads((PSEUDO_DIR / "manifest.json").read_text())["potentials"],
    }
    write_provenance(prov, run_dir / "provenance.json")

    print(json.dumps({"run_dir": str(run_dir), "ok": parsed.ok,
                      "E_tot_Ry": parsed.total_energy_ry,
                      "E_F_eV": parsed.fermi_ev,
                      "wall_s": sum(r["wall_time_s"] for r in runs)}, indent=2))
    return 0 if parsed.ok and n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
