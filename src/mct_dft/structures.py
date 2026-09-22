"""Crystal structure construction and validation (ASE-based).

Zincblende primitive cell: FCC primitive lattice, two-atom basis
   cation at (0,0,0), anion at (1/4,1/4,1/4)  (fractional, primitive cell).
Space group of zincblende: F-43m (No. 216).
"""
from __future__ import annotations

import math

import numpy as np
from ase import Atoms
from ase.build import bulk as ase_bulk
from ase.neighborlist import neighbor_list

# Literature starting lattice constants (Angstrom), used as *starting points*
# for DFT optimisation, notvera as results.
LATTICE_START = {"CdTe": 6.482, "HgTe": 6.453}


def zincblende_primitive(species_a: str, species_b: str, a_ang: float) -> Atoms:
    """Primitive (2-atom) zincblende cell via ASE, species_a = cation.

    Layout: cation at (0,0,0), anion at (1/4,1/4,1/4) of the primitive cell.
    """
    atoms = ase_bulk(f"{species_a}{species_b}", "zincblende",
                     a=a_ang, cubic=False)
    # ASE zincblende basis: cation at (0,0,0), anion at (1/4,1/4,1/4) of the
    # primitive cell — exactly the required basis.
    atoms.set_pbc(True)
    return atoms


def zincblende_conventional(species_a: str, species_b: str, a_ang: float) -> Atoms:
    try:
        atoms = ase_bulk(f"{species_a}{species_b}", "zincblende",
                         a=a_ang, cubic=True)
    except ValueError:
        # some ASE versions require per-species names; fall back to manual
        a = a_ang
        pos = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        pos_b = [(0.25, 0.25, 0.25), (0.75, 0.75, 0.25), (0.75, 0.25, 0.75),
                 (0.25, 0.75, 0.75)]
        atoms = Atoms([species_a] * 4 + [species_b] * 4,
                      scaled_positions=pos + pos_b,
                      cell=[(a, 0, 0), (0, a, 0), (0, 0, a)], pbc=True)
    atoms.set_pbc(True)
    return atoms


def validate_structure(atoms: Atoms, nn_tol_ang: float = 0.05) -> dict:
    """Basic structural sanity checks. Returns a dict of diagnostics.

    For zincblende the nearest-neighbour (cation-anion) distance must equal
    sqrt(3)/4 * a_conv and every atom must have exactly 4 nearest neighbours.
    """
    vol = float(atoms.get_volume())
    # primitive zincblende volume = a^3/4 -> a_conv = (4V)^(1/3);
    # NN distance = sqrt(3)/4 a; use a 25% margin for the neighbour search
    # conventional zb cell holds 8 atoms: a_conv^3 = 8 V/N (exact for any
    # supercell built from the same per-atom volume)
    a_conv = (8.0 * vol / len(atoms)) ** (1.0 / 3.0)
    d_nn = math.sqrt(3.0) / 4.0 * a_conv
    cutoff = 1.25 * d_nn
    i, j, d = neighbor_list("ijd", atoms, cutoff=cutoff)
    # nearest-neighbour distance = minimum bond length found
    dmin = float(d.min())
    counts: dict[int, int] = {}
    for idx in i[np.isclose(d, dmin, atol=nn_tol_ang)]:
        counts[int(idx)] = counts.get(int(idx), 0) + 1
    coordination = sorted(set(counts.values()))
    return {
        "n_atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "volume_ang3": vol,
        "volume_per_atom_ang3": vol / len(atoms),
        "nn_distance_ang": dmin,
        "nn_coordination": coordination,
        # expected a-zb: nn = sqrt(3)/4 * a; primitive-cell volume = a^3/4
        "implied_a_conventional_ang": dmin * 4.0 / math.sqrt(3.0),
    }


def _half_cell_dist(atoms: Atoms) -> float:  # kept for backward compat
    lat = atoms.cell.lengths()
    return float(0.60 * max(lat))
