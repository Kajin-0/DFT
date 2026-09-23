"""Band-edge bookkeeping: VBM/CBM, gap, direct/indirect.

Never force a positive gap: if band ordering is inverted (VBM>CBM) we report
the signed gap and mark the system 'inverted'/'overlap'.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import HC_EV_UM


@dataclass
class BandEdges:
    vbm_ev: float
    cbm_ev: float
    gap_ev: float
    k_vbm: np.ndarray
    k_cbm: np.ndarray
    i_vbm: int
    i_cbm: int
    direct: bool
    metallic: bool
    status: str  # 'semiconductor', 'semimetal/overlap', 'metallic', ...
    occupations_reference_ev: float  # energy used as occupancy threshold


def band_edges(
    kpts_cryst: np.ndarray,
    eigenvalues_ev: np.ndarray,
    reference_ev: float | None = None,
    tol_ev: float = 1e-4,
    n_valence: int | None = None,
) -> BandEdges:
    """Locate VBM/CBM from a dense eigenvalue grid.

    Two modes:

    * **n_valence given (preferred)** — the band-count definition: VBM =
      max_k E_{n_valence}(k), CBM = min_k E_{n_valence+1}(k). This can produce
      a *negative* gap, signalling an indirect band overlap / band inversion
      (e.g. semimetallic HgTe). Honesty for narrow/negative gap systems lives
      in this mode.
    * **reference_ev fallback** — occupancy threshold; only valid for true
      semiconductors (Fermi lies in a real gap), can never return negative.
    """
    E = np.asarray(eigenvalues_ev, dtype=float)      # (nk, nb)
    nk, nb = E.shape

    if n_valence is not None and 0 < n_valence < nb:
        vb_col = E[:, n_valence - 1]
        cb_col = E[:, n_valence]
        i_vbm = int(np.argmax(vb_col))
        i_cbm = int(np.argmin(cb_col))
        vbm = float(vb_col[i_vbm])
        cbm = float(cb_col[i_cbm])
        gap = cbm - vbm
        ref = reference_ev if reference_ev is not None else (vbm + cbm) / 2
        # Γ can legitimately occur at multiple path indices: directness is a
        # property of the k-VECTOR, not the list index
        kvbm = np.asarray(kpts_cryst)[i_vbm]
        kcbm = np.asarray(kpts_cryst)[i_cbm]
        direct = bool(np.allclose(kvbm % 1.0, kcbm % 1.0, atol=1e-4))
        if gap < -tol_ev:
            status = "inverted/overlap (band-count gap < 0): not a normal gap"
        elif abs(gap) <= tol_ev:
            status = "gapless touching bands"
        else:
            status = ("direct-gap semiconductor" if direct
                      else "indirect-gap semiconductor")
        return BandEdges(
            vbm_ev=vbm, cbm_ev=cbm, gap_ev=float(gap),
            k_vbm=kvbm, k_cbm=kcbm,
            i_vbm=i_vbm, i_cbm=i_cbm, direct=direct, metallic=False,
            status=status, occupations_reference_ev=float(ref))

    if reference_ev is None:
        reference_ev = float(np.median(E))

    occ_mask = E < reference_ev - tol_ev
    unocc_mask = E > reference_ev + tol_ev

    vbm, cbm = -np.inf, np.inf
    i_vbm = i_cbm = -1
    if occ_mask.any():
        i_flat = int(np.argmax(np.where(occ_mask, E, -np.inf)))
        i_vbm = i_flat // nb
        vbm = float(E.flat[i_flat])
    if unocc_mask.any():
        i_flat = int(np.argmin(np.where(unocc_mask, E, np.inf)))
        i_cbm = i_flat // nb
        cbm = float(E.flat[i_flat])

    gap = cbm - vbm
    metallic = not (occ_mask.any() and unocc_mask.any())
    direct = (not metallic and i_vbm >= 0 and i_cbm >= 0 and bool(
        np.allclose(np.asarray(kpts_cryst)[i_vbm] % 1.0,
                    np.asarray(kpts_cryst)[i_cbm] % 1.0, atol=1e-4)))

    if metallic:
        status = "metallic (no occupied/unoccupied split at reference)"
    elif gap < 0:
        status = "inverted/overlap (indirect band overlap): not a normal gap"
    elif abs(gap) < 1e-6 and not direct:
        status = "zero-gap / touching bands"
    else:
        status = ("direct-gap semiconductor" if direct
                  else "indirect-gap semiconductor")

    return BandEdges(
        vbm_ev=vbm, cbm_ev=cbm, gap_ev=float(gap),
        k_vbm=np.asarray(kpts_cryst)[i_vbm] if i_vbm >= 0 else np.full(3, np.nan),
        k_cbm=np.asarray(kpts_cryst)[i_cbm] if i_cbm >= 0 else np.full(3, np.nan),
        i_vbm=i_vbm, i_cbm=i_cbm, direct=direct, metallic=metallic,
        status=status, occupations_reference_ev=reference_ev,
    )


def n_valence_from_electrons(n_electrons: float, noncollinear: bool) -> int:
    """Number of valence bands: ceil(ne/2) collinear, ceil(ne) spinor."""
    import math
    return math.ceil(n_electrons) if noncollinear else math.ceil(n_electrons / 2)


def cutoff_um_from_gap(gap_ev: float) -> float | None:
    """Cutoff wavelength from a *positive, physically meaningful* gap."""
    if gap_ev <= 0 or not np.isfinite(gap_ev):
        return None
    return HC_EV_UM / gap_ev
