#!/usr/bin/env python3
"""Download pseudopotentials and build pseudopotentials/manifest.json.

Source: official Quantum ESPRESSO pseudopotential library (PSlibrary 1.0.0,
fully relativistic, PBE, rrkjus ultrasoft) — same generator family and XC for
all three elements. Files are verified by parsing their UPF headers, never by
filenames; the manifest records sha256 of every file actually used.

Fallbacks to older PSlibrary releases are tried in order if the primary file
is unavailable. The UPF files themselves are git-ignored; this script plus the
manifest (sha256 pins) makes the set reproducible.
"""
from __future__ import annotations

import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from mct_dft.pseudopotentials import (  # noqa: E402
    parse_upf_header, sha256_file, validate_soc_set,
)

BASE = "https://pseudopotentials.quantum-espresso.org/upf_files/"

CANDIDATES = {
    # role "soc": fully-relativistic rrkj ultrasoft, PBE (has_so=true)
    # role "nosoc": scalar-relativistic twin — SAME generator family + XC,
    #               used as the SOC-off control. (QE >= 7 refuses to run a
    #               fully-relativistic PP with lspinorb=.false.)
    "Hg": {
        "soc": [
            "Hg.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
            "Hg.rel-pbe-spnl-rrkjus_psl.1.0.0.UPF",
            "Hg.rel-pbe-n-rrkjus_psl.1.0.0.UPF",
        ],
        "nosoc": [
            "Hg.pbe-spn-rrkjus_psl.1.0.0.UPF",
            "Hg.pbe-spnl-rrkjus_psl.1.0.0.UPF",
            "Hg.pbe-n-rrkjus_psl.1.0.0.UPF",
        ],
    },
    "Cd": {
        "soc": [
            "Cd.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
            "Cd.rel-pbe-n-rrkjus_psl.1.0.0.UPF",
        ],
        "nosoc": [
            "Cd.pbe-spn-rrkjus_psl.1.0.0.UPF",
            "Cd.pbe-n-rrkjus_psl.1.0.0.UPF",
        ],
    },
    "Te": {
        "soc": [
            "Te.rel-pbe-dn-rrkjus_psl.1.0.0.UPF",
            "Te.rel-pbe-n-rrkjus_psl.1.0.0.UPF",
        ],
        "nosoc": [
            "Te.pbe-dn-rrkjus_psl.1.0.0.UPF",
            "Te.pbe-n-rrkjus_psl.1.0.0.UPF",
        ],
    },
}
# role "nc": norm-conserving scalar-relativistic — used ONLY for optics,
# because QE's epsilon.x does not support ultrasoft pseudopotentials.
# role "pawopt": PAW kjpaw scalar-relativistic (tested: epsilon.x rejects PAW
# too) -- kept for completeness but not used for production optics.
PAW_CANDIDATES = {
    "Hg": ["Hg.pbe-spn-kjpaw_psl.1.0.0.UPF", "Hg.pbe-n-kjpaw_psl.1.0.0.UPF"],
    "Cd": ["Cd.pbe-n-kjpaw_psl.1.0.0.UPF", "Cd.pbe-spn-kjpaw_psl.1.0.0.UPF"],
    "Te": ["Te.pbe-n-kjpaw_psl.1.0.0.UPF", "Te.pbe-dn-kjpaw_psl.1.0.0.UPF"],
}
for el, files in PAW_CANDIDATES.items():
    CANDIDATES[el]["pawopt"] = files

# role "nc": norm-conserving scalar-relativistic ONCV (SG15 set). Used ONLY
# for the epsilon.x optics branch (epsilon.x supports neither USPP nor PAW).
SG15_BASE = "http://www.quantum-simulation.org/potentials/sg15_oncv/upf/"
SG15_VERSION = "SG15 ONCVPSP v1.2 (Hamann, Phys. Rev. B 88 (2013) 085117)"

SRC_VERSION = {
    "psl.1.0.0": "PSlibrary 1.0.0 (Dal Corso, Comput. Mater. Sci. 95 (2014) 337)",
}


def source_version(fname: str) -> str:
    for key, tag in SRC_VERSION.items():
        if key in fname:
            return tag
    return "Quantum ESPRESSO upf_files repository"


def fetch(url: str, dest: Path) -> bool:
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = r.read()
        if len(data) < 2000 or b"<PP_HEADER" not in data:
            return False
        dest.write_bytes(data)
        return True
    except Exception as exc:
        print(f"  fetch failed: {url} ({exc})")
        return False


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ppdir = root / "pseudopotentials"
    ppdir.mkdir(exist_ok=True)
    manifest_entries = []
    soc_metas = []
    # SG15 NC set (optics branch only)
    sg15 = []
    for element in ("Hg", "Cd", "Te"):
        name = f"{element}_ONCV_PBE-1.2.upf"
        dest = ppdir / name
        if not dest.exists():
            print(f"[{element}/nc] fetching SG15 {name}")
            if not fetch(SG15_BASE + name, dest):
                print(f"ERROR: SG15 NC PP for {element} unavailable")
                return 1
        meta = parse_upf_header(str(dest))
        if meta.element.lower() != element.lower() or meta.is_ultrasoft \
                or meta.relativistic.lower() != "scalar":
            print(f"ERROR: SG15 {name} rejected: el={meta.element} "
                  f"us={meta.is_ultrasoft} rel={meta.relativistic}")
            return 1
        sg15.append((name, meta, dest))

    for element, roles in CANDIDATES.items():
        for role, candidates in roles.items():
            chosen = None
            for name in candidates:
                dest = ppdir / name
                if not dest.exists():
                    print(f"[{element}/{role}] trying {name}")
                    if not fetch(BASE + name, dest):
                        continue
                meta = parse_upf_header(str(dest))
                ok_el = meta.element.lower() == element.lower()
                ok_soc = meta.soc_capable() if role == "soc" else True
                ok_scalar = meta.relativistic.lower() == "scalar" \
                    if role == "nosoc" else True
                if role == "nc":
                    ok_scalar = (meta.relativistic.lower() == "scalar"
                                 and not meta.is_ultrasoft)
                if role == "pawopt":
                    ok_scalar = (meta.relativistic.lower() == "scalar"
                                 and meta.is_paw)
                if not (ok_el and ok_soc and ok_scalar):
                    print(f"  rejected {name}: element={meta.element!r} "
                          f"rel={meta.relativistic!r} has_so={meta.has_so}")
                    dest.unlink(missing_ok=True)
                    continue
                chosen = (name, meta, dest)
                break
            if chosen is None:
                print(f"ERROR: no usable {role} PBE UPF for {element}")
                return 1
            name, meta, dest = chosen
            if role == "soc":
                soc_metas.append(meta)
            manifest_entries.append({
                "original_filename": name,
                "local_path": str(dest.relative_to(root)),
                "element": element,
                "role": role,  # soc = production SOC; nosoc = SOC-off control
                "source": BASE + name,
                "source_release": source_version(name),
                "xc_functional": meta.functional,
                "pseudo_type": "ultrasoft" if meta.is_ultrasoft else meta.pseudo_type,
                "relativistic": meta.relativistic,
                "has_so": meta.has_so,
                "z_valence": meta.z_valence,
                "l_max": meta.l_max,
                "core_correction": meta.core_correction,
                "suggested_ecutwfc_ry": meta.suggested_ecutwfc_ry,
                "suggested_ecutrho_ry": meta.suggested_ecutrho_ry,
                "sha256": sha256_file(str(dest)),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            })
            print(f"[{element}/{role}] OK {name}  Zv={meta.z_valence} "
                  f"rel={meta.relativistic} so={meta.has_so}")

    for name, meta, dest in sg15:
        manifest_entries.append({
            "original_filename": name,
            "local_path": str(dest.relative_to(root)),
            "element": meta.element,
            "role": "nc",
            "source": SG15_BASE + name,
            "source_release": SG15_VERSION,
            "xc_functional": meta.functional,
            "pseudo_type": "norm-conserving",
            "relativistic": meta.relativistic,
            "has_so": meta.has_so,
            "z_valence": meta.z_valence,
            "l_max": meta.l_max,
            "core_correction": meta.core_correction,
            "suggested_ecutwfc_ry": meta.suggested_ecutwfc_ry,
            "suggested_ecutrho_ry": meta.suggested_ecutrho_ry,
            "sha256": sha256_file(str(dest)),
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        })
        print(f"[{meta.element}/nc] OK {name} Zv={meta.z_valence} "
              f"rel={meta.relativistic} NC")

    problems = validate_soc_set(soc_metas, xc_expected="PBE")
    if problems:
        print("SOC validation FAILED:")
        for p in problems:
            print("  -", p)
        return 2

    manifest = {
        "comment": ("All potentials: PSlibrary rrkjus family, PBE. role=soc are "
                    "fully relativistic with has_so (production); role=nosoc are "
                    "the scalar-relativistic twins used as the lspinorb=.false. "
                    "control (QE refuses FR PPs without lspinorb). UPF files are "
                    "git-ignored; sha256 pins + this script = reproducibility."),
        "validated_for_soc": True,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "potentials": manifest_entries,
    }
    out = ppdir / "manifest.json"
    out.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
