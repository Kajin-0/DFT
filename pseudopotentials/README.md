# Pseudopotentials

**Family:** PSlibrary 1.0.0 (A. Dal Corso, Comput. Mater. Sci. 95 (2014) 337)
**XC:** PBE for every element — no functional mixing.
**Type:** projector-augmented rrkj ultrasoft (`rrkjus`).

Two roles, both pinned by SHA-256 in `manifest.json`:

| role   | files | use |
|--------|-------|-----|
| `soc`   | `*.rel-pbe-*` fully relativistic, `has_so="T"` | production SOC runs (`noncolin=.true.`, `lspinorb=.true.`) |
| `nosoc` | `*.pbe-*` scalar-relativistic twins, same generator | SOC-off control (QE refuses FR files with `lspinorb=.false.`) |

UPF files themselves are **git-ignored**. Regenerate with a bit-identical
result via:

```bash
python scripts/fetch_pseudopotentials.py
```

which downloads from the official QE pseudopotential server, parses the UPF
headers (element, functional, relativistic character, `has_so`, Z_valence),
rejects anything that fails SOC validation, and rewrites the manifest with
SHA-256 hashes and retrieval timestamps. Headers — not filenames — decide
fitness for spin-orbit production runs.
