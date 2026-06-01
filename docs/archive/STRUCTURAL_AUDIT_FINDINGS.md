# Structural Audit Findings — 2026-05-01

## Critical Issue: WT Structure File Mismatch

The `5i6z_A.pdb` file and all derivative mutant PDB files are **not** chain A
from the 5I6Z crystal structure. They originate from a different source.

### Evidence

| File | Res 277 | Res 348 | Res 439 | TM6 Cα RMSD vs raw |
|---|---|---|---|---|
| `5i6z.pdb` (true 5I6Z) | SER | **ALA** | SER | — |
| `5i6z_A.pdb` | TYR | **THR** | VAL | 37.2 Å |
| `5i6Z_S348T_3.pdb` | TYR | **THR** | VAL | 37.2 Å |
| `5i6Z_S348_A.pdb` | TYR | **THR** | VAL | 37.2 Å |

A 37 Å RMSD means these files share no structural overlap with 5I6Z —
they are a completely different protein structure (or a significantly different
coordinate frame from a homology model).

The earlier RMSD of 0.030 Å between "WT" and "S348T" reflected that the two
files were nearly identical, not that the mutation was well-tolerated.

### Root Cause

The `5i6z_A.pdb` was likely generated through homology modeling or processed
by a different tool that:
1. Used a different template or applied its own numbering
2. Already incorporated the S348T substitution before saving as the "WT"
3. Did not apply the S348A substitution correctly (file still shows THR at 348)

### What Was Done

1. True chain A extracted from `5i6z.pdb` → `data/structures/5i6z_A_true.pdb`
   - 544 ATOM residues, ALA at position 348 (canonical per 5I6Z SEQADV)

2. ChimeraX script written to generate correct mutants:
   - `scripts/generate_mutants.cxc` — uses `swapaa` on `5i6z_A_true.pdb`
   - Output: `5i6z_A_S348T_correct.pdb` (ALA348 → THR)
   - Output: `5i6z_A_S348A_correct.pdb` (ALA348 remains ALA — same as WT in
     this crystal structure context)

3. Docking pipeline prepared:
   - `scripts/docking_prep.py` — receptor PDBQT generation and Vina configs
   - `output/receptor_wt.pdbqt` and `output/vina_wt.conf` ready
   - `output/LIGAND_PREP_INSTRUCTIONS.txt` — how to generate ligand PDBQT

---

## 5I6Z Engineered Mutations (SEQADV)

The 5I6Z crystal structure itself has several engineered mutations vs canonical
human SERT (UniProt P31645):

| PDB Residue | Canonical |
|---|---|
| ALA 291 | ILE 291 |
| SER 439 | THR 439 |
| ALA 554 | CYS 554 |
| ALA 580 | CYS 580 |

Position 348 is **not** listed in SEQADV, meaning ALA348 is canonical for
this construct.

---

## Next Steps

### Immediate (required before docking)
1. Open ChimeraX, load `data/structures/5i6z_A_true.pdb`, run `generate_mutants.cxc`
2. Prepare escitalopram PDBQT (see `output/LIGAND_PREP_INSTRUCTIONS.txt`)
3. Place `vina.exe` in `scripts/` or on PATH
4. Run: `python scripts/docking_prep.py`

### Interpretation note on "S348T"
In the 5I6Z context, the mutation at position 348 is A348T (Alanine → Threonine).
If the biological intent is S348T (Serine → Threonine based on another species or
numbering scheme), a model starting from canonical human SERT sequence and
homology modeling is required, not the crystal structure directly.
