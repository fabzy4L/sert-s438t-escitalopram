# Session Summary — 2026-05-01

## What Was Accomplished

### Critical Discovery: Structural Integrity Issue

Running the audit revealed that `5i6z_A.pdb` and all three mutant files are **not from the 5I6Z crystal structure** — they share a 37 Å RMSD with the raw PDB at the TM6 helix. Every file had THR at position 348, which is why the prior "WT vs S348T" RMSD was 0.030 Å (comparing essentially identical structures). The prior analysis was comparing the same structure against itself.

### Files Built and Committed

| File | Description |
|---|---|
| `scripts/structural_audit.py` | Full audit of all PDB files; runs pairwise RMSD, identifies binding pocket centroid |
| `data/structures/5i6z_A_true.pdb` | The actual 5I6Z chain A (ALA at pos 348, clean baseline) |
| `scripts/generate_mutants.cxc` | ChimeraX script to produce correct S348T/S348A via `swapaa` from the true WT |
| `scripts/docking_prep.py` | Generates receptor PDBQT files + Vina configs; auto-runs docking if vina.exe is present |
| `output/receptor_wt.pdbqt` | WT receptor ready for docking |
| `output/vina_wt.conf` / `vina_s348t.conf` | AutoDock Vina config files with correct binding box |
| `docs/STRUCTURAL_AUDIT_FINDINGS.md` | Full documentation of the discrepancy and corrective steps |

---

## Next 3 Actions (in order)

### 1. Generate Correct Mutants in ChimeraX

Open ChimeraX → File > Open → select `scripts/generate_mutants.cxc`.

This runs `swapaa` on the true WT to produce `5i6z_A_S348T_correct.pdb`.

### 2. Prepare the Ligand PDBQT

See `output/LIGAND_PREP_INSTRUCTIONS.txt` for full options. Easiest path:

1. Download OpenBabel from [openbabel.org](https://openbabel.org)
2. Run:

```bash
obabel output/escitalopram.pdb -O output/escitalopram.pdbqt --partialcharge gasteiger -h
```

### 3. Install AutoDock Vina and Run Docking

1. Download `vina.exe` from [github.com/ccsb-scripps/AutoDock-Vina/releases](https://github.com/ccsb-scripps/AutoDock-Vina/releases)
2. Place it in `scripts/`
3. Run:

```bash
python scripts/docking_prep.py
```

The script will auto-run both WT and S348T docking and print the binding affinity comparison (kcal/mol).

---

## Key Technical Details

- **Binding pocket centroid** (S1 site, from Coleman et al. 2016 residues): x=33.1, y=187.3, z=141.0
- **Docking box size**: 25 x 25 x 25 Å
- **True WT residue 348**: ALA (canonical in 5I6Z construct; no SEQADV entry)
- **S348T in 5I6Z context**: technically A348T (Ala → Thr)
- **5I6Z engineered mutations**: A291 (canonical ILE), S439 (canonical THR), A554 (canonical CYS), A580 (canonical CYS)

For full details on the structural discrepancy, see `docs/STRUCTURAL_AUDIT_FINDINGS.md`.
