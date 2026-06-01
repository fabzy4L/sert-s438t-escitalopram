# RECEPTOR_DESIGN — Project Gaps (Updated 2026-05-02)

The core computational pipeline investigates the **S438T** mutation in SERT and its effect on escitalopram binding affinity.

---

## Completed

- Structural audit and provenance correction (scripts/structural_audit.py)
- True WT extraction: data/structures/5i6z_A_true.pdb
- In-silico mutagenesis for S438T: data/structures/5i6z_A_S438T.pdb
- Ligand and receptor PDBQT generation for S438T
- AutoDock Vina 1.2.5 docking — WT and S438T, 10 poses each
- Research article: docs/research/S438T_RESEARCH_ARTICLE.md / .html
- Validation script: scripts/validation_redock.py (68P re-dock from 5I71)
- S1 pocket minimization script: scripts/S438T_Minimization_Workflow.cxc
- Multi-seed docking script: scripts/docking_multiseed.py (exhaustiveness=32, 5 seeds)
- Notebook paths verified — no hardcoded paths present
- 5I6Z_edited.fasta confirmed as WT crystal reference sequence (pipeline uses PDB directly)

---

## Remaining: Run the Scripts

### Step 1 — Validate box placement
```
python scripts/validation_redock.py
```
Pass criterion: 68P centroid inside grid; RMSD < 2.0 Å if OpenBabel available.
Output: output/68P_crystal_reference.pdb, output/validation_report.txt

### Step 2 — Minimize S1 pocket (ChimeraX)
Open ChimeraX → File > Open → scripts/S438T_Minimization_Workflow.cxc
Output: data/structures/5i6z_A_S438T_minimized.pdb
Then re-run docking_prep.py pointing at the minimized structure before Step 3.

### Step 3 — Re-generate receptor PDBQT with polar hydrogens
```
obabel output/receptor_s438t.pdbqt -O output/receptor_s438t_H.pdbqt -h
obabel output/receptor_wt.pdbqt    -O output/receptor_wt_H.pdbqt    -h
```
Update vina_wt.conf and vina_s438t.conf to reference the _H variants.

### Step 4 — Multi-seed docking (publication-grade)
```
python scripts/docking_multiseed.py
```
Runs 5 seeds × 2 receptors at exhaustiveness=32. Reports mean ± SD.
Output: output/multiseed_docking_report.txt

---

## Phase 2 (Future Work)

### Canonical Biological WT
The 5I6Z construct has engineered mutations (A291I, S439T, A554C, A580C).
Build a canonical SERT model from UniProt P31645 via homology modeling
for a clinically relevant baseline.

### Extended Analysis
- MM-GBSA re-scoring of top docked poses
- 100 ns MD simulation (WT vs S438T) in explicit POPC lipid bilayer to capture
  the desolvation and kinetic effects — the only method that can reproduce
  Andersen's (2009) 320-fold Ki increase.
