# RECEPTOR_DESIGN — Project Review (Updated 2026-05-01)

## Project Overview

**Goal:** Structural investigation of the Serotonin Transporter (SERT) and how the **S438T** point mutation affects Escitalopram binding affinity. This study pivoted from an initial focus on A348T after identifying a digit-swap error in the site numbering.

**Authoritative output:** `docs/research_article.md` / `docs/research_article.html`

---

## Current Status

| Area | Status | Output |
|---|---|---|
| Wild-type structures (5I6Z, 5I71, 2A65) | Done | `data/structures/` |
| Structural audit & provenance correction | Done | `scripts/structural_audit.py`, `data/structures/5i6z_A_true.pdb` |
| **Mutation Correction (Pivot)** | **Done** | **A348T → S438T** |
| In-silico mutagenesis (S438T, g+ rotamer) | Done | `data/structures/5i6z_A_S438T.pdb` |
| Escitalopram 3D structure | Done | `output/escitalopram.pdb` |
| Ligand PDBQT preparation | Done | `output/escitalopram.pdbqt` |
| Receptor PDBQT preparation | Done | `output/receptor_wt.pdbqt`, `output/receptor_s438t.pdbqt` |
| Molecular docking (AutoDock Vina 1.2.5) | Done | `output/docked_wt.pdbqt`, `output/docked_s438t.pdbqt` |
| Research article | Done | `docs/research_article.md/.html` |

---

## Key Results (Pivot to S438T)

| Metric | Wild-Type (SER438) | Mutant (THR438) |
|---|---|---|
| Best docking affinity | −8.1 kcal/mol | −8.3 kcal/mol |
| ΔΔG | — | −0.2 kcal/mol |

**Conclusion:** The S438T mutation in the S1 pocket maintains high affinity for escitalopram in the rigid-receptor docking model. The predicted shift is minimal (−0.2 kcal/mol), suggesting that escitalopram may be less sensitive to the S438T steric clash than citalopram.

---

## Important Correction

This project originally investigated **A348T** due to a digit-swap error in the mutation site numbering (misreading 438 as 348). Position 348 is a peripheral site in TM6, whereas 438 is a primary binding pocket residue. All results for 348 have been superseded by the S438T analysis.

---

## Remaining Work

1. Re-run docking at exhaustiveness = 32 with multiple seeds
2. Regenerate ligand PDBQT with OpenBabel Gasteiger charges
3. Localized energy minimization of the S1 pocket with THR438
4. MD simulation (100 ns, explicit POPC bilayer) to capture dynamic effects
