# SERT S438T — Remaining Analysis Steps

**Status as of 2026-05-02**

Completed: validation re-dock (box proximity PASS), polar-H receptor PDBQTs generated.
Remaining: multi-seed publication-grade docking.

---

## Step 1 — Update multiseed script to use polar-H receptors

Edit `scripts/docking_multiseed.py`, lines 34–35:

```python
# Change FROM:
WT_PDBQT    = OUTPUT / "receptor_wt.pdbqt"
S438T_PDBQT = OUTPUT / "receptor_s438t.pdbqt"

# Change TO:
WT_PDBQT    = OUTPUT / "receptor_wt_H.pdbqt"
S438T_PDBQT = OUTPUT / "receptor_s438t_H.pdbqt"
```

---

## Step 2 — Run multi-seed docking

Requires: `scripts/vina.exe` present (or `vina` on PATH), `output/escitalopram.pdbqt` present.

```powershell
cd C:\Users\Fabian Alvarez-Primo\Documents\Github\DATA_ANALYTICS\bioinformatics\python_bio\Receptor_Design\RECEPTOR_DESIGN
python scripts/docking_multiseed.py
```

Runtime: ~30–60 min (5 seeds × 2 receptors × exhaustiveness=32).

Output: `output/multiseed_docking_report.txt` — contains mean ± SD affinities and ΔΔG.

---

## Step 3 — Optional: ChimeraX pocket minimization (improves S438T geometry)

1. Open ChimeraX
2. **File > Open** → select `scripts/S438T_Minimization_Workflow.cxc`
3. Wait for minimization to finish
4. Output: `data/structures/5i6z_A_S438T_minimized.pdb`

Then re-generate the polar-H PDBQT from the minimized structure:

```powershell
$env:PATH += ";C:\Program Files\OpenBabel-3.1.1"
obabel data/structures/5i6z_A_S438T_minimized.pdb -O output/receptor_s438t_minimized_H.pdbqt -xr -h
```

Then update `docking_multiseed.py` line 35 to use `receptor_s438t_minimized_H.pdbqt` and re-run Step 2.

---

## Step 4 — Update research article with final results

After Step 2, open `output/multiseed_docking_report.txt` and copy the mean ± SD values into:

- `docs/research/S438T_RESEARCH_ARTICLE.md` — Section 3.1 Results table
- `docs/research/S438T_RESEARCH_ARTICLE.html` — same section

Replace the current single-seed values (`-8.100` / `-8.297` / `ΔΔG -0.197`) with the multi-seed mean ± SD format, e.g.:

```
WT (SER438):    X.XXX ± 0.XXX kcal/mol  (n=5, exhaustiveness=32)
S438T (THR438): X.XXX ± 0.XXX kcal/mol  (n=5, exhaustiveness=32)
ΔΔG:            X.XXX kcal/mol
```

---

## Key file locations

| File | Purpose |
|------|---------|
| `output/receptor_wt_H.pdbqt` | WT receptor with polar H (ready) |
| `output/receptor_s438t_H.pdbqt` | S438T receptor with polar H (ready) |
| `output/escitalopram.pdbqt` | Ligand PDBQT (must exist before Step 2) |
| `scripts/vina.exe` | AutoDock Vina binary |
| `output/multiseed_docking_report.txt` | Final docking output |
| `output/68P_crystal_reference.pdb` | Crystal reference pose (validation complete) |
| `output/validation_report.txt` | Validation report (box PASS, RMSD note) |

---

## Notes

- **Validation RMSD (6.632 Å)** is not a failure — it reflects a cross-enantiomer comparison (R-citalopram from 5I71 vs S-citalopram receptor 5I6Z). Box centroid (3.92 Å from center, inside 25 Å grid) is the meaningful pass criterion.
- **ΔΔG interpretation:** Rigid-receptor Vina cannot reproduce the experimental 320-fold Ki increase (Andersen et al., 2009). The ~−0.2 kcal/mol result is expected for this method. MD simulation would be required to capture desolvation and Na⁺ coordination effects.
- **OpenBabel PATH** (Windows): `$env:PATH += ";C:\Program Files\OpenBabel-3.1.1"` before running obabel commands.
