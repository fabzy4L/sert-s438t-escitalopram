# Step 3 — ChimeraX Pocket Minimization (S438T)

Improves S438T receptor geometry before final publication-grade docking.

---

## Manual Steps in ChimeraX

1. Open **ChimeraX**
2. **File > Open** → navigate to and select:
   ```
   scripts/S438T_Minimization_Workflow.cxc
   ```
3. Wait for minimization to complete
4. Confirm output file exists:
   ```
   data/structures/5i6z_A_S438T_minimized.pdb
   ```

---

## After Minimization — Convert to Polar-H PDBQT

Run in PowerShell from the project root:

```powershell
$env:PATH += ";C:\Program Files\OpenBabel-3.1.1"
obabel data/structures/5i6z_A_S438T_minimized.pdb -O output/receptor_s438t_minimized_H.pdbqt -xr -h
```

---

## Update Docking Script

Edit `scripts/docking_multiseed.py` line 35:

```python
# Change FROM:
S438T_PDBQT = OUTPUT / "receptor_s438t_H.pdbqt"

# Change TO:
S438T_PDBQT = OUTPUT / "receptor_s438t_minimized_H.pdbqt"
```

---

## Re-run Multi-seed Docking

```powershell
python scripts/docking_multiseed.py
```

Output: `output/multiseed_docking_report.txt`

---

## Notes

- This step is optional but recommended — the minimized structure better captures the THR438 rotamer geometry
- Runtime: ~30–60 min (5 seeds × exhaustiveness=32)
- After docking completes, update `docs/research/S438T_RESEARCH_ARTICLE.md` Section 3.1 with the new mean ± SD values
