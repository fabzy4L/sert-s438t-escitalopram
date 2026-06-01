# SERT S438T — Publication Roadmap

**Current baseline:** Multi-seed rigid docking (AutoDock Vina 1.2.5, n=5, exhaustiveness=32)
- WT: −8.575 ± 0.013 kcal/mol
- S438T (minimized): −8.306 ± 0.015 kcal/mol
- ΔΔG: +0.269 kcal/mol (within Vina noise; cannot reproduce experimental 320-fold Ki shift)

Each phase below builds directly on the previous one.

---

## Phase 1 — Flexible Receptor Docking

**Goal:** Allow key S1 pocket residues to flex during docking, capturing induced-fit at THR438.

**Tools:** AutoDock Vina 1.2.5 (already installed), AutoDockTools or MGLTools (free)

**Steps:**
1. Define flexible residues: Y95, D98, T438, S439 (core pharmacophoric contacts)
2. Split receptor into rigid + flexible PDBQT using `prepare_flexreceptor4.py` (MGLTools)
3. Run multi-seed docking with `--flex receptor_flex.pdbqt` flag
4. Compare ΔΔG to rigid baseline

**Expected improvement:** Pocket can accommodate THR methyl group by minor side-chain adjustment — expected to increase ΔΔG signal modestly.

**Output:** `output/flexdock_report.txt`, updated Section 3.1 in article

---

## Phase 2 — Ensemble Receptor Docking

**Goal:** Replace single crystal structure with a conformational ensemble to account for receptor dynamics.

**Tools:** GROMACS or OpenMM (free), existing Python pipeline

**Steps:**
1. Run short (~50 ns) MD on apo WT and apo S438T receptors (explicit solvent, Na⁺/Cl⁻)
2. Cluster trajectory → extract 5–10 representative conformations per receptor
3. Dock escitalopram into each cluster representative (multi-seed)
4. Report ensemble-averaged ΔΔG ± SD

**Expected improvement:** Captures pocket breathing and Na⁺ coordination geometry — the primary mechanism per Andersen (2009).

**Output:** Ensemble docking report, RMSD cluster figure, Na⁺ coordination distance analysis

---

## Phase 3 — MM-GBSA Rescoring

**Goal:** Rescore best docked poses with a more physically rigorous free energy model including implicit solvation.

**Tools:** AmberTools 23 (free), `MMPBSA.py`

**Steps:**
1. Take top-ranked docked pose from Phase 1 for each receptor
2. Run short restrained MD (10 ns) on each complex in AmberTools
3. Compute MM-GBSA binding free energies with per-residue decomposition
4. Decompose ΔΔG by residue to identify which contacts drive the affinity change

**Expected improvement:** Captures desolvation penalty of THR438 methyl group; per-residue decomposition is a publication-standard figure.

**Output:** ΔG_bind table, per-residue decomposition heatmap (key figure for paper)

---

## Phase 4 — MD Simulation + MM-PBSA (Primary Publication Result)

**Goal:** Full dynamic treatment of both complexes; the minimum standard for a computational ΔΔG paper.

**Tools:** GROMACS 2024 (free, GPU-accelerated), CHARMM36m force field, CGenFF for escitalopram

**Steps:**
1. Parameterize escitalopram with CGenFF (online server, free)
2. Embed each receptor in POPC lipid bilayer (CHARMM-GUI membrane builder)
3. Add explicit TIP3P water, 150 mM NaCl
4. Run 100–200 ns production MD on WT+ESC and S438T+ESC complexes
5. MM-PBSA over last 50 ns trajectory
6. Analyze: Na⁺ coordination at S/T438, water occupancy, RMSF of pocket residues

**Expected improvement:** Directly captures the thermodynamic basis of the 320-fold Ki increase. This phase produces the primary result for publication.

**Key analyses for paper:**
- ΔΔG_bind (MM-PBSA) vs experimental ΔΔG from Ki ratio
- Na⁺ coordination distance time series (S438 vs T438)
- Pocket volume over trajectory (MDpocket or fpocket)
- RMSF comparison WT vs S438T

---

## Phase 5 — Alchemical Free Energy Perturbation (Optional, Gold Standard)

**Goal:** Compute ΔΔG via thermodynamic integration to directly compare with experimental Ki values.

**Tools:** GROMACS + `pmx` (free), or OpenFE (Python, free)

**Steps:**
1. Set up S438→T438 alchemical perturbation in the apo and holo (+ escitalopram) states
2. Run λ-windows (typically 11–20 windows, 5–10 ns each)
3. Compute ΔΔG via BAR or MBAR estimator
4. Compare to experimental ΔΔG = RT·ln(320) ≈ +3.4 kcal/mol (37°C)

**Expected improvement:** Quantitatively reproduce the experimental 320-fold Ki increase. Strongest possible computational result.

**Note:** Requires ~200–400 ns total simulation time; GPU cluster recommended.

---

## Figures Required for Submission

| Figure | Phase | Description |
|--------|-------|-------------|
| Fig 1 | Baseline | S1 pocket: WT vs S438T overlay, escitalopram pose |
| Fig 2 | Baseline | Docking affinity comparison (bar ± SD, n=5) |
| Fig 3 | Phase 2/4 | Na⁺ coordination distance time series |
| Fig 4 | Phase 3/4 | Per-residue MM-GBSA decomposition heatmap |
| Fig 5 | Phase 4 | Pocket volume / RMSF comparison |
| Fig 6 | Phase 5 | FEP ΔΔG vs experimental (correlation plot) |

Figures 1–2 are already producible from the current baseline.

---

## Target Journals

| Journal | IF | Scope fit |
|---|---|---|
| *ACS Chemical Neuroscience* | ~4.5 | SERT pharmacology + computation |
| *Journal of Chemical Information and Modeling* | ~5.6 | Computational methods focus |
| *PLOS Computational Biology* | ~4.3 | Open access, methods-forward |
| *Biochemistry* | ~3.5 | Structural/biophysical baseline study |

**Minimum viable manuscript:** Phases 1–3 completed (flexible docking + MM-GBSA). Targets *Biochemistry* or *PLOS Comp Bio*.

**Full study:** Phase 4 completed (MD + MM-PBSA). Targets *JCIM* or *ACS Chem Neuro*.

---

## Projected Timeline

Assumes part-time work (~10–15 hrs/week), Windows machine, consumer GPU (no HPC access).

| Phase | Work | Compute | Wall Time |
|-------|------|---------|-----------|
| **Phase 1** — Flexible docking | MGLTools setup, PDBQT prep, analysis | ~1 hr (Vina) | **1–2 weeks** |
| **Phase 2** — Ensemble docking | GROMACS/WSL2 setup, MD config | 2–4 days/system (consumer GPU) | **6–10 weeks** |
| **Phase 3** — MM-GBSA | AmberTools conda install, rescoring | 1–2 days | **2–3 weeks** |
| **Phase 4** — MD + MM-PBSA | CHARMM-GUI bilayer, equilibration, production | 3–6 weeks/system (consumer GPU) | **2–4 months** |
| **Manuscript draft** | Writing, figures, revision | — | **4–6 weeks** |

**Minimum viable paper (Phases 1–3 + writing): ~4–5 months**

**Full study (Phases 1–4 + writing): ~8–12 months**

**Primary bottleneck: compute.** SERT is a membrane protein — a full bilayer system (~150k atoms) runs ~10–20 ns/day on a consumer GPU. 200 ns production MD = 10–20 days of uninterrupted GPU time per complex (40+ days across both systems). HPC access would collapse Phase 4 from months to 2–3 weeks. Cloud HPC options (AWS, Google Cloud) can be cost-effective for burst compute when needed.

---

## Dependencies / Environment Notes

- GROMACS: install via WSL2 (Ubuntu) on Windows — significantly easier than native Windows build
- AmberTools: free via conda (`conda install -c conda-forge ambertools`)
- CGenFF: free web server (https://cgenff.umaryland.edu)
- CHARMM-GUI: free web server for membrane builder
- OpenFE: `pip install openfe` (Python FEP framework, well-documented)
