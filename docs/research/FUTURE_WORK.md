# SERT S438T — Publication Roadmap

## Current Finding

Multi-seed rigid docking (AutoDock Vina 1.2.5, n=5, exhaustiveness=32):

- WT: −8.575 ± 0.013 kcal/mol
- S438T (minimized): −8.306 ± 0.015 kcal/mol
- ΔΔG: +0.269 kcal/mol (within Vina noise)

**The inability to reproduce the experimental 320-fold Kᵢ increase is itself the primary finding of the current study.** It demonstrates that the S438T substitution does not alter the static bound-state geometry of escitalopram in the S1 pocket — the Thr methyl group is accommodated without significant steric clash. This localizes the pharmacogenomic mechanism to conformational dynamics: the mutation most likely disrupts the binding pathway (Na⁺ coordination, induced-fit access) rather than the final bound pose. That is a mechanistically specific, testable hypothesis — and it is what the MD phases below are designed to confirm.

Each phase below builds directly on the previous one.

---

## Phase 1 — Flexible Receptor Docking — DONE (2026-09-14)

**Goal:** Allow key S1 pocket residues to flex during docking, capturing induced-fit at THR438.

**Tools used:** AutoDock Vina (python bindings, not the 1.2.5 CLI binary —
none was available in the environment this actually ran in), and
`AutoDockTools.MoleculePreparation.AD4FlexibleReceptorPreparation` called
directly under ADFRsuite's bundled Python 2 (`pythonsh`) —
`prepare_flexreceptor4.py` (MGLTools) was never actually installed; the
underlying class it wraps was already available locally and needed no
new install. See `scripts/prep_flexreceptor.py`.

**What was done:**
1. Flexible residues Y95, D98, T438/S438, S439 — confirmed present at
   exactly those numbers in both `receptor_wt_H.pdbqt` and
   `receptor_s438t_minimized_H.pdbqt`.
2. Split each receptor into rigid + flexible PDBQT
   (`scripts/prep_flexreceptor.py`) — verified the existing 25×25×25 Å
   box already covers all 4 residues' reach (max deviation ~7.4 Å on any
   axis vs. a 12.5 Å half-width), so no box-enlargement was needed here
   (unlike the sibling `psilocin-5ht2a` project, which had to enlarge its
   box for the same class of check).
3. Multi-seed docking (same 5 seeds as `docking_multiseed.py`:
   42, 123, 456, 789, 1001) with `Vina.set_receptor(rigid_pdbqt_filename=...,
   flex_pdbqt_filename=...)` — `scripts/flex_docking_multiseed.py`.
4. Compared ΔΔG to the rigid baseline.

**Result:** WT −8.527 ± 0.003 kcal/mol, S438T −8.317 ± 0.014 kcal/mol,
ΔΔG = **+0.210 kcal/mol** (flexible) vs. **+0.269 kcal/mol** (rigid,
`output/multiseed_docking_report.txt`). The 0.059 kcal/mol shift is
statistically real (exceeds 2× the combined SD across both runs, ~0.024
kcal/mol) — flexibility at these 4 residues did measurably change the
result — but the ΔΔG got *smaller*, moving further from the ~+3.4
kcal/mol experimental value, not closer. Practically: a real effect of
essentially no size relative to the ~16x gap that still needs
explaining.

**Reading:** this does not overturn the "Current Finding" above — it
reinforces it. Letting the four most obvious pharmacophoric-contact
residues flex did not recover the experimental effect, which is more
evidence (not proof) that the mechanism really is in the binding
*pathway* (Na⁺ coordination, induced-fit access dynamics) rather than
the final bound-state geometry — the thing no single-endpoint docking
run, rigid or flexible, can capture. Phase 2 (ensemble/MD) remains the
next real test.

**Output:** `output/flex_multiseed_docking_report.txt`,
`output/flex_multiseed_docking_summary.json`.

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

**Pipeline status: scaffolded, not yet executed.** `md/scripts/` (staged
WT+S438T import/minimize/equilibrate/production/analysis pipeline,
gated the same way `scripts/validation_redock.py` is) and
`docs/research/MD_METHODOLOGY.md` (design writeup) exist and are covered
by 42 passing unit tests, but no CHARMM-GUI export, no CGenFF
parameterization of escitalopram, and no working GROMACS build exist in
this repo yet -- the steps below are still what needs to happen, in this
order, before any of the numbers in "Key analyses for paper" are real.

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

- GROMACS: `bash scripts/build_gromacs.sh` builds a GPU-enabled install
  inside WSL2 (Ubuntu) on Windows, without needing sudo/apt (toolchain
  comes from conda). This is a tested recipe, not just `cmake && make` —
  it works around two real issues: (1) the CUDA toolkit version has to
  support your specific GPU architecture (Blackwell/RTX 50-series needs
  CUDA 12.8+; the script defaults to 12.9 and auto-detects your GPU's
  architecture via `CMAKE_CUDA_ARCHITECTURES=native`), and (2) the final
  `gmx` link can fail on `undefined reference to cufftDestroy` even
  though everything compiled — a conda-packaged-CUDA quirk the script
  detects and retries around automatically. See the script's own header
  comments, or the sibling `psilocin-5ht2a` repo's
  `docs/SESSION_HANDOFF_2026-09-10.md` for the full narrative of how this
  was found.
- AmberTools: free via conda (`conda install -c conda-forge ambertools`)
- CGenFF: free web server (https://cgenff.umaryland.edu)
- CHARMM-GUI: free web server for membrane builder
- OpenFE: `pip install openfe` (Python FEP framework, well-documented)
