# Phase 4 (MD) Methodology

**Scope:** GROMACS system building, equilibration, production, and
comparative analysis for the WT+escitalopram and S438T+escitalopram
complexes. Docking (Phases 1-3) methodology and the "Current Finding"
this phase exists to explain live in `docs/research/FUTURE_WORK.md`; this
document picks up where that one leaves off, once CHARMM-GUI/CGenFF
exports exist for both systems.

**Status: scaffolded, not yet executed.** No CHARMM-GUI export, no CGenFF
parameterization of escitalopram, and no working GROMACS build exist in
this repo yet. Everything below describes what `md/scripts/` does once
those exist, not results that exist today. This mirrors the sibling
project `psilocin-5ht2a`'s `md/scripts/` scaffold and
`docs/METHODOLOGY.md`, adapted for two parallel systems instead of one.

## 1. Why two systems change the pipeline shape

`psilocin-5ht2a` runs one receptor through one MD pipeline. This repo's
whole point is a WT-vs-S438T *comparison* (Na+ coordination distance,
RMSF, pocket volume — `FUTURE_WORK.md`'s Phase 4 "Key analyses"), so every
staged script takes `--system wt|s438t` and writes into a
system-specific subdirectory (`md/system_prep/<system>/`,
`md/analysis/<system>/`) rather than one shared path. A final stage,
`06_compare_systems.py`, reads both systems' `05_analysis.py` output and
produces the actual comparison — it's the only stage that touches both
systems at once.

```
CHARMM-GUI export (WT)         CHARMM-GUI export (S438T)
       |                              |
       v                              v
01_import_charmm_gui.py       01_import_charmm_gui.py
  --system wt                   --system s438t
       |                              |
       v                              v
02_minimize.py --system wt    02_minimize.py --system s438t
       |                              |
       v                              v
03_equilibrate.py (nvt, npt)  03_equilibrate.py (nvt, npt)
  --system wt                   --system s438t
       |                              |
       v                              v
04_production.py --system wt  04_production.py --system s438t
       |                              |
       v                              v
05_analysis.py --system wt    05_analysis.py --system s438t
       |                              |
       +-------------+  +------------+
                     v  v
              06_compare_systems.py
```

## 2. Design pattern: every automatable stage is a gate

Same convention as `scripts/validation_redock.py` (this repo's docking
validation, `--strict`, machine-readable `gate_status.json`) and
`psilocin-5ht2a/md/scripts/`: each stage writes a human-readable
`*_report.txt`, a `gate_status.json` (`{"gate": "PASS"|"FAIL", ...
criteria}`), and accepts `--strict` to exit non-zero on FAIL so a stage
can't silently feed a bad result to the next one — for BOTH systems
independently, since a bad WT run shouldn't be masked by a good S438T run
or vice versa.

That convention exists *because* of a real incident, not preemptively:
`scripts/validation_redock.py` reported a redock FAIL for a while that
turned out to be a parsing bug (it graded the worst of 5 Vina poses
instead of the best one, fixed in commit `3d18fcd` — the same bug and fix
as `psilocin-5ht2a/docs/REDOCK_BUG_HANDOFF.md`, this repo's script was
ported from that one). Direct consequences here:

1. **The report/gate plumbing is factored into `md/scripts/gate_utils.py`**
   rather than reimplemented per script or per system.
2. **Every gate has a real unit test** (`tests/`) against a synthetic
   fixture — a fake `gmx` log, a fake `.xvg`, a fake fpocket `_info.txt`,
   a synthetic RMSF/coordination series — not just "the script ran
   without an exception."
3. **`05_analysis.py`'s residue-438 atom name is system-aware, not
   shared.** Ser438 (WT) has a hydroxyl oxygen named `OG`; Thr438
   (S438T) has the analogous oxygen named `OG1` (CHARMM36 naming — Thr
   has an extra methyl, `CG2`, that shifts the numbering). Defaulting to
   one shared atom name for both systems would silently find nothing for
   whichever system it's wrong for — exactly the "plausible-looking
   default silently selects the wrong thing" failure class the redock
   bug was. `RESID_438_ATOM_DEFAULT` is a dict keyed by system, and
   `tests/test_05_analysis.py` pins the two values as required to differ.
4. **`06_compare_systems.py` matches RMSF by resid, not array position**
   (`compare_rmsf()`), since the two systems' `05_analysis.py` runs
   aren't guaranteed to enumerate residues in the same order — matching
   by position would risk comparing the wrong pair of residues.

## 3. Stage details

### 3.1 Import (`01_import_charmm_gui.py`)

Validates one system's CHARMM-GUI Membrane Builder download plus its
CGenFF ligand parameters for escitalopram (`FUTURE_WORK.md` Phase 4 step
1: "Parameterize escitalopram with CGenFF"). Beyond the base
topology/index/coordinates/`.mdp` checks (same as `psilocin-5ht2a`'s
version), this also requires a CGenFF `.str` stream file (or split
`.rtf`+`.prm`) — a missing one resolves everything else fine and then
fails deep inside `gmx grompp` on an unknown residue, which is a much
worse place to find out than here.

### 3.2 Minimization (`02_minimize.py`)

Identical gate logic to `psilocin-5ht2a`'s version (Fmax convergence
parsed from the `mdrun` log, default threshold 1000 kJ/mol/nm) — this
part of the pipeline doesn't differ between a GPCR and a transporter, or
between WT and a point mutant. `--system` only changes the default output
directory and the report/gate_status.json's `system` field.

### 3.3 Equilibration (`03_equilibrate.py`)

Same NVT-then-NPT, tail-stability-gated design as `psilocin-5ht2a`
(temperature 2 K, pressure 50 bar, density 5 kg/m³ default tolerances
over the last 100 ps). Run once per system per phase.

### 3.4 Production (`04_production.py`)

Long GPU-offloaded `gmx mdrun`, 100-200 ns per `FUTURE_WORK.md`. Not
pass/fail gated (same reasoning as `psilocin-5ht2a`: there's no single
"production succeeded" number). SERT is a full membrane transporter
(~150k atoms per `FUTURE_WORK.md`'s own estimate, larger than
`psilocin-5ht2a`'s GPCR system) — expect roughly 10-20 ns/day on the
scoped RTX 5060 8GB, meaningfully slower than the sibling project.
`--smoke-test` exists for the same GPU-utilization sanity check before
committing both systems to a multi-day run each.

### 3.5 Per-system analysis (`05_analysis.py`)

Uses MDAnalysis to compute, per system:

- **Pocket-residue RMSF** (default resids 95, 98, 438, 439 — Y95, D98,
  [S/T]438, S439, the same four residues `FUTURE_WORK.md` Phase 1 already
  flags as core pharmacophoric contacts), aligned to the pocket backbone
  first so whole-protein drift over a 100-200 ns membrane-protein run
  doesn't inflate the numbers.
- **Na+ coordination distance** — nearest `SOD` ion (CHARMM36/CGenFF's
  sodium residue name, not the generic PDB `NA`) to residue 438's
  reference hydroxyl atom, per frame. This is the direct structural test
  of `FUTURE_WORK.md`'s central hypothesis: "the mutation most likely
  disrupts the binding pathway (Na⁺ coordination, induced-fit access)
  rather than the final bound pose."
- **Pocket volume** (fpocket, best-effort — `SKIPPED` in the report if
  fpocket isn't installed, not a hard failure) on a supplied single-frame
  PDB. `select_largest_pocket()` picks the largest pocket by volume as a
  default S1-site proxy; its docstring flags this needs verifying against
  residue proximity once a real trajectory exists, not trusted on rank
  alone.

Gated on trajectory completeness (`trajectory_sanity()`: enough frames,
non-zero), **not** a physical pass/fail threshold — RMSF and Na+
coordination distance are this project's actual publication results, not
a correctness check on the pipeline (contrast with `psilocin-5ht2a`'s
`05_analysis.py`, where the ligand-RMSD-plateau gate IS a correctness
check on pose stability).

### 3.6 System comparison (`06_compare_systems.py`)

The only stage in this phase without an equivalent in `psilocin-5ht2a`
(that project has one system). Reads both systems' `05_analysis.py`
output (`rmsf_per_resid.csv`, `na_coordination_per_frame.csv`) and
produces:

- Per-residue RMSF delta (S438T − WT), matched by resid.
- Na+ coordination distance summary stats for each system plus the delta
  in means (not a time-aligned per-frame comparison — the two systems'
  trajectories are independent runs, not synchronized).

Gated on `inputs_compatible()` — both systems must have RMSF data over
the *same* residue set, or the comparison doesn't mean anything. This is
a compatibility gate, not a physical-result gate, same spirit as 3.5.

## 4. What's still a scaffold vs. what's real

Everything in `md/scripts/` is real, runnable code with real parsing/gate
logic (42 tests, `tests/`), but none of it has been run against real
GROMACS output — there isn't any yet. Specifically not done in this pass:

- No `gmx` build exists in this repo's tracked state.
- No CHARMM-GUI export or CGenFF parameterization has been generated or
  imported for either system.
- Water occupancy analysis (`FUTURE_WORK.md`'s Phase 4 list) has no
  scaffold yet.
- MM-PBSA (`FUTURE_WORK.md` step 5, AmberTools) is out of scope for this
  GROMACS-focused pass — a separate phase, per `FUTURE_WORK.md`'s own
  Phase 3/4 split.
- Phase 5 (alchemical FEP) has no scaffold — stretch goal per
  `FUTURE_WORK.md`.
