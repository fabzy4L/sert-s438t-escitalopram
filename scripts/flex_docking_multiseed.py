"""
Flexible-receptor multi-seed docking for SERT S438T mutation study --
docs/research/FUTURE_WORK.md Phase 1 ("Flexible Receptor Docking"),
finally run for real.

Rigid docking (docking_multiseed.py) could not reproduce the
experimental 320-fold Ki increase from S438T: DeltaDeltaG = +0.269
kcal/mol, within Vina's ~0.5 kcal/mol noise floor. This reruns the same
WT vs S438T comparison with the 4 S1 pocket residues (Y95, D98,
T438/S438, S439) treated as flexible, to test whether letting the
pocket move captures any of the induced-fit/Na+-coordination effect
rigid docking misses.

Uses the AutoDock Vina Python bindings (rigid_pdbqt_filename +
flex_pdbqt_filename) rather than this repo's usual CLI-binary +
config-file approach (docking_multiseed.py) -- no standalone `vina`
binary is available in the environment this was actually run in, only
the vina Python package. Same box, same exhaustiveness (32), same 5
seeds as docking_multiseed.py, for direct comparability.

Requires the rigid/flex receptor splits from prep_flexreceptor.py
(gitignored build artifacts, regenerate if missing):

    ~/ADFRsuite-1.1dev/bin/pythonsh scripts/prep_flexreceptor.py \
        --receptor output/receptor_wt_H.pdbqt \
        --residues 95,98,438,439 --out-dir output
    ~/ADFRsuite-1.1dev/bin/pythonsh scripts/prep_flexreceptor.py \
        --receptor output/receptor_s438t_minimized_H.pdbqt \
        --residues 95,98,438,439 --out-dir output

Usage (run in the vina micromamba env / any env with the vina package):
    python scripts/flex_docking_multiseed.py

Output:
  output/flex_multiseed_docking_report.txt
  output/flex_multiseed_docking_summary.json
"""
import json
from pathlib import Path

import numpy as np
from vina import Vina

ROOT   = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output"

WT_RIGID    = OUTPUT / "receptor_wt_H_rigid.pdbqt"
WT_FLEX     = OUTPUT / "receptor_wt_H_flex.pdbqt"
S438T_RIGID = OUTPUT / "receptor_s438t_minimized_H_rigid.pdbqt"
S438T_FLEX  = OUTPUT / "receptor_s438t_minimized_H_flex.pdbqt"
LIG_PDBQT   = OUTPUT / "escitalopram.pdbqt"
REPORT      = OUTPUT / "flex_multiseed_docking_report.txt"
SUMMARY     = OUTPUT / "flex_multiseed_docking_summary.json"

BOX_CENTER     = (33.06, 187.25, 141.04)
BOX_SIZE       = (25, 25, 25)  # unchanged from the rigid screen -- verified
                                # the 4 flex residues (Y95, D98, T438/S438,
                                # S439) all stay within this box's reach
                                # (max deviation ~7.4 A on any axis vs a
                                # 12.5 A half-width), unlike psilocin-5ht2a's
                                # tighter box which had to be enlarged.
EXHAUSTIVENESS = 32
NUM_MODES      = 10
SEEDS          = [42, 123, 456, 789, 1001]


def run_seed(rigid_pdbqt: str, flex_pdbqt: str, label: str, seed: int) -> float:
    v = Vina(sf_name="vina", seed=seed)
    v.set_receptor(rigid_pdbqt_filename=rigid_pdbqt, flex_pdbqt_filename=flex_pdbqt)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)
    v.set_ligand_from_file(str(LIG_PDBQT))
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=NUM_MODES)
    affinity = float(v.energies(n_poses=NUM_MODES)[0][0])
    print(f"    [{label}] seed={seed}: {affinity:.3f} kcal/mol")
    return affinity


def run_receptor(rigid_pdbqt: Path, flex_pdbqt: Path, label: str) -> tuple:
    print(f"\n  {label}  ({rigid_pdbqt.name} + {flex_pdbqt.name})")
    affinities = [run_seed(str(rigid_pdbqt), str(flex_pdbqt), label, seed) for seed in SEEDS]
    mean = float(np.mean(affinities))
    std = float(np.std(affinities, ddof=1)) if len(affinities) > 1 else 0.0
    print(f"    --> mean = {mean:.3f}  SD = {std:.3f} kcal/mol  (n={len(affinities)})")
    return mean, std, affinities



# Rigid-screen reference values, from output/multiseed_docking_report.txt
# (docking_multiseed.py has no JSON summary to parse programmatically --
# these are the exact numbers reported there, WT -8.575 +/- 0.013,
# S438T -8.306 +/- 0.015, n=5 seeds, same box/exhaustiveness/seed set).
RIGID_WT_MEAN, RIGID_WT_STD = -8.575, 0.013
RIGID_S438T_MEAN, RIGID_S438T_STD = -8.306, 0.015
RIGID_DDG = RIGID_S438T_MEAN - RIGID_WT_MEAN  # +0.269


def write_report(wt: tuple, s438t: tuple) -> str:
    wt_mean, wt_std, wt_scores = wt
    s4_mean, s4_std, s4_scores = s438t
    ddg = s4_mean - wt_mean
    rigid_ddg = RIGID_DDG

    lines = []
    lines.append("SERT + Escitalopram -- Flexible-Receptor Multi-Seed Docking Results")
    lines.append(f"AutoDock Vina (python bindings) | exhaustiveness={EXHAUSTIVENESS} | "
                  f"seeds={SEEDS} | flex residues: Y95, D98, T438/S438, S439")
    lines.append("=" * 70)
    lines.append("")
    lines.append("SETUP")
    lines.append(f"  Receptor WT   : {WT_RIGID.name} + {WT_FLEX.name}")
    lines.append(f"  Receptor S438T: {S438T_RIGID.name} + {S438T_FLEX.name}")
    lines.append(f"  Ligand        : {LIG_PDBQT.name}")
    lines.append(f"  Box centre    : {BOX_CENTER}")
    lines.append(f"  Box size      : {BOX_SIZE[0]} x {BOX_SIZE[1]} x {BOX_SIZE[2]} A "
                  f"(unchanged -- verified to already cover the flex residues' reach)")
    lines.append("")
    lines.append("-" * 70)
    lines.append("PER-SEED RESULTS")
    lines.append("")
    for label, scores in [("WT (SER438)", wt_scores), ("S438T (THR438)", s4_scores)]:
        lines.append(f"  {label}")
        lines.append(f"  {'Seed':<8} {'Best Affinity (kcal/mol)'}")
        lines.append(f"  {'-'*8} {'-'*24}")
        for seed, aff in zip(SEEDS, scores):
            lines.append(f"  {seed:<8} {aff:.3f}")
        lines.append("")

    lines.append("-" * 70)
    lines.append("SUMMARY")
    lines.append("")
    lines.append(f"  WT (SER438)    : {wt_mean:.3f} +/- {wt_std:.3f} kcal/mol")
    lines.append(f"  S438T (THR438) : {s4_mean:.3f} +/- {s4_std:.3f} kcal/mol")
    lines.append(f"  Delta-Delta-G  : {ddg:+.3f} kcal/mol  (S438T - WT, flexible)")
    lines.append(f"  Rigid Delta-Delta-G (docking_multiseed.py, for comparison): "
                  f"{rigid_ddg:+.3f} kcal/mol")
    lines.append(f"  Experimental Delta-Delta-G (320-fold Ki increase, Andersen et al. "
                  f"2009): ~+3.4 kcal/mol (RT*ln(320) at 37C)")
    lines.append("")
    lines.append("  Interpretation:")
    ddg_shift = ddg - rigid_ddg
    # combined SD across all four legs (WT rigid, S438T rigid, WT flex,
    # S438T flex) -- the right comparison for "did flexibility change
    # DeltaDeltaG", not just the flexible run's own internal noise.
    shift_combined_sd = float(np.sqrt(
        wt_std**2 + s4_std**2 + RIGID_WT_STD**2 + RIGID_S438T_STD**2))
    lines.append(f"  DeltaDeltaG shift (flexible - rigid): {ddg_shift:+.3f} kcal/mol; "
                  f"combined SD across both runs: {shift_combined_sd:.3f} kcal/mol.")
    if abs(ddg_shift) < 2 * shift_combined_sd:
        lines.append(f"  That shift is within ~2x the combined SD -- flexibility at these "
                      f"4 residues did not produce a result distinguishable from the rigid "
                      f"screen's noise. Same conclusion either way: {abs(ddg):.3f} kcal/mol "
                      f"is far short of the ~3.4 kcal/mol experimental value.")
    else:
        direction = "smaller (moving away from the experimental gap)" if ddg_shift < 0 else \
                    "larger (moving toward the experimental gap)"
        lines.append(f"  That shift ({abs(ddg_shift):.3f}) exceeds 2x the combined SD -- a "
                      f"real, reproducible change, not noise. DeltaDeltaG got {direction}. "
                      f"Either way, {abs(ddg):.3f} kcal/mol remains far short of the ~3.4 "
                      f"kcal/mol experimental value -- this is a real effect of essentially "
                      f"no practical size.")
    lines.append("")
    lines.append("-" * 70)
    lines.append("NOTE: flexibility at these 4 residues alone does not model explicit")
    lines.append("Na+ coordination chemistry or solvation/desolvation -- both flagged in")
    lines.append("FUTURE_WORK.md as the more likely mechanism (conformational dynamics")
    lines.append("during the binding pathway, not the final bound-state geometry this")
    lines.append("docking approach -- rigid or flexible -- can ever capture). MD (Phase")
    lines.append("2/4) remains the next real test, not further docking variants.")

    report_text = "\n".join(lines)
    REPORT.write_text(report_text)
    return report_text


def main() -> None:
    print("=" * 60)
    print("FLEXIBLE-RECEPTOR MULTI-SEED DOCKING  |  exhaustiveness=32  |  5 seeds")
    print("=" * 60)

    missing = [p for p in [WT_RIGID, WT_FLEX, S438T_RIGID, S438T_FLEX, LIG_PDBQT] if not p.exists()]
    if missing:
        raise SystemExit("ERROR: missing files:\n" + "\n".join(f"  {p}" for p in missing))

    wt = run_receptor(WT_RIGID, WT_FLEX, "WT (SER438)")
    s438t = run_receptor(S438T_RIGID, S438T_FLEX, "S438T (THR438)")

    SUMMARY.write_text(json.dumps({
        "seeds": SEEDS, "exhaustiveness": EXHAUSTIVENESS, "box_size_angstrom": BOX_SIZE,
        "wt": {"mean_kcal_mol": wt[0], "std_kcal_mol": wt[1], "scores": wt[2]},
        "s438t": {"mean_kcal_mol": s438t[0], "std_kcal_mol": s438t[1], "scores": s438t[2]},
        "ddg_kcal_mol": s438t[0] - wt[0],
    }, indent=2))

    print("\n" + "=" * 60)
    report = write_report(wt, s438t)
    print(report)
    print(f"\nReport:  {REPORT}")
    print(f"Summary: {SUMMARY}")


if __name__ == "__main__":
    main()
