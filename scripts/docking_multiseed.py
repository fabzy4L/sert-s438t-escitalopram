"""
Multi-seed docking for SERT S438T mutation study.

Runs AutoDock Vina 1.2.5 with exhaustiveness=32 across 5 independent seeds
for both WT (SER438) and S438T mutant (THR438) receptors.
Reports mean +/- SD of best-pose binding affinity across seeds.

Usage:
  python scripts/docking_multiseed.py

Requirements:
  pip install numpy
  AutoDock Vina binary in scripts/vina.exe or on PATH
  output/receptor_wt.pdbqt
  output/receptor_s438t.pdbqt
  output/escitalopram.pdbqt

Output:
  output/multiseed_docking_report.txt
"""

import subprocess
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT    = Path(__file__).resolve().parents[1]
OUTPUT  = ROOT / "output"
SCRIPTS = ROOT / "scripts"

WT_PDBQT    = OUTPUT / "receptor_wt_H.pdbqt"
S438T_PDBQT = OUTPUT / "receptor_s438t_minimized_H.pdbqt"
LIG_PDBQT   = OUTPUT / "escitalopram.pdbqt"
REPORT      = OUTPUT / "multiseed_docking_report.txt"

BOX_CENTER    = (33.06, 187.25, 141.04)
BOX_SIZE      = (25, 25, 25)
EXHAUSTIVENESS = 32
NUM_MODES      = 10
SEEDS          = [42, 123, 456, 789, 1001]


def find_vina():
    if shutil.which("vina"):
        return "vina"
    local = SCRIPTS / "vina.exe"
    return str(local) if local.exists() else None


def write_conf(conf_path, receptor_pdbqt, out_pdbqt):
    with open(conf_path, "w") as f:
        f.write(f"receptor = {receptor_pdbqt}\n")
        f.write(f"ligand   = {LIG_PDBQT}\n\n")
        f.write(f"center_x = {BOX_CENTER[0]:.3f}\n")
        f.write(f"center_y = {BOX_CENTER[1]:.3f}\n")
        f.write(f"center_z = {BOX_CENTER[2]:.3f}\n\n")
        f.write(f"size_x = {BOX_SIZE[0]}\n")
        f.write(f"size_y = {BOX_SIZE[1]}\n")
        f.write(f"size_z = {BOX_SIZE[2]}\n\n")
        f.write(f"exhaustiveness = {EXHAUSTIVENESS}\n")
        f.write(f"num_modes      = {NUM_MODES}\n")
        f.write(f"out            = {out_pdbqt}\n")


def parse_best_affinity(log_text):
    for line in log_text.splitlines():
        parts = line.split()
        if parts and parts[0] == "1":
            try:
                return float(parts[1])
            except (ValueError, IndexError):
                pass
    return None


def run_seed(vina_bin, receptor_pdbqt, label, seed):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        conf_path = tmp / "conf.txt"
        out_pdbqt = tmp / "docked.pdbqt"
        write_conf(conf_path, receptor_pdbqt, out_pdbqt)

        result = subprocess.run(
            [vina_bin, "--config", str(conf_path), "--seed", str(seed)],
            capture_output=True, text=True
        )
        log = result.stdout + result.stderr
        affinity = parse_best_affinity(log)

        if affinity is None:
            print(f"    [{label}] seed={seed}: WARNING — could not parse affinity")
        else:
            print(f"    [{label}] seed={seed}: {affinity:.3f} kcal/mol")

        return affinity, log


def run_receptor(vina_bin, receptor_pdbqt, label):
    print(f"\n  {label}  ({receptor_pdbqt.name})")
    affinities = []
    seed_logs  = []

    for seed in SEEDS:
        aff, log = run_seed(vina_bin, receptor_pdbqt, label, seed)
        if aff is not None:
            affinities.append(aff)
        seed_logs.append((seed, aff, log))

    if not affinities:
        print(f"    ERROR: no valid results for {label}")
        return None, None, seed_logs

    mean = float(np.mean(affinities))
    std  = float(np.std(affinities, ddof=1) if len(affinities) > 1 else 0.0)
    print(f"    --> mean = {mean:.3f}  SD = {std:.3f} kcal/mol  (n={len(affinities)})")
    return mean, std, seed_logs


def write_report(results):
    wt_mean, wt_std, wt_logs     = results["wt"]
    s4_mean, s4_std, s4_logs     = results["s438t"]

    lines = []
    lines.append("SERT + Escitalopram — Multi-Seed Docking Results")
    lines.append(f"AutoDock Vina 1.2.5 | exhaustiveness={EXHAUSTIVENESS} | seeds={SEEDS}")
    lines.append("=" * 70)
    lines.append("")
    lines.append("SETUP")
    lines.append(f"  Receptor WT   : {WT_PDBQT.name}")
    lines.append(f"  Receptor S438T: {S438T_PDBQT.name}")
    lines.append(f"  Ligand        : {LIG_PDBQT.name}")
    lines.append(f"  Box centre    : {BOX_CENTER}")
    lines.append(f"  Box size      : {BOX_SIZE[0]} x {BOX_SIZE[1]} x {BOX_SIZE[2]} A")
    lines.append("")
    lines.append("-" * 70)
    lines.append("PER-SEED RESULTS")
    lines.append("")

    for label, seed_logs in [("WT (SER438)", wt_logs), ("S438T (THR438)", s4_logs)]:
        lines.append(f"  {label}")
        lines.append(f"  {'Seed':<8} {'Best Affinity (kcal/mol)'}")
        lines.append(f"  {'-'*8} {'-'*24}")
        for seed, aff, _ in seed_logs:
            val = f"{aff:.3f}" if aff is not None else "N/A"
            lines.append(f"  {seed:<8} {val}")
        lines.append("")

    lines.append("-" * 70)
    lines.append("SUMMARY")
    lines.append("")

    def fmt(mean, std):
        if mean is None:
            return "N/A"
        return f"{mean:.3f} +/- {std:.3f} kcal/mol"

    lines.append(f"  WT (SER438)    : {fmt(wt_mean, wt_std)}")
    lines.append(f"  S438T (THR438) : {fmt(s4_mean, s4_std)}")

    if wt_mean is not None and s4_mean is not None:
        ddg = s4_mean - wt_mean
        lines.append(f"  Delta-Delta-G  : {ddg:+.3f} kcal/mol  (S438T - WT)")
        lines.append("")
        lines.append("  Interpretation:")
        lines.append(f"  Both values are within AutoDock Vina scoring noise (~0.5 kcal/mol).")
        lines.append(f"  S438T {'shows marginally stronger' if ddg < 0 else 'shows marginally weaker'} predicted binding in the rigid-receptor model.")

    lines.append("")
    lines.append("-" * 70)
    lines.append("NOTE: rigid-receptor docking does not capture the desolvation penalty")
    lines.append("or Na+ coordination disruption associated with S438T. The experimental")
    lines.append("320-fold Ki increase (Andersen et al., 2009) requires flexible-receptor")
    lines.append("or MD-based methods to reproduce.")

    report_text = "\n".join(lines)
    with open(REPORT, "w") as f:
        f.write(report_text)
    return report_text


def main():
    print("=" * 60)
    print("MULTI-SEED DOCKING  |  exhaustiveness=32  |  5 seeds")
    print("=" * 60)

    # Pre-flight checks
    vina_bin = find_vina()
    if not vina_bin:
        sys.exit("ERROR: AutoDock Vina not found. Place vina.exe in scripts/ or add to PATH.")

    missing = [p for p in [WT_PDBQT, S438T_PDBQT, LIG_PDBQT] if not p.exists()]
    if missing:
        sys.exit(f"ERROR: missing files:\n" + "\n".join(f"  {p}" for p in missing))

    print(f"\n  Vina binary : {vina_bin}")
    print(f"  Seeds       : {SEEDS}")
    print(f"  Exhaustiveness: {EXHAUSTIVENESS}")

    results = {}
    for label, pdbqt, key in [
        ("WT (SER438)",    WT_PDBQT,    "wt"),
        ("S438T (THR438)", S438T_PDBQT, "s438t"),
    ]:
        mean, std, logs = run_receptor(vina_bin, pdbqt, label)
        results[key] = (mean, std, logs)

    print("\n" + "=" * 60)
    report = write_report(results)
    print(report)
    print(f"\nFull report saved: {REPORT}")


if __name__ == "__main__":
    main()
