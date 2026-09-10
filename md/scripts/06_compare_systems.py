"""
WT-vs-S438T comparison: the actual headline result docs/research/
FUTURE_WORK.md Phase 4 asks for (Fig 3: Na+ coordination distance time
series S438 vs T438; Fig 5: RMSF comparison). Reads both systems' output
from 05_analysis.py (md/analysis/wt/ and md/analysis/s438t/) rather than
loading trajectories again -- this script is pure post-processing over
already-computed per-system results.

Gated on the two systems' outputs actually being comparable (same pocket
residues, both have usable Na+ coordination series), NOT on any physical
threshold -- like 05_analysis.py, "gate" here means "is this comparison
meaningful," not "is the biology good."

Usage:
    python 06_compare_systems.py \
        --wt-dir ../analysis/wt --s438t-dir ../analysis/s438t --strict

Requires: numpy.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from gate_utils import GateLogger, write_gate_status, enforce_strict


def inputs_compatible(wt_resids: list[int], mut_resids: list[int]) -> tuple[bool, str]:
    """Pure gate: the RMSF comparison only means something if both systems
    were analyzed over the same pocket-residue set. Unit-tested."""
    if not wt_resids or not mut_resids:
        return False, "one or both systems have no RMSF data"
    if set(wt_resids) != set(mut_resids):
        return False, (f"pocket-residue sets differ: WT={sorted(wt_resids)} "
                        f"vs S438T={sorted(mut_resids)} -- re-run "
                        f"05_analysis.py with matching --pocket-resids")
    return True, f"{len(wt_resids)} matched residues"


def compare_rmsf(wt_resids: list[int], wt_rmsf: np.ndarray,
                  mut_resids: list[int], mut_rmsf: np.ndarray) -> list[dict]:
    """Per-residue RMSF delta (S438T - WT), matched by resid rather than
    by array position -- the two systems' 05_analysis.py runs aren't
    guaranteed to enumerate residues in the same order. Pure function,
    unit-tested."""
    wt_by_resid = dict(zip(wt_resids, wt_rmsf))
    mut_by_resid = dict(zip(mut_resids, mut_rmsf))
    rows = []
    for resid in sorted(set(wt_by_resid) & set(mut_by_resid)):
        wt_v, mut_v = wt_by_resid[resid], mut_by_resid[resid]
        rows.append({
            "resid": resid,
            "rmsf_wt_angstrom": float(wt_v),
            "rmsf_s438t_angstrom": float(mut_v),
            "delta_angstrom": float(mut_v - wt_v),
        })
    return rows


def compare_na_coordination(wt_series: np.ndarray, mut_series: np.ndarray) -> dict:
    """Summary stats for each system's Na+-coordination-distance series
    plus the delta in means -- not a time-aligned per-frame comparison,
    since the two systems' trajectories aren't synchronized frame-for-
    frame (independent runs). Pure function, unit-tested."""
    if len(wt_series) == 0 or len(mut_series) == 0:
        return {"wt_mean_angstrom": None, "wt_std_angstrom": None,
                "s438t_mean_angstrom": None, "s438t_std_angstrom": None,
                "delta_mean_angstrom": None,
                "wt_n": len(wt_series), "s438t_n": len(mut_series)}
    wt_mean, mut_mean = float(np.mean(wt_series)), float(np.mean(mut_series))
    return {
        "wt_mean_angstrom": wt_mean, "wt_std_angstrom": float(np.std(wt_series)),
        "s438t_mean_angstrom": mut_mean, "s438t_std_angstrom": float(np.std(mut_series)),
        "delta_mean_angstrom": mut_mean - wt_mean,
        "wt_n": len(wt_series), "s438t_n": len(mut_series),
    }


def _read_rmsf_csv(path: Path) -> tuple[list[int], np.ndarray]:
    if not path.exists():
        return [], np.array([])
    resids, values = [], []
    with open(path) as f:
        reader = csv.reader(f)
        next(reader, None)  # header
        for row in reader:
            resids.append(int(float(row[0])))
            values.append(float(row[1]))
    return resids, np.array(values)


def _read_na_csv(path: Path) -> np.ndarray:
    if not path.exists():
        return np.array([])
    values = []
    with open(path) as f:
        reader = csv.reader(f)
        next(reader, None)  # header
        for row in reader:
            values.append(float(row[0]))
    return np.array(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wt-dir", default="../analysis/wt")
    parser.add_argument("--s438t-dir", default="../analysis/s438t")
    parser.add_argument("--out-dir", default="../analysis/comparison")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    wt_dir, mut_dir = Path(args.wt_dir), Path(args.s438t_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gl = GateLogger()
    gl.log("=" * 60)
    gl.log("WT vs S438T COMPARISON")
    gl.log("=" * 60)

    gl.log(f"\n[1] Reading {wt_dir} and {mut_dir} ...")
    wt_resids, wt_rmsf = _read_rmsf_csv(wt_dir / "rmsf_per_resid.csv")
    mut_resids, mut_rmsf = _read_rmsf_csv(mut_dir / "rmsf_per_resid.csv")
    wt_na = _read_na_csv(wt_dir / "na_coordination_per_frame.csv")
    mut_na = _read_na_csv(mut_dir / "na_coordination_per_frame.csv")

    ok, reason = inputs_compatible(wt_resids, mut_resids)
    gl.log(f"  RMSF inputs: {reason} -> {'OK' if ok else 'NOT COMPARABLE'}")

    rmsf_rows = compare_rmsf(wt_resids, wt_rmsf, mut_resids, mut_rmsf) if ok else []
    gl.log("\n[2] Per-residue RMSF (S438T - WT)...")
    for row in rmsf_rows:
        gl.log(f"  resid {row['resid']}: WT={row['rmsf_wt_angstrom']:.3f} A, "
               f"S438T={row['rmsf_s438t_angstrom']:.3f} A, "
               f"delta={row['delta_angstrom']:+.3f} A")

    na_summary = compare_na_coordination(wt_na, mut_na)
    gl.log("\n[3] Na+ coordination distance (S438/T438) ...")
    if na_summary["wt_mean_angstrom"] is not None:
        gl.log(f"  WT:    {na_summary['wt_mean_angstrom']:.2f} +/- "
               f"{na_summary['wt_std_angstrom']:.2f} A (n={na_summary['wt_n']})")
        gl.log(f"  S438T: {na_summary['s438t_mean_angstrom']:.2f} +/- "
               f"{na_summary['s438t_std_angstrom']:.2f} A (n={na_summary['s438t_n']})")
        gl.log(f"  Delta (S438T - WT): {na_summary['delta_mean_angstrom']:+.2f} A")
    else:
        gl.log(f"  SKIPPED -- missing data (WT n={na_summary['wt_n']}, "
               f"S438T n={na_summary['s438t_n']})")

    verdict = "PASS" if ok else "FAIL"
    gl.log(f"\n  GATE (comparison validity, not a physical result): {verdict}")

    gl.write_report(out_dir / "comparison_report.txt")
    write_gate_status(out_dir / "gate_status.json", verdict, {
        "rmsf_compatible": ok,
        "rmsf_compatibility_reason": reason,
        "n_rmsf_residues_compared": len(rmsf_rows),
        "na_coordination": na_summary,
    })
    if rmsf_rows:
        (out_dir / "rmsf_comparison.json").write_text(json.dumps(rmsf_rows, indent=2))

    enforce_strict(args.strict, verdict,
        f"GATE FAILED: {reason}. Fix the mismatch (or re-run 05_analysis.py "
        f"for both systems with the same --pocket-resids) before treating "
        f"any WT-vs-S438T numbers here as meaningful.")


if __name__ == "__main__":
    main()
