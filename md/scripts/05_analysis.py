"""
Per-system post-production analysis (docs/research/FUTURE_WORK.md Phase 4
"Key analyses for paper"): RMSF of the pocket-adjacent residues (Y95, D98,
438, S439 -- Phase 1's flexible-residue set), the Na+ coordination distance
time series at residue 438, and pocket volume (fpocket, best-effort). Run
once per system (--system wt|s438t); 06_compare_systems.py combines both
systems' output into the actual WT-vs-S438T comparison FUTURE_WORK.md
calls for.

Gated on trajectory completeness/sanity (enough frames, no NaNs), NOT on a
physical pass/fail threshold -- unlike psilocin-5ht2a's redock gate, RMSF
and Na+ coordination distance ARE the publication result here, not a
correctness check on the pipeline. A "gate" in this file means "is this
trajectory even usable for the comparison downstream," not "is the
biology good."

Residue-438 atom-name trap: Ser438 (WT) has a hydroxyl oxygen named OG;
Thr438 (S438T) has the analogous oxygen named OG1, not OG -- CHARMM36
residue templates name them differently because Thr has an extra methyl
(CG2). Silently reusing "OG" for both systems would find nothing (and,
worse, would find nothing silently) for S438T. --resid438-atom defaults
per --system rather than sharing one default, specifically to not repeat
the class of bug in docs/REDOCK_BUG_HANDOFF.md (this repo's own
validation_redock.py had it too, fixed in commit 3d18fcd) where a
plausible-looking default silently pointed at the wrong atom/frame/pose.

Uses MDAnalysis (see requirements.txt) for trajectory/RMSF handling rather
than hand-parsing XTC frames, for the same reason psilocin-5ht2a's
05_analysis.py does: multi-frame binary trajectory parsing is exactly the
class of problem a maintained library should own. The MDAnalysis- and
fpocket-dependent code is isolated in run_analysis()/run_fpocket() so the
pure gate/geometry/parsing functions below them stay importable and
testable without either installed (see tests/test_05_analysis.py).

Usage:
    python 05_analysis.py --system wt \
        --topology ../system_prep/wt/md.tpr --trajectory ../system_prep/wt/md.xtc \
        --strict
    python 05_analysis.py --system s438t \
        --topology ../system_prep/s438t/md.tpr --trajectory ../system_prep/s438t/md.xtc \
        --strict

Requires: MDAnalysis, numpy. fpocket on PATH for pocket volume (optional --
SKIPPED if absent, not a hard failure).
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from gate_utils import GateLogger, write_gate_status, enforce_strict

SYSTEMS = ["wt", "s438t"]
RESID_438_ATOM_DEFAULT = {"wt": "OG", "s438t": "OG1"}  # Ser vs Thr, see module docstring
POCKET_RESIDS_DEFAULT = "95 98 438 439"  # Y95, D98, [S/T]438, S439 -- FUTURE_WORK.md Phase 1
# CHARMM36/CGenFF (CHARMM-GUI's force field) names the sodium ion residue
# SOD, not the generic PDB "NA" -- get this wrong and the Na+ selection
# silently matches nothing.
NA_RESNAME_DEFAULT = "SOD"
MIN_FRAMES_DEFAULT = 10


def trajectory_sanity(n_frames: int, min_frames: int = MIN_FRAMES_DEFAULT) -> tuple[bool, str]:
    """Completeness gate, not a physical-result gate: is there enough
    trajectory to say anything at all. Pure function, unit-tested."""
    if n_frames == 0:
        return False, "trajectory has zero frames"
    if n_frames < min_frames:
        return False, f"only {n_frames} frames (< {min_frames} minimum)"
    return True, f"{n_frames} frames"


def na_coordination_distance(na_pos: np.ndarray, target_pos: np.ndarray) -> float:
    """Distance from one Na+ ion position to the residue-438 reference
    atom. Pure function -- the "nearest ion" reduction happens in
    run_analysis() across all Na+ positions in a frame, this just does the
    one-pair geometry so it's testable without a trajectory."""
    return float(np.linalg.norm(na_pos - target_pos))


FPOCKET_POCKET_HEADER = re.compile(r"^Pocket\s+(\d+)\s*:")
FPOCKET_VOLUME_LINE = re.compile(r"Volume\s*:\s*([\d.eE+-]+)")
FPOCKET_SCORE_LINE = re.compile(r"^\s*Score\s*:\s*([\d.eE+-]+)")


def parse_fpocket_info(text: str) -> list[dict]:
    """Parses fpocket's <name>_info.txt output into a list of
    {pocket_id, score, volume_A3} dicts, in file order. Pure function --
    unit-tested against a synthetic fpocket-format fixture in
    tests/test_05_analysis.py, no fpocket binary required."""
    pockets = []
    current = None
    for line in text.splitlines():
        header = FPOCKET_POCKET_HEADER.match(line)
        if header:
            if current is not None:
                pockets.append(current)
            current = {"pocket_id": int(header.group(1)), "score": None, "volume_A3": None}
            continue
        if current is None:
            continue
        vol = FPOCKET_VOLUME_LINE.search(line)
        if vol:
            current["volume_A3"] = float(vol.group(1))
            continue
        score = FPOCKET_SCORE_LINE.match(line)
        if score:
            current["score"] = float(score.group(1))
    if current is not None:
        pockets.append(current)
    return pockets


def select_largest_pocket(pockets: list[dict]) -> dict | None:
    """Picks the largest pocket by volume as a default proxy for the S1
    binding pocket. Pure function -- unit-tested. Caveat (real-run note,
    not enforced here): fpocket ranks pockets across the WHOLE structure,
    so "largest" isn't guaranteed to be the orthosteric S1 site once a
    real trajectory exists -- verify against residue proximity (e.g. does
    the pocket's alpha-sphere centroid sit near Y95/D98/438/S439) before
    trusting this for the paper, don't just take pocket 1 on faith."""
    if not pockets:
        return None
    return max(pockets, key=lambda p: p["volume_A3"] or 0.0)


def run_fpocket(pdb_path: Path) -> list[dict]:
    """Runs fpocket on a single-frame PDB and returns parsed pockets.
    SKIPPED (empty list) if fpocket isn't installed -- this is a
    best-effort analysis, not a hard pipeline requirement."""
    if not shutil.which("fpocket"):
        return []
    subprocess.run(["fpocket", "-f", str(pdb_path)], check=True,
                    cwd=pdb_path.parent)
    info_path = pdb_path.parent / f"{pdb_path.stem}_out" / f"{pdb_path.stem}_info.txt"
    if not info_path.exists():
        return []
    return parse_fpocket_info(info_path.read_text())


def run_analysis(topology: Path, trajectory: Path, system: str,
                  pocket_resids: list[int], na_resname: str, resid438_atom: str):
    import MDAnalysis as mda
    from MDAnalysis.analysis import align, rms

    u = mda.Universe(str(topology), str(trajectory))
    n_frames = len(u.trajectory)

    # Align to the first frame on the pocket-residue backbone before RMSF
    # -- unaligned RMSF conflates whole-protein drift/rotation (real, over
    # a 100-200 ns membrane-protein run) with true local fluctuation.
    pocket_sel = " or ".join(f"resid {r}" for r in pocket_resids)
    ref = u.copy()
    align.AlignTraj(u, ref, select=f"({pocket_sel}) and backbone", in_memory=True).run()

    pocket_atoms = u.select_atoms(f"({pocket_sel}) and name CA")
    rmsf_values = rms.RMSF(pocket_atoms).run().results.rmsf if len(pocket_atoms) else np.array([])

    target = u.select_atoms(f"resid 438 and name {resid438_atom}")
    na_ions = u.select_atoms(f"resname {na_resname}")
    na_series = []
    if len(target) == 1 and len(na_ions) > 0:
        for _ts in u.trajectory:
            dists = [na_coordination_distance(na.position, target.positions[0])
                     for na in na_ions]
            na_series.append(min(dists))

    return n_frames, pocket_atoms.resids.tolist() if len(pocket_atoms) else [], \
        rmsf_values, np.array(na_series)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", choices=SYSTEMS, required=True)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--pocket-resids", default=POCKET_RESIDS_DEFAULT,
                         help=f"space-separated resids for RMSF, default "
                              f"'{POCKET_RESIDS_DEFAULT}' (Y95, D98, 438, "
                              f"S439 per FUTURE_WORK.md Phase 1)")
    parser.add_argument("--na-resname", default=NA_RESNAME_DEFAULT,
                         help=f"sodium ion residue name in this topology, "
                              f"default '{NA_RESNAME_DEFAULT}' (CHARMM36/"
                              f"CGenFF convention -- override if a "
                              f"different force field was used)")
    parser.add_argument("--resid438-atom", default=None,
                         help="reference atom on residue 438 for Na+ "
                              "coordination distance -- defaults to OG "
                              "(Ser, --system wt) or OG1 (Thr, --system "
                              "s438t), see module docstring for why these "
                              "differ")
    parser.add_argument("--min-frames", type=int, default=MIN_FRAMES_DEFAULT)
    parser.add_argument("--pocket-pdb", default=None,
                         help="single-frame PDB for fpocket volume "
                              "analysis (e.g. the first production frame); "
                              "omit to skip pocket-volume analysis")
    parser.add_argument("--out-dir", default=None,
                         help="defaults to ../analysis/<system>")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    resid438_atom = args.resid438_atom or RESID_438_ATOM_DEFAULT[args.system]
    pocket_resids = [int(r) for r in args.pocket_resids.split()]
    out_dir = Path(args.out_dir) if args.out_dir else Path(f"../analysis/{args.system}")
    out_dir.mkdir(parents=True, exist_ok=True)

    gl = GateLogger()
    gl.log("=" * 60)
    gl.log(f"POST-PRODUCTION ANALYSIS -- {args.system.upper()}")
    gl.log("=" * 60)

    gl.log(f"\n[1] Loading {args.trajectory} ...")
    n_frames, rmsf_resids, rmsf_values, na_series = run_analysis(
        Path(args.topology), Path(args.trajectory), args.system,
        pocket_resids, args.na_resname, resid438_atom)

    ok, reason = trajectory_sanity(n_frames, args.min_frames)
    gl.log(f"  {reason} -> {'OK' if ok else 'INSUFFICIENT'}")

    gl.log(f"\n[2] Pocket-residue RMSF (resids {pocket_resids}, "
           f"{len(rmsf_resids)} CA atoms matched)...")
    for resid, val in zip(rmsf_resids, rmsf_values):
        gl.log(f"  resid {resid}: {val:.3f} A")

    gl.log(f"\n[3] Na+ ({args.na_resname}) coordination distance to "
           f"resid 438 {resid438_atom} ...")
    if len(na_series):
        gl.log(f"  Mean: {na_series.mean():.2f} A, "
               f"range: {na_series.min():.2f}-{na_series.max():.2f} A "
               f"(n={len(na_series)} frames)")
    else:
        gl.log(f"  SKIPPED -- couldn't find resid 438 atom '{resid438_atom}' "
               f"and/or any {args.na_resname} ions in this topology")

    fpocket_pockets = []
    largest_pocket = None
    if args.pocket_pdb:
        gl.log(f"\n[4] Pocket volume (fpocket on {args.pocket_pdb})...")
        if not shutil.which("fpocket"):
            gl.log("  SKIPPED -- fpocket not on PATH")
        else:
            fpocket_pockets = run_fpocket(Path(args.pocket_pdb))
            largest_pocket = select_largest_pocket(fpocket_pockets)
            if largest_pocket:
                gl.log(f"  Largest pocket: #{largest_pocket['pocket_id']}, "
                       f"volume {largest_pocket['volume_A3']:.1f} A^3 "
                       f"(verify this is actually the S1 site, see "
                       f"select_largest_pocket()'s docstring)")
            else:
                gl.log("  fpocket ran but found no pockets, or output "
                       "couldn't be parsed")
    else:
        gl.log("\n[4] Pocket volume: skipped (--pocket-pdb not given)")

    verdict = "PASS" if ok else "FAIL"
    gl.log(f"\n  GATE (trajectory usability, not a physical result): {verdict}")

    gl.write_report(out_dir / "analysis_report.txt")
    write_gate_status(out_dir / "gate_status.json", verdict, {
        "system": args.system,
        "n_frames": n_frames,
        "min_frames": args.min_frames,
        "pocket_resids": pocket_resids,
        "rmsf_resids": rmsf_resids,
        "na_resname": args.na_resname,
        "resid438_atom": resid438_atom,
        "na_coordination_mean_angstrom": float(na_series.mean()) if len(na_series) else None,
        "n_fpocket_pockets": len(fpocket_pockets),
        "largest_pocket_volume_A3": largest_pocket["volume_A3"] if largest_pocket else None,
    })
    if len(rmsf_values):
        np.savetxt(out_dir / "rmsf_per_resid.csv",
                   np.column_stack([rmsf_resids, rmsf_values]),
                   delimiter=",", header="resid,rmsf_angstrom", comments="")
    if len(na_series):
        np.savetxt(out_dir / "na_coordination_per_frame.csv", na_series,
                   delimiter=",", header="distance_angstrom", comments="")

    enforce_strict(args.strict, verdict,
        f"GATE FAILED ({args.system}): {reason}. Not enough trajectory to "
        f"feed into 06_compare_systems.py yet.")


if __name__ == "__main__":
    main()
