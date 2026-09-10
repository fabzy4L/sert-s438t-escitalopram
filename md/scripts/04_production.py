"""
Production MD run driver: long GPU-offloaded gmx mdrun for one system
(--system wt|s438t). Not gated pass/fail like the earlier stages -- there's
no single "did production work" numeric criterion at this point, that's
what 05_analysis.py's trajectory-completeness checks are for. This script
just runs it (or a short --smoke-test slice) and reminds you to check GPU
utilization before committing to the full 100-200 ns
(docs/research/FUTURE_WORK.md Phase 4), since this repo needs to run BOTH
systems to completion for the WT-vs-S438T comparison to mean anything.

Usage:
    # short GPU-utilization sanity check first, per FUTURE_WORK.md's own
    # "primary bottleneck: compute" warning:
    python 04_production.py --system wt --gro ../system_prep/wt/npt.gro \
        --top ../system_prep/wt/topol.top --mdp ../system_prep/wt/production.mdp \
        --smoke-test

    # then the real run (repeat for --system s438t):
    python 04_production.py --system wt --gro ../system_prep/wt/npt.gro \
        --top ../system_prep/wt/topol.top --mdp ../system_prep/wt/production.mdp

Requires: gmx on PATH, built with CUDA support.
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

SYSTEMS = ["wt", "s438t"]
SMOKE_TEST_STEPS = 5000  # a few ps at a typical 2 fs timestep


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", choices=SYSTEMS, required=True)
    parser.add_argument("--gro", required=True)
    parser.add_argument("--top", required=True)
    parser.add_argument("--mdp", required=True)
    parser.add_argument("--out-dir", default=None,
                         help="defaults to ../system_prep/<system>")
    parser.add_argument("--deffnm", default="md")
    parser.add_argument(
        "--smoke-test", action="store_true",
        help=f"override the .mdp's nsteps down to {SMOKE_TEST_STEPS} (a "
             f"few ps) and run that instead of full production -- use "
             f"this first with nvidia-smi open in another terminal to "
             f"confirm GPU utilization before committing to a multi-day "
             f"run per-system (both WT and S438T need a full run)")
    parser.add_argument("--ntomp", type=int, default=8)
    args = parser.parse_args()

    if not shutil.which("gmx"):
        sys.exit("ERROR: gmx not on PATH")

    out_dir = Path(args.out_dir) if args.out_dir else Path(f"../system_prep/{args.system}")
    out_dir.mkdir(parents=True, exist_ok=True)
    tpr = out_dir / f"{args.deffnm}.tpr"

    grompp_cmd = ["gmx", "grompp", "-f", str(args.mdp), "-c", str(args.gro),
                  "-p", str(args.top), "-o", str(tpr),
                  "-po", str(out_dir / f"{args.deffnm}_mdout.mdp")]
    if args.smoke_test:
        print(f"[SMOKE TEST -- {args.system}] overriding nsteps to "
              f"{SMOKE_TEST_STEPS} -- watch `nvidia-smi` in another "
              f"terminal now")
    subprocess.run(grompp_cmd, cwd=out_dir, check=True)

    mdrun_cmd = ["gmx", "mdrun", "-deffnm", args.deffnm,
                 "-nb", "gpu", "-pme", "gpu", "-bonded", "gpu",
                 "-ntmpi", "1", "-ntomp", str(args.ntomp)]
    if args.smoke_test:
        mdrun_cmd += ["-nsteps", str(SMOKE_TEST_STEPS)]
    subprocess.run(mdrun_cmd, cwd=out_dir, check=True)

    print(f"\nDone ({args.system}). "
          f"{'Smoke test' if args.smoke_test else 'Production'} "
          f"run output: {out_dir / (args.deffnm + '.xtc')}")
    if not args.smoke_test:
        print("Reminder: SERT is a membrane protein (~150k atoms per "
              "FUTURE_WORK.md) -- expect ~10-20 ns/day on a consumer GPU. "
              "If this run's rate looks far off that, check nvidia-smi "
              "before assuming it's just slow. Both --system wt and "
              "--system s438t need to reach comparable trajectory length "
              "before 06_compare_systems.py's comparison is meaningful.")


if __name__ == "__main__":
    main()
