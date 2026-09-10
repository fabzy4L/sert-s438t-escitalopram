"""
Import and validate a CHARMM-GUI Membrane Builder output archive into
md/system_prep/<system>/, for one of the two parallel systems this repo
tracks (WT or S438T).

CHARMM-GUI's membrane builder download is a directory (or a .tar.gz/.zip
of one) containing GROMACS-ready topology, coordinates, staged .mdp files,
and -- since escitalopram is parameterized via CGenFF, not a standard force
field residue (docs/research/FUTURE_WORK.md, Phase 4 step 1) -- a CGenFF
stream file (.str) and/or ligand-specific .itp that topol.top #includes.
This script does NOT talk to charmm-gui.org or cgenff.umaryland.edu --
those steps are interactive (upload the receptor, place escitalopram, pick
a POPC bilayer; separately submit escitalopram's structure to CGenFF) and
have to happen in a browser first. This script only validates and stages
what comes back from those downloads, so a missing or partial export fails
loudly here instead of silently breaking 02_minimize.py three steps later
-- and specifically catches the common CGenFF mistake of downloading the
.mol2/parameter files but forgetting the .str stream file topol.top
actually needs.

Usage:
    python 01_import_charmm_gui.py --system wt \
        --input path/to/wt-charmm-gui-download --strict
    python 01_import_charmm_gui.py --system s438t \
        --input path/to/s438t-charmm-gui-download --strict

Requires: nothing beyond stdlib.
"""
import argparse
import shutil
import sys
import tarfile
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from gate_utils import GateLogger, write_gate_status, enforce_strict

SYSTEMS = ["wt", "s438t"]

# CHARMM-GUI's gromacs/ output directory uses these canonical names.
REQUIRED_TOPOLOGY = ["topol.top"]
REQUIRED_INDEX = ["index.ndx"]
# The initial coordinate file name varies by CHARMM-GUI version/step
# (step5_input.gro is typical for Membrane Builder); accept common variants
# rather than hardcoding one.
COORDINATE_CANDIDATES = ["step5_input.gro", "step5_input.pdb", "system.gro"]
# CHARMM-GUI bundles a staged sequence of equilibration + production .mdp
# files; require at least one exists rather than hardcoding exact step
# numbers/counts, since that varies with membrane complexity.
MDP_GLOB = "*.mdp"
# CGenFF ligand parameters: the stream file (.str, has both RTF and PRM
# sections) is what topol.top actually #includes for escitalopram. Some
# CHARMM-GUI Ligand Reader flows instead emit split .rtf/.prm -- accept
# either so this isn't over-strict about which CGenFF workflow was used.
CGENFF_STREAM_GLOB = "*.str"
CGENFF_SPLIT_GLOBS = ["*.rtf", "*.prm"]


def find_archive_root(path: Path) -> Path:
    """If path is a CHARMM-GUI tar.gz/zip, extract it and return the
    directory to validate; if it's already a directory, return it as-is."""
    if path.is_dir():
        return path
    if not path.exists():
        sys.exit(f"ERROR: {path} does not exist")
    name = path.name
    if name.endswith((".tar.gz", ".tgz")):
        extract_to = path.parent / name.split(".tar")[0].removesuffix(".tgz")
        with tarfile.open(path) as tf:
            tf.extractall(extract_to)
        return extract_to
    if path.suffix == ".zip":
        extract_to = path.parent / path.stem
        with zipfile.ZipFile(path) as zf:
            zf.extractall(extract_to)
        return extract_to
    sys.exit(f"ERROR: {path} is neither a directory nor a .tar.gz/.zip archive")


def validate_charmm_gui_dir(root: Path) -> dict:
    """Returns a dict describing what was found/missing. Pure function, no
    side effects beyond reading the filesystem -- kept separate from
    main() so it's directly unit-testable against a synthetic fixture
    directory (see tests/test_01_import_charmm_gui.py) without needing a
    real CHARMM-GUI download."""
    topol = locate(root, REQUIRED_TOPOLOGY)
    index = locate(root, REQUIRED_INDEX)
    coords = locate(root, COORDINATE_CANDIDATES)
    mdp_files = sorted(root.rglob(MDP_GLOB))
    cgenff_stream = sorted(root.rglob(CGENFF_STREAM_GLOB))
    cgenff_split = [p for glob in CGENFF_SPLIT_GLOBS for p in sorted(root.rglob(glob))]

    missing = []
    if topol is None:
        missing.append(f"topology ({'/'.join(REQUIRED_TOPOLOGY)})")
    if index is None:
        missing.append(f"index ({'/'.join(REQUIRED_INDEX)})")
    if coords is None:
        missing.append(f"initial coordinates ({'/'.join(COORDINATE_CANDIDATES)})")
    if not mdp_files:
        missing.append("no .mdp files found (need at least equilibration + production)")
    if not cgenff_stream and not cgenff_split:
        missing.append(
            "no CGenFF ligand parameters found (need a .str stream file, "
            "or split .rtf + .prm) -- escitalopram won't resolve in "
            "topol.top without one, see docs/research/FUTURE_WORK.md Phase 4"
        )

    return {
        "topol": str(topol) if topol else None,
        "index": str(index) if index else None,
        "coordinates": str(coords) if coords else None,
        "mdp_files": [str(p) for p in mdp_files],
        "cgenff_files": [str(p) for p in (cgenff_stream + cgenff_split)],
        "missing": missing,
    }


def locate(root: Path, names) -> Path | None:
    for name in names:
        hits = list(root.rglob(name))
        if hits:
            return hits[0]
    return None


def stage_into(out_dir: Path, found: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for key in ("topol", "index", "coordinates"):
        src = found.get(key)
        if src:
            shutil.copy2(src, out_dir / Path(src).name)
    for mdp in found["mdp_files"]:
        shutil.copy2(mdp, out_dir / Path(mdp).name)
    for cgenff in found["cgenff_files"]:
        shutil.copy2(cgenff, out_dir / Path(cgenff).name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", choices=SYSTEMS, required=True,
                         help="which of the two parallel systems this "
                              "CHARMM-GUI export is for -- determines the "
                              "default staging directory")
    parser.add_argument("--input", required=True,
                         help="CHARMM-GUI Membrane Builder download for "
                              "this system: a directory, or a .tar.gz/.zip "
                              "archive of one")
    parser.add_argument("--out-dir", default=None,
                         help="defaults to ../system_prep/<system>")
    parser.add_argument("--strict", action="store_true",
                         help="exit 1 if required files are missing, so "
                              "02_minimize.py can't run against an "
                              "incomplete/miscopied CHARMM-GUI export")
    args = parser.parse_args()
    out_dir = Path(args.out_dir) if args.out_dir else Path(f"../system_prep/{args.system}")

    gl = GateLogger()
    gl.log("=" * 60)
    gl.log(f"CHARMM-GUI IMPORT VALIDATION -- {args.system.upper()}")
    gl.log("=" * 60)

    root = find_archive_root(Path(args.input))
    gl.log(f"\n[1] Validating {root} ...")
    found = validate_charmm_gui_dir(root)

    gl.log(f"  Topology:      {found['topol'] or 'MISSING'}")
    gl.log(f"  Index:         {found['index'] or 'MISSING'}")
    gl.log(f"  Coordinates:   {found['coordinates'] or 'MISSING'}")
    gl.log(f"  .mdp files:    {len(found['mdp_files'])} found")
    for mdp in found["mdp_files"]:
        gl.log(f"    - {mdp}")
    gl.log(f"  CGenFF params: {len(found['cgenff_files'])} found")
    for f in found["cgenff_files"]:
        gl.log(f"    - {f}")

    verdict = "PASS" if not found["missing"] else "FAIL"
    if found["missing"]:
        gl.log("\n  Missing:")
        for m in found["missing"]:
            gl.log(f"    - {m}")

    if verdict == "PASS":
        gl.log(f"\n[2] Staging into {out_dir} ...")
        stage_into(out_dir, found)
        gl.log("  Done.")
    else:
        gl.log("\n[2] Skipped staging -- required files missing.")

    gl.log(f"\n  GATE: {verdict}")

    out_dir.mkdir(parents=True, exist_ok=True)
    gl.write_report(out_dir / "import_report.txt")
    write_gate_status(out_dir / "gate_status.json", verdict, {
        "system": args.system,
        "topol": found["topol"],
        "index": found["index"],
        "coordinates": found["coordinates"],
        "n_mdp_files": len(found["mdp_files"]),
        "n_cgenff_files": len(found["cgenff_files"]),
        "missing": found["missing"],
        "source": str(root),
    })

    enforce_strict(args.strict, verdict,
        f"GATE FAILED ({args.system}): {root} is missing required CHARMM-GUI "
        f"output ({', '.join(found['missing'])}). Re-download from "
        f"charmm-gui.org Membrane Builder / cgenff.umaryland.edu before "
        f"running 02_minimize.py.")


if __name__ == "__main__":
    main()
