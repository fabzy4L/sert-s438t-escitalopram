"""
Docking box validation: re-dock co-crystallized (R)-citalopram (68P) from 5I71
into the 5I6Z WT receptor.

Method:
  1. Align 5I71 chain A to 5I6Z chain A on shared CA atoms (BioPython Superimposer)
  2. Apply the transformation to 68P atoms -> 5I6Z coordinate frame
  3. Save transformed reference pose: output/68P_crystal_reference.pdb
  4. Check box centre proximity (centroid within grid = quick pass)
  5. If OpenBabel is on PATH: prepare 68P.pdbqt and run Vina re-dock
     Pass criterion: best-pose heavy-atom RMSD < 2.0 A vs crystal pose

Usage:
  python scripts/validation_redock.py

Requirements:
  pip install biopython numpy
  Optional: openbabel (obabel) on PATH for full re-dock
"""

import subprocess
import shutil
import sys
from pathlib import Path

import numpy as np

try:
    from Bio import PDB
    from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
except ImportError:
    sys.exit("ERROR: pip install biopython")

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT       = Path(__file__).resolve().parents[1]
STRUCTURES = ROOT / "data" / "structures"
OUTPUT     = ROOT / "output"
SCRIPTS    = ROOT / "scripts"
OUTPUT.mkdir(exist_ok=True)

REF_PDB    = STRUCTURES / "5i6z_A_true.pdb"   # 5I6Z chain A (receptor)
MOB_PDB    = STRUCTURES / "5i71.pdb"           # 5I71 full structure
WT_PDBQT   = OUTPUT / "receptor_wt.pdbqt"
LIG_PDBQT  = OUTPUT / "escitalopram.pdbqt"

REF_POSE   = OUTPUT / "68P_crystal_reference.pdb"
DOCK_PDBQT = OUTPUT / "68P_redock.pdbqt"
DOCK_CONF  = OUTPUT / "vina_68P_validation.conf"
DOCK_LOG   = OUTPUT / "vina_68P_validation.log"
REPORT     = OUTPUT / "validation_report.txt"

# S1 box from existing docking setup
BOX_CENTER = np.array([33.06, 187.25, 141.04])
BOX_SIZE   = np.array([25.0,  25.0,   25.0])

LIGAND_RESNAME = "68P"
LIGAND_CHAIN   = "A"
RMSD_THRESHOLD = 2.0  # Angstrom


# ── 1. Align 5I71 to 5I6Z on CA atoms ────────────────────────────────────
def align_5i71_to_5i6z():
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", str(REF_PDB))
    mob_struct = parser.get_structure("mob", str(MOB_PDB))

    ref_ca, mob_ca = {}, {}
    for res in ref_struct[0]["A"]:
        if res.id[0] == " " and "CA" in res:
            ref_ca[res.id[1]] = res["CA"]
    for res in mob_struct[0]["A"]:
        if res.id[0] == " " and "CA" in res:
            mob_ca[res.id[1]] = res["CA"]

    common = sorted(set(ref_ca) & set(mob_ca))
    if len(common) < 50:
        sys.exit(f"ERROR: only {len(common)} common CA atoms — check chain IDs.")

    sup = Superimposer()
    sup.set_atoms([ref_ca[r] for r in common], [mob_ca[r] for r in common])
    sup.apply(list(mob_struct.get_atoms()))

    print(f"  Aligned 5I71 -> 5I6Z: {len(common)} CA atoms, RMSD = {sup.rms:.3f} A")
    return mob_struct


# ── 2. Extract transformed 68P and save reference pose ────────────────────
class LigandSelect(Select):
    def accept_residue(self, residue):
        return residue.get_resname().strip() == LIGAND_RESNAME

def extract_reference_pose(aligned_mob):
    io = PDBIO()
    io.set_structure(aligned_mob)
    io.save(str(REF_POSE), LigandSelect())
    print(f"  Reference pose saved: {REF_POSE.name}")

def get_ligand_coords(pdb_path):
    """Return Nx3 array of heavy-atom coords from a PDB."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("lig", str(pdb_path))
    coords = []
    for atom in struct.get_atoms():
        if atom.element != "H":
            coords.append(atom.get_vector().get_array())
    return np.array(coords)

def get_ligand_centroid(pdb_path):
    return get_ligand_coords(pdb_path).mean(axis=0)


# ── 3. Box proximity check ─────────────────────────────────────────────────
def check_box_proximity(centroid):
    half = BOX_SIZE / 2
    inside = np.all(np.abs(centroid - BOX_CENTER) <= half)
    dist   = np.linalg.norm(centroid - BOX_CENTER)
    return inside, dist


# ── 4. Prepare 68P PDBQT via OpenBabel (if available) ────────────────────
def prepare_ligand_pdbqt():
    if not shutil.which("obabel"):
        return None
    cmd = [
        "obabel", str(REF_POSE),
        "-O", str(DOCK_PDBQT),
        "--partialcharge", "gasteiger",
        "-h", "--gen3d"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if DOCK_PDBQT.exists():
        print(f"  68P PDBQT prepared via OpenBabel: {DOCK_PDBQT.name}")
        return DOCK_PDBQT
    print(f"  OpenBabel conversion failed: {result.stderr.strip()}")
    return None


# ── 5. Vina re-dock ────────────────────────────────────────────────────────
def find_vina():
    if shutil.which("vina"):
        return "vina"
    local = SCRIPTS / "vina.exe"
    return str(local) if local.exists() else None

def write_validation_conf():
    with open(DOCK_CONF, "w") as f:
        f.write(f"receptor = {WT_PDBQT}\n")
        f.write(f"ligand   = {DOCK_PDBQT}\n\n")
        f.write(f"center_x = {BOX_CENTER[0]:.3f}\n")
        f.write(f"center_y = {BOX_CENTER[1]:.3f}\n")
        f.write(f"center_z = {BOX_CENTER[2]:.3f}\n\n")
        f.write(f"size_x = {int(BOX_SIZE[0])}\n")
        f.write(f"size_y = {int(BOX_SIZE[1])}\n")
        f.write(f"size_z = {int(BOX_SIZE[2])}\n\n")
        f.write("exhaustiveness = 16\n")
        f.write("num_modes      = 5\n")
        f.write(f"out            = {DOCK_PDBQT}\n")

def run_vina_validation(vina_bin):
    result = subprocess.run(
        [vina_bin, "--config", str(DOCK_CONF)],
        capture_output=True, text=True
    )
    log = result.stdout + result.stderr
    with open(DOCK_LOG, "w") as f:
        f.write(log)
    return log

def parse_best_affinity(log):
    for line in log.splitlines():
        parts = line.split()
        if parts and parts[0] == "1":
            try:
                return float(parts[1])
            except (ValueError, IndexError):
                pass
    return None


# ── 6. RMSD between two PDB files (heavy atoms) ───────────────────────────
def get_ligand_coords_by_name(pdb_path):
    """Return dict of {atom_name: coords} for heavy atoms only."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("lig", str(pdb_path))
    coords = {}
    for atom in struct.get_atoms():
        if atom.element not in ("H", None):
            coords[atom.get_name().strip()] = atom.get_vector().get_array()
    return coords

def compute_rmsd(ref_pdb, mob_pdb):
    """RMSD matched by atom name to avoid ordering artefacts."""
    ref = get_ligand_coords_by_name(ref_pdb)
    mob = get_ligand_coords_by_name(mob_pdb)
    common = sorted(set(ref) & set(mob))
    if not common:
        return None
    ref_arr = np.array([ref[k] for k in common])
    mob_arr = np.array([mob[k] for k in common])
    diff = ref_arr - mob_arr
    return float(np.sqrt((diff ** 2).sum(axis=1).mean()))


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    lines = []
    def log(msg=""):
        print(msg)
        lines.append(msg)

    log("=" * 60)
    log("DOCKING BOX VALIDATION — 68P RE-DOCK")
    log("=" * 60)

    # Step 1 & 2: align and extract
    log("\n[1] Aligning 5I71 to 5I6Z...")
    aligned = align_5i71_to_5i6z()
    log("\n[2] Extracting transformed 68P reference pose...")
    extract_reference_pose(aligned)

    # Step 3: box check
    log("\n[3] Box proximity check...")
    centroid = get_ligand_centroid(REF_POSE)
    inside, dist = check_box_proximity(centroid)
    log(f"  Crystal pose centroid (5I6Z frame): {centroid.round(2)}")
    log(f"  Box centre:                          {BOX_CENTER}")
    log(f"  Distance centroid -> box centre:     {dist:.2f} A")
    log(f"  Centroid inside grid:                {'YES' if inside else 'NO'}")

    proximity_pass = dist <= max(BOX_SIZE) / 2
    log(f"  Proximity result: {'PASS' if proximity_pass else 'FAIL'}")

    # Step 4 & 5: full re-dock if tools available
    full_rmsd = None
    vina_affinity = None

    log("\n[4] Full re-dock (requires OpenBabel + Vina)...")
    lig_pdbqt = prepare_ligand_pdbqt()
    vina_bin  = find_vina()

    if lig_pdbqt and vina_bin and WT_PDBQT.exists():
        write_validation_conf()
        log(f"  Running Vina ({vina_bin})...")
        vina_log = run_vina_validation(vina_bin)
        vina_affinity = parse_best_affinity(vina_log)
        if vina_affinity:
            log(f"  Best re-dock affinity: {vina_affinity:.3f} kcal/mol")
        log("  Computing heavy-atom RMSD vs crystal pose...")
        full_rmsd = compute_rmsd(REF_POSE, DOCK_PDBQT)
        if full_rmsd is not None:
            verdict = "PASS" if full_rmsd < RMSD_THRESHOLD else "FAIL"
            log(f"  RMSD best pose vs crystal: {full_rmsd:.3f} A  [{verdict} (<{RMSD_THRESHOLD} A)]")
        else:
            log("  RMSD could not be computed from docked output.")
    else:
        missing = []
        if not lig_pdbqt:  missing.append("OpenBabel (obabel) not on PATH")
        if not vina_bin:   missing.append("AutoDock Vina binary not found")
        if not WT_PDBQT.exists(): missing.append(f"{WT_PDBQT.name} not found")
        log(f"  Skipped — missing: {', '.join(missing)}")
        log("  Proximity check above is the available validation.")

    # Summary
    log("\n" + "=" * 60)
    log("SUMMARY")
    log("=" * 60)
    log(f"  Box proximity (centroid within grid): {'PASS' if proximity_pass else 'FAIL'}")
    if full_rmsd is not None:
        verdict = "PASS" if full_rmsd < RMSD_THRESHOLD else "FAIL"
        log(f"  Re-dock RMSD < {RMSD_THRESHOLD} A:                 {verdict}  ({full_rmsd:.3f} A)")
    if vina_affinity:
        log(f"  Re-dock best affinity:                {vina_affinity:.3f} kcal/mol")
    log(f"\n  Reference pose: {REF_POSE}")
    log(f"  Full report:    {REPORT}")

    with open(REPORT, "w") as f:
        f.write("\n".join(lines))

if __name__ == "__main__":
    main()
