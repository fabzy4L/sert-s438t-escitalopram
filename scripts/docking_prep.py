"""
Docking preparation pipeline for SERT + Escitalopram (S438T study).

Generates PDBQT files for WT and S438T receptors and writes AutoDock Vina
config files using the S1 binding pocket centroid from structural_audit.py.

Requirements:
  pip install biopython numpy
  AutoDock Vina binary: https://github.com/ccsb-scripps/AutoDock-Vina/releases
  Place vina.exe (Windows) in scripts/ or add it to PATH.

Usage:
  python scripts/docking_prep.py
  vina --config output/vina_wt.conf
  vina --config output/vina_s438t.conf

For multi-seed publication-grade runs (exhaustiveness=32):
  python scripts/docking_multiseed.py
"""

import sys
import os
from pathlib import Path

try:
    from Bio import PDB
    from Bio.PDB import PDBParser, PDBIO, Select
    import numpy as np
except ImportError:
    sys.exit("ERROR: pip install biopython numpy")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
STRUCTURES   = PROJECT_ROOT / "data" / "structures"
OUTPUT       = PROJECT_ROOT / "output"
OUTPUT.mkdir(exist_ok=True)

# Binding pocket centroid from structural_audit.py (5I6Z S1 site)
BOX_CENTER = (33.06, 187.25, 141.04)
BOX_SIZE   = (25, 25, 25)

# AutoDock4 atom type map: (element, is_aromatic, is_hydrogen_acceptor, is_donor)
# Simplified mapping sufficient for protein receptor preparation
AA_AUTODOCK_TYPES = {
    "C":  "C",  "CA": "C",  "CB": "C",  "CD": "C",  "CD1": "C",  "CD2": "C",
    "CE": "C",  "CE1": "C", "CE2": "C", "CE3": "C", "CG": "C",   "CG1": "C",
    "CG2": "C", "CH2": "C", "CZ": "C",  "CZ2": "C", "CZ3": "C",
    "N":  "NA", "ND1": "NA","ND2": "N", "NE": "N",  "NE1": "N",  "NE2": "N",
    "NH1": "N", "NH2": "N", "NZ": "N",
    "O":  "OA", "OD1": "OA","OD2": "OA","OE1": "OA","OE2": "OA","OG": "OA",
    "OG1": "OA","OH": "OA",
    "S":  "SA", "SD": "SA", "SG": "SA",
}

def get_autodock_type(atom_name, element):
    """Return AutoDock atom type for a protein atom."""
    upper = atom_name.strip().upper()
    if upper in AA_AUTODOCK_TYPES:
        return AA_AUTODOCK_TYPES[upper]
    el = element.strip().upper() if element else upper[0]
    return {"C": "C", "N": "N", "O": "OA", "S": "SA", "H": "H"}.get(el, "C")

def pdb_to_pdbqt_receptor(pdb_path: Path, pdbqt_path: Path):
    """Convert a receptor PDB to PDBQT format (rigid, no torsions, zero charges)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", str(pdb_path))

    lines = []
    serial = 1
    for model in structure.get_models():
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if residue.get_id()[0] != " ":
                    continue  # skip HETATM
                for atom in residue.get_atoms():
                    coord = atom.get_vector().get_array()
                    occ   = atom.get_occupancy() or 1.0
                    bfac  = atom.get_bfactor()   or 0.0
                    el    = atom.element or ""
                    atype = get_autodock_type(atom.get_name(), el)
                    name  = atom.get_name().strip()
                    resn  = residue.get_resname().strip()
                    chain_id = chain.get_id()
                    resseq   = residue.get_id()[1]

                    line = (
                        f"ATOM  {serial:5d} {name:<4s} {resn:3s} {chain_id}"
                        f"{resseq:4d}    "
                        f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                        f"{occ:6.2f}{bfac:6.2f}    "
                        f"{'0.000':>6s} {atype:<2s}\n"
                    )
                    lines.append(line)
                    serial += 1
        break  # first model only

    with open(pdbqt_path, "w") as f:
        f.writelines(lines)
    print(f"  Receptor PDBQT: {pdbqt_path.name}  ({serial-1} atoms)")

def write_vina_config(config_path: Path, receptor_pdbqt: Path,
                      ligand_pdbqt: Path, output_pdbqt: Path,
                      center, size, exhaustiveness=8, n_poses=10):
    with open(config_path, "w") as f:
        f.write(f"receptor = {receptor_pdbqt}\n")
        f.write(f"ligand   = {ligand_pdbqt}\n\n")
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n\n")
        f.write(f"size_x = {size[0]}\n")
        f.write(f"size_y = {size[1]}\n")
        f.write(f"size_z = {size[2]}\n\n")
        f.write(f"exhaustiveness = {exhaustiveness}\n")
        f.write(f"num_modes      = {n_poses}\n")
        f.write(f"out            = {output_pdbqt}\n")
    print(f"  Vina config: {config_path.name}")

def escitalopram_pdbqt_note(ligand_pdbqt_path: Path):
    """
    Write a placeholder PDBQT for escitalopram and instructions for proper prep.
    Full PDBQT (with Gasteiger charges and torsions) requires openbabel or meeko+rdkit.
    """
    note_path = ligand_pdbqt_path.with_name("LIGAND_PREP_INSTRUCTIONS.txt")
    with open(note_path, "w") as f:
        f.write("ESCITALOPRAM LIGAND PDBQT PREPARATION\n")
        f.write("=" * 50 + "\n\n")
        f.write("Option A — OpenBabel (recommended, free):\n")
        f.write("  Install: https://openbabel.org/wiki/Get_Open_Babel\n")
        f.write("  Run:\n")
        f.write("    obabel output/escitalopram.pdb -O output/escitalopram.pdbqt\n")
        f.write("           --partialcharge gasteiger -h\n\n")
        f.write("Option B — conda + meeko (most accurate):\n")
        f.write("  conda install -c conda-forge rdkit\n")
        f.write("  pip install meeko\n")
        f.write("  mk_prepare_ligand.py -i output/escitalopram.pdb\n")
        f.write("                       -o output/escitalopram.pdbqt\n\n")
        f.write("Option C — Use pre-prepared PDBQT from ZINC:\n")
        f.write("  Search 'escitalopram' at https://zinc.docking.org\n")
        f.write("  Download PDBQT directly.\n\n")
        f.write("Once prepared, place the file at:\n")
        f.write(f"  {ligand_pdbqt_path}\n")
    print(f"  Ligand prep instructions: {note_path.name}")
    return note_path

def check_autodock_vina():
    """Check if the vina binary is available."""
    import shutil
    if shutil.which("vina"):
        return "vina"
    local_vina = Path(__file__).parent / "vina.exe"
    if local_vina.exists():
        return str(local_vina)
    return None

def run_docking(config_path: Path, vina_binary: str):
    import subprocess
    log_path = config_path.with_suffix(".log")
    print(f"\nRunning: {vina_binary} --config {config_path}")
    result = subprocess.run(
        [vina_binary, "--config", str(config_path)],
        capture_output=True, text=True
    )
    combined = result.stdout + result.stderr
    print(combined)
    with open(log_path, "w") as f:
        f.write(combined)
    if result.returncode != 0:
        print(f"Vina exit code: {result.returncode}")
    return result.returncode == 0

def parse_vina_results(log_path: Path):
    """Extract binding affinities from vina log."""
    if not log_path.exists():
        return []
    results = []
    with open(log_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit():
                try:
                    results.append({
                        "mode": int(parts[0]),
                        "affinity_kcal_mol": float(parts[1]),
                        "rmsd_lb": float(parts[2]),
                        "rmsd_ub": float(parts[3]),
                    })
                except ValueError:
                    continue
    return results

def print_docking_comparison(wt_results, s438t_results):
    print("\n" + "=" * 60)
    print("DOCKING RESULTS: ESCITALOPRAM BINDING AFFINITY")
    print("=" * 60)
    for label, results in [("WT (SER438)", wt_results), ("S438T (THR438)", s438t_results)]:
        if results:
            best = results[0]["affinity_kcal_mol"]
            print(f"  {label}: best pose = {best:.2f} kcal/mol")
        else:
            print(f"  {label}: no results")
    if wt_results and s438t_results:
        delta = s438t_results[0]["affinity_kcal_mol"] - wt_results[0]["affinity_kcal_mol"]
        print(f"\n  Delta-Delta-G (S438T - WT): {delta:+.3f} kcal/mol")
        print(f"  Note: differences < 0.5 kcal/mol are within Vina scoring noise.")
        print(f"  For publication-grade results run: python scripts/docking_multiseed.py")

if __name__ == "__main__":
    print("SERT + Escitalopram Docking Preparation")
    print("=" * 50)

    # Receptor files — use correctly generated mutants if available,
    # fall back to old files for now
    wt_pdb    = STRUCTURES / "5i6z_A_true.pdb"
    s438t_pdb = STRUCTURES / "5i6z_A_S438T.pdb"

    if not wt_pdb.exists():
        print("ERROR: Run structural_audit.py first to extract 5i6z_A_true.pdb")
        sys.exit(1)

    if not s438t_pdb.exists():
        print("ERROR: 5i6z_A_S438T.pdb not found.")
        print("  Run scripts/generate_s438t.py first to generate the S438T mutant.")
        sys.exit(1)

    # Convert receptors to PDBQT
    print("\n1. Preparing receptor PDBQT files...")
    wt_pdbqt    = OUTPUT / "receptor_wt.pdbqt"
    s438t_pdbqt = OUTPUT / "receptor_s438t.pdbqt"
    pdb_to_pdbqt_receptor(wt_pdb,    wt_pdbqt)
    pdb_to_pdbqt_receptor(s438t_pdb, s438t_pdbqt)

    # Ligand
    print("\n2. Ligand preparation...")
    ligand_pdbqt = OUTPUT / "escitalopram.pdbqt"
    if not ligand_pdbqt.exists():
        print("  escitalopram.pdbqt not found — generating prep instructions.")
        escitalopram_pdbqt_note(ligand_pdbqt)
    else:
        print(f"  Found: {ligand_pdbqt.name}")

    # Vina config files
    print("\n3. Writing AutoDock Vina config files...")
    wt_conf    = OUTPUT / "vina_wt.conf"
    s438t_conf = OUTPUT / "vina_s438t.conf"
    out_wt     = OUTPUT / "docked_wt.pdbqt"
    out_s438t  = OUTPUT / "docked_s438t.pdbqt"

    write_vina_config(wt_conf,    wt_pdbqt,    ligand_pdbqt, out_wt,    BOX_CENTER, BOX_SIZE)
    write_vina_config(s438t_conf, s438t_pdbqt, ligand_pdbqt, out_s438t, BOX_CENTER, BOX_SIZE)

    # Try to run docking if vina is available and ligand is ready
    print("\n4. Checking for AutoDock Vina binary...")
    vina_bin = check_autodock_vina()
    if vina_bin and ligand_pdbqt.exists():
        print(f"  Found: {vina_bin}")
        run_docking(wt_conf,    vina_bin)
        run_docking(s438t_conf, vina_bin)
        wt_results    = parse_vina_results(wt_conf.with_suffix(".log"))
        s438t_results = parse_vina_results(s438t_conf.with_suffix(".log"))
        print_docking_comparison(wt_results, s438t_results)
    else:
        if not vina_bin:
            print("  AutoDock Vina not found.")
        if not ligand_pdbqt.exists():
            print("  Ligand PDBQT not prepared yet.")
        print("\n  TO RUN DOCKING:")
        print("  1. Install AutoDock Vina:")
        print("     https://github.com/ccsb-scripps/AutoDock-Vina/releases")
        print("     Place vina.exe in scripts/ or add to PATH")
        print("  2. Prepare ligand PDBQT (see output/LIGAND_PREP_INSTRUCTIONS.txt)")
        print("  3. Re-run: python scripts/docking_prep.py")
        print("\n  OR run docking manually:")
        print(f"     vina --config output/vina_wt.conf")
        print(f"     vina --config output/vina_s438t.conf")
        print("\n  For publication-grade runs (exhaustiveness=32, 5 seeds):")
        print(f"     python scripts/docking_multiseed.py")

    print("\nDone. Output files in:", OUTPUT)
