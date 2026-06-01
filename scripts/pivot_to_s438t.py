"""
Pivot Script: S438T Analysis Pipeline
1. Generates S438T PDBQT.
2. Generates WT PDBQT (from true WT).
3. Writes Vina configs.
4. Runs Vina docking.
5. Prints comparison.
"""

import sys
import os
from pathlib import Path
import subprocess

# Add current dir to path to import from docking_prep
sys.path.append(str(Path(__file__).parent))
try:
    import docking_prep
except ImportError:
    pass # We will define functions locally if needed

PROJECT_ROOT = Path(__file__).resolve().parents[1]
STRUCTURES   = PROJECT_ROOT / "data" / "structures"
OUTPUT       = PROJECT_ROOT / "output"
SCRIPTS      = PROJECT_ROOT / "scripts"

# Redefine necessary functions from docking_prep to be self-contained
def pdb_to_pdbqt_receptor(pdb_path: Path, pdbqt_path: Path):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", str(pdb_path))

    # Simplified mapping
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

    lines = []
    serial = 1
    for model in structure.get_models():
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if residue.get_id()[0] != " ": continue
                for atom in residue.get_atoms():
                    coord = atom.get_vector().get_array()
                    atype = AA_AUTODOCK_TYPES.get(atom.get_name().strip().upper(), "C")
                    line = (
                        f"ATOM  {serial:5d} {atom.get_name().strip():<4s} {residue.get_resname():3s} {chain.get_id()}"
                        f"{residue.get_id()[1]:4d}    "
                        f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                        f"  1.00  0.00    "
                        f" 0.000 {atype:<2s}\n"
                    )
                    lines.append(line)
                    serial += 1
        break
    with open(pdbqt_path, "w") as f:
        f.writelines(lines)

def write_vina_config(config_path, receptor, ligand, out, center, size=(25,25,25)):
    with open(config_path, "w") as f:
        f.write(f"receptor = {receptor}\n")
        f.write(f"ligand   = {ligand}\n\n")
        f.write(f"center_x = {center[0]}\n")
        f.write(f"center_y = {center[1]}\n")
        f.write(f"center_z = {center[2]}\n\n")
        f.write(f"size_x = {size[0]}\n")
        f.write(f"size_y = {size[1]}\n")
        f.write(f"size_z = {size[2]}\n\n")
        f.write(f"exhaustiveness = 8\n")
        f.write(f"num_modes      = 10\n")
        f.write(f"out            = {out}\n")

def run_vina(config, vina_bin):
    log = config.with_suffix(".log")
    cmd = [vina_bin, "--config", str(config)]
    print(f"  Running: {' '.join(map(str, cmd))}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    with open(log, "w") as f:
        f.write(res.stdout + res.stderr)
    return log

def parse_best(log):
    with open(log) as f:
        for line in f:
            if "   1 " in line:
                return float(line.split()[1])
    return None

if __name__ == "__main__":
    wt_pdb = STRUCTURES / "5i6z_A_true.pdb"
    mut_pdb = STRUCTURES / "5i6z_A_S438T.pdb"
    ligand_pdbqt = OUTPUT / "escitalopram.pdbqt"
    vina_bin = SCRIPTS / "vina.exe"

    print("Pivot: WT vs S438T Docking")
    
    # 1. PDBQT
    wt_pdbqt = OUTPUT / "receptor_wt.pdbqt"
    mut_pdbqt = OUTPUT / "receptor_s438t.pdbqt"
    pdb_to_pdbqt_receptor(wt_pdb, wt_pdbqt)
    pdb_to_pdbqt_receptor(mut_pdb, mut_pdbqt)
    
    # 2. Configs
    center = (33.06, 187.25, 141.04) # S1 pocket
    wt_conf = OUTPUT / "vina_wt.conf"
    mut_conf = OUTPUT / "vina_s438t.conf"
    write_vina_config(wt_conf, wt_pdbqt, ligand_pdbqt, OUTPUT / "docked_wt.pdbqt", center)
    write_vina_config(mut_conf, mut_pdbqt, ligand_pdbqt, OUTPUT / "docked_s438t.pdbqt", center)
    
    # 3. Run
    wt_log = run_vina(wt_conf, vina_bin)
    mut_log = run_vina(mut_conf, vina_bin)
    
    # 4. Results
    wt_best = parse_best(wt_log)
    mut_best = parse_best(mut_log)
    
    print("\n" + "="*40)
    print("FINAL RESULTS (Pivot to S438T)")
    print("="*40)
    print(f"WT Affinity:    {wt_best} kcal/mol")
    print(f"S438T Affinity: {mut_best} kcal/mol")
    if wt_best and mut_best:
        print(f"Delta G:        {mut_best - wt_best:+.3f} kcal/mol")
