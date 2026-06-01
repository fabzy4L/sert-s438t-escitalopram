"""
In-silico mutagenesis: SER438 -> THR in 5i6z_A_true.pdb
Corrects the previous study of A348T.
Places OG1 and CG2 using NERF geometry with the g+ rotamer (chi1 = +60 deg).
"""

import sys
import numpy as np
from pathlib import Path

try:
    from Bio import PDB
    from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom
    from Bio.PDB.vectors import Vector
except ImportError:
    sys.exit("ERROR: pip install biopython")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
STRUCTURES   = PROJECT_ROOT / "data" / "structures"

# ── Geometry ───────────────────────────────────────────────────────────────────
def place_atom(p1, p2, p3, bond_length, angle_deg, dihedral_deg):
    p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)
    bc = p3 - p2
    bc /= np.linalg.norm(bc)
    ab = p2 - p1
    ab /= np.linalg.norm(ab)
    n = np.cross(ab, bc)
    norm_n = np.linalg.norm(n)
    if norm_n < 1e-8:
        n = np.array([0.0, 0.0, 1.0])
    else:
        n /= norm_n
    m = np.cross(n, bc)
    angle = np.radians(angle_deg)
    dihed = np.radians(dihedral_deg)
    d = bond_length * (
        -np.cos(angle) * bc
        + np.sin(angle) * np.cos(dihed) * m
        + np.sin(angle) * np.sin(dihed) * n
    )
    return p3 + d

def get_coord(residue, atom_name):
    return np.array(residue[atom_name].get_vector().get_array())

# ── Mutagenesis ────────────────────────────────────────────────────────────────
def ser_to_thr(structure, chain_id, resseq, chi1_deg=60.0):
    """
    Replace SER at resseq with THR by removing OG and adding OG1, CG2.
    """
    model = list(structure.get_models())[0]
    chain = model[chain_id]

    target = None
    for res in chain:
        if res.get_id()[1] == resseq and res.get_resname() == "SER":
            target = res
            break
    if target is None:
        print(f"  ERROR: SER at {chain_id}{resseq} not found.")
        return structure

    # Remove OG
    if "OG" in target:
        target.detach_child("OG")

    # Backbone coords
    N  = get_coord(target, "N")
    CA = get_coord(target, "CA")
    CB = get_coord(target, "CB")

    # Place OG1: N-CA-CB-OG1 dihedral = chi1, CB-OG1 = 1.430 A, CA-CB-OG1 = 109.5
    og1_coord = place_atom(N, CA, CB, bond_length=1.430, angle_deg=109.5, dihedral_deg=chi1_deg)

    # Place CG2: N-CA-CB-CG2 dihedral = chi1 + 120, CB-CG2 = 1.521 A, CA-CB-CG2 = 110.5
    cg2_coord = place_atom(N, CA, CB, bond_length=1.521, angle_deg=110.5, dihedral_deg=chi1_deg + 120.0)

    # Mutate residue name
    target.resname = "THR"

    # Add OG1
    og1 = Atom.Atom("OG1", og1_coord, 1.0, 0.0, " ", " OG1", None, "O")
    target.add(og1)

    # Add CG2
    cg2 = Atom.Atom("CG2", cg2_coord, 1.0, 0.0, " ", " CG2", None, "C")
    target.add(cg2)

    print(f"  Mutated SER{resseq} -> THR{resseq}  |  "
          f"OG1 @ ({og1_coord[0]:.2f}, {og1_coord[1]:.2f}, {og1_coord[2]:.2f})  |  "
          f"CG2 @ ({cg2_coord[0]:.2f}, {cg2_coord[1]:.2f}, {cg2_coord[2]:.2f})")
    return structure

# ── Main ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = PDBParser(QUIET=True)

    wt_path = STRUCTURES / "5i6z_A_true.pdb"
    if not wt_path.exists():
        sys.exit("ERROR: Run structural_audit.py first to produce 5i6z_A_true.pdb")

    print("=" * 60)
    print("IN-SILICO MUTAGENESIS: SER438 -> THR (Correct Pivot)")
    print("=" * 60)

    print("\n1. Generating S438T mutant (SER438 -> THR, chi1=+60 deg g+ rotamer)...")
    s_thr = parser.get_structure("wt", str(wt_path))
    s_thr = ser_to_thr(s_thr, "A", 438, chi1_deg=60.0)
    out_thr = STRUCTURES / "5i6z_A_S438T.pdb"
    io = PDBIO()
    io.set_structure(s_thr)
    io.save(str(out_thr))
    print(f"  Saved: {out_thr.name}")

    print("\nDone. Mutant file:")
    print(f"  {out_thr}")
