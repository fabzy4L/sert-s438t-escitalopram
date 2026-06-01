"""
In-silico mutagenesis: ALA348 -> THR in 5i6z_A_true.pdb
Places OG1 and CG2 using NERF geometry with the g+ rotamer (chi1 = +60 deg),
which is the most common THR rotamer in transmembrane helices.
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
    """
    NERF atom placement.
    Places D bonded to p3, given the chain p1-p2-p3.
    angle_deg  = angle at p3 (p2-p3-D)
    dihedral_deg = dihedral p1-p2-p3-D
    """
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
def ala_to_thr(structure, chain_id, resseq, chi1_deg=60.0):
    """
    Replace ALA at resseq with THR by adding OG1 and CG2.
    Uses g+ rotamer (chi1=+60) by default, typical for TM helices.
    Returns the modified structure.
    """
    model = list(structure.get_models())[0]
    chain = model[chain_id]

    # Find the target residue
    target = None
    for res in chain:
        if res.get_id()[1] == resseq and res.get_resname() == "ALA":
            target = res
            break
    if target is None:
        print(f"  WARNING: ALA at {chain_id}{resseq} not found. Skipping.")
        return structure

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

    print(f"  Mutated ALA{resseq} -> THR{resseq}  |  "
          f"OG1 @ ({og1_coord[0]:.2f}, {og1_coord[1]:.2f}, {og1_coord[2]:.2f})  |  "
          f"CG2 @ ({cg2_coord[0]:.2f}, {cg2_coord[1]:.2f}, {cg2_coord[2]:.2f})")
    return structure

def ala_to_ser(structure, chain_id, resseq, chi1_deg=60.0):
    """
    Replace ALA at resseq with SER by adding OG.
    """
    model = list(structure.get_models())[0]
    chain = model[chain_id]

    target = None
    for res in chain:
        if res.get_id()[1] == resseq and res.get_resname() == "ALA":
            target = res
            break
    if target is None:
        print(f"  WARNING: ALA at {chain_id}{resseq} not found. Skipping.")
        return structure

    N  = get_coord(target, "N")
    CA = get_coord(target, "CA")
    CB = get_coord(target, "CB")

    og_coord = place_atom(N, CA, CB, bond_length=1.417, angle_deg=109.5, dihedral_deg=chi1_deg)

    target.resname = "SER"
    og = Atom.Atom("OG", og_coord, 1.0, 0.0, " ", " OG ", None, "O")
    target.add(og)

    print(f"  Mutated ALA{resseq} -> SER{resseq}  |  "
          f"OG @ ({og_coord[0]:.2f}, {og_coord[1]:.2f}, {og_coord[2]:.2f})")
    return structure

# ── RMSD validation ────────────────────────────────────────────────────────────
def backbone_rmsd(s1, s2, chain_id, start, end):
    parser = PDBParser(QUIET=True)
    coords1, coords2 = [], []
    model1 = list(s1.get_models())[0]
    model2 = list(s2.get_models())[0]
    for resseq in range(start, end + 1):
        for atom_name in ["N", "CA", "C", "O"]:
            try:
                c1 = model1[chain_id][(" ", resseq, " ")][atom_name].get_vector().get_array()
                c2 = model2[chain_id][(" ", resseq, " ")][atom_name].get_vector().get_array()
                coords1.append(c1)
                coords2.append(c2)
            except KeyError:
                pass
    if not coords1:
        return None
    diff = np.array(coords1) - np.array(coords2)
    return float(np.sqrt(np.mean(diff ** 2)))

# ── Main ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = PDBParser(QUIET=True)

    wt_path = STRUCTURES / "5i6z_A_true.pdb"
    if not wt_path.exists():
        sys.exit("ERROR: Run structural_audit.py first to produce 5i6z_A_true.pdb")

    print("=" * 60)
    print("IN-SILICO MUTAGENESIS: 5I6Z CHAIN A")
    print("=" * 60)

    # --- S348T ---
    print("\n1. Generating S348T mutant (ALA348 -> THR, chi1=+60 deg g+ rotamer)...")
    s_thr = parser.get_structure("wt", str(wt_path))
    s_thr = ala_to_thr(s_thr, "A", 348, chi1_deg=60.0)
    out_thr = STRUCTURES / "5i6z_A_S348T_correct.pdb"
    io = PDBIO()
    io.set_structure(s_thr)
    io.save(str(out_thr))
    print(f"  Saved: {out_thr.name}")

    # Verify
    s_thr_check = parser.get_structure("thr", str(out_thr))
    model = list(s_thr_check.get_models())[0]
    res348 = None
    for r in model["A"]:
        if r.get_id()[1] == 348:
            res348 = r
            break
    if res348:
        atoms = [a.get_name() for a in res348]
        print(f"  Residue 348 atoms: {atoms}")

    # --- S348A (WT = ALA, so this is identical to WT in 5I6Z context) ---
    print("\n2. Generating S348A file (ALA348 already ALA in 5I6Z WT)...")
    import shutil
    out_ala = STRUCTURES / "5i6z_A_S348A_correct.pdb"
    shutil.copy(wt_path, out_ala)
    print(f"  Copied WT directly -> {out_ala.name}")
    print("  NOTE: In the 5I6Z crystal structure, position 348 is ALA.")
    print("  S348A in 5I6Z context = WT. No structural change needed.")

    # --- RMSD between WT and S348T (backbone 340-360) ---
    print("\n3. Backbone RMSD validation (TM6 residues 340-360)...")
    s_wt  = parser.get_structure("wt",  str(wt_path))
    s_thr2 = parser.get_structure("thr", str(out_thr))
    rmsd = backbone_rmsd(s_wt, s_thr2, "A", 340, 360)
    if rmsd is not None:
        print(f"  WT vs S348T backbone RMSD (340-360): {rmsd:.4f} A")
        if rmsd < 0.5:
            print("  Expected: very low (<0.5 A) because only side-chain atoms differ.")
        else:
            print("  WARNING: unexpectedly high RMSD — check atom placement.")
    else:
        print("  Could not compute RMSD.")

    print("\nDone. Mutant files:")
    print(f"  {out_thr}")
    print(f"  {out_ala}")
