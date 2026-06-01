"""
Generate AutoDock Vina-compatible PDBQT for escitalopram.
Uses the RDKit-generated PDB (output/escitalopram.pdb) and applies:
  - AutoDock4 atom type assignment based on element and connectivity
  - Simplified Gasteiger-like partial charges
  - Torsion tree (ROOT / BRANCH / ENDBRANCH / TORSDOF) from CONECT records

No external dependencies beyond biopython and numpy.
"""

import sys
import numpy as np
from pathlib import Path
from collections import defaultdict

try:
    from Bio.PDB import PDBParser
except ImportError:
    sys.exit("ERROR: pip install biopython")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT       = PROJECT_ROOT / "output"

PDB_IN   = OUTPUT / "escitalopram.pdb"
PDBQT_OUT = OUTPUT / "escitalopram.pdbqt"

# ── Connectivity from CONECT records ──────────────────────────────────────────
def parse_conect(pdb_path):
    """Return adjacency dict {serial: set of bonded serials} (1-indexed)."""
    adj = defaultdict(set)
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("CONECT"):
                parts = line.split()
                src = int(parts[1])
                for dest in parts[2:]:
                    d = int(dest)
                    adj[src].add(d)
                    adj[d].add(src)
    return adj

def bond_count(adj, serial):
    return len(adj[serial])

# ── AutoDock4 atom type assignment ────────────────────────────────────────────
def get_ad4_type(serial, atom_name, element, adj):
    """
    Assign AutoDock4 atom type based on element, atom name, and connectivity.
    """
    el = element.strip().upper() if element else atom_name.strip()[0].upper()
    neighbors = adj[serial]

    if el == "F":
        return "F"
    if el == "H":
        # HD if bonded to O or N (H-bond donor H)
        for nb in neighbors:
            pass  # We don't have element of neighbor easily here
        return "H"
    if el == "C":
        # Aromatic carbons: those in 6-membered rings (heuristic: exactly 2 C neighbors in ring)
        # Simple heuristic: atom name starts with 'C' and has exactly 3 heavy-atom neighbors
        heavy_nb = [n for n in neighbors if n <= 24]  # atoms 1-24 are heavy atoms
        if len(heavy_nb) == 3:
            return "A"  # likely aromatic
        return "C"
    if el == "N":
        # NA if can be H-bond acceptor (e.g., amine, cyano)
        return "NA"
    if el == "O":
        return "OA"  # all oxygens are H-bond acceptors

    return el

# ── Partial charges (simplified Gasteiger-like) ───────────────────────────────
# These are approximations sufficient for comparative docking (WT vs mutant).
# For publication-quality docking, use openbabel or meeko.
CHARGE_BY_TYPE = {
    "F":  -0.218,
    "A":  -0.073,   # aromatic C
    "C":   0.000,   # aliphatic C
    "OA": -0.251,   # ether O in isobenzofuran
    "NA": -0.312,   # amine N; cyano N gets adjusted below
    "H":   0.057,   # non-polar H
}

# Per-atom charge overrides (serial: charge)
CHARGE_OVERRIDES = {
    # N1 = cyano N (atom 15): more negative
    15: -0.413,
    # N2 = tertiary amine N (atom 22): protonated at pH 7.4
    # Using neutral form for simplicity; protonated would add +1 total charge
    22: -0.312,
    # O1 = ether O (atom 9)
    9:  -0.270,
    # F1 (atom 1)
    1:  -0.218,
    # C7 = quaternary C (atom 8)
    8:   0.180,
    # C bonded to CN (nitrile, atom 13=C11): slight positive
    13:  0.100,
}

# ── Torsion tree for escitalopram ─────────────────────────────────────────────
# Rigid core: isobenzofuran bicycle (atoms 2-18, the two fused rings + O-CH2)
# Flexible parts:
#   1. Fluorophenyl rotation: C4(5)-C7(8)  branch atoms: 1-7
#   2. Propyl chain: C7(8)-C16(19)  chain: 19-20-21-22-23-24
#   3. N-methyl 1: N2(22)-C19(23)  branch: 23
#   4. N-methyl 2: N2(22)-C20(24)  branch: 24
#
# ROOT = the rigid bicyclic core (atoms 9-18 = O1, C8-C17)
# TORSDOF = 4 (4 active rotatable bonds: C4-C7, C7-C16, C16-C17, C17-C18)
# Terminal methyls on N are counted but contribute little to binding pose.

RIGID_CORE = set(range(9, 19))          # atoms 9-18 (O1 + isobenzofuran ring)
FLUOROPHENYL = set(range(1, 8))         # atoms 1-7 (F + benzene ring)
PROPYL_CHAIN = {19, 20, 21, 22, 23, 24} # C16, C17, C18, N2, C19, C20
N_METHYL_1   = {23}
N_METHYL_2   = {24}

# ── Build PDBQT lines ─────────────────────────────────────────────────────────
def read_hetatm(pdb_path):
    """Parse HETATM records: returns list of dicts."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                serial   = int(line[6:11])
                name     = line[12:16].strip()
                resname  = line[17:20].strip()
                chain    = line[21].strip() or " "
                resseq   = int(line[22:26])
                x        = float(line[30:38])
                y        = float(line[38:46])
                z        = float(line[46:54])
                occ      = float(line[54:60]) if len(line) > 54 else 1.0
                bfac     = float(line[60:66]) if len(line) > 60 else 0.0
                element  = line[76:78].strip() if len(line) > 76 else name[0]
                atoms.append({
                    "serial": serial, "name": name, "resname": resname,
                    "chain": chain, "resseq": resseq,
                    "x": x, "y": y, "z": z,
                    "occ": occ, "bfac": bfac, "element": element,
                })
    return atoms

def atom_line(serial, name, resname, chain, resseq, x, y, z, charge, ad4type, record="ATOM"):
    padded_name = f" {name:<3s}" if len(name) < 4 else name
    return (
        f"{record:<6s}{serial:5d} {padded_name:<4s} {resname:3s} {chain}"
        f"{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.0:6.2f}{0.00:6.2f}    "
        f"{charge:6.3f} {ad4type:<2s}\n"
    )

def h_on(serial, atoms, adj):
    """Return list of H atom serials bonded to this heavy atom."""
    return [a["serial"] for a in atoms
            if a["serial"] > 24 and serial in adj.get(a["serial"], set())]

def write_pdbqt(atoms, adj, out_path):
    """
    Write PDBQT with torsion tree rooted at C7 (atom 8, quaternary carbon).

    Tree layout:
      ROOT: isobenzofuran bicyclic (atoms 8-18) + their H atoms
        BRANCH 8  5  : fluorophenyl rotation (atoms 1-7 + H)
        BRANCH 8 19  : propyl chain
          BRANCH 19 20
            BRANCH 20 21
              BRANCH 21 22  : N2 amine
                BRANCH 22 23: N-methyl 1
                BRANCH 22 24: N-methyl 2
    TORSDOF 5
    """
    by_serial = {a["serial"]: a for a in atoms}

    def atom_ln(serial):
        a = by_serial[serial]
        ad4type = get_ad4_type(serial, a["name"], a["element"], adj)
        charge  = CHARGE_OVERRIDES.get(serial, CHARGE_BY_TYPE.get(ad4type, 0.0))
        return atom_line(serial, a["name"], a["resname"], a["chain"] or " ",
                         a["resseq"], a["x"], a["y"], a["z"],
                         charge, ad4type)

    # Rigid core: C7(8) + isobenzofuran ring (9-18) + their H atoms (29-33)
    ROOT_HEAVY = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
    ROOT_H     = [h for s in ROOT_HEAVY for h in h_on(s, atoms, adj)]

    # Fluorophenyl heavy atoms in ring order
    FPHENYL_HEAVY = [5, 4, 3, 2, 1, 6, 7]
    FPHENYL_H     = [h for s in FPHENYL_HEAVY for h in h_on(s, atoms, adj)]

    lines = [
        "REMARK Escitalopram PDBQT\n",
        "REMARK Charges: simplified Gasteiger approximation\n",
        "REMARK For accurate results, regenerate with openbabel --partialcharge gasteiger\n",
    ]

    # ROOT
    lines.append("ROOT\n")
    for s in ROOT_HEAVY:
        lines.append(atom_ln(s))
    for s in ROOT_H:
        lines.append(atom_ln(s))
    lines.append("ENDROOT\n")

    # BRANCH 1: fluorophenyl rotation around C7(8)-C4(5) bond
    lines.append("BRANCH  8  5\n")
    for s in FPHENYL_HEAVY:
        lines.append(atom_ln(s))
    for s in FPHENYL_H:
        lines.append(atom_ln(s))
    lines.append("ENDBRANCH  8  5\n")

    # BRANCH 2: propyl chain, C7(8)-C16(19)
    lines.append("BRANCH  8 19\n")
    lines.append(atom_ln(19))
    for s in h_on(19, atoms, adj):
        lines.append(atom_ln(s))

    lines.append("BRANCH 19 20\n")
    lines.append(atom_ln(20))
    for s in h_on(20, atoms, adj):
        lines.append(atom_ln(s))

    lines.append("BRANCH 20 21\n")
    lines.append(atom_ln(21))
    for s in h_on(21, atoms, adj):
        lines.append(atom_ln(s))

    lines.append("BRANCH 21 22\n")
    lines.append(atom_ln(22))

    lines.append("BRANCH 22 23\n")
    lines.append(atom_ln(23))
    for s in h_on(23, atoms, adj):
        lines.append(atom_ln(s))
    lines.append("ENDBRANCH 22 23\n")

    lines.append("BRANCH 22 24\n")
    lines.append(atom_ln(24))
    for s in h_on(24, atoms, adj):
        lines.append(atom_ln(s))
    lines.append("ENDBRANCH 22 24\n")

    lines.append("ENDBRANCH 21 22\n")
    lines.append("ENDBRANCH 20 21\n")
    lines.append("ENDBRANCH 19 20\n")
    lines.append("ENDBRANCH  8 19\n")

    # Active torsions: C7-C4, C7-C16, C16-C17, C17-C18, C18-N2 = 5
    lines.append("TORSDOF 5\n")

    with open(out_path, "w") as f:
        f.writelines(lines)

# ── Main ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("ESCITALOPRAM PDBQT GENERATION")
    print("=" * 60)

    if not PDB_IN.exists():
        sys.exit(f"ERROR: {PDB_IN} not found")

    atoms = read_hetatm(PDB_IN)
    adj   = parse_conect(PDB_IN)

    print(f"\nInput:  {PDB_IN.name}")
    print(f"  Heavy atoms: {sum(1 for a in atoms if a['serial'] <= 24)}")
    print(f"  H atoms:     {sum(1 for a in atoms if a['serial'] > 24)}")
    print(f"  Total bonds: {sum(len(v) for v in adj.values()) // 2}")

    print("\nAtom type assignments:")
    for a in atoms:
        if a["serial"] <= 24:
            ad4 = get_ad4_type(a["serial"], a["name"], a["element"], adj)
            chg = CHARGE_OVERRIDES.get(a["serial"], CHARGE_BY_TYPE.get(ad4, 0.0))
            print(f"  {a['serial']:3d}  {a['name']:<4s}  type={ad4:<3s}  q={chg:+.3f}")

    write_pdbqt(atoms, adj, PDBQT_OUT)
    print(f"\nSaved: {PDBQT_OUT}")
    print(f"Lines: {sum(1 for _ in open(PDBQT_OUT))}")
    print("\nNOTE: Charges are simplified approximations.")
    print("For peer-reviewed docking, regenerate with openbabel --partialcharge gasteiger")
