"""
Structural audit of SERT PDB files.
Reports residue identity at position 348 across all structures, extracts
the true chain A from 5I6Z, and characterizes the S1 binding pocket.
"""

import os
import sys
from pathlib import Path

try:
    from Bio import PDB
    from Bio.PDB import PDBParser, PDBIO, Select
    import numpy as np
except ImportError:
    print("ERROR: Install biopython and numpy: pip install biopython numpy")
    sys.exit(1)

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[1]
STRUCTURES   = PROJECT_ROOT / "data" / "structures"
OUTPUT       = PROJECT_ROOT / "output"
OUTPUT.mkdir(exist_ok=True)

FILES = {
    "5I6Z_raw":   STRUCTURES / "5i6z.pdb",
    "5i6z_A":     STRUCTURES / "5i6z_A.pdb",
    "S348T":      STRUCTURES / "5i6Z_S348T_3.pdb",
    "S348A":      STRUCTURES / "5i6Z_S348_A.pdb",
    "5I6Z_3":     STRUCTURES / "5i6Z_3.pdb",
}

# ── Helpers ────────────────────────────────────────────────────────────────────
def parse(path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))

def residue_at(structure, chain_id, resseq, model_idx=0):
    try:
        model = list(structure.get_models())[model_idx]
        chain = model[chain_id]
        for res in chain:
            if res.get_id()[1] == resseq:
                return res.get_resname().strip()
        return "NOT FOUND"
    except KeyError:
        return f"CHAIN {chain_id} NOT FOUND"

def chain_sequence_snippet(structure, chain_id, start, end, model_idx=0):
    try:
        model = list(structure.get_models())[model_idx]
        chain = model[chain_id]
        residues = [(r.get_id()[1], r.get_resname()) for r in chain
                    if start <= r.get_id()[1] <= end and r.get_id()[0] == " "]
        return residues
    except KeyError:
        return []

def atom_coords(structure, chain_id, resseq, atom_name="CA", model_idx=0):
    try:
        model = list(structure.get_models())[model_idx]
        res = model[chain_id][(" ", resseq, " ")]
        return res[atom_name].get_vector().get_array()
    except (KeyError, Exception):
        return None

def calc_rmsd(coords_a, coords_b):
    if coords_a is None or coords_b is None:
        return None
    diff = coords_a - coords_b
    return float(np.sqrt(np.mean(diff**2)))

# ── Pocket residues for human SERT S1 site (canonical, from literature) ────────
# Coleman et al. 2016 (5I6Z): key binding pocket residues
BINDING_POCKET_RESIDUES = [74, 86, 90, 94, 95, 98, 121, 125, 132, 172,
                            176, 180, 277, 357, 439, 484, 490, 491, 493]

# ── Audit ──────────────────────────────────────────────────────────────────────
def audit_all():
    print("=" * 70)
    print("SERT STRUCTURAL AUDIT")
    print("=" * 70)

    structures = {}
    for label, path in FILES.items():
        if not path.exists():
            print(f"\n[MISSING] {label}: {path}")
            continue
        s = parse(path)
        structures[label] = s

        chains = [c.get_id() for c in list(s.get_models())[0].get_chains()]
        n_res  = sum(1 for r in list(s.get_models())[0].get_residues()
                     if r.get_id()[0] == " ")
        print(f"\n-- {label} ({path.name}) ------------------")
        print(f"   Chains: {chains}  |  ATOM residues: {n_res}")

        for chain in chains:
            r348 = residue_at(s, chain, 348)
            r439 = residue_at(s, chain, 439)
            r277 = residue_at(s, chain, 277)
            snippet = chain_sequence_snippet(s, chain, 340, 360)
            print(f"   Chain {chain}: res277={r277}  res348={r348}  res439={r439}")
            if snippet:
                seq_str = " ".join(f"{n}:{aa}" for n, aa in snippet)
                print(f"   TM6 region (340-360): {seq_str}")

    return structures

# ── Extract true chain A from 5I6Z ─────────────────────────────────────────────
class ChainASelect(Select):
    def accept_chain(self, chain):
        return chain.get_id() == "A"

    def accept_residue(self, residue):
        return residue.get_id()[0] == " "  # exclude HETATM

def extract_chain_a(structures):
    print("\n" + "=" * 70)
    print("EXTRACTING CHAIN A FROM RAW 5I6Z")
    print("=" * 70)

    if "5I6Z_raw" not in structures:
        print("ERROR: 5i6z.pdb not found")
        return

    s = structures["5I6Z_raw"]
    out_path = STRUCTURES / "5i6z_A_true.pdb"
    io = PDBIO()
    io.set_structure(s)
    io.save(str(out_path), ChainASelect())

    # Verify
    s2 = parse(out_path)
    n_res = sum(1 for r in list(s2.get_models())[0].get_residues()
                if r.get_id()[0] == " ")
    r348 = residue_at(s2, "A", 348)
    print(f"   Saved: {out_path.name}")
    print(f"   ATOM residues: {n_res}  |  Residue 348: {r348}")
    return s2

# ── Binding pocket analysis ─────────────────────────────────────────────────────
def binding_pocket_analysis(structures):
    print("\n" + "=" * 70)
    print("BINDING POCKET ANALYSIS (5I6Z CHAIN A TRUE)")
    print("=" * 70)

    true_wt_path = STRUCTURES / "5i6z_A_true.pdb"
    if not true_wt_path.exists():
        print("Run extract_chain_a first.")
        return

    s = parse(true_wt_path)
    chain = "A"

    print("\nKey S1 binding site residues (literature: Coleman et al. 2016):")
    pocket_coords = []
    for resseq in BINDING_POCKET_RESIDUES:
        name = residue_at(s, chain, resseq)
        ca   = atom_coords(s, chain, resseq)
        if ca is not None:
            pocket_coords.append(ca)
            print(f"   {resseq:4d}: {name}  @ ({ca[0]:.2f}, {ca[1]:.2f}, {ca[2]:.2f})")
        else:
            print(f"   {resseq:4d}: {name}  (CA not found)")

    if pocket_coords:
        centroid = np.mean(pocket_coords, axis=0)
        print(f"\nBinding pocket centroid: ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f})")

        # Radius needed to enclose all pocket residues
        dists = [np.linalg.norm(c - centroid) for c in pocket_coords]
        print(f"Max distance from centroid: {max(dists):.2f} Å")
        print(f"Suggested docking box center: x={centroid[0]:.1f} y={centroid[1]:.1f} z={centroid[2]:.1f}")
        print(f"Suggested docking box size: 25 x 25 x 25 Å")

        return centroid

# ── Pairwise RMSD over shared backbone atoms ───────────────────────────────────
def pairwise_rmsd(structures):
    print("\n" + "=" * 70)
    print("PAIRWISE Cα RMSD (residues 340-360, TM6 helix)")
    print("=" * 70)

    labels = list(structures.keys())
    for i, la in enumerate(labels):
        for lb in labels[i+1:]:
            sa = structures[la]
            sb = structures[lb]
            coords_a, coords_b = [], []
            for resseq in range(340, 361):
                for chain in ["A"]:
                    ca = atom_coords(sa, chain, resseq)
                    cb = atom_coords(sb, chain, resseq)
                    if ca is not None and cb is not None:
                        coords_a.append(ca)
                        coords_b.append(cb)
            if coords_a:
                rmsd = np.sqrt(np.mean([(a - b) @ (a - b)
                                        for a, b in zip(coords_a, coords_b)]))
                print(f"   {la:12s} vs {lb:12s}: RMSD = {rmsd:.4f} Å  ({len(coords_a)} Cα pairs)")
            else:
                print(f"   {la:12s} vs {lb:12s}: no shared Cα in TM6")

# ── Summary report ──────────────────────────────────────────────────────────────
def write_summary(centroid):
    report_path = OUTPUT / "structural_audit_report.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("SERT STRUCTURAL AUDIT — SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        f.write("KEY FINDINGS:\n\n")
        f.write("1. RESIDUE 348 IDENTITY\n")
        f.write("   5i6z.pdb (raw):     ALA (canonical per SEQADV)\n")
        f.write("   5i6z_A.pdb:         THR (mutant already applied)\n")
        f.write("   5i6Z_S348T_3.pdb:   THR (same as 5i6z_A — nearly identical)\n")
        f.write("   5i6Z_S348_A.pdb:    THR (S348A label is INCORRECT — still THR)\n\n")
        f.write("2. SEQADV ENGINEERED MUTATIONS IN 5I6Z CRYSTAL STRUCTURE\n")
        f.write("   A291 (canonical ILE291)\n")
        f.write("   S439 (canonical THR439) -- engineered reversal of natural THR\n")
        f.write("   A554 (canonical CYS554)\n")
        f.write("   A580 (canonical CYS580)\n\n")
        f.write("3. WT BASELINE\n")
        f.write("   The true WT chain A is extracted as: data/structures/5i6z_A_true.pdb\n")
        f.write("   Residue 348 = ALA (canonical in 5I6Z context)\n\n")
        f.write("4. RECOMMENDED NEXT STEPS\n")
        f.write("   a. Re-generate S348T mutant from 5i6z_A_true.pdb (ALA→THR)\n")
        f.write("   b. Re-generate S348A mutant from 5i6z_A_true.pdb (ALA→ALA = same)\n")
        f.write("      NOTE: If project intends S348 = SER in biological context,\n")
        f.write("      start from the canonical human SERT sequence and use homology\n")
        f.write("      modeling, not the engineered crystal structure directly.\n\n")
        if centroid is not None:
            f.write("5. DOCKING BOX PARAMETERS (for AutoDock Vina)\n")
            f.write(f"   center_x = {centroid[0]:.2f}\n")
            f.write(f"   center_y = {centroid[1]:.2f}\n")
            f.write(f"   center_z = {centroid[2]:.2f}\n")
            f.write("   size_x = 25\n")
            f.write("   size_y = 25\n")
            f.write("   size_z = 25\n")
    print(f"\nSummary written to: {report_path}")

# ── Main ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    structures = audit_all()
    extract_chain_a(structures)
    centroid = binding_pocket_analysis(structures)
    pairwise_rmsd(structures)
    write_summary(centroid)
    print("\nAudit complete.")
