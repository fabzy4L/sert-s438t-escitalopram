"""
Split a prepared receptor PDBQT into rigid + flexible PDBQT for the S1
pocket residues Y95, D98, T438/S438, S439 (docs/research/FUTURE_WORK.md
Phase 1: "Allow key S1 pocket residues to flex during docking, capturing
induced-fit at THR438"). Flexible side chains, not backbone.

FUTURE_WORK.md documented this as needing `prepare_flexreceptor4.py`
(MGLTools) -- not installed here. Ported from psilocin-5ht2a's identical
situation (docs/setup.md there, "Flexible-residue docking: tooling gap
closed"): the underlying library class
(AutoDockTools.MoleculePreparation.AD4FlexibleReceptorPreparation) is
bundled inside a local ADFRsuite-1.1dev install and importable directly
under its Python 2 interpreter (pythonsh), bypassing the missing CLI
wrapper entirely. Confirmed working here too, on this repo's hand-rolled
receptor_wt_H.pdbqt / receptor_s438t_minimized_H.pdbqt (built by
docking_prep.py's own biopython-based PDB->PDBQT converter, not
ADFRsuite's prepare_receptor -- AD4FlexibleReceptorPreparation only needs
MolKit to parse the file and residue numbers to match, so the receptor's
origin doesn't matter).

Must run under ADFRsuite's bundled Python 2 interpreter (pythonsh), not
this repo's normal python3 env:

    ~/ADFRsuite-1.1dev/bin/pythonsh scripts/prep_flexreceptor.py \
        --receptor output/receptor_wt_H.pdbqt \
        --residues 95,98,438,439 \
        --out-dir output

Run once per receptor (WT and S438T each need their own split, since
they're different structures with residue 438 as SER vs THR).
"""
import getopt
import os
import sys

from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4FlexibleReceptorPreparation


def main():
    opts, _ = getopt.getopt(sys.argv[1:], "", ["receptor=", "residues=", "out-dir="])
    opt = dict(opts)
    receptor_filename = opt["--receptor"]
    residue_numbers = [s.strip() for s in opt["--residues"].split(",")]
    out_dir = opt.get("--out-dir", "output")

    rec = Read(receptor_filename)[0]
    rec.buildBondsByDistance()

    all_res = []
    for num in residue_numbers:
        matched = [r for r in rec.chains.residues if r.number == num]
        if not matched:
            print "WARNING: no residue numbered", num, "found"
            continue
        all_res.extend(matched)
    print "matched", len(all_res), "of", len(residue_numbers), "requested residues:"
    for r in all_res:
        print "  ", r.full_name()

    base = os.path.splitext(os.path.basename(receptor_filename))[0]
    rigid_out = os.path.join(out_dir, base + "_rigid.pdbqt")
    flex_out = os.path.join(out_dir, base + "_flex.pdbqt")

    AD4FlexibleReceptorPreparation(
        rec,
        mode="automatic",
        rigid_filename=rigid_out,
        flexres_filename=flex_out,
        residues=all_res,
    )
    print "wrote", rigid_out
    print "wrote", flex_out


if __name__ == "__main__":
    main()
