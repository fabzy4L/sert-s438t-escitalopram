# Project Workflow: SERT Receptor Design

## Purpose
This project aims to study the structural biology of the Serotonin Transporter (SERT) and its interactions with SSRIs like Escitalopram. The goal is to understand how specific mutations (e.g., S348T) affect ligand binding and transporter function.

## Workflow Steps

### Phase 1: Data Acquisition
- Download wild-type SERT structures from the PDB (e.g., 5I6Z, 5I71, 2A65).
- Obtain the primary amino acid sequence in FASTA format.

### Phase 2: Ligand Preparation
- Use `notebooks/Escitalopram_PythonPDB.ipynb` to generate a 3D conformation of Escitalopram.
- Input: SMILES string (`Fc1ccc(cc1)[C@@]3(OCc2cc(C#N)ccc23)CCCN(C)C`).
- Output: `output/escitalopram.pdb`.

### Phase 3: Sequence Engineering
- Modify the FASTA sequence to introduce mutations.
- Use `notebooks/REVERSE TRANSLATION.ipynb` to determine the corresponding DNA sequence for synthetic biology applications.

### Phase 4: Structural Modeling (In-Silico Mutagenesis)
- Use tools like Expasy ProtMod to apply mutations to the sequence.
- Model the 3D structure of the mutant using the wild-type as a template (Homology Modeling or Refinement).
- *Relevant Files:* `5i6Z_S348T_3.pdb`, `5i6Z_S348_A.pdb`.

### Phase 5: Visualization and Analysis
- Load wild-type and mutant structures into ChimeraX.
- Overlay structures to identify conformational shifts in the binding pocket.
- Analyze ligand-residue interactions.
- *Saved Sessions:* `Receptor_Analysis.cxs`, `2a65_v_516z.cxs`.

## Future Goals
- Perform Molecular Dynamics (MD) simulations to study the stability of the mutant-ligand complexes.
- Quantitative docking studies to estimate binding affinities.
