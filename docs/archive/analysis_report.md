# ChimeraX Structural Analysis Report: SERT S348T Mutation

## 1. Executive Summary
This report summarizes the structural comparison between the Wild-Type (WT) Human Serotonin Transporter (SERT) and the S348T mutant models. The analysis focuses on global structural stability (RMSD) and local interactions at the mutation site (Residue 348).

## 2. Methodology
The session was initialized by closing previous structures and loading the following models:
- **#1:** `5i6z_A.pdb` (WT Reference)
- **#2:** `5i6Z_S348T_3.pdb` (Mutant Model A)
- **#3:** `5i6Z_S348_A.pdb` (Mutant Model B)
- **#4:** `escitalopram.pdb` (Ligand)

## 3. Structural Alignment Results
Structural alignment was performed using the Matchmaker tool to calculate the Root Mean Square Deviation (RMSD).

| Comparison | RMSD (Å) | Pruned Atom Pairs | Alignment Score |
|:---|:---|:---|:---|
| **WT (#1) vs. Mutant (#2)** | 0.030 | 544 | 2890.1 |
| **WT (#1) vs. Mutant (#3)** | 0.005 | 542 | 2901.7 |

**Observation:** Both mutant models show extremely high structural similarity to the WT, with RMSD values < 0.05 Å, indicating that the S348T mutation does not significantly disrupt the protein's fold.

## 4. Local Interaction Analysis (Residue 348)
Detailed analysis was performed on the environment surrounding the mutation site.

### 4.1 Hydrogen Bonding
- **Scope:** Residues 340-360
- **Result:** **102 Hydrogen Bonds** identified.
- **Significance:** High density of H-bonds suggests a stable local helix/loop configuration in the region surrounding the mutation.

### 4.2 Atomic Contacts
- **Scope:** Residue 348 (4.0 Å radius)
- **Result:** **154 Contacts** detected.
- **Significance:** The mutation site is highly packed, with significant interaction potential with neighboring residues and potentially the ligand.

## 5. Visualization Configuration
The following visual styles were applied for final inspection:
- **Backbone:** Cartoon representation.
- **Region 340-360:** Stick representation.
- **Residue 348:** Sphere representation (Color-coded: White for #1, Yellow for #2, Magenta for #3).
- **Ligand (#4):** Stick representation.

---
*Report generated from ChimeraX Log Analysis (output/1.html)*
