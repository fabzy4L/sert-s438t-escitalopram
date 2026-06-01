# Computational Analysis of the S438T Mutation in the Human Serotonin Transporter: Effect on Escitalopram Binding Affinity

**Fabian A. Alvarez-Primo, Ph.D.**

*Computational Structural Biology | May 1, 2026*

---

## Abstract

The human Serotonin Transporter (SERT, *SLC6A4*) is the primary molecular target for Selective Serotonin Reuptake Inhibitors (SSRIs). Understanding the structural determinants of ligand affinity is critical for drug design and explaining clinical resistance. This study investigates the **S438T** mutation — a site within the S1 binding pocket identified as a key determinant of SSRI sensitivity (Plenge *et al.*, 2020; Andersen *et al.*, 2009). Utilizing the 5I6Z crystal structure, we performed *in-silico* mutagenesis and molecular docking with AutoDock Vina 1.2.5. Our results indicate that escitalopram binds to the S438T mutant with a best-pose affinity of **−8.3 kcal/mol**, compared to **−8.1 kcal/mol** for the wild-type, a $\Delta\Delta G$ of **−0.2 kcal/mol**. While the Serine-to-Threonine substitution introduces a methyl group into the S1 pocket, escitalopram maintains a high-affinity binding mode in the rigid-receptor model. These results provide a computational baseline for SERT–SSRI interaction studies at the S438 site.

---

## 1. Introduction

Selective serotonin reuptake inhibitors (SSRIs) such as escitalopram are the first-line treatment for major depressive disorder and anxiety. These drugs function by binding to the central orthosteric site (S1) of the Serotonin Transporter (SERT), preventing the reuptake of serotonin into the presynaptic neuron.

The structural biology of SERT was significantly advanced by the resolution of the human SERT construct in complex with antidepressants (PDB: 5I6Z) [1]. Among the residues defining the S1 pocket, **Serine 438** is of particular interest. Previous experimental work established that the S438T mutation — introducing a single γ-methyl group via Serine-to-Threonine substitution — creates a steric clash that dramatically reduces the affinity of citalopram and escitalopram [2, 3]. This study presents a direct computational assessment of S438T using *in-silico* mutagenesis and molecular docking to quantify the effect on escitalopram binding affinity.

---

## 2. Methods

### 2.1 Receptor Preparation
The wild-type receptor was derived from chain A of PDB 5I6Z. The **S438T** mutation was introduced *in-silico* using a Python-based implementation of the Natural Extension Reference Frame (NERF) algorithm. The threonine side chain was placed in the **g+ rotamer** ($\chi1 = +60^\circ$), the most statistically favored conformation for threonine in transmembrane helices.

### 2.2 Ligand Preparation
The 3D structure of **escitalopram** was geometry-optimized and converted to PDBQT format. Torsional degrees of freedom were assigned to the propyl-amine chain and the fluorophenyl ring to allow conformational flexibility during docking.

### 2.3 Molecular Docking
Docking was performed using **AutoDock Vina 1.2.5**. The grid box was centered on the S1 binding pocket (33.06, 187.25, 141.04) with dimensions of 25 × 25 × 25 Å. To obtain statistically robust affinity estimates, docking was run with exhaustiveness=32 across 5 independent random seeds (42, 123, 456, 789, 1001) for both the wild-type (SER438) and the mutant (THR438) receptors. Mean best-pose affinity ± standard deviation was reported across seeds.

---

## 3. Results

### 3.1 Comparative Binding Affinity

The docking results for escitalopram are summarized in the table below.

| Receptor State | Residue 438 | Mean Affinity ± SD (kcal/mol) | n |
| :--- | :---: | :---: | :---: |
| Wild-Type | SER | −8.575 ± 0.013 | 5 |
| **S438T Mutant** | **THR** | **−8.324 ± 0.012** | 5 |

(exhaustiveness=32; seeds=[42, 123, 456, 789, 1001])

The S438T mutation resulted in a marginal decrease in predicted binding affinity ($\Delta\Delta G = +0.251$ kcal/mol). In the 5I6Z crystal conformation, the additional methyl group of Threonine 438 does not sterically ablate escitalopram binding — the ligand accommodates the substitution within the S1 pocket geometry — but the multi-seed analysis reveals a consistent, reproducible weakening of the predicted affinity.

---

## 4. Discussion

The +0.251 kcal/mol shift suggests that **escitalopram** (the (S)-enantiomer) shows marginally weaker predicted binding at the S438T mutant relative to wild-type. In the rigid-receptor model, the pocket accommodates the Threonine side chain without fully disrupting the primary pharmacophoric contacts with Y95 and D98.

However, rigid docking does not capture the dynamic "locking" mechanism or the desolvation penalties associated with the S438T substitution. The +0.251 kcal/mol difference falls within Vina's scoring noise (~0.5 kcal/mol), making the binding essentially equivalent to wild-type in this structural snapshot. This contrasts with the experimentally observed 320-fold increase in $K_i$ for escitalopram reported by Andersen *et al.* (2009), which is attributed to desolvation costs, Na⁺ coordination disruption at S438, and induced-fit conformational changes invisible to a static-receptor model. These effects must be addressed in future flexible-receptor or molecular dynamics studies.

---

## 5. Conclusion

This study provides a computational assessment of the SERT S438T mutation using *in-silico* mutagenesis and multi-seed molecular docking. Escitalopram retains high affinity for the S438T mutant (−8.324 ± 0.012 kcal/mol vs. −8.575 ± 0.013 kcal/mol wild-type; $\Delta\Delta G = +0.251$ kcal/mol; n=5, exhaustiveness=32). These results establish a robust rigid-docking baseline and set the stage for molecular dynamics simulations to resolve the kinetic and thermodynamic consequences of the S438T substitution.

---

## 6. References

1. Coleman, J.A., *et al.* (2016). X-ray structures of the human serotonin transporter in complex with an antidepressant. *Science*, 352(6290), 1478–1480.
2. Andersen, J., *et al.* (2009). Location of the antidepressant binding site in the serotonin transporter. *Journal of Biological Chemistry*, 284(16), 10276–10284.
3. Plenge, P., *et al.* (2020). The mechanism of a high-affinity allosteric inhibitor of the serotonin transporter. *Nature Communications*, 11, 1491.
4. Trott, O., & Olson, A.J. (2010). AutoDock Vina: Improving the speed and accuracy of docking. *Journal of Computational Chemistry*, 31(2), 455–461.
