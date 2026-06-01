# Structural and Computational Analysis of the S438T Mutation in the Human Serotonin Transporter and Its Effect on Escitalopram Binding

**Fabian A. Alvarez-Primo, Ph.D.**

*Computational Structural Biology | May 1, 2026*

---

## Abstract

The human Serotonin Transporter (SERT, *SLC6A4*) is the principal molecular target of selective serotonin reuptake inhibitors (SSRIs). Understanding how residue-level mutations within the SERT binding pocket alter drug affinity is fundamental to both mechanistic pharmacology and the rational design of next-generation antidepressants. This study investigates the structural and computational consequences of the **S438T** mutation in SERT, a site identified in the literature (Plenge *et al.*, 2020) as a critical determinant of SSRI sensitivity. A significant portion of this study documents a course-correction from an initial investigation of the **A348T** position — identified here as a digit-swap error in the prior research workflow. Using the 5I6Z crystal structure as a reference, in-silico mutagenesis was performed to replace Serine 438 with Threonine using NERF geometry and the favored g+ rotamer. Molecular docking of escitalopram against both the wild-type (WT) and S438T models using AutoDock Vina 1.2.5 yielded best-pose binding free energies of **−8.1 kcal/mol** and **−8.3 kcal/mol**, respectively, corresponding to a $\Delta\Delta G$ of **−0.2 kcal/mol**. While the Plenge paper suggests a steric clash for certain ligands, our docking indicates that escitalopram maintains a high-affinity binding mode in the S438T variant, with the predicted affinity difference falling within the scoring function's noise floor. These findings underscore the importance of precision in site-specific mutagenesis workflows and provide a corrected baseline for future SERT–SSRI interaction studies.

---

## 1. Introduction

Major depressive disorder (MDD) affects an estimated 280 million people worldwide, representing one of the leading causes of disability globally [1]. Selective serotonin reuptake inhibitors (SSRIs) remain the most widely prescribed antidepressants, with escitalopram — the (S)-enantiomer of citalopram — consistently demonstrating superior tolerability and clinical efficacy within the class [2]. SSRIs act by competitively blocking the serotonin reuptake transporter (SERT, *SLC6A4*), thereby prolonging serotonergic neurotransmission in the synaptic cleft.

The molecular architecture of human SERT was resolved in 2016 by Coleman *et al.* through X-ray crystallography (PDB: 5I6Z) [3]. This structural data revealed the S1 (primary) and S2 (allosteric) binding sites and defined key pharmacophoric interactions. The substitution of serine to threonine at residue **438** — located in the S1 binding pocket — has been shown to induce a "steric clash" that significantly reduces the affinity of citalopram, serving as a hallmark tool for distinguishing binding sites [11].

Crucially, this project initially focused on the **348** position (A348T). A rigorous comparative audit and literature review revealed that this was a **digit-swap error** (348 vs 438). Position 348 is a peripheral residue in TM6 with no direct contact with the S1 pocket, whereas 438 is a direct ligand-contact residue. This study pivot documents the correction of this error and the subsequent analysis of the biologically relevant S438T site.

---

## 2. Methods

### 2.1 Structural Inputs and Data Sources

The human SERT structure was obtained from the Protein Data Bank (PDB entry 5I6Z) [3]. Chain A was extracted to produce the canonical wild-type receptor (`5i6z_A_true.pdb`). 

### 2.2 In-Silico Mutagenesis (S438T)

The **S438T** substitution (Serine to Threonine) was introduced computationally using the NERF (Natural Extension Reference Frame) algorithm [5]. Unlike the previous A348T model, which targeted an Alanine, this mutation replaced the SER438 residue. Threonine side-chain atoms (OG1, CG2) were placed on the original SER CB carbon using the g+ rotamer (χ1 = +60°), which is statistically favored in transmembrane helices.

### 2.3 Molecular Docking

Molecular docking was performed with AutoDock Vina 1.2.5. The box center (33.06, 187.25, 141.04) was centered on the S1 binding pocket. Docking was performed independently against the WT and S438T receptor models using an exhaustiveness of 8.

---

## 3. Results

### 3.1 Correction of the 348 vs 438 Discrepancy

A structural audit revealed that residue 348 (ALA in 5I6Z) is located on the periphery of the transporter, while residue 438 (SER in 5I6Z) is a primary component of the S1 binding pocket. The initial investigation into "S348T" was based on a numerical error. Re-analysis of the source literature (Plenge *et al.*, 2020) confirmed that **S438T** is the mutation of interest for studying SSRI resistance and steric hindrance.

### 3.2 S438T Binding Affinity

AutoDock Vina results for escitalopram binding are summarized in Table 3.

**Table 3. Escitalopram docking results (S438T vs. WT)**

| Receptor Model | Best Affinity (kcal/mol) | $\Delta\Delta G$ (kcal/mol) |
|:---|:---:|:---:|
| Wild-Type (SER438) | −8.1 | — |
| **S438T Mutant (THR438)** | **−8.3** | **−0.2** |

The results show a minor shift in binding affinity (−0.2 kcal/mol) for the S438T mutant. While the Plenge paper reports a significant affinity loss for citalopram at this site, our docking of escitalopram suggests that the (S)-enantiomer may accommodate the additional methyl group of Threonine 438 with minimal thermodynamic penalty under the rigid-receptor conditions of the docking run.

---

## 4. Discussion

The pivot from A348T to S438T transformed the study from an analysis of a peripheral, functionally silent site to a high-impact pocket residue. The computed $\Delta\Delta G$ of −0.2 kcal/mol is below the standard error of the Vina scoring function (~0.5 kcal/mol), suggesting that the S438T mutation does not ablate escitalopram binding in the 5I6Z structural context. 

However, the "steric clash" hypothesis often involves subtle backbone rearrangements or desolvation effects not captured in rigid-receptor docking. Future work utilizing molecular dynamics (MD) will be required to determine if the S438T mutation alters the *residence time* or *binding kinetics* of escitalopram, even if the equilibrium binding free energy remains similar.

---

## 5. Conclusion

This study successfully pivoted to the biologically relevant **S438T** mutation in SERT after identifying a critical digit-swap error in the initial research plan. The computational results indicate that escitalopram binds to the S438T mutant with an affinity comparable to the wild-type (−8.3 vs −8.1 kcal/mol). This finding provides a corrected foundation for understanding how SERT pocket mutations influence SSRI efficacy.

---

## 6. References

1. World Health Organization. (2023). *Depressive disorder (depression)*. WHO Fact Sheet.
2. Cipriani, A., *et al.* (2018). Comparative efficacy and acceptability of 21 antidepressant drugs. *Lancet*, 391(10128), 1357–1366.
3. Coleman, J.A., *et al.* (2016). X-ray structures of the human serotonin transporter. *Science*, 352(6290), 1478–1480.
4. Andersen, J., *et al.* (2011). Molecular basis for SSRI binding to SERT. *Journal of Biological Chemistry*, 286(32), 27834–27843.
5. Parsons, H.M. (1999). NERF algorithm. *Journal of Computational Chemistry*, 20(11), 1127–1137.
6. Trott, O., & Olson, A.J. (2010). AutoDock Vina. *Journal of Computational Chemistry*, 31(2), 455–461.
11. Plenge, P., *et al.* (2020). The mechanism of a high-affinity allosteric inhibitor of the serotonin transporter. *Nature Communications*, 11, 1491.

---

*Correspondence: fpalvarez23@gmail.com*
*Data and scripts: github.com/fabzy4L/DATA_ANALYTICS — bioinformatics/python_bio/Receptor_Design/RECEPTOR_DESIGN*
