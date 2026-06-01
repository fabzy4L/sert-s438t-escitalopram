# Structural Comparison Report: LeuT (2A65) vs. hSERT (5I6Z)

## 1. Context and Significance
This report details the relationship between two foundational structures in the Neurotransmitter:Sodium Symporter (NSS) family. This comparison is critical for understanding the structural biology of the human Serotonin Transporter (hSERT) using its bacterial homolog, LeuT, as a high-resolution template.

## 2. Model Overview

| Feature | PDB 2A65 (LeuT) | PDB 5I6Z (hSERT) |
| :--- | :--- | :--- |
| **Protein** | Leucine Transporter (LeuT) | Human Serotonin Transporter (hSERT) |
| **Organism** | *Aquifex aeolicus* (Bacterium) | *Homo sapiens* (Human) |
| **Resolution** | **1.65 Å** (High Fidelity) | **4.53 Å** (Low Fidelity) |
| **Ligand** | Leucine (Natural Substrate) | (S)-Citalopram (SSRI Antidepressant) |
| **Status** | Historical "Gold Standard" | Modern Therapeutic Target |

## 3. Key Structural Insights

### 3.1 The "LeuT-Fold" Architecture
Both proteins share the conserved **LeuT-fold**, characterized by:
- **Scaffold Domain:** TMs 3, 4, 8, and 9 (relatively static).
- **Bundle Domain:** TMs 1, 2, 6, and 7 (undergo significant rotation/movement during transport).
- **Inverted Repeat:** The first 5 TMs are related to the second 5 TMs by a pseudo-twofold rotation axis in the membrane plane.

### 3.2 Conformational States
The comparison in the `2a65_v_516z.cxs` session likely highlights the transition between different steps of the transport cycle:
- **2A65 (Occluded State):** The protein is closed to both the extracellular and intracellular environments, trapping the substrate in the core.
- **5I6Z (Outward-Open State):** The extracellular vestibule is open to the synaptic cleft, allowing drugs like Escitalopram to enter and block the transporter.

### 3.3 The Allosteric Site (S2)
A major distinction discovered in the 5I6Z structure is the **allosteric site**.
- **Central Site (S1):** Where the primary drug molecule binds to block transport.
- **Allosteric Site (S2):** Located in the extracellular vestibule. In 5I6Z, a second citalopram molecule was found here, which "locks" the first molecule into the central site, explaining why certain SSRIs stay bound to the protein for a long time.

## 4. Why This Comparison Matters for Receptor Design
1. **Template Validation:** Because 2A65 has much higher resolution, it is used to validate the side-chain orientations in lower-resolution human models.
2. **Drug Mechanism:** Comparing these structures shows how an inhibitor (like citalopram in 5I6Z) prevents the protein from reaching the occluded state (seen in 2A65), thereby halting serotonin reuptake.
3. **Mutation Mapping:** Understanding the conserved residues between LeuT and hSERT allows us to predict the impact of mutations like **S438T** based on their equivalent positions in the high-resolution bacterial model.

---
*Theoretical supplement to the ChimeraX session: 2a65_v_516z.cxs*
