# ChimeraX Results Generation — SERT Receptor Design

This guide walks through generating quantitative results and publication-quality figures for the SERT S348 mutation project. All commands are entered in the ChimeraX command bar (bottom of the screen).

---

## Prerequisites

- ChimeraX installed and open
- All `.pdb` files available in `data/structures/`
- `output/escitalopram.pdb` available

---

## Step 1 — Load All Structures

Open the wild-type (chain A only), both mutants, and the generated ligand.

```
open ../data/structures/5i6z_A.pdb
open ../data/structures/5i6Z_S348T_3.pdb
open ../data/structures/5i6Z_S348_A.pdb
open ../output/escitalopram.pdb
```

Assign labels so you can track models by name. In the **Models panel** (right side), rename them:
- `#1` → WT (5i6z_A)
- `#2` → S348T
- `#3` → S348A
- `#4` → Escitalopram

---

## Step 2 — Superpose Structures (Structure Alignment)

Align both mutants onto the wild-type using MatchMaker. This corrects for any positional differences before comparison.

```
matchmaker #2 to #1
matchmaker #3 to #1
```

ChimeraX will report the **RMSD score** in the log for each alignment.

**Record in `Annotations.txt`:**
- WT vs S348T RMSD: ___ Å
- WT vs S348A RMSD: ___ Å

> Low RMSD (< 0.5 Å) = global fold preserved; differences are local to the mutation site.

---

## Step 3 — Calculate Local RMSD at the Binding Pocket

Measure RMSD specifically around residue 348 and its neighbors (binding pocket region, roughly residues 340–360).

```
rmsd #2:340-360 to #1:340-360
rmsd #3:340-360 to #1:340-360
```

**Record in `Annotations.txt`:**
- Binding pocket RMSD (WT vs S348T): ___ Å
- Binding pocket RMSD (WT vs S348A): ___ Å

---

## Step 4 — Visualize the Binding Pocket

Show ribbon for full structure, then switch to stick representation for residues near position 348.

```
# Ribbon for all models
style #1,2,3 ribbon

# Show sticks for binding pocket residues (adjust range as needed)
style #1:340-360 stick
style #2:340-360 stick
style #3:340-360 stick

# Color each model distinctly
color #1 cornflower blue
color #2 tomato
color #3 forest green
```

Zoom into the binding pocket:

```
view #1:348
```

---

## Step 5 — Place Escitalopram in the Binding Site

5I6Z is co-crystallized with Escitalopram, so the native binding pose is embedded in the structure. To show it:

```
# Open full 5i6z (contains the co-crystallized ligand ESC)
open ../data/structures/5i6z.pdb

# Select and display only the ligand (residue name ESC or LIG — check Model panel)
show #5:ESC
style #5:ESC stick
color #5:ESC yellow
```

If the residue name is different, use the **Select > Residues** panel to find the ligand by browsing the model.

To use your generated escitalopram instead:

```
# Show the RDKit-generated ligand
show #4
style #4 stick
color #4 orange
```

---

## Step 6 — Analyze Hydrogen Bonds

Find hydrogen bonds between the ligand and SERT residues in each structure.

```
# H-bonds between escitalopram and WT SERT
hbonds #5:ESC intersubmodel true reveal true

# H-bonds between escitalopram and S348T
hbonds #2 restrict #5:ESC intersubmodel true reveal true
```

ChimeraX will draw dashed lines for each H-bond and list them in the log.

**Record in `Annotations.txt`:**
- H-bond partners in WT: ___
- H-bond partners in S348T: ___
- H-bond partners in S348A: ___
- Any bonds gained or lost due to mutation: ___

---

## Step 7 — Analyze Contacts / Van der Waals Interactions

```
# Contacts within 4 Å of residue 348 in WT
contacts #1:348 radius 4.0 reveal true

# Repeat for mutants
contacts #2:348 radius 4.0 reveal true
contacts #3:348 radius 4.0 reveal true
```

Note which residues are within contact distance in each variant and whether the mutation opens or closes the pocket.

---

## Step 8 — Generate Molecular Surface of Binding Pocket

Visualize the shape and electrostatics of the binding pocket.

```
# Surface for WT only (cleaner view)
surface #1

# Color surface by electrostatic potential (requires APBS or Coulombic coloring)
coulombic #1
```

To show just the pocket region, hide the full surface and show a clipping plane:

```
clip front 10 fromCenter #1:348
```

---

## Step 9 — Export Figures

Save images for each key view. Use consistent orientation across all figures.

```
# Set view, then save
save ../output/fig1_overlay_WT_S348T_S348A.png width 1920 height 1080 supersample 3

# Binding pocket close-up
view #1:348
save ../output/fig2_binding_pocket_WT.png width 1920 height 1080 supersample 3

# WT with ligand
save ../output/fig3_WT_escitalopram.png width 1920 height 1080 supersample 3

# S348T with ligand
hide #1,3
show #2
save ../output/fig4_S348T_escitalopram.png width 1920 height 1080 supersample 3

# S348A with ligand
hide #2
show #3
save ../output/fig5_S348A_escitalopram.png width 1920 height 1080 supersample 3
```

> `supersample 3` increases image quality (3× antialiasing). Use `supersample 4` for publication.

---

## Step 10 — Save ChimeraX Session

Save the final session so the full analysis state is preserved.

```
save ../data/structures/Final_Analysis.cxs
```

---

## Step 11 — Record Findings in Annotations.txt

Open `docs/guides/Annotations.txt` and fill in:

```
=== SERT Mutation Analysis — ChimeraX Results ===

Date: [date]

--- RMSD (Global) ---
WT vs S348T: ___ Å
WT vs S348A: ___ Å

--- RMSD (Binding Pocket, residues 340-360) ---
WT vs S348T: ___ Å
WT vs S348A: ___ Å

--- Hydrogen Bonds with Escitalopram ---
WT:    [list residues]
S348T: [list residues — note gained/lost bonds]
S348A: [list residues — note gained/lost bonds]

--- Contacts within 4 Å of position 348 ---
WT:    [list]
S348T: [list]
S348A: [list]

--- Observations ---
[Describe conformational shifts, pocket opening/closing, sidechain reorientation]

--- Figures Generated ---
fig1 — Structural overlay (WT + S348T + S348A)
fig2 — Binding pocket close-up (WT)
fig3 — WT with Escitalopram
fig4 — S348T with Escitalopram
fig5 — S348A with Escitalopram
```

---

## Output Checklist

- [ ] RMSD values recorded (global + local)
- [ ] H-bond analysis complete for all 3 structures
- [ ] Contact analysis complete for all 3 structures
- [ ] 5 figures exported to `output/`
- [ ] `Annotations.txt` filled in
- [ ] `Final_Analysis.cxs` saved
