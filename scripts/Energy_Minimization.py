# 1. Preparation: Add Hydrogens and Charges
addh
addcharge

# 2. Select the Region of Interest (TM6 Helix 340-360)
# This includes the S348T mutation site and surrounding residues
select :340-360

# 3. Define the Constraints
# We fix the rest of the protein to maintain global structural integrity (RMSD ~0.030A)
# while allowing the local environment to relax
# Note: Ensure the model ID matches your current session (defaulting to #1)
freeze ~:340-360

# 4. Energy Minimization (Steepest Descent + Conjugate Gradient)
# nsteps: 100 iterations of Steepest Descent to remove major clashes
# cgsteps: 100 iterations of Conjugate Gradient for refinement
# stepsize: 0.02 (Standard for protein systems)
minimize nsteps 100 cgsteps 100 stepsize 0.02 interval 10

# 5. Post-Minimization Analysis
# Calculate the displacement of the Threonine methyl group
select :348
measure rmsd #1:340-360 to_structure #1:340-360