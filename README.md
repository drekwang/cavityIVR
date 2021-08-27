# cavityIVR
Simple notebook for computing survival curves of a cavity-modified unimolecular dissociation reaction.
The code in this repo form the foundation of the paper by D. S. Wang et al., "Cavity-modified unimolecular dissociation reactions via intramolecular vibrational energy redistribution."

Download all the files into a single folder.
Run 'main.m' in MATLAB.
This script should compute the survival probability vs. time curve for N sets of given conditions.
These conditions including the cavity field strength and angle, how many initial states to consider, the dipole moment of the molecule, etc.
Each set of conditions are fed into 'dissociation_cavity_function' (a function that computes unimolecular dissociations in a cavity).
This function 1) computes all the initial states with an identical energy given some constraints, and 2) propagates each and notes when the molecule dissociates.
This information is used to plot survival vs. time curves.
In the paper and while doing the resesarch, we do a few other things that can be easily implemented based on the code here.
1) We plot Fourier transforms of the trajectories. Just apply fft() to the time-dependent trajectories.
2) We explore what happens with varying initial energies. Just change the hard-coded value of 34 kcal/mol.
3) We explored what happens with and without the ultrastrong coupling terms in the equations of motion. The terms are labelled in 'odefun_Gmattrix_theta.m'. Just comment them out.

If you have further questions about using this code, please feel free to contact me.
