# Beta-Helix-Design
Program for remodelling and design of tau beta-helix protein scaffold using PyRosetta.

This program takes as input a 5-repeat $\beta$-helix AlphaFold model ("bhp_model.pdb") and applies to the AlphaFold model a blueprint file "bpfile", wherein the last 15 residues are remodelled to a $\alpha$-helical conformation, using Rosetta's Remodel mover. Amino acids in the $\alpha$-helix and 3-residue linkers are designed using the FastDesign protocol in Rosetta. Backbone, bond angles, bond lengths, torsional angles, are all free to move during the design process. 10,000 decoys (i.e. designed protein candidates) are generated. 
