# Beta-Barrel-Identifier

beta_barrel_V5.py
Takes a list of PDB identifiers (e.g. "1RVZ, 1gfl, 4fxF") or a text file of comma-separated PDB identifiers and returns two text files as output in your current working directory, one of them giving the sequence of each strand within all beta barrels in the given pdbs, and the other giving the sequence of each strand within all beta sheet in the given pdbs

betabarrels_pymol_cmds.py
same as beta_barrel_V5.py except that it outputs commands recognized by pymol into a text file called pymol_barrel_test.txt to color barrels and sheets differently

betabarrels_pymol_test.py
PyMOL script to color the pdbs entered in the betabarrels_pymol_cmds.py script to visualize the barrels and sheets
