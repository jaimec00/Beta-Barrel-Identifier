# Beta-Barrel-Identifier

betabarrels_V5.py:
Takes a list of PDB identifiers (e.g. "1RVZ, 1gfl, 4fxF") or a text file of comma-separated PDB identifiers and returns two text files as output in your current working directory (or a specified absolute path), one of them giving the sequence of each strand within all beta barrels in the given PDBs, and the other giving the sequence of each strand within all beta sheet in the given PDBs.

betabarrels_pymol_cmds.py:
Same as betabarrels_V5.py except that it outputs commands recognized by PyMOL into a text file called pymol_barrel_test.txt to color barrels and sheets differently

betabarrels_pymol_test.py:
PyMOL script to color the PDBs entered in the betabarrels_pymol_cmds.py script to visualize the barrels and sheets. Sheets are blue, barrels are red, and everything else is green.
