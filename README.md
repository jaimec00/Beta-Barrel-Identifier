# Beta-Barrel-Identifier


betabarrels_V6.py:
Takes a list of PDB identifiers (e.g. "1RVZ 1gfl 4fxF") or a text file of PDB identifiers and returns three text files as output in your current working directory, one of them giving the sequence of each strand within all beta barrels in the given PDBs, another giving the sequence of each strand within all beta sheets in the given PDBs, and the final text file contains commands that are recognized by PyMOL (through betabarrels_pymol_test.py) to visualize the output.

example use: 

python3 betabarrels_V6.py 1rvz 1gfl 3rbm
python3 betabarrels_V6.py ./input.txt

betabarrels_pymol_test.py:
PyMOL script to color the PDBs entered in betabarrels_V6.py to visualize the barrels and sheets. Sheets are blue, barrels are red, and everything else is green.

example use (within PyMOL console):

run betabarrels_pymol_test.py

