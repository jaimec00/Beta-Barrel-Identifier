# Beta-Barrel-Identifier

Takes a list of PDB identifiers (e.g. "1RVZ, 1gfl, 4fxF") or a text file of comma-separated PDB identifiers and returns two text files as output in your current working directory, one of them giving information about all beta-sheets in the list (sheets.txt), and the other of all beta barrels in the list (barrels.txt).

Each file contains the relevant information of each beta sheet/barrel, such as
1) its unique identifier, which in this case is the pdb ID and the sheet ID, separated by an underscore (e.g. 1GFL_AA)
2) whether the sheet/barrel is made up of parallel or antiparallel strands
3) the information of each strand that makes up the beta sheet
    a) amino acid sequence
    b) the chain where that strand is located
    c) the residues which that strand spans.

Example output:\n
5DPJ_AA1
    ANTI-PARALLEL # This beta barrel is made up of antiparallel strands
    ['VVPILVELDGDV', 'A', '11-22'] #[Amino acid sequence, chain letter, residue range]
    ['HKFSVRGEGGD', 'A', '25-36']
    ['KLTLKFIC', 'A', '41-48']
    ['HMVLLEFVTAA', 'A', '217-227']
    ...
