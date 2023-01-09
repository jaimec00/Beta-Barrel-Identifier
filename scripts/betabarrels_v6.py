# beta_barrelsV6.py by Jaime Cardenas
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# imports
import sys
import os
import requests
from collections import defaultdict
import pandas as pd
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
def main():
    # Get a list of PDB files and clean it
    pdb_list = sys.argv[1:] #input("Enter a list of PDB Identifiers (pdb1 pdb2 pdb3 ...) or a text file with a list of PDBs: ")
    if ".txt" in pdb_list[0]:
        with open(pdb_list[0], 'r') as p:
            pdb_list = p.read()
            pdb_list = pdb_list.replace(',', ' ').replace('\n', ' ').replace('  ', ' ').upper().split(' ')
    
    # Create empty list to store sheets and barrels
    beta_sheets = []
    beta_barrels = []
    
    # Iterate over the PDB's
    for pdb in pdb_list:
        # Get a list of all the beta structures
        beta_structures = get_beta_structures(pdb)

        # Seperate into sheets and barrels and append to corresponding list
        for beta_structure in beta_structures:
            if beta_structure.type == "BARREL":
                beta_barrels.append(beta_structure)
            else:
                beta_sheets.append(beta_structure)
    
    # Write all beta barrels to barrels.txt and all beta sheets to sheets.txt, and a third file (pymol_commands.txt) contains
    # commands that PyMOL can read so you can visualize the output
    if not os.path.isdir("out"): os.mkdir("out")
    with (
        open(os.path.join("out", "barrels.txt"), "w") as b, 
        open(os.path.join("out", "sheets.txt"), "w") as s, 
        open(os.path.join("out","pymol_commands.txt"), "w") as p
        ):

        p.write(",".join(pdb_list) + "\n")

        for beta_barrel in beta_barrels:            
            for strand in beta_barrel.strands:
                b.write(strand.sequence + "\n")
                p.write(f"red, model {strand.pdbID} and chain {strand.start_chain} and resi {strand.start_position}-{strand.end_position}" + "\n")
            b.write("\n")

        for beta_sheet in beta_sheets:
            for strand in beta_sheet.strands:
                s.write(strand.sequence + "\n")
                p.write(f"blue, model {strand.pdbID} and chain {strand.start_chain} and resi {strand.start_position}-{strand.end_position}" + "\n")
            s.write("\n")
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Define function to extract all of the beta structures within a PDB
def get_beta_structures(pdb):
    # Define the AA code to convert three-letter codes to one-letter codes
    aa_code = pd.DataFrame({"three_letter":['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 
                                            'ASN', 'PRO','GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR'], 
                            "one_letter":['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 
                                            'S', 'T', 'V', 'W', 'Y']
                            })

    # Fetch the PDB info from the RCSB
    url = f"https://files.rcsb.org/view/{pdb}.pdb"
    pdb_file = requests.get(url).content.decode("UTF-8").split("\n")
    # ------------------------------------------------------------------------------------------------------------------
    # Create empty lists to store the PDB sequence and the list of beta strands in the PDB
    pdb_seq = []
    beta_strands = []

    # Iterate over each line
    for line in pdb_file:
        # Each line containing "SHEET" is a single beta strand. Extract the info from that line
        if "SHEET" == line[:len("SHEET")]:
            pdb_id = pdb
            sheet_id = line[11:14].replace(' ', '')
            strand_id = int(line[7:10].replace(' ', ''))
            start_chain = (line[21].replace(' ', ''))
            start_position = int(line[22:26].replace(' ', ''))
            end_chain = (line[32].replace(' ', ''))
            end_position = int(line[33:37].replace(' ', ''))
            sense = line[38:40].replace(' ', '')

            # Make a BetaStrand object with the appropriate info and append it to the beta_strand list
            beta_strands.append( BetaStrand(pdb_id, sheet_id, strand_id, start_chain, start_position, end_chain, end_position, sense) )

        pdb_atom = line[13:16].replace(' ', '')
        # Each line containing "ATOM" and whose atom is "CA" (the alpha carbon) represents a residue in the protein sequence. 
        # Extract the info
        if ("ATOM" == line[:len("ATOM")] or "HETATM" == line[:len("HETATM")]) and "CA" == pdb_atom:
            chain = line[21].replace(' ', '')
            position = int(line[22:26].replace(' ', ''))
            
            # This is the three letter code for the amino acid
            aa_three = line[17:20].replace(' ', '')
            # If the amino acid is canonical, replace it with the corresponding one letter code
            if aa_three in aa_code.three_letter.values: 
                aa_one = aa_code.loc[aa_code.three_letter == aa_three, "one_letter"].values[0]
            # If the amino acid is non-canonical, replace it with "X"
            else: 
                aa_one = "X"
            
            # Make an Amino_Acid object and append it to the pdb_seq list
            pdb_seq.append( Amino_Acid(chain, position, aa_one) )
    # ------------------------------------------------------------------------------------------------------------------
    # Create an empty dictionary to store the beta structures (the keys) made up of beta strands (the list within the key)
    beta_structures_pre = defaultdict(list)
    for strand in beta_strands:
        
        # For each strand, find the starting and ending index in the PDB sequence 
        strand.start_idx = [idx for idx, aa in enumerate(pdb_seq) if aa.position == strand.start_position and aa.chain == strand.start_chain][0]
        strand.end_idx = [idx for idx, aa in enumerate(pdb_seq) if aa.position == strand.end_position and aa.chain == strand.end_chain][0]
        strand.idx_positions = [*range(strand.start_idx, strand.end_idx+1)]
        
        # Find the sequence of each strand using the idx positions found in the last few lines
        strand.sequence = "".join([pdb_seq[idx].aa for idx in strand.idx_positions])
        
        # Add each strand to its corresponding beta structure by grouping strand with identical sheetID's
        beta_structures_pre[f"{strand.pdbID}_{strand.sheetID}"].append(strand)
        
    # Create a list to store the output of this function (the list of beta structures)
    beta_structures = []
    for structure in beta_structures_pre:
        beta_structures.append( BetaStructure( beta_structures_pre[structure] ) )

    return beta_structures
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Define the Amino_Acid class
class Amino_Acid:
    def __init__(self, chain, position, aa):
        self.chain = chain
        self.position = position
        self.aa = aa
        
    def __str__(self):
        return self.aa

    def __repr__(self):
        return self.aa
# ----------------------------------------------------------------------------------------------------------------------
# Define BetaStrand class
class BetaStrand:
    def __init__(self, pdbID, sheetID, strandID, start_chain, start_position, end_chain, end_position, sense):
        self.pdbID = pdbID
        self.sheetID = sheetID
        self.strandID = strandID
        self.start_chain = start_chain
        self.start_position = start_position
        self.end_chain = end_chain
        self.end_position = end_position
        if sense == "-1":
            self.sense = "ANTI-PARALLEL"
        elif sense == "1":
            self.sense = "PARALLEL"
        else:
            self.sense = "N/A"

    def __str__(self):
        return f"{self.pdbID}_{self.sheetID}_{self.strandID}"

    def __repr__(self):
        return f"{self.pdbID}_{self.sheetID}_{self.strandID}"
# ----------------------------------------------------------------------------------------------------------------------
# Define BetaStructure Class
class BetaStructure:
    def __init__(self, strands):
        self.strands = strands
        self.pdbID = strands[0].pdbID
        self.sheetID = strands[0].sheetID

        # Create a list that contains every beta strand's index values in the PDB sequence
        all_idx = []
        for strand in self.strands:
            all_idx += strand.idx_positions

        # If the set(removes all repeats) is shorter than the original list, it is a barrel. Else, it is a sheet
        if len( set(all_idx) ) < len( all_idx ):
            self.type = "BARREL"
        else:
            self.type = "SHEET"

    def __str__(self):
        return f"{self.pdbID}_{self.sheetID}"

    def __repr__(self):
        return f"{self.pdbID}_{self.sheetID}"
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Run the main function
if __name__ == "__main__":
    main()

