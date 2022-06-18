# beta_barrelsV5.py by Jaime Cardenas
# ----------------------------------------------------------------------------------------------------------------------
import requests
from collections import defaultdict


# ----------------------------------------------------------------------------------------------------------------------
# Define main function
def main():
    # Get list of pdb files and clean
    pdb_list = input("Enter a list of PDB Identifiers (pdb1, pdb2, pdb3, ...) or a text file with a list of PDBs: ")
    if ".txt" in pdb_list:
        with open(pdb_list, 'r') as p:
            pdb_list = p.read()
    pdb_list = pdb_list.replace(' ', '').replace('\n', ',').split(',')

    # Setup empty list for beta sheets and beta barrels
    beta_sheets = list()
    beta_barrels = list()

    # Append the beta barrels and beta sheets from each pdb in the pdb list into their corresponding list from above
    for pdb in pdb_list:
        pdb = pdb.upper()
        pdb_info = get_pdb_info(pdb)

        # Get the sequence of the protein from the PDB
        pdb_seq = pdb_info[1]
        # Get the strand information of the protein from the PDB
        beta_strands = pdb_info[0]

        # Define the strands that make up a beta structure as either a sheet or a barrel
        beta_structures = is_barrel(beta_strands, pdb_seq)
        beta_sheets_pre = beta_structures[0]
        beta_barrels_pre = beta_structures[1]

        # Get the sequence from the strand object for each strand in the lists defined above
        for sheet in beta_sheets_pre:
            for strand in sheet.strands:
                seq = strand.get_strand_sequence(pdb_seq)[0]
                if seq is not None:
                    beta_sheets.append(seq)
        for barrel in beta_barrels_pre:
            for strand in barrel.strands:
                seq = strand.get_strand_sequence(pdb_seq)[0]
                if seq is not None:
                    beta_barrels.append(seq)

    # Write the sequences to their appropriate files
    write_file(beta_sheets, beta_barrels)
# ----------------------------------------------------------------------------------------------------------------------
# Define BetaStrand class
class BetaStrand:
    def __init__(self, pdbID, sheetID, strandID, start_chain, start_resi, end_chain, end_resi, sense):
        self.pdbID = pdbID
        self.sheetID = sheetID
        self.strandID = strandID
        self.start_chain = start_chain
        self.start_resi = start_resi
        self.end_chain = end_chain
        self.end_resi = end_resi
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

    # Define a function of the BetaStrand class to get its sequence
    def get_strand_sequence(self, pdb_sequence):
        # Create placeholders
        start_idx = "PlaceHolder"
        end_idx = "PlaceHolder"

        for idx, residue in enumerate(pdb_sequence):
            chain = residue[0]
            pos = residue[1]

            # Find the index value of the starting residue and of the ending residue
            if self.start_chain == chain and self.start_resi == pos:
                start_idx = idx
            if self.end_chain == chain and self.end_resi == pos:
                end_idx = idx
                break

        # If the above code failed to find the starting and ending residues (e.g. the line is HETATM), then the
        # strand sequence and the strand sequence index are set to none, else the strand sequence and its index are
        # assigned their values
        if start_idx != "PlaceHolder" and end_idx != "PlaceHolder":
            strand_sequence_idx = [*range(start_idx, end_idx+1)]
            strand_sequence = "".join([pdb_sequence[i][2] for i in strand_sequence_idx])
        else:
            strand_sequence_idx = None
            strand_sequence = None

        return strand_sequence, strand_sequence_idx
# ----------------------------------------------------------------------------------------------------------------------
# Define BetaStructure Class
class BetaStructure:
    def __init__(self, pdb_sequence, strands):
        self.strands = strands
        self.pdbID = strands[0].pdbID
        self.sheetID = strands[0].sheetID

        # This is done so that if a structure is made up of just one strand, it will not be assigned as a beta barrel
        if len(self.strands) == 1:
            self.type = "SHEET"
        else:
            # Get a list of lists, with each inner list containing the index positions of a single strand's sequence.
            # This is done for each strand in the beta structure
            strand_sequence_idx_list = list()
            for strand in self.strands:
                strand_sequence_idx_list.append( strand.get_strand_sequence(pdb_sequence)[1] )

            # Create a counter to check if any of the index positions of the strands within a beta structure overlap
            betas_temp = 0
            # for each strand in the structure, check if any of that strand's index positions are in another strand
            # within the same structure
            for n in range( len(strand_sequence_idx_list) ):
                # Get the index positions of the n strand
                strand_sequence_idx_n = strand_sequence_idx_list[n]
                # Get the index positions of all the strands that are not n
                strand_sequence_idx_list_not_n = [strand_sequence_idx_list[i] for i in range(len(strand_sequence_idx_list))
                                             if i != n]


                for position_idx_n in strand_sequence_idx_n:
                    for strand_sequence_idx_not_n in strand_sequence_idx_list_not_n:
                        if position_idx_n in strand_sequence_idx_not_n:
                            betas_temp += 1            

            if betas_temp > 0:
                self.type = "BARREL"
            else:
                self.type = "SHEET"

    def __str__(self):
        return f"{self.pdbID}_{self.sheetID}"

    def __repr__(self):
        return f"{self.pdbID}_{self.sheetID}"
# ----------------------------------------------------------------------------------------------------------------------
# Define function to get a list of beta strands of the input PDB and a list of tuples, with each tuple containing
# the chain, residue number, and one-letter code of the amino acid at that position
def get_pdb_info(pdb):
    aa_code = [['ALA', 'A'], ['CYS', 'C'], ['ASP', 'D'], ['GLU', 'E'], ['PHE', 'F'], ['GLY', 'G'], ['HIS', 'H'],
               ['ILE', 'I'], ['LYS', 'K'], ['LEU', 'L'], ['MET', 'M'], ['ASN', 'N'], ['PRO', 'P'], ['GLN', 'Q'],
               ['ARG', 'R'], ['SER', 'S'], ['THR', 'T'], ['VAL', 'V'], ['TRP', 'W'], ['TYR', 'Y']]
    url = "https://files.rcsb.org/view/%s.pdb" % pdb
    pdb_file = requests.get(url).content.decode("UTF-8").split("\n")
    # ------------------------------------------------------------------------------------------------------------------
    pdb_seq = list()
    beta_strands = list()

    for line in pdb_file:
        if "SHEET" == line[:len("SHEET")]:
            pdb_id = pdb
            sheet_id = line[11:14].replace(' ', '')
            strand_id = int(line[7:10].replace(' ', ''))
            start_chain = (line[21].replace(' ', ''))
            start_resi = int(line[22:26].replace(' ', ''))
            end_chain = (line[32].replace(' ', ''))
            end_resi = int(line[33:37].replace(' ', ''))
            sense = int(line[38:40].replace(' ', ''))

            beta_strands.append(
                BetaStrand(pdb_id, sheet_id, strand_id, start_chain, start_resi, end_chain, end_resi, sense) )

        pdb_atom = line[13:16].replace(' ', '')
        if "ATOM" == line[:len("ATOM")] and "CA" == pdb_atom:
            pdb_chain = line[21].replace(' ', '')
            pdb_resi = int(line[22:26].replace(' ', ''))
            pdb_aa_three_letter = line[17:20].replace(' ', '')
            for idx, threeletter in enumerate(aa_code):
                if threeletter[0] == pdb_aa_three_letter:
                    pdb_aa_one_letter = threeletter[1]
                    pdb_seq.append((pdb_chain, pdb_resi, pdb_aa_one_letter))
                    break
                elif idx == 19:
                    pdb_seq.append((pdb_chain, pdb_resi, "X"))

    return beta_strands, pdb_seq
# ----------------------------------------------------------------------------------------------------------------------
# Define function to take a list of beta strands and categorize them as either part of a beta barrel or part
# of a beta sheet
def is_barrel(beta_strands, pdb_seq):
    beta_structures_pre = defaultdict(list)
    beta_structures = list()
    beta_sheets = list()
    beta_barrels = list()

    for strand in beta_strands:
        key = f"{strand.pdbID}_{strand.sheetID}"
        beta_structures_pre[key].append(strand)

    for structure in beta_structures_pre:
        beta_structures.append(BetaStructure(pdb_seq, beta_structures_pre[structure]))

    for structure in beta_structures:
        if structure.type == "SHEET":
            beta_sheets.append(structure)
        elif structure.type == "BARREL":
            beta_barrels.append(structure)

    return beta_sheets, beta_barrels
# ----------------------------------------------------------------------------------------------------------------------
# Define a function to write all the sequences of beta barrels to a file called barrels.txt, and all
# the sequences of beta sheets to a file called sheets.txt
def write_file(beta_sheets, beta_barrels):
    with open("sheets.txt", "w") as s:
        for strand_seq in beta_sheets:
            s.write(strand_seq + '\n')

    with open("barrels.txt", 'w') as b:
        for strand_seq in beta_barrels:
            b.write(strand_seq + '\n')
# ----------------------------------------------------------------------------------------------------------------------
# Run the main function at the top of this file
main()
