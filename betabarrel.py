# ----------------------------------------------------------------------------------------------------------------------
# betabarrel.py by Jaime Cardenas
# ----------------------------------------------------------------------------------------------------------------------
import requests
# ----------------------------------------------------------------------------------------------------------------------
def get_sheets():
    '''
    Takes a list of PDB identifiers (e.g. "1RVZ, 1gfl, 4fxF") or a text file of comma-separated PDB identifiers and
    returns two text files as output in your current working directory, one of them giving information about all
    beta-sheets in the list (sheets.txt), and the other of all beta barrels in the list (barrels.txt).

    Each file contains the relevant information of each beta sheet/barrel, such as
    1) its unique identifier, which in this case is the pdb ID and the sheet ID, separated by an underscore (e.g. 1GFL_AA)
    2) whether the sheet/barrel is made up of parallel or antiparallel strands
    3) the information of each strand that makes up the beta sheet
        a) amino acid sequence
        b) the chain where that strand is located
        c) the residues which that strand spans.

    Example output:
    5DPJ_AA1
        ANTI-PARALLEL
        ['VVPILVELDGDV', 'A', '11-22']
        ['HKFSVRGEGGD', 'A', '25-36']
        ['KLTLKFIC', 'A', '41-48']
        ['HMVLLEFVTAA', 'A', '217-227']
        ...
    '''
    # ------------------------------------------------------------------------------------------------------------------
    print("\nTakes a list of PDB identifiers (e.g. '1RVZ, 1gfl, 4fxF') or a text file of comma-separated PDB identifiers\n"
          "and returns two text files as output in your current working directory, one of them giving information about \n"
          "all beta-sheets in the list (sheets.txt), and the other of all beta barrels in the list (barrels.txt).\n\n"
        
        "Each file contains the relevant information of each beta sheet/barrel, such as \n"
        "1) its unique identifier, which in this case is the pdb ID and the sheet ID, separated by an underscore (e.g. 1GFL_AA)\n"
        "2) whether the sheet/barrel is made up of parallel or antiparralel strands\n"
        "3) the information of each strand that makes up the beta sheet\n"
        "   a) amino acid sequence\n"
        "   b) the chain where that strand is located\n"
        "   c) the residues which that strand spans\n\n"
   
        "Example output:\n"
        "5DPJ_AA1\n"
        "    ANTI-PARALLEL\n"
        "    ['VVPILVELDGDV', 'A', '11-22']\n"
        "    ['HKFSVRGEGGD', 'A', '25-36']\n"
        "    ['KLTLKFIC', 'A', '41-48']\n"
        "    ['HMVLLEFVTAA', 'A', '217-227']\n"
        "    ...")

    pdb_list = input("Enter a list of PDB Identifiers (pdb1, pdb2, pdb3, ...) or a text file with a list of PDBs: ")
    if ".txt" in pdb_list:
        with open(pdb_list, 'r') as p:
            pdb_list = p.read()
    pdb_list = pdb_list.replace(' ', '').split(',')
    # ------------------------------------------------------------------------------------------------------------------
    dics = get_sheets_and_barrels(pdb_list)
    write_file(dics[0], dics[1])
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
def get_pdb_info(pdb):
    aa_code = [['ALA', 'A'], ['CYS', 'C'], ['ASP', 'D'], ['GLU', 'E'], ['PHE', 'F'], ['GLY', 'G'], ['HIS', 'H'],
               ['ILE', 'I'], ['LYS', 'K'], ['LEU', 'L'], ['MET', 'M'], ['ASN', 'N'], ['PRO', 'P'], ['GLN', 'Q'],
               ['ARG', 'R'], ['SER', 'S'], ['THR', 'T'], ['VAL', 'V'], ['TRP', 'W'], ['TYR', 'Y']]
    url = "https://files.rcsb.org/view/%s.pdb" % pdb
    pdb_file = requests.get(url).content.decode("UTF-8").split("\n")
    # ------------------------------------------------------------------------------------------------------------------
    pdb_seq = []
    strands = []
    # ------------------------------------------------------------------------------------------------------------------
    for line in pdb_file:
        pdb_atom = line[13:17].replace(' ', '')
        if "SHEET" == line[:len("SHEET")]:
            sheet_id = f"{pdb.upper()}_{line[11:14].replace(' ', '')}"
            strand_id = int(line[7:10].replace(' ', ''))
            start = (line[21].replace(' ', ''), int(line[22:26].replace(' ', '')))
            end = (line[32].replace(' ', ''), int(line[33:37].replace(' ', '')))
            sense = int(line[38:40].replace(' ', ''))
            strands.append([sheet_id,strand_id, start, end, sense])
        # --------------------------------------------------------------------------------------------------------------
        if "ATOM" == line[:len("ATOM")] and "CA" == pdb_atom:
            pdb_chain = line[21].replace(' ', '')
            pdb_resi = int(line[22:26].replace(' ', ''))
            pdb_resi_3 = line[17:20].replace(' ', '')
            for threeletter in aa_code:
                if threeletter[0] == pdb_resi_3:
                    pdb_resi_1 = threeletter[1]
            pdb_seq.append((pdb_chain, pdb_resi, pdb_resi_1))
    # ------------------------------------------------------------------------------------------------------------------
    return strands, pdb_seq
# ----------------------------------------------------------------------------------------------------------------------
def get_sheets_and_barrels(pdb_list):
    sheets = {}
    barrels = {}
    # ------------------------------------------------------------------------------------------------------------------
    for pdb in pdb_list:
        strands = get_pdb_info(pdb)[0]
        pdb_seq = get_pdb_info(pdb)[1]
        sheets_pre = {}
        barrels_pre = {}
        delete = []
        # --------------------------------------------------------------------------------------------------------------
        for strand in strands:
            sheet_id_ = strand[0]
            sheet_info = strand[1:]
            if sheet_id_ not in sheets_pre:
                sheets_pre[sheet_id_] = [sheet_info]
            else:
                sheets_pre[sheet_id_].append(sheet_info)
        # ------------------------------------------------------------------------------------------------------------------
        for sheet in sheets_pre:
            start_resi_first = sheets_pre[sheet][0][1][1]
            end_resi_first = sheets_pre[sheet][0][2][1]
            start_resi_last = sheets_pre[sheet][-1][1][1]
            end_resi_last = sheets_pre[sheet][-1][2][1]
            for i in range(start_resi_first, end_resi_first):
                if i in range(start_resi_last, end_resi_last):
                    barrels_pre[sheet] = sheets_pre[sheet][:]
                    delete.append(sheet)
                    break
        # ------------------------------------------------------------------------------------------------------------------
        for i in delete:
            del sheets_pre[i]
        # ------------------------------------------------------------------------------------------------------------------
        for sheet in sheets_pre:
            if sheets_pre[sheet][1][3] == -1:
                sense_ = 'ANTI-PARALLEL'
            else:
                sense_ = 'PARALLEL'
            sheets[sheet] = [sense_]
            for strand in sheets_pre[sheet]:
                seq = ''
                start_resi = strand[1][1]
                end_resi = strand[2][1]
                chain = strand[1][0]
                for resi in range(start_resi, end_resi+1):
                    for aa in pdb_seq:
                        pdb_chain_ = aa[0]
                        pdb_resi_ = aa[1]
                        pdb_one_letter = aa[2]
                        if resi == pdb_resi_ and chain == pdb_chain_:
                            seq += pdb_one_letter
                sheets[sheet].append([seq, chain, f'{str(start_resi)}-{str(end_resi)}'])
        # ------------------------------------------------------------------------------------------------------------------
        for barrel in barrels_pre:
            try:
                if barrels_pre[barrel][1][3] == -1:
                    sense_ = 'ANTI-PARALLEL'
                else:
                    sense_ = 'PARALLEL'
            except IndexError:
                sense_ = 'N/A'
            barrels[barrel] = [sense_]
            for strand in barrels_pre[barrel]:
                seq = ''
                start_resi = strand[1][1]
                end_resi = strand[2][1]
                chain = strand[1][0]
                for resi in range(start_resi, end_resi+1):
                    for aa in pdb_seq:
                        pdb_chain_ = aa[0]
                        pdb_resi_ = aa[1]
                        pdb_one_letter = aa[2]
                        if resi == pdb_resi_ and chain == pdb_chain_:
                            seq += pdb_one_letter
                # ------------------------------------------------------------------------------------------------------
                barrels[barrel].append([seq, chain, f'{str(start_resi)}-{str(end_resi)}'])
        # --------------------------------------------------------------------------------------------------------------
    return sheets, barrels
    # ------------------------------------------------------------------------------------------------------------------
def write_file(sheets_dic, barrels_dic):
    with open("sheets.txt", "w") as s:
        for sheet in sheets_dic:
            s.write(f"{sheet}\n")
            strands = sheets_dic[sheet]
            for strand in strands:
                s.write("    %s\n" % str(strand))

    with open("barrels.txt", "w") as b:
        for barrel in barrels_dic:
            b.write(f"{barrel}\n")
            strands = barrels_dic[barrel]
            for strand in strands:
                b.write("    %s\n" % str(strand))


get_sheets()

