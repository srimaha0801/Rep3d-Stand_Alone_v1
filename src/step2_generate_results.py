import os,re,time,warnings,sqlite3,string,openpyxl,json
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from collections import defaultdict

## Code to use 'residue_num.db' for residue number

# Dictionary mapping three-letter amino acid codes to one-letter codes
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Suppress PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)

###############################################################################
# GLOBAL DICTIONARY to store actual PDB residue numbers:
# residues_map[file][chain] = list of real residue numbers in order
###############################################################################


def get_user_coordinates(usr_coord):
    try:
        coordinates = [tuple(map(float, coord.split())) for coord in re.findall(r'\((.*?)\)', usr_coord)]
        if len(coordinates) < 5:
            print("Please enter at least 5 residues.")
        return coordinates
    except ValueError:
        print("Invalid input! Please enter valid coordinates.")


def generate_sequences(ca_atoms):
    sequences = {chain_id: "".join(three_to_one.get(resname, '?') for _, resname in ca_atoms[chain_id])
                 for chain_id in ca_atoms}
    unique_sequences = {}
    for chain_id, sequence in sequences.items():
        if sequence not in unique_sequences.values():
            unique_sequences[chain_id] = sequence
    return unique_sequences


def calculate_distances(ca_atoms):
    distances = defaultdict(list)
    for chain_id, residues in ca_atoms.items():
        for ca, residue_name in residues:
            if ca != '?':
                ca_coord = np.array(ca.get_coord())
                distance = np.linalg.norm(ca_coord)
                distances[chain_id].append((distance, residue_name))
            else:
                distances[chain_id].append(('?', residue_name))
    return distances


def distance_representation(distance):
    num_intervals_A = 84
    num_intervals_B = 10

    max_distance_A = 1000
    min_distance_B = 1000
    max_distance_B = 3000

    interval_length_A = max_distance_A / num_intervals_A
    interval_length_B = (max_distance_B - min_distance_B) / num_intervals_B

    if isinstance(distance, str):
        return '?'
    elif 0 <= distance <= max_distance_A:
        interval_index = int(distance / interval_length_A)
        symbols = (string.printable * (num_intervals_A // len(string.printable) + 1))[:num_intervals_A]
        return symbols[interval_index]
    elif min_distance_B <= distance <= max_distance_B:
        interval_index = int((distance - min_distance_B) / interval_length_B)
        symbols = (string.printable * (num_intervals_A + num_intervals_B // len(string.printable) + 1)
                   )[num_intervals_A:num_intervals_A + num_intervals_B]
        return symbols[interval_index]
    else:
        return '?'
    
    
def extract_ca_atoms(structure):
    ca_atoms = defaultdict(list)
    model = structure[0]
    for chain in model:
        chain_id = chain.get_id()
        for residue in chain:
            resname = residue.get_resname()
            if resname in three_to_one:
                if "CA" in residue:
                    ca_atoms[chain_id].append((residue["CA"], resname))
                else:
                    print(f"Missing CA atom in residue {residue.get_id()} in chain {chain_id} (Residue Name: {resname})")
                    ca_atoms[chain_id].append(('?', resname))
    return ca_atoms


def extract_one_symbolsuence(structure):
    sequences = defaultdict(str)
    model = structure[0]
    for chain in model:
        chain_id = chain.get_id()
        for residue in chain:
            resname = residue.get_resname()
            if resname in three_to_one:
                sequences[chain_id] += three_to_one.get(resname, '?')
    return sequences


def build_filter_query(filters):
    """Build WHERE clause dynamically based on filters."""
    conditions = []
    params = []
    valid_methods = [
        "x-ray diffraction",
        "electron microscopy",
        "solution nmr",
        "solid-state nmr"
    ]

    if "organism" in filters and filters["organism"].strip() != "":
        conditions.append("LOWER(pdb_header.organism) = ?")
        params.append(filters["organism"])

    if "method" in filters and filters["method"].strip() != "":
        method_value = filters["method"].strip()
        if method_value in valid_methods:
            conditions.append("LOWER(pdb_header.method) = ?")
            params.append(method_value)
        elif method_value == "hybrid":
            conditions.append("(LOWER(pdb_header.method) LIKE '%;%' OR LOWER(pdb_header.method) LIKE '%hybrid%')")
        elif method_value == "others":
            placeholders = ",".join("?" * len(valid_methods))
            conditions.append(f"(LOWER(pdb_header.method) NOT IN ({placeholders}) AND LOWER(pdb_header.method) NOT LIKE '%;%')")
            params.extend(valid_methods)

    if "resolution" in filters and filters["resolution"].strip() != "":
        parts = filters["resolution"].split()
        if parts[0] == "between":
            conditions.append("pdb_header.resolution BETWEEN ? AND ?")
            params.extend([float(parts[1]), float(parts[2])])
        elif parts[0] == "greater_than":
            conditions.append("pdb_header.resolution >= ?")
            params.append(float(parts[1]))
        elif parts[0] == "less_than":
            conditions.append("pdb_header.resolution <= ?")
            params.append(float(parts[1]))
        elif parts[0] == "equal":
            conditions.append("pdb_header.resolution = ?")
            params.append(float(parts[1]))

    if "rfac" in filters and filters["rfac"].strip() != "":
        parts = filters["rfac"].split()
        if parts[0] == "between":
            conditions.append("pdb_header.r_value BETWEEN ? AND ?")
            params.extend([float(parts[1]), float(parts[2])])
        elif parts[0] == "greater_than":
            conditions.append("pdb_header.r_value > ?")
            params.append(float(parts[1]))
        elif parts[0] == "less_than":
            conditions.append("pdb_header.r_value < ?")
            params.append(float(parts[1]))
        elif parts[0] == "equal":
            conditions.append("pdb_header.r_value = ?")
            params.append(float(parts[1]))
            
    
    where_clause = " AND ".join(conditions) if conditions else "1=1"
    # print(where_clause,params)
    return where_clause, params


def get_sequences_with_filters(symbols_db_cursor, filters, user_sequence_str):
    cursor = symbols_db_cursor

    where_clause, params = build_filter_query(filters)

    query = f"""
        SELECT sequences.file, sequences.chain, sequences.sequence
        FROM sequences
        JOIN pdb_header ON sequences.file = pdb_header.file
        WHERE {where_clause}
    """

    if user_sequence_str:
        query += " AND sequences.sequence LIKE ?"
        params.append(f"%{user_sequence_str}%")

    cursor.execute(query, params)
    
    results = cursor.fetchall()

    return results

def get_sequences_without_filters(symbols_db_cursor, user_sequence_str):
    cursor = symbols_db_cursor

    query = f"""
        SELECT sequences.file, sequences.chain, sequences.sequence
        FROM sequences WHERE sequences.sequence GLOB '*{user_sequence_str}*';
    """

    cursor.execute(query)
    
    results = cursor.fetchall()
    return results

def get_uniprot_id(uniprot,chain):
    uniprot_list = uniprot.split(',')

    for item in uniprot_list:
        if item.startswith(chain + "-"):
            return item.split("-", 1)[1]
    return "N/A"

def add_tsv_header(usr_coord,user_sequence_str,file):
    coord = usr_coord
    coord_len = len(user_sequence_str)
    file = file
    
    with open(file, "w") as f:
        f.truncate()
        f.write("# Title: Rep3D : An algorithm to identify structurally similar repeats\n")
        f.write("# Authors: Gurleen Kaur, Madhumathi Sanjeevi, Srimaha Gandhi and Kanagaraj Sekar\n")
        f.write("# Institution: IISc, Bangalore.\n\n")

        f.write(f"# Input coordinates: {coord}\n")
        f.write(f"# Number of coordinates: {coord_len}\n\n")
        
        f.write("# ===================== Results =====================\n\n")

def append_tsv_footer(output_path, execution_time,total_results,unique_pdb_ids):
    """
    Appends execution time information at the end of a TSV file.
    """
    with open(output_path, "a") as f:
        f.write("\n")
        f.write("# ===================== Summary =====================\n")
        f.write(f"# No. of. Results from Rep3d : {total_results}\n")
        f.write(f"# No. of PDB IDs : {len(unique_pdb_ids)}\n")
        f.write(f"# Execution time (seconds): {execution_time:.3f}\n")


def main(usr_coord,results_path,symbols_db_path,amino_acid_seq_db_path,filters):
    
    start_time = time.time()

    try:
        os.makedirs(results_path, exist_ok=True)
    except Exception as e:
        print(f"❌ Failed to create directory '{results_path}': {e}")
        return
        
    if not os.path.isfile(symbols_db_path):
        print(f"❌ The file '{symbols_db_path}' does not exist")
        return
    
    if not os.path.isfile(amino_acid_seq_db_path):
        print(f"❌ The file '{amino_acid_seq_db_path}' does not exist")
        return
      
    write_header_match = True
    
    symbols_db_conn = sqlite3.connect(symbols_db_path)
    symbols_db_cursor = symbols_db_conn.cursor()
    amino_acid_db_conn = sqlite3.connect(amino_acid_seq_db_path)
    amino_acid_seq_cursor = amino_acid_db_conn.cursor()
    
    # distance calculation and symbol assignment for user-provided coords
    user_coordinates = get_user_coordinates(usr_coord)
    user_sequence = []
    for i, coord in enumerate(user_coordinates, start=1):
        distance = np.linalg.norm(coord)
        representation = distance_representation(distance)
        user_sequence.append(representation)

    user_sequence_str = ''.join(user_sequence)
    print("User-Provided Residue Symbol Sequence:", user_sequence_str, "Length : ",len(user_sequence_str))
    
    match_pos_out = os.path.join(results_path,"A_Results.tsv")
    if os.path.exists(match_pos_out):
        add_tsv_header(usr_coord,user_sequence_str,match_pos_out)
        
    else:
        add_tsv_header(usr_coord,user_sequence_str,match_pos_out)
                   
    headers = os.path.join(results_path,"B_Headers.tsv")
    if os.path.exists(headers):
        add_tsv_header(usr_coord,user_sequence_str,headers)
    else:
        add_tsv_header(usr_coord,user_sequence_str,headers)    
    
    if filters:
        matched_symbols = get_sequences_with_filters(symbols_db_cursor, filters, user_sequence_str)
    else:
        matched_symbols = get_sequences_without_filters(symbols_db_cursor,  user_sequence_str)
    written_pdbchains = set()
    total_results = 0
    unique_pdb_ids = set()

    for file, chain, stored_sequence in matched_symbols:
        start_idx = 0
        while True:
            pos = stored_sequence.find(user_sequence_str, start_idx)

            if pos == -1:
                break

            start_idx = pos
            end_idx = pos+len(user_sequence_str)
            matching_letters = stored_sequence[start_idx:end_idx]
            amino_acid_seq_cursor.execute("SELECT sequence FROM sequences WHERE file=? AND chain=?", (file, chain))
            result = amino_acid_seq_cursor.fetchone()
            if result:
                amino_acid_seq = result[0]
                matching_amino_acids = amino_acid_seq[start_idx:end_idx]
                real_start = '?'
                real_end = '?'
                symbols_db_cursor.execute("SELECT real_resnum FROM residues_map WHERE file =? AND chain=?", (file, chain))
                row = symbols_db_cursor.fetchone()
                symbols_db_cursor.execute("select classification,organism,method,resolution,r_value,uniprot from pdb_header where file=?",(file,))
                header_row = symbols_db_cursor.fetchone()
                if header_row:
                    classification, organism, method, resolution, r_value, uniprot = header_row
                else:
                    classification = organism = method = uniprot  = "N/A"
                    resolution = r_value = 0
                
                if row:
                    res_num = row[0]
                    res_nums = res_num.split(',')
                    if end_idx - 1 < len(res_nums):
                        real_start = res_nums[start_idx]
                        real_end = res_nums[end_idx - 1]
                    else:
                        start_idx = pos +1
                        continue
                else:
                    start_idx = pos +1
                    continue
                unique_pdb_ids.add(file)
                print(f"{file}\t{chain}\t{real_start}\t{real_end}\t{matching_letters}\t{matching_amino_acids}")
                df_row = pd.DataFrame({
                    "PdbId_Chain": [f"{file}_{chain}"],
                    "Start_Position": [real_start],
                    "End_Position": [real_end],
                    "Sequence": [matching_amino_acids],
                    "Length": len(matching_amino_acids)
                })
                df_row.to_csv(match_pos_out,sep="\t",mode='a', index=False, header=write_header_match)
                total_results += 1

                pdb_chain = f"{file}_{chain}"
                
                if pdb_chain not in written_pdbchains:
                    df_header = pd.DataFrame({
                        "PdbId_Chain" : [f"{file}_{chain}"],
                        "UniProt_ID": [get_uniprot_id(uniprot,chain)],
                        "Method": [method],
                        "Classification": [classification],
                        "Organism": [organism.upper()],
                        "Resolution": [resolution],
                        "R-Value": [r_value],
                    })
                    df_header.to_csv(headers,sep="\t", mode='a', index=False, header=write_header_match)
                    written_pdbchains.add(pdb_chain)
                write_header_match = False

            else:
                print(f"Amino acid sequence not found for file: {file}, Chain: {chain}")
            start_idx = pos +1
            
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {end_time - start_time:.3f} seconds")
    append_tsv_footer(match_pos_out, execution_time,total_results,unique_pdb_ids)
    
    print("✅ Script execution completed successfully!\n")
    print(f"✅ Your output has been saved to: {match_pos_out}")
    
    
if __name__ == "__main__":
    
    usr_coord = "(93.700  58.178  27.351),(96.744  59.681  25.631),(96.149  57.652  22.469),(92.397  57.209  22.083),(91.176  60.109  24.295),(88.548  57.909  25.994),(87.960  57.747  29.736),(88.675  54.128  30.605),(89.402  52.270  33.815),(91.867  51.285  34.904),(94.807  52.999  33.234),(98.584  52.733  33.095),(100.343  55.919  34.205),(103.970  56.265  33.062),(106.697  57.575  35.378),(107.753  60.290  32.975),(107.819  63.973  33.842),(104.239  64.889  32.893),(103.153  61.417  34.056),(101.195  60.494  30.919),(98.140  58.252  31.152),(97.292  55.437  28.726),(94.398  52.959  28.763),(95.173  49.211  28.798),(93.911  48.581  25.243),(96.173  51.251  23.742),(99.026  49.981  25.858),(98.518  46.504  24.405),(98.534  47.885  20.872),(101.599  49.980  21.714),(103.693  47.053  22.986),(102.590  44.680  20.240),(103.405  47.251  17.555),(106.805  47.853  19.121),(107.379  44.091  19.244),(106.577  43.663  15.558),(108.863  46.537  14.556),(111.637  45.272  16.828),(111.688  48.673  18.542),(112.210  47.708  22.177),(114.818  49.646  24.114),(117.680  47.684  25.584),(118.211  48.797  29.151),(121.797  47.599  29.612),(121.223  47.109  33.331),(117.958  45.122  33.195),(116.510  43.618  30.002),(114.360  45.112  27.249),(111.630  47.694  27.850),(109.594  50.597  26.501),(109.390  54.252  27.575),(106.356  56.570  27.624),(104.432  56.484  24.327),(105.032  60.231  23.903),(108.781  59.671  23.880),(109.620  61.375  27.161)"
    # usr_coord = "(78.899  -3.971  9.342),(82.405  -5.330  8.754),(83.729  -1.907  7.755),(81.814  0.690  9.759),(80.533  -1.527  12.628),(77.017  -0.039  12.415),(73.811  -2.063  12.391),(72.227  -1.205  9.056),(69.524  -2.769  6.929),(69.686  -4.376  4.556),(73.176  -5.815  4.194),(75.243  -7.725  1.658),(76.524  -11.054  2.984),(79.492  -12.533  1.090),(79.706  -16.245  0.265),(83.085  -16.649  1.903),(83.808  -19.112  4.684),(82.905  -16.933  7.682),(80.171  -15.322  5.556),(81.135  -11.704  6.265),(78.477  -8.997  6.386),(78.822  -5.513  4.878),(76.313  -2.649  4.727),(75.345  -1.215  1.315),(77.112  2.135  1.849),(80.450  0.485  2.644),(79.947  -1.865  -0.263),(79.500  1.126  -2.567),(82.680  2.724  -1.250),(84.478  -0.617  -1.567),(83.594  -1.151  -5.241),(84.202  2.463  -6.223),(87.609  2.452  -4.529),(88.501  -0.794  -6.273),(87.350  0.702  -9.579),(89.541  3.775  -9.160),(92.591 1.684  -8.257),(91.927  -0.726  -11.123),(91.899  -3.591  -8.613),(89.088  -5.768  -9.944),(89.754  -9.485  -10.178),(89.490  -11.070  -13.588),(87.492  -14.257  -13.294),(88.664  -15.980  -16.483),(85.364  -17.833  -16.771),(82.963  -14.894  -16.295),(84.036  -11.233  -16.352),(85.486  -8.980  -13.661),(84.483  -9.132  -9.993),(85.449  -8.569  -6.367),(85.858  -11.001  -3.456),(85.009  -10.490  0.231),(86.418  -7.209  1.585),(88.430  -9.165  4.175),(90.287  -10.987  1.417),(88.889  -14.461  2.002)"
    # # 📌 Paste your symbols database path (symbols.db)
    # symbols_db_path = "/path/to/your/databases/symbols.db"

    # # 📌 Paste your amino acid sequence database (amino_acid_seq.db)
    # amino_acid_seq_db_path = "/path/to/your/databases/amino_acid_seq.db"

    # # 📌 Output directory where results will be stored
    # results_path = "/path/to/your/output/"
    
    # 📌 Paste your symbols database path (symbols.db)
    symbols_db_path = "/home/srimaha/Pictures/scripts/data/databases/coord_52/coord_52_symbols.db"

    # 📌 Paste your amino acid sequence database (amino_acid_seq.db)
    amino_acid_seq_db_path = "/home/srimaha/Pictures/scripts/data/databases/coord_52/coord_52_amino_acid_seq.db"

    # 📌 Output directory where results will be stored
    results_path = "/home/srimaha/Pictures/scripts/data/output/coord_52/"
    
    # Filters
    # If you want filters → put True, else False
    filter_choice = True

    if filter_choice:
        filters = {
            # filter: organism name or leave it ""
            "organism": "homo sapiens",
            # filter: experiment method (or "")
            "method": "x-ray diffraction",
            # filter: resolution condition (e.g., greater_than 2.0, between 1.0 2.5)
            "resolution": "greater_than 2.5",
            # filter: R-factor condition (or "")
            "rfac": "",
            # filter: list of PDB IDs(e.g., ["9rm2", "9otp", "9ofx", "9r96", "9j4p"]) or leave it []
            "pdb_ids": []
        }
    else:
        filters = {}  # no filters applied
            
    print("🚀 Script started...\n")
    main(usr_coord,results_path,symbols_db_path,amino_acid_seq_db_path,filters)