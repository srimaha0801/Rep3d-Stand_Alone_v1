import os,sys 
import time
import warnings
import numpy as np
import pandas as pd
import sqlite3,re
import string,gzip,shutil
from Bio.PDB import PDBParser, MMCIFParser, PDBExceptions
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB import parse_pdb_header
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import tempfile
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from get_fasta import fasta_seq_main

files = []
seen_ids = set()

VALID_EXT = ('.pdb', '.cif', '.ent')

def extract_id(filename):
    return filename.split('.')[0].lower()

def add_file_or_unzip(path):
    filename = os.path.basename(path)
    pdb_id = extract_id(filename)
    # Accept .ent
    if filename.endswith('.ent'):
        files.append(path)
        seen_ids.add(pdb_id)
        return

    # Accept .cif
    if filename.endswith('.cif'):
        files.append(path)
        return

    # Accept .pdb only if .ent is not present
    if filename.endswith('.pdb'):
        if pdb_id not in seen_ids:
            files.append(path)
        return

    # Accept compressed formats
    if filename.endswith(tuple(ext + '.gz' for ext in VALID_EXT)):
        extracted = os.path.join(tempfile.gettempdir(), filename[:-3])
        if not os.path.exists(extracted):
            with gzip.open(path, "rb") as f_in, open(extracted, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        add_file_or_unzip(extracted)
        return

# Dictionary mapping three-letter amino acid codes to one-letter codes
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Suppress PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)

residues_map = defaultdict(dict)

def parse_structure(file):
    parser = None
    if file.lower().endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(PERMISSIVE=True, QUIET=True)

    try:
        structure = parser.get_structure("protein", file)
        if len(structure) == 0:
            print(f"Skipping empty structure: {file}")
            return None
        return structure
    except (PDBExceptions.PDBConstructionException, ValueError) as e:
        print(f"Skipping file due to parse error: {file} -> {e}")
        return None


def normalize_pdb_name(filepath):
    name = os.path.basename(filepath).lower()
    if name.startswith("pdb"):
        name = name[3:]
    if name.endswith(".ent") or name.endswith(".cif") or name.endswith(".pdb"):
        name = name[:-4]
    return name


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


def extract_one_letter_sequence(structure):
    sequences = defaultdict(str)
    model = structure[0]
    for chain in model:
        chain_id = chain.get_id()
        for residue in chain:
            resname = residue.get_resname()
            if resname in three_to_one:
                sequences[chain_id] += three_to_one.get(resname, '?')
    return sequences
  

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
        try:
            return symbols[interval_index]
        except IndexError:
            return '?'
    elif min_distance_B <= distance <= max_distance_B:
        interval_index = int((distance - min_distance_B) / interval_length_B)
        symbols = (string.printable * (num_intervals_A + num_intervals_B // len(string.printable) + 1)
                   )[num_intervals_A:num_intervals_A + num_intervals_B]
        try:
            return symbols[interval_index]
        except IndexError:
            return '?'
    else:
        return '?'


def safe_float(value):
    try:
        if isinstance(value, list):
            value = value[0]
        if value == "?" or value is None:
            return 0.0
        return float(value)
    except:
        return 0.0

def get_pdb_header_info_cif(file_path):
    organism = "N/A"
    r_value = 0.0
    chains = set()
    classification = "N/A"
    method = "N/A"
    resolution = 0.0
    uniprot_ids = "N/A"
    # print(file_path)
    try:
        cif_dict = MMCIF2Dict(file_path)

        # Experimental method
        method = ",".join(cif_dict.get("_exptl.method", ["?"])).upper()

        # Resolution
        if "_refine.ls_d_res_high" in cif_dict:
            resolution = safe_float(cif_dict.get("_refine.ls_d_res_high", "?"))

        # Classification / Keywords
        if "_struct_keywords.pdbx_keywords" in cif_dict:
            classification = cif_dict["_struct_keywords.pdbx_keywords"][0]
        elif "_struct_keywords.text" in cif_dict:
            classification = cif_dict["_struct_keywords.text"][0]

        # Organism
        
        organism_keys = [
            "_entity_src_gen.pdbx_gene_src_scientific_name",
            "_entity_src_nat.pdbx_organism_scientific",
            "_entity_src_gen.pdbx_organism_scientific"
        ]

        organisms = []

        for key in organism_keys:
            if key in cif_dict:
                for org in cif_dict[key]:
                    if org and org != "?":
                        # clean extra text
                        org = org.strip("'").strip('"')

                        # remove subsp.* etc
                        org = re.sub(r'\s+subsp.*', '', org, flags=re.IGNORECASE)

                        # remove anything inside parentheses
                        if "(" in org:
                            org = org.split("(")[0].strip()

                        organisms.append(org)

        # Unique list
        organisms = list(set(organisms))
        organism = "; ".join(organisms)
        # print(organisms)
        # organism_keys = [
        #     "_entity_src_gen.pdbx_gene_src_scientific_name",
        #     "_entity_src_nat.pdbx_organism_scientific",
        #     "_entity_src_gen.pdbx_organism_scientific"
        # ]
        # for key in organism_keys:
        #     if key in cif_dict:
        #         organism = cif_dict[key][0]
        #         # clean up subsp. or extra
        #         organism = re.sub(r'\s+subsp.*', '', organism, flags=re.IGNORECASE)
        #         if "(" in organism:
        #             organism = organism.split("(")[0].strip()
        #             break
                    
        # Chains
        if "_atom_site.label_asym_id" in cif_dict:
            chains = set(cif_dict["_atom_site.label_asym_id"])
        elif "_atom_site.auth_asym_id" in cif_dict:
            chains = set(cif_dict["_atom_site.auth_asym_id"])

        # R-value
        if "_refine.ls_R_factor_R_work" in cif_dict:
            r_value = safe_float(cif_dict.get("_refine.ls_R_factor_R_work", "?"))

        # UniProt IDs (if present)
        if "_struct_ref_seq.pdbx_db_accession" in cif_dict and "_struct_ref_seq.pdbx_strand_id" in cif_dict:
            uniprot_ids_list = []
            accessions = cif_dict["_struct_ref_seq.pdbx_db_accession"]
            strands = cif_dict["_struct_ref_seq.pdbx_strand_id"]
            for chain, acc in zip(strands, accessions):
                uniprot_ids_list.append(f"{chain}-{acc}")
            uniprot_ids = ",".join(uniprot_ids_list)

        chains_str = ",".join(sorted(chains)) if chains else "?"

        return {
            "file": os.path.basename(file_path),
            "organism": organism,
            "method": method,
            "resolution": resolution,
            "r_value": round(r_value * 100, 3),
            "chains": chains_str,
            "classification": classification,
            "uniprot_ids": uniprot_ids,
            # "organisms":organisms
        }

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return {
            "file": os.path.basename(file_path),
            "organism": organism,
            "method": method,
            "resolution": 0.0,
            "r_value": 0.0,
            "chains": chains_str,
            "classification": classification,
            "uniprot_ids": uniprot_ids
        }


def get_pdb_header_info(file_path):
    organism = "N/A"
    r_value = 0
    chains = set()
    classification = "N/A"
    method = "N/A"
    resolution = 0

    try:
        with open(file_path) as f:
            lines = f.readlines()
 
        # Extract classification from HEADER line
        for line in lines:
            if line.startswith("HEADER"):
                classification = line[10:50].strip()
                break

        # Extract method and resolution from REMARK 2 or other lines
        with open(file_path) as f:
            header = parse_pdb_header(f)
        method = header.get("structure_method") or "N/A"
        resolution = header.get("resolution") or 0

        # Extract organism from SOURCE lines
        organisms = []
        pattern_subsp = re.compile(r'\s+subsp.*', re.IGNORECASE)
        for line in lines:
            if line.startswith("SOURCE") and "ORGANISM_SCIENTIFIC" in line:
                
                # Extract text after ORGANISM_SCIENTIFIC:
                organism = line.split("ORGANISM_SCIENTIFIC:")[1].strip(" ;\n")

                # Remove parentheses content
                if "(" in organism:
                    organism = organism.split("(")[0].strip()

                # Remove subsp. or variants
                organism = pattern_subsp.sub('', organism).strip()
                organism = organism.upper()
                organisms.append(organism)
                
        organisms = list(set(organisms))
        organism = "; ".join(organisms)
        # Extract R-value from REMARK 3 lines
        for line in lines:
            if line.startswith("REMARK   3   R VALUE"):
                match = re.search(r"R VALUE *\(WORKING SET\) *: *([\d\.]+)", line)
                if match:
                    r_value = float(match.group(1))
                    break
                else:
                    match = re.search(r"R VALUE *: *([\d\.]+)", line)
                    if match:
                        r_value = float(match.group(1))
                        break

        # Extract chains from ATOM/HETATM lines
        for line in lines:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                chains.add(line[21].strip())
        chains_str = ",".join(sorted(chains)) if chains else "N/A"
        
        return {
            "file": os.path.basename(file_path),
            "organism": organism,
            "method": method.upper(),
            "resolution": resolution,
            "r_value":  round(r_value * 100, 3),
            "chains": chains_str,
            "classification": classification,
            "uniprot_ids": format_uniprot_ids(extract_uniprot_ids_from_pdb(file_path))
        }
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return {
            "file": os.path.basename(file_path),
            "organism": organism,
            "method": method,
            "resolution": resolution,
            "r_value":  round(r_value * 100, 2),
            "chains": "N/A",
            "classification": classification,
            "uniprot_ids": []
        }

def extract_uniprot_ids_from_pdb(file_path):
    uniprot_mapping = {}
    try:
        with open(file_path) as f:
            lines = f.readlines()

        for line in lines:
            if line.startswith("DBREF"):
                parts = line.split()
                if len(parts) >= 8:
                    chain = parts[2]
                    uniprot_id = parts[6]
                    if chain not in uniprot_mapping:
                        uniprot_mapping[chain] = set()
                    uniprot_mapping[chain].add(uniprot_id)
    except Exception as e:
        print(f"Error extracting UniProt IDs from {file_path}: {e}")
    return uniprot_mapping


def format_uniprot_ids(uniprot_mapping):
    uniprot_list = []
    for chain, ids in uniprot_mapping.items():
        for uniprot_id in ids:
            uniprot_list.append(f"{chain}-{uniprot_id}")
    return ",".join(uniprot_list)


def create_database(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE IF NOT EXISTS sequences
                  (file TEXT, chain TEXT, sequence TEXT)''')
    
    cursor.execute('''CREATE TABLE IF NOT EXISTS residues_map
                  (file TEXT, chain TEXT, idx_in_chain INTEGER, real_resnum TEXT)''')
    
    conn.commit()
    return conn


def store_sequence(conn, file, chain, sequence):
    cursor = conn.cursor()
    strip_file = file.strip()
    file_name = normalize_pdb_name(strip_file)
    cursor.execute("INSERT INTO sequences (file, chain, sequence) VALUES (?, ?, ?)",
                   (file_name, chain, sequence))
    conn.commit()


def update_database(files, symbols_db_conn, amino_acid_seq_conn):
    cursor = symbols_db_conn.cursor()
    cursor.execute("SELECT file FROM sequences")
    processed_files = set(row[0] for row in cursor.fetchall())
    
    resmap_cursor = symbols_db_conn.cursor()
    for pfile in processed_files:
        resmap_cursor.execute("SELECT chain, real_resnum FROM residues_map WHERE file=?", (pfile,))
        rows = resmap_cursor.fetchall()
        
        if rows:
            chain_map = defaultdict(list)
            for chain_id, rnum in rows:
                chain_map[chain_id].append(rnum)
            residues_map[pfile] = dict(chain_map)
        else:
            chain_map = {}
            for chain_id, rnum_str in rows:
                chain_map[chain_id] = [int(x) for x in rnum_str.split(",")]
            residues_map[pfile] = chain_map
            
    ########################################################################
    normalized_files = [normalize_pdb_name(f) for f in files]
    
    new_files = [files[i] for i, nf in enumerate(normalized_files) if nf not in processed_files]
    
    # We expect process_file to return local residue maps as well
    with ProcessPoolExecutor() as executor:
        results = list(filter(None, executor.map(process_file, new_files)))
        print(results)
    # Merge data from each subprocess back into the main DB + global residues_map
    for file, letter_sequences, amino_acid_sequences, local_map in results:
        print(file)
        for chain, letter_sequence in letter_sequences.items():
            store_sequence(symbols_db_conn, file, chain, letter_sequence)
        
        for chain, amino_acid_sequence in amino_acid_sequences.items():
            store_sequence(amino_acid_seq_conn, file, chain, amino_acid_sequence)

        # Merge the local residue-number map into the global residues_map
        residues_map[file] = local_map

        ins_cursor = symbols_db_conn.cursor()
        
        for chain_id, num_list in local_map.items():
            
            strip_file = file.strip()
            file_name = normalize_pdb_name(strip_file)
            real_nums = ",".join(str(n) for n in num_list)
            
            ins_cursor.execute("INSERT INTO residues_map (file, chain, real_resnum) VALUES (?, ?, ?)",
                                (file_name, chain_id, str(real_nums)))
        
        
        try:
            if ".cif" in file:
                headers = get_pdb_header_info_cif(file)
                
            else:
                headers = get_pdb_header_info(file)
            
            ins_cursor.execute("""
                CREATE TABLE IF NOT EXISTS pdb_header (
                    file TEXT PRIMARY KEY,
                    classification TEXT,
                    organism TEXT,
                    uniprot TEXT,
                    method TEXT,
                    resolution REAL,
                    r_value REAL,
                    chains TEXT
                )
            """)
            
            ins_cursor.execute("""
                INSERT INTO pdb_header (file, classification, organism, uniprot, method, resolution, r_value, chains)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                file_name, headers['classification'], headers['organism'], headers["uniprot_ids"],
                headers['method'], headers['resolution'], headers['r_value'], headers['chains']
            ))
        except Exception as e:
            error = e

        symbols_db_conn.commit()


def process_file(file):
    """
    Returns a 4-tuple:
      (file, letter_sequences, all_sequences, local_res_map)
    where local_res_map[chain] = [list_of_real_residue_numbers].
    """
    print(file)
    structure = parse_structure(file)
    if structure is None:
        return None
    # Calculate amino acid sequence (one-letter) for each chain
    all_sequences = extract_one_letter_sequence(structure)

    # Extract CA atoms, compute distances -> distance-based letter sequence
    ca_atoms = extract_ca_atoms(structure)
    distances = calculate_distances(ca_atoms)
    letter_sequences = {
        chain_id: "".join(distance_representation(dist) for dist, _ in distances[chain_id])
        for chain_id in ca_atoms
    }
    # print(letter_sequences)
    # Build local map of real residue numbers in the same order as letter_sequences
    local_res_map = {}
    for chain_id, residue_list in ca_atoms.items():
        chain_nums = []
        for ca_atom, _ in residue_list:
            if ca_atom != '?':
                parent_res = ca_atom.get_parent()
                _, real_resnum, _ = parent_res.id
                if real_resnum in chain_nums:
                    continue
                else:
                    chain_nums.append(real_resnum)
            else:
                chain_nums.append('?')
        local_res_map[chain_id] = chain_nums
    # print(local_res_map)
    return (file, letter_sequences, all_sequences, local_res_map)


def main(input_folder, output_dir):
    start_time = time.time()
    folder_name = os.path.basename(input_folder)
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    symbols_db_path = os.path.join(output_dir, f'{folder_name}_symbols.db')
    amino_acid_seq_db_path = os.path.join(output_dir, f'{folder_name}_amino_acid_seq.db')

    # Create databases
    symbols_db = create_database(symbols_db_path)
    amino_acid_seq_db = create_database(amino_acid_seq_db_path)
    fasta_seq_db = fasta_seq_main(input_folder, amino_acid_seq_db_path)
    for filename in os.listdir(input_folder):
        if filename.endswith(".gz"):
            gz_path = os.path.join(input_folder, filename)
            out_path = os.path.join(input_folder, filename[:-3])
            
            if os.path.exists(out_path):
                continue

            with gzip.open(gz_path, 'rb') as f_in:
                with open(out_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

    
    files = []
    folder_files = [os.path.join(input_folder, file) for file in os.listdir(input_folder) if file.endswith(('.pdb', '.cif', '.ent'))]
    files.extend(folder_files)
    
    if not files:
        print("❌ No valid structure files found. Please provide files in correct format (.pdb, .ent, .cif)")
        return
    
    update_database(files, symbols_db, amino_acid_seq_db)

    end_time = time.time()
    print(f"Execution time: {end_time - start_time:.2f} seconds")
    print("✅ Script execution completed successfully!\n")
    print(f"✅ Letter Sequence Database has been saved to: {symbols_db_path}")
    print(f"✅ Amino Acid Sequence Database has been saved to: {amino_acid_seq_db_path}")


if __name__ == "__main__":
       
    # 📌 Paste your input folder path containing ('.pdb', '.cif', '.ent') files
    input_folder = r"path/to/input"
    
    # 📌 Paste your output folder path where database files will be saved
    output_db_dir = r"path/to/out"
    
    print("🚀 Script started...")
    # for folder in input_folder:
    if not os.path.isdir(input_folder):
        print(f"❌ Input folder does not exist: {input_folder}")
        sys.exit(1)
            
    main(input_folder, output_db_dir)