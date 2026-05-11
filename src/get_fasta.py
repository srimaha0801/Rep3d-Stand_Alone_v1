import os
import requests
import sqlite3


# ---------------------------
# Extract PDB ID from filename
# ---------------------------
def extract_pdb_id(filename):

    name = os.path.basename(filename).lower()

    if name.startswith("pdb"):
        name = name[3:]

    return name.split(".")[0].lower()


# ---------------------------
# Fetch FASTA from RCSB
# ---------------------------
def get_pdb_fasta(pdb_id):

    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"

    try:
        response = requests.get(url)

        if response.status_code == 200:
            return response.text

        return None

    except:
        return None


# ---------------------------
# Parse FASTA → REAL chain IDs
# ---------------------------
def parse_fasta_real_chains(fasta_text):

    sequences = []

    if not fasta_text:
        return sequences

    entries = fasta_text.strip().split(">")[1:]

    for entry in entries:

        lines = entry.split("\n")

        header = lines[0]

        sequence = "".join(lines[1:])

        parts = header.split("|")

        if len(parts) < 2:
            continue

        chain_part = parts[1]

        # Example:
        # "Chains A, B"
        # "Chain J"

        if "Chains" in chain_part:
            chains = chain_part.replace("Chains", "").strip()

        elif "Chain" in chain_part:
            chains = chain_part.replace("Chain", "").strip()

        else:
            chains = chain_part.strip()

        # chain_ids = [c.strip().lower() for c in chains.split(",")]

        # for chain_id in chain_ids:
        #     sequences.append((chain_id, sequence))'
        chain_ids = []

        for c in chains.split(","):

            c = c.strip()

            # Handle auth chain IDs
            if "[auth" in c:

                auth_part = c.split("[auth")[-1]
                real_chain = auth_part.replace("]", "").strip()

                chain_ids.append(real_chain)

            else:
                chain_ids.append(c)

        for chain_id in chain_ids:
            sequences.append((chain_id, sequence))

    return sequences


# ---------------------------
# Create DB
# ---------------------------
def create_database(db_name):

    conn = sqlite3.connect(db_name)

    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS fasta_sequences (
            pdb_id TEXT,
            chain TEXT,
            sequence TEXT
        )
    """)

    conn.commit()

    return conn


# ---------------------------
# Store Data
# ---------------------------
def store_sequences(conn, pdb_id, sequences):

    cursor = conn.cursor()

    for chain, seq in sequences:

        cursor.execute(
            """
            INSERT INTO fasta_sequences
            (pdb_id, chain, sequence)
            VALUES (?, ?, ?)
            """,
            (
                pdb_id.lower(),
                chain,
                seq
            )
        )

    conn.commit()


# ---------------------------
# MAIN
# ---------------------------
def fasta_seq_main(input_folder, db_name):

    if not os.path.exists(input_folder):
        return False

    conn = create_database(db_name)

    files = [
        f for f in os.listdir(input_folder)
        if f.endswith((".pdb", ".ent", ".cif"))
    ]

    if not files:
        conn.close()
        return False

    total_processed = 0

    for file in files:

        pdb_id = extract_pdb_id(file)

        fasta_text = get_pdb_fasta(pdb_id)

        if fasta_text:

            sequences = parse_fasta_real_chains(fasta_text)

            if sequences:

                store_sequences(conn, pdb_id, sequences)

                total_processed += 1

    conn.close()

    return True

# if __name__ == "__main__":

#     input_folder = r"path/to/input_folder"

#     db_name = r"path/to/fasta_sequences.db"

#     result = fasta_seq_main(input_folder, db_name)

#     if result:
#         print("FASTA extraction completed successfully.")
#     else:
#         print("FASTA extraction failed.")