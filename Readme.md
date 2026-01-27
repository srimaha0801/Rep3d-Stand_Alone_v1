# Rep3D: An Algorithm to Identify Structurally Similar Repeats

## Authors
Gurleen Kaur, Madhumathi Sanjeevi, Srimaha Gandhi and Kanagaraj Sekar 

## Institution
Indian Institute of Science (IISc), Bangalore

---

## Overview
Rep3D is a standalone algorithm designed to identify **structurally similar repeats** in protein structures using 3D coordinate information and sequence-based symbolic representations.

This repository contains scripts to:
- Generate required SQLite databases
- Identify and extract matching structural repeats

---

## Requirements
- Python **3.8 or higher**
- SQLite3

---

## Setup

### 1. Create a Virtual Environment

```

python -m venv env

To activate the environment:
      • Linux/macOS:   source env/bin/activate
      • Windows (PowerShell):   .\env\Scripts\Activate.ps1
      • Windows (CMD):          .\env\Scripts\activate.bat

```

### 2. Install required dependencies:
```
   pip install wheel requests numpy pandas bio biopython openpyxl
```
### 3. Navigate to the source directory
``` 
    Rep3D-stand_alone_v1/src
```
### 4. Run step1_generate_databases.py 
```
Before running, open: step1_generate_databases.py

At the bottom of the file, inside the "if __name__ == '__main__':" block,

paste your input folder path (containing pdb/ent/cif files) and output folder path.

Example:

    if __name__ == "__main__":
        # 📌 Paste your input folder path containing ('.pdb', '.cif', '.ent') files
        input_folder = r"/path/to/input/Bacillus_Thermoproteolyticus"
        # 📌 Paste your output folder path where database files will be saved
        output_dir = r"/path/to/databases/Bacillus_Thermoproteolyticus"
        main(input_folders, output_dir)

Run Step 1 using:
    python step1_generate_databases.py

```
### 5. Run step2_generate_results.py

```
Before running, open:
    step2_generate_results.py

Inside the "if __name__ == '__main__':" block, provide:

    • XYZ coordinates of the C-alpha (Cα) atoms
    • Path to letter sequence DB
    • Path to amino acid sequence DB
    • Output directory path

Notes for XYZ coordinates:
    - Provide the coordinates for **Cα atoms only**.
    - Coordinates must be in the following format:
          (x1 y1 z1),(x2 y2 z2),(x3 y3 z3),...
    - You must provide a **minimum of 5 Cα coordinates** for a valid search.

    - Example:
          usr_coord = "(32.732  36.723   6.405),(33.094  34.136   9.151),(30.430  31.813   7.657),(32.172  31.710   4.261),(35.248  30.628   6.169),(33.312  27.796   7.786)"
    - You can include more than 5 coordinates as needed for longer fragments.

```
### Optional Filtering

```
You can optionally filter your structure search results using the variable:
    filter_choice
    - `filter_choice = True`  → apply filters

    - `filter_choice = False` → run without filters

Filter Rules:

    - Filters are optional — leave a field as `""` or `[]` if not needed.

    - **Organism** must be written in **lowercase**.

          Example: "homo sapiens", "bacillus thermoproteolyticus"

    - **Experimental Method** must be written in **lowercase**.

      Supported method options:
        • "x-ray diffraction"
        • "electron microscopy"
        • "solution nmr"
        • "solid-state nmr"
        • "hybrid"
        • "others"

    - **Resolution** and **R-factor** formats supported:
        • "greater_than X"
        • "less_than X"
        • "equal X"
        • "between X,Y"

      Example:
          "greater_than 2.0"
          "between 1.5,2.5"
    - **PDB ID filtering**:
      Provide a list of PDB IDs in lowercase like:
          ["9rm2", "9otp", "9ofx", "9r96"]

Example filter block:

filters = {
    "organism": "homo sapiens",
    "method": "x-ray diffraction",
    "resolution": "greater_than 0.1",
    "rfac": "",
    "pdb_ids": []
}

Example:

    if __name__ == "__main__":

        # 📌 Paste your XYZ coordinates here
        usr_coord = "(32.732  36.723   6.405),(33.094  34.136   9.151),(30.430  31.813   7.657),(32.172  31.710   4.261),(35.248  30.628   6.169),(33.312  27.796   7.786)"

        # 📌 Paste your letter sequence database (symbols.db)
        symbols_db_path = r"/path/to/databases/Bacillus_Thermoproteolyticus/Bacillus_Thermoproteolyticus_symbols.db"

        # 📌 Paste your amino acid sequence database (amino_acid_seq.db)
        amino_acid_seq_db_path = r"/path/to/databases/Bacillus_Thermoproteolyticus/Bacillus_Thermoproteolyticus_amino_acid_seq.db"

        # 📌 Output directory where results will be stored
        results_path = r"/path/to/output/Bacillus_Thermoproteolyticus/"

        main(usr_coord, results_path, symbols_db_path, amino_acid_seq_db_path)

Run Step 2 using:
    python step2_generate_results.py

```

### Output

After execution, a TSV result file will be generated in the specified results folder.