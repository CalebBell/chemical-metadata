import os
from chemicals.identifiers import (
    ChemicalMetadataDB,
    folder,
    PUBCHEM_LARGE_DB_NAME,
    PUBCHEM_SMALL_DB_NAME,
    PUBCHEM_CATION_DB_NAME,
    PUBCHEM_ANION_DB_NAME,
    PUBCHEM_IONORGANIC_DB_NAME,
    PUBCHEM_EXAMPLE_DB_NAME
)

def load_db_cas_numbers(filepath):
    """Load a database file and return a set of its CAS numbers"""
    db = ChemicalMetadataDB(elements=False, main_db=None, user_dbs=[filepath])
    return set(db.CAS_index.keys())

def read_db_lines(filepath):
    """Read all lines from a database file"""
    with open(filepath, 'r', encoding='utf-8') as f:
        return f.readlines()

def write_filtered_db(input_filepath, output_filepath, cas_to_remove):
    """Write a new database file excluding lines with specified CAS numbers"""
    lines = read_db_lines(input_filepath)
    filtered_lines = []
    
    for line in lines:
        values = line.rstrip('\n').split('\t')
        if len(values) >= 2:  # Ensure line has enough columns
            cas = int(values[1].replace('-', ''))  # CAS is second column
            if cas not in cas_to_remove:
                filtered_lines.append(line)
    
    with open(output_filepath, 'w', encoding='utf-8') as f:
        f.writelines(filtered_lines)

def main():
    # Define file paths using constants
    cation_db = os.path.join(folder, PUBCHEM_CATION_DB_NAME)
    anion_db = os.path.join(folder, PUBCHEM_ANION_DB_NAME)
    inorganic_db = os.path.join(folder, PUBCHEM_IONORGANIC_DB_NAME)
    example_db = os.path.join(folder, PUBCHEM_EXAMPLE_DB_NAME)
    small_db = os.path.join(folder, PUBCHEM_SMALL_DB_NAME)
    large_db = os.path.join(folder, PUBCHEM_LARGE_DB_NAME)
    
    # Collect CAS numbers from specialized databases
    cas_to_remove = set()
    specialized_dbs = [cation_db, anion_db, inorganic_db, example_db]
    
    for db_path in specialized_dbs:
        if os.path.exists(db_path):
            db_cas = load_db_cas_numbers(db_path)
            cas_to_remove.update(db_cas)
    
    # Filter and overwrite the small and large databases
    write_filtered_db(small_db, small_db, cas_to_remove)
    write_filtered_db(large_db, large_db, cas_to_remove)

if __name__ == '__main__':
    main()