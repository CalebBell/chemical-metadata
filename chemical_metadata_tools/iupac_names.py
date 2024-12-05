import os
import sys

def is_debugger_active():
    return hasattr(sys, 'gettrace') and sys.gettrace() is not None

def standardize_chemical_names(tsv_file):
    """
    Load chemical names from TSV and return standardized SMILES-to-name mapping.
    Selects shortest lowercase name, with fewest non-letter characters as tiebreaker.
    Names are first sorted alphabetically to ensure consistent selection.
    
    Args:
        tsv_file: Path to TSV file with Name and SMILES columns
        
    Returns:
        dict: Maps SMILES strings to their standardized names
    """
    # Initialize dict to collect names for each SMILES
    smiles_to_names = {}
    
    # Read TSV file directly
    with open(tsv_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            # Split on tab, strip whitespace
            try:
                name, smiles = line.strip().split('\t')
                # Initialize list if SMILES not seen before
                if smiles not in smiles_to_names:
                    smiles_to_names[smiles] = []
                smiles_to_names[smiles].append(name)
            except ValueError:
                continue  # Skip malformed lines
    
    def count_non_letters(name):
        return sum(not c.isalpha() for c in name)
    
    # Process each SMILES to get standardized name
    standardized = {}
    for smiles, names in smiles_to_names.items():
        # Convert to lowercase and sort alphabetically first
        names = sorted(name.lower() for name in names)
        # Then sort by length and non-letter count
        best_name = min(names, key=lambda x: (len(x), count_non_letters(x)))
        standardized[smiles] = best_name
    
    return standardized

file_path = os.path.join(os.path.dirname(__file__), '..', 'opsin_2.7.0_name_to_smiles_mapping_from_chemicals_metadata.tsv')
if is_debugger_active() or 1:
    iupac_standard_names = {}
else:
    iupac_standard_names = standardize_chemical_names(file_path)
