from pathlib import Path
import pandas as pd
from typing import Dict, Set, List, Tuple
import logging
from dataclasses import dataclass

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class PropertyChange:
    """Container for tracking changes in a chemical property"""
    old_value: str
    new_value: str
    property_name: str

@dataclass
class SynonymChanges:
    """Container for tracking changes in synonyms"""
    added: Set[str]
    removed: Set[str]

class ChemicalDiffAnalyzer:
    def __init__(self, old_file: Path, new_file: Path):
        """Initialize with paths to old and new TSV files"""
        # Column names for the TSV files
        self.columns = [
            'cid', 'CAS', 'formula', 'molecular_weight', 
            'smiles', 'inchi', 'inchikey', 'iupac_name', 'common_name'
        ]
        
        # Properties to compare (excluding CAS and synonyms)
        self.compare_properties = [
            'cid', 'formula', 'molecular_weight', 
            'smiles', 'inchi', 'inchikey', 'iupac_name', 'common_name',
        ]
        
        # Load data
        self.old_data = self._load_tsv(old_file)
        self.new_data = self._load_tsv(new_file)
        
        # Get sets of CAS numbers
        self.old_cas_numbers = set(self.old_data['CAS'])
        self.new_cas_numbers = set(self.new_data['CAS'])

    def _load_tsv(self, file_path: Path) -> pd.DataFrame:
        """Load TSV file with variable number of synonym columns"""
        # Read all lines and process manually
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Split lines and create list of dictionaries
        data = []
        for line in lines:
            parts = line.strip().split('\t')
            row = dict(zip(self.columns, parts[:len(self.columns)]))
            row['synonyms'] = parts[len(self.columns):] if len(parts) > len(self.columns) else []
            data.append(row)
        
        return pd.DataFrame(data)

    def find_removed_compounds(self) -> Set[str]:
        """Find CAS numbers that were in old data but not in new"""
        return self.old_cas_numbers - self.new_cas_numbers

    def find_added_compounds(self) -> Set[str]:
        """Find CAS numbers that are in new data but not in old"""
        return self.new_cas_numbers - self.old_cas_numbers

    def find_synonym_changes(self) -> Dict[str, SynonymChanges]:
        common_cas = self.old_cas_numbers & self.new_cas_numbers
        changes: Dict[str, SynonymChanges] = {}

        for cas in common_cas:
            old_row = self.old_data[self.old_data['CAS'] == cas].iloc[0]
            new_row = self.new_data[self.new_data['CAS'] == cas].iloc[0]
            
            # Filter out IUPAC and common names from synonyms
            # remove only those from the new row
            old_synonyms = set(old_row['synonyms']) - {new_row['iupac_name'], new_row['common_name']}
            new_synonyms = set(new_row['synonyms']) - {new_row['iupac_name'], new_row['common_name']}

            added_synonyms = new_synonyms - old_synonyms
            removed_synonyms = old_synonyms - new_synonyms
            
            if added_synonyms or removed_synonyms:
                changes[cas] = SynonymChanges(added_synonyms, removed_synonyms)
        
        return changes

    def find_property_changes(self) -> Dict[str, List[PropertyChange]]:
        """Find changes in properties for compounds present in both datasets"""
        common_cas = self.old_cas_numbers & self.new_cas_numbers
        changes: Dict[str, List[PropertyChange]] = {}

        for cas in common_cas:
            old_row = self.old_data[self.old_data['CAS'] == cas].iloc[0]
            new_row = self.new_data[self.new_data['CAS'] == cas].iloc[0]
            
            compound_changes = []
            
            for prop in self.compare_properties:
                old_val = str(old_row[prop])
                new_val = str(new_row[prop])
                
                # Handle molecular weight comparison with tolerance
                if prop == 'molecular_weight':
                    try:
                        old_float = float(old_val)
                        new_float = float(new_val)
                        if abs(old_float - new_float) > 0.001:  # Tolerance of 0.001
                            compound_changes.append(PropertyChange(old_val, new_val, prop))
                    except ValueError:
                        if old_val != new_val:
                            compound_changes.append(PropertyChange(old_val, new_val, prop))
                else:
                    if old_val != new_val:
                        compound_changes.append(PropertyChange(old_val, new_val, prop))
            
            if compound_changes:
                changes[cas] = compound_changes
        
        return changes

    def generate_report(self, output: Path):
        """Generate a comprehensive report of all changes"""
        removed = self.find_removed_compounds()
        added = self.find_added_compounds()
        changes = self.find_property_changes()
        synonym_changes = self.find_synonym_changes()
        if isinstance(output, (str, Path)):
            f = open(output, 'w', encoding='utf-8')
            should_close = True
        else:
            f = output
            should_close = False

        # Report header
        f.write("Chemical Metadata Difference Report\n")
        f.write("=================================\n\n")
        
        # Summary
        f.write("Summary:\n")
        f.write(f"- Compounds removed: {len(removed)}\n")
        f.write(f"- Compounds added: {len(added)}\n")
        f.write(f"- Compounds with property changes: {len(changes)}\n")
        f.write(f"- Compounds with synonym changes: {len(synonym_changes)}\n\n")
        
        # Removed compounds
        if removed:
            f.write("\nRemoved Compounds:\n")
            f.write("------------------\n")
            for cas in sorted(removed):
                old_data = self.old_data[self.old_data['CAS'] == cas].iloc[0]
                f.write(f"CAS: {cas}\n")
                f.write(f"  IUPAC Name: {old_data['iupac_name']}\n")
                f.write(f"  Common Name: {old_data['common_name']}\n")
                f.write(f"  Formula: {old_data['formula']}\n\n")
                f.write(f"  CID: {old_data['cid']}\n")
        
        # Added compounds
        if added:
            f.write("\nAdded Compounds:\n")
            f.write("----------------\n")
            for cas in sorted(added):
                new_data = self.new_data[self.new_data['CAS'] == cas].iloc[0]
                f.write(f"CAS: {cas}\n")
                f.write(f"  IUPAC Name: {new_data['iupac_name']}\n")
                f.write(f"  Common Name: {new_data['common_name']}\n")
                f.write(f"  Formula: {new_data['formula']}\n")
                f.write(f"  Molecular Weight: {new_data['molecular_weight']}\n")
                f.write(f"  SMILES: {new_data['smiles']}\n")
                f.write(f"  InChI: {new_data['inchi']}\n")
                f.write(f"  InChIKey: {new_data['inchikey']}\n")
                f.write(f"  CID: {new_data['cid']}\n")
                if 'synonyms' in new_data and new_data['synonyms']:
                    f.write("  Synonyms:\n")
                    for synonym in new_data['synonyms']:
                        f.write(f"    - {synonym}\n")
                f.write("\n")        
        # Property changes
        if changes:
            f.write("\nProperty Changes:\n")
            f.write("----------------\n")
            for cas in sorted(changes.keys()):
                f.write(f"CAS: {cas}\n")
                for change in changes[cas]:
                    f.write(f"  {change.property_name}:\n")
                    f.write(f"    - Old: {change.old_value}\n")
                    f.write(f"    + New: {change.new_value}\n")
                f.write("\n")
        
        # Synonym changes
        if synonym_changes:
            f.write("\nSynonym Changes:\n")
            f.write("---------------\n")
            for cas in sorted(synonym_changes.keys()):
                changes = synonym_changes[cas]
                compound_name = self.new_data[self.new_data['CAS'] == cas].iloc[0]['common_name']
                f.write(f"CAS: {cas} ({compound_name})\n")
                
                # Sort changes alphabetically
                for removed in sorted(changes.removed):
                    f.write(f"  - {removed}\n")
                for added in sorted(changes.added):
                    f.write(f"  + {added}\n")
                f.write("\n")
        if should_close:
            f.close()

def main():
    import argparse
    from io import StringIO
    import sys
    
    parser = argparse.ArgumentParser(description='Compare two chemical metadata TSV files')
    parser.add_argument('old_file', type=Path, help='Path to the old TSV file')
    parser.add_argument('new_file', type=Path, help='Path to the new TSV file')
    parser.add_argument('output_file', type=Path, nargs='?', help='Optional path to write the report')
    
    args = parser.parse_args()
    
    # Verify files exist
    if not args.old_file.exists():
        logger.error(f"Old file not found: {args.old_file}")
        return
    if not args.new_file.exists():
        logger.error(f"New file not found: {args.new_file}")
        return
    
    analyzer = ChemicalDiffAnalyzer(args.old_file, args.new_file)
    
    if args.output_file:
        analyzer.generate_report(args.output_file)
        logger.info(f"Report generated: {args.output_file}")
    else:
        output = StringIO()
        analyzer.generate_report(output)
        print(output.getvalue())
        output.close()

if __name__ == '__main__':
    main()