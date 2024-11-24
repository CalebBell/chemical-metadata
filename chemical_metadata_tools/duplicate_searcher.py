import json
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict

from os.path import join, dirname

import chemicals
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from chemicals.identifiers import ChemicalMetadata

FOLDER = chemicals.identifiers.folder

PUBCHEM_CATION_DB_NAME = 'Cation db.tsv'
PUBCHEM_ANION_DB_NAME = 'Anion db.tsv'
PUBCHEM_IONORGANIC_DB_NAME = 'Inorganic db.tsv'
PUBCHEM_EXAMPLE_DB_NAME = 'chemical identifiers example user db.tsv'


def load_chemical_file(filepath: str) -> List[ChemicalMetadata]:
    """Load chemicals from a TSV file"""
    chemicals = []
    with open(filepath, encoding='utf-8') as f:
        for line in f:
            try:
                values = line.rstrip('\n').split('\t')
                (pubchemid, CAS, formula, MW, smiles, InChI, InChI_key, 
                 iupac_name, common_name) = values[0:9]
                
                CAS = int(CAS.replace('-', ''))
                synonyms = values[7:]
                pubchemid = int(pubchemid)
                
                obj = ChemicalMetadata(
                    pubchemid, CAS, formula, float(MW), smiles,
                    InChI, InChI_key, iupac_name, common_name,
                    synonyms
                )
                chemicals.append(obj)
                
            except Exception as e:
                print(f"Error processing line in {filepath}:")
                print(f"  {str(e)}")
                print(f"  Line: {line.strip()}")
                continue
                
    return chemicals


@dataclass
class ChemicalEntry:
    """Simplified chemical entry for duplicate tracking"""
    pubchemid: int
    CAS: int
    formula: str
    MW: float
    smiles: str
    InChI: str
    InChI_key: str
    iupac_name: str
    common_name: str

class ChemicalDatabaseValidator:
    """Validates chemical database entries and manages duplicates"""
    
    def __init__(self, duplicate_dir: str = "duplicate_records"):
        self.duplicate_dir = Path(duplicate_dir)
        self.duplicate_dir.mkdir(exist_ok=True)
        
        # Initialize duplicate tracking files
        self.duplicate_files = {
            'pubchem': self.duplicate_dir / 'pubchem_duplicates.json',
            'cas': self.duplicate_dir / 'cas_duplicates.json',
            'smiles': self.duplicate_dir / 'smiles_duplicates.json',
            'inchi': self.duplicate_dir / 'inchi_duplicates.json',
            'inchi_key': self.duplicate_dir / 'inchi_key_duplicates.json',
            'formula': self.duplicate_dir / 'formula_duplicates.json'
        }
        
        # Track current values for each identifier type
        self.current_values = {
            'pubchem': {},
            'cas': {},
            'smiles': {},
            'inchi': {},
            'inchi_key': {},
            'formula': {}
        }
        
        # Load existing duplicates if any
        self.duplicates = self._load_duplicate_records()
    
    def _load_duplicate_records(self) -> Dict[str, Dict]:
        """Load existing duplicate records or create new ones"""
        records = {}
        for key, file_path in self.duplicate_files.items():
            if file_path.exists():
                with open(file_path, 'r') as f:
                    records[key] = json.load(f)
            else:
                records[key] = {}
        return records
    
    def _save_duplicate_record(self, identifier_type: str):
        """Save a specific duplicate record file"""
        with open(self.duplicate_files[identifier_type], 'w') as f:
            json.dump(self.duplicates[identifier_type], f, indent=2)

    def _create_entry(self, obj: 'ChemicalMetadata') -> ChemicalEntry:
        """Create a simplified entry from a ChemicalMetadata object"""
        return ChemicalEntry(
            pubchemid=obj.pubchemid,
            CAS=obj.CAS,
            formula=obj.formula,
            MW=obj.MW,
            smiles=obj.smiles,
            InChI=obj.InChI,
            InChI_key=obj.InChI_key,
            iupac_name=obj.iupac_name,
            common_name=obj.common_name
        )

    def check_and_record_duplicate(self, obj: 'ChemicalMetadata', source_file: str) -> List[str]:
        """Check for duplicates and record them if found"""
        conflicts = []
        entry = self._create_entry(obj)
        entry_dict = asdict(entry)
        entry_with_source = {
            'entry': entry_dict,
            'source_file': source_file
        }
        is_example_db = "example user db.tsv" in source_file
        # Check each identifier
        checks = {
            'pubchem': (obj.pubchemid, entry.pubchemid != -1),
            'cas': (obj.CAS, True),
            'smiles': (obj.smiles, obj.smiles != ''),
            'inchi': (obj.InChI, obj.InChI != ''),
            'inchi_key': (obj.InChI_key, obj.InChI_key != ''),
            'formula': (obj.formula, not is_example_db)
        }
        
        for key, (value, should_check) in checks.items():
            if should_check:
                str_value = str(value)
                
                if str_value in self.current_values[key]:
                    # Found a duplicate
                    if str_value not in self.duplicates[key]:
                        # First duplicate found - add both original and current entry
                        self.duplicates[key][str_value] = [
                            self.current_values[key][str_value],
                            entry_with_source
                        ]
                        self._save_duplicate_record(key)
                    else:
                        # Additional duplicate - append to existing list
                        self.duplicates[key][str_value].append(entry_with_source)
                        self._save_duplicate_record(key)
                    
                    conflicts.append(f"Found duplicate {key}: {value} in {source_file}")
                else:
                    # First occurrence - just store it
                    self.current_values[key][str_value] = entry_with_source
                    
        return conflicts

    def get_duplicate_summary(self) -> Dict[str, int]:
        """Get summary of duplicates found for each identifier type"""
        return {
            key: len(records)
            for key, records in self.duplicates.items()
        }

    def get_duplicates_for_identifier(self, identifier_type: str, value: str) -> List[Dict]:
        """Get all duplicate entries for a specific identifier value"""
        return self.duplicates.get(identifier_type, {}).get(str(value), [])

    def export_duplicate_report(self, output_file: str = "duplicate_report.txt"):
        """Export a human-readable report of all duplicates"""
        with open(output_file, 'w') as f:
            f.write("Chemical Database Duplicate Report\n")
            f.write("=================================\n\n")
            
            for identifier_type, records in self.duplicates.items():
                if records:  # Only show identifier types that have duplicates
                    f.write(f"\n{identifier_type.upper()} Duplicates:\n")
                    f.write("-" * (len(identifier_type) + 11) + "\n")
                    
                    for value, entries in records.items():
                        f.write(f"\nValue: {value}\n")
                        for i, entry in enumerate(entries, 1):
                            f.write(f"  Entry {i} from {entry['source_file']}:\n")
                            for k, v in entry['entry'].items():
                                f.write(f"    {k}: {v}\n")

                        
def main():
    # Setup output directory with timestamp
    output_dir = Path(__file__).parent.parent
    
    # Initialize validator
    validator = ChemicalDatabaseValidator(str(output_dir / "duplicates"))
    
    # Define database files to check
    db_files = [
        PUBCHEM_CATION_DB_NAME,
        PUBCHEM_ANION_DB_NAME,
        PUBCHEM_IONORGANIC_DB_NAME,
        PUBCHEM_EXAMPLE_DB_NAME
    ]
    
    total_chemicals = 0
    total_conflicts = 0
    
    # Process each database
    for db_name in db_files:
        db_path = os.path.join(FOLDER, db_name)
        if not os.path.exists(db_path):
            print(f"Warning: Database file not found: {db_path}")
            continue
            
        print(f"\nProcessing {db_name}...")
        chemicals = load_chemical_file(db_path)
        total_chemicals += len(chemicals)
        
        file_conflicts = 0
        for chemical in chemicals:
            conflicts = validator.check_and_record_duplicate(chemical, db_name)
            if conflicts:
                file_conflicts += len(conflicts)
                
        total_conflicts += file_conflicts
        print(f"Found {file_conflicts} conflicts in {len(chemicals)} entries")
    
    # Generate summary report
    # summary_path = output_dir / "analysis_summary.txt"
    # with open(summary_path, 'w') as f:
    #     f.write("Chemical Database Duplicate Analysis Summary\n")
    #     f.write("=========================================\n\n")
    #     f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
    #     f.write("Files Processed:\n")
    #     for db_name in db_files:
    #         f.write(f"  - {db_name}\n")
        
    #     f.write(f"\nTotal Chemicals Processed: {total_chemicals}\n")
    #     f.write(f"Total Conflicts Found: {total_conflicts}\n\n")
        
    #     f.write("Duplicate Summary by Identifier:\n")
    #     for id_type, count in validator.get_duplicate_summary().items():
    #         f.write(f"  {id_type}: {count} duplicates\n")
    
    # Generate detailed report
    # validator.export_duplicate_report(str(output_dir / "detailed_duplicates.txt"))
    
    print(f"\nAnalysis complete!")
    print(f"Total chemicals processed: {total_chemicals}")
    print(f"Total conflicts found: {total_conflicts}")
    print(f"\nReports generated in: {output_dir}")
    print(f"  - analysis_summary.txt")
    print(f"  - detailed_duplicates.txt")
    print(f"  - duplicates/ (JSON files)")

if __name__ == "__main__":
    main()
