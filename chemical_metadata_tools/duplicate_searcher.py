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
from collections import OrderedDict
from chemicals.identifiers import (ChemicalMetadata, PUBCHEM_LARGE_DB_NAME, PUBCHEM_SMALL_DB_NAME, PUBCHEM_EXAMPLE_DB_NAME, PUBCHEM_CATION_DB_NAME, PUBCHEM_ANION_DB_NAME, PUBCHEM_IONORGANIC_DB_NAME)

FOLDER = chemicals.identifiers.folder



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

def load_all_preferences() -> Dict:
    """Load and combine all preference files from various locations"""
    preferences = {
        "preferred_cas": set(),
        "unpreferred_cas": set()
    }
    
    # Locations to check for preference files
    preference_locations = [
        Path(FOLDER) / "inorganic_preferences.json",
        Path(FOLDER) / "anion_preferences.json",
        Path(FOLDER) / "cation_preferences.json",
        # Path(FOLDER) / "organic_preferences.json",
        # Add other potential locations
    ]
    
    for pref_file in preference_locations:
        if pref_file.exists():
            try:
                with open(pref_file) as f:
                    data = json.load(f)
                    # Convert CAS numbers to standard format (no hyphens) for comparison
                    preferences["preferred_cas"].update(
                        int(cas.replace('-', '')) for cas in data.get("preferred_cas", [])
                    )
                    preferences["unpreferred_cas"].update(
                        int(cas.replace('-', '')) for cas in data.get("unpreferred_cas", [])
                    )
            except Exception as e:
                print(f"Error loading preferences from {pref_file}: {e}")
    
    return preferences

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
            'formula': self.duplicate_dir / 'formula_duplicates.json',
            'common_name': self.duplicate_dir / 'common_name_duplicates.json',
        }
        
        # Track current values for each identifier type
        self.current_values = {
            'pubchem': {},
            'cas': {},
            'smiles': {},
            'inchi': {},
            'inchi_key': {},
            'formula': {},
            'common_name': {},
        }
        
        # Load existing duplicates if any
        self.duplicates = self._load_duplicate_records()
    
    def _load_duplicate_records(self) -> Dict[str, Dict]:
        """Load empty records dicts"""
        records = {}
        for key in self.duplicate_files.keys():
            records[key] = {}
        return records
    
    def _save_duplicate_record(self, identifier_type: str):
        """Save a specific duplicate record file"""
        preferences = load_all_preferences()
        with open(self.duplicate_files[identifier_type], 'w') as f:
            # Sort dictionary by keys before saving
            to_save = sorted(self.duplicates[identifier_type].items())
            sorted_data = OrderedDict(to_save)

            # Create new dict excluding preference-handled duplicates
            filtered_data = OrderedDict()
            for entry, entry_data in sorted_data.items():
                # Count how many preferred/unpreferred CAS numbers we have
                preferred_count = sum(1 for item in entry_data 
                                if item['entry']['CAS'] in preferences['preferred_cas'])
                unpreferred_count = sum(1 for item in entry_data 
                                    if item['entry']['CAS'] in preferences['unpreferred_cas'])
                # Keep entry if it's not explained by preferences
                # (should have exactly one preferred and rest unpreferred)
                if not (preferred_count == 1 and 
                    unpreferred_count == len(entry_data) - 1):
                    filtered_data[entry] = entry_data

            # Ensure consistent JSON formatting
            json.dump(filtered_data, f, indent=2, sort_keys=True)

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
        is_example_db = PUBCHEM_EXAMPLE_DB_NAME in source_file or PUBCHEM_SMALL_DB_NAME in source_file or PUBCHEM_LARGE_DB_NAME in source_file
        # Check each identifier
        # Note that formulas are only checked for inorganic/ion db as organics match all the time for obvious reasons
        checks = {
            'pubchem': (obj.pubchemid, entry.pubchemid != -1),
            'cas': (obj.CAS, True),
            'smiles': (obj.smiles, obj.smiles != ''),
            'inchi': (obj.InChI, obj.InChI != ''),
            'inchi_key': (obj.InChI_key, obj.InChI_key != ''),
            'formula': (obj.formula, not is_example_db),
            'common_name': (obj.common_name, obj.common_name != ''),
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
                    else:
                        # Additional duplicate - append to existing list
                        self.duplicates[key][str_value].append(entry_with_source)
                    
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
        PUBCHEM_EXAMPLE_DB_NAME,
        PUBCHEM_SMALL_DB_NAME,
        PUBCHEM_LARGE_DB_NAME # 400-800 conflicts to sort out, also failing to remove some redundant data automatically
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
        # print(f"Found {file_conflicts} conflicts in {len(chemicals)} entries")
    for key in validator.duplicate_files.keys():
        validator._save_duplicate_record(key)

    
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
