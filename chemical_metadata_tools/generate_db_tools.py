 
from __future__ import annotations
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
import logging
from chemicals.elements import molecular_weight, nested_formula_parser
from dataclasses import dataclass

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from pubchempy import get_compounds, Compound, get_json
import json
from collections import Counter
from thermo import serialize_formula
from fcache.cache import FileCache
from chemical_metadata_tools.iupac_names import iupac_standard_names
from chemical_metadata_tools.parse_CAS_data import lower_case_first_letter_name
from chemical_metadata_tools.synonym_utils import fix_synonym_case

from collections import defaultdict

import re
from typing import Optional

def is_inchikey(s: str) -> bool:
    """
    Validates if a string is a properly formatted InChIKey.
    
    InChIKey format: XXXXXXXXXXXXXX-YYYYYYYYFV-P where:
    - First block: 14 uppercase characters (connectivity hash)
    - Second block: 8 characters (remaining layers hash) + 
                   1 character ('S' for standard or 'N' for nonstandard) +
                   1 character (version, currently 'A' for version 1)
    - Third block: 1 character for protonation state 
                  ('N' for none, 'O'/'P'/etc for added, 'M'/'L'/etc for removed)
    
    Args:
        s: String to validate
        
    Returns:
        bool: True if string is valid InChIKey format, False otherwise
    """
    if not isinstance(s, str):
        return False
    s = s.upper()  # also looking to remove lower case inchi keys if any noncompliant ones exist    
    # Check total length (14 + 10 + 1 + 2 hyphens = 27)
    if len(s) != 27:
        return False
        
    # Split into blocks
    blocks = s.split('-')
    if len(blocks) != 3:
        return False
        
    block1, block2, block3 = blocks
    
    # Validate first block (14 uppercase letters)
    if len(block1) != 14 or not block1.isalpha() or not block1.isupper():
        return False
        
    # Validate second block (8 chars + S/N + A)
    if len(block2) != 10:
        return False
    
    # First 8 chars should be uppercase letters
    if not block2[:8].isalpha() or not block2[:8].isupper():
        return False
        
    # 9th char should be S or N
    if block2[8] not in {'S', 'N'}:
        return False
        
    # 10th char should be version identifier (currently only A)
    if block2[9] != 'A':
        return False
        
    # Validate third block (single char for protonation state)
    if len(block3) != 1 or not block3.isalpha() or not block3.isupper():
        return False
        
    return True
def remove_inchikeys(strings: list[str]) -> list[str]:
    """
    Returns a new list with all InChIKey strings removed.
    
    Args:
        strings: List of strings to filter
        
    Returns:
        List of strings with InChIKeys removed
    """
    return [s for s in strings if not is_inchikey(s)]
def remove_inchikeys(strings: list[str]) -> list[str]:
    """
    Returns a new list excluding strings that contain InChIKeys as substrings.
    Uses regex for fast initial filtering before strict validation.
    
    Args:
        strings: List of strings to filter
        
    Returns:
        List of strings without any InChIKeys in them
    """
    # Pattern matches potential InChIKey format: 14 chars - 10 chars - 1 char
    # Using \S (non-whitespace) instead of [A-Z] for initial fast scan
    potential_pattern = re.compile(r'\S{14}-\S{10}-\S{1}')
    
    def contains_inchikey(s: str) -> bool:
        # First quick check - any potential matches?
        potential_matches = potential_pattern.findall(s)
        if not potential_matches:
            return False
            
        # Verify any potential matches with strict validator
        return any(is_inchikey(match) for match in potential_matches)
    
    return [s for s in strings if not contains_inchikey(s)]
def deduplicate_names(names):
    """
    Deduplicate strings that differ only in capitalization, preferring versions with more lowercase letters.
    
    Args:
        names: List of strings to deduplicate
        
    Returns:
        List of strings with case-sensitive duplicates removed, preferring lowercase versions
    """
    # Group names by their lowercase equivalent
    case_groups = defaultdict(list)
    for name in names:
        case_groups[name.lower()].append(name)
    
    # Choose the preferred version (with most lowercase letters) from each group
    preferred_names = {}
    for lower_name, variations in case_groups.items():
        # Select the version with the most lowercase letters, breaking ties lexicographically
        preferred_names[lower_name] = max(
            variations, key=lambda var: (sum(1 for c in var if c.islower()), -ord(var[0]))
        )
    
    # Maintain original order, only taking the first occurrence of each deduplicated name
    result = []
    seen = set()
    for name in names:
        lower_name = name.lower()
        if lower_name not in seen:
            seen.add(lower_name)
            result.append(preferred_names[lower_name])
    
    return result

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class ChemicalData:
    """Container for chemical compound data"""
    cid: int
    iupac_name: str
    name: str
    names: List[str]
    molecular_weight: float
    smiles: str
    inchi_val: str
    inchikey: str
    formula: str

class ChemicalMetadataProcessor:
    def __init__(self, cache_dir: Path, workdir_path: Path, require_charge=None, require_not_charge=None, ignore_CASs=set()):
        self.cache = FileCache('myapp', flag='cs', serialize=True, 
                             app_cache_dir=str(cache_dir))
        
        # List of CAS numbers that are problematic
        self.ignored_CASs = ignore_CASs
        self.workdir_path = workdir_path
        
        # Load cached data
        self.syn_data = self._load_json(workdir_path/'Good synoynms by CAS.json')
        self.pdf_data = self._load_json(workdir_path/'Parsed CAS metadata.json')
        
        # Process names
        self.all_user_names = self._collect_user_names()
        self.all_names, self.dup_names = self._collect_all_names()

        self.require_charge=require_charge
        self.require_not_charge=require_not_charge

    def _load_json(self, filename: str) -> Dict:
        """Load and parse a JSON file"""
        with open(filename, 'r', encoding='utf-8') as f:
            return json.load(f)

    def _collect_user_names(self) -> Set[str]:
        """Collect all user-defined synonyms"""
        names = []
        for d in self.syn_data.values():
            if 'synonyms' in d:
                names.extend(d['synonyms'])
        return set(names)

    def _collect_all_names(self) -> Tuple[Set[str], Set[str]]:
        """Collect all names and find duplicates"""
        all_names = []
        for CAS, d in self.pdf_data.items():
            if 'Other Names' in d:
                all_names.extend(d['Other Names'])
            if 'Name' in d:
                all_names.append(d['Name'])
            if 'Deleted CAS' in d:
                all_names.extend(d['Deleted CAS'])
            if 'Alternate CAS' in d:
                all_names.extend(d['Alternate CAS'])
            all_names.append(CAS)
        
        # Find duplicates
        name_counts = Counter(all_names)
        dup_names = {item for item, count in name_counts.items() if count > 1}
        
        return set(all_names), dup_names

    def _get_mol_data(self, filepath: Path, CAS: str) -> Optional[ChemicalData]:
        """Get molecular data from file or fallback to stored data"""
        try:
            # Try to read from mol file first; maybe from scifinder/common chemistry or pubchem
            if "PUBCHEM_COMPOUND_CID" in open(str(filepath)).read():
                orig_mol = mol = next(Chem.SDMolSupplier(str(filepath)))
                cid = int(mol.GetProp('PUBCHEM_COMPOUND_CID'))
            else:
                orig_mol = mol = Chem.MolFromMolFile(str(filepath))
                cid = -1
            if mol is None:
                raise ValueError(f"Cannot read mol file {filepath}")
            
            # Convert to InChI and back for better ion handling. If this is removed, many ions get ruined
            # inchi_val = inchi.MolToInchi(mol, options='/FixedH')
            inchi_val = inchi.MolToInchi(mol)
            if not inchi_val:
                logger.debug(f"Could not serialize mol to inchi {CAS}")
                raise ValueError(f"Could not serialize mol to inchi {CAS}")
            mol = inchi.MolFromInchi(inchi_val)
            # If we were supposed to have a charge and we don't, inchi messed up and we can only take MW
            found_charge = Chem.GetFormalCharge(mol)
            if (self.require_charge is not None and found_charge != self.require_charge) or (self.require_not_charge is not None and self.require_not_charge == found_charge):
                # try to go back to the first path for the case of broken ions
                # mol = Chem.MolFromMolFile(str(filepath))
                # found_charge = Chem.GetFormalCharge(mol)
                # if (self.require_charge is not None and found_charge != self.require_charge) or (self.require_not_charge is not None and self.require_not_charge == found_charge):
                logger.debug(f"Rdkit did not handle mole charge correctly for {CAS}")
                raise ValueError(f"Rdkit did not handle mole charge correctly for {CAS}")
            
            smiles = Chem.MolToSmiles(mol, True)
            inchikey = (inchi.InchiToInchiKey(inchi_val) if inchi_val else '')
            mw = Descriptors.MolWt(mol)
            formula = CalcMolFormula(mol, True, True)

            if cid != -1:
                if orig_mol.GetProp('PUBCHEM_MOLECULAR_FORMULA') and  orig_mol.GetProp('PUBCHEM_MOLECULAR_FORMULA') != formula:
                    # rdkit likely added hs in a bad way
                    formula = orig_mol.GetProp('PUBCHEM_MOLECULAR_FORMULA')
                    mw = float(orig_mol.GetProp('PUBCHEM_MOLECULAR_WEIGHT')) # has to be PUBCHEM_MOLECULAR_WEIGHT, see 7784-13-6 cid 24564
                    inchi_val = orig_mol.GetProp('PUBCHEM_IUPAC_INCHI')
                    inchikey = orig_mol.GetProp('PUBCHEM_IUPAC_INCHIKEY')
                    smiles = orig_mol.GetProp('PUBCHEM_OPENEYE_CAN_SMILES')


            return ChemicalData(
                cid=cid,
                iupac_name='',
                name='',
                names=[''],
                molecular_weight=mw,
                smiles=smiles,
                inchi_val=inchi_val,
                inchikey=inchikey,
                formula=formula
            )
            
        except Exception as e:
            logger.debug(f"Mol file processing failed for {CAS}: {str(e)}")
            return self._get_fallback_data(CAS)

    def _get_fallback_data(self, CAS: str) -> Optional[ChemicalData]:
        """Get data from stored synonyms or PubChem when mol file fails"""
        if CAS not in self.syn_data:
            logger.error(f"No fallback data available for {CAS}")
            return None
            
        d = self.syn_data[CAS]
        if 'pubchem' in d:
            return self._get_pubchem_data(d['pubchem'])
        else:
            return ChemicalData(
                cid=-1,
                iupac_name='',
                name='',
                names=d.get('synonyms', ['']),
                molecular_weight=(float(d['MW']) if ('MW' in d and d['MW'] is not None ) else None),
                smiles=d.get('smiles', ''),
                inchi_val=d.get('inchi', ''),
                inchikey=d.get('inchikey', ''),
                formula=d.get('formula', '')
            )

    def _get_pubchem_data(self, compound_id: int) -> Optional[ChemicalData]:
        """Fetch and cache PubChem data"""
        cache_key = str(compound_id)
        
        if cache_key in self.cache:
            data = self.cache[cache_key]
            return ChemicalData(*data)
        
        try:
            pc = Compound.from_cid(compound_id)
            data = ChemicalData(
                cid=pc.cid,
                iupac_name=pc.iupac_name,
                name=pc.iupac_name,
                names=pc.synonyms,
                molecular_weight=float(pc.molecular_weight),
                smiles=pc.canonical_smiles,
                inchi_val=pc.inchi,
                inchikey=pc.inchikey,
                formula=pc.molecular_formula
            )
            self.cache[cache_key] = (
                data.cid, data.iupac_name, data.iupac_name, data.names, float(data.molecular_weight),
                data.smiles, data.inchi_val, data.inchikey, data.formula
            )
            return data
        except Exception as e:
            logger.error(f"PubChem lookup failed for {compound_id}: {str(e)}")
            return None

    def _get_pubchem_additional_data(self, mol_data: ChemicalData) -> Optional[Tuple[int, str, List[str]]]:
        """Get additional data from PubChem using InChIKey.
        IUPAC name is the only name offered by pubchempy."""
        try:
            if mol_data.cid > 0:
                cache_key = 'cid_lookup' + str(mol_data.cid)
                known_cid = True
            else:
                cache_key = str(mol_data.inchikey)
                known_cid = False
            if cache_key in self.cache:
                return self.cache[cache_key]
            
            if known_cid:
                compounds = get_compounds(mol_data.cid, 'cid')
            else:
                compounds = get_compounds(mol_data.inchikey, 'inchikey')
            if compounds:
                pc = compounds[0]
                data = (pc.cid, pc.iupac_name, pc.synonyms)
                self.cache[cache_key] = data
                return data
            else:
                self.cache[cache_key] = (-1, '', [''])
                return (-1, '', [''])
        except Exception as e:
            logger.debug(f"PubChem additional data lookup failed: {str(e)}")
            return (-1, '', [''])

    def _combine_with_scifinder(self, CAS: str, mol_data: ChemicalData, 
                              pubchem_data: Optional[Tuple[int, str, List[str]]]) -> ChemicalData:
        """Combine molecular data with SciFinder data"""

        syn_dict = self.syn_data.get(CAS, {})
        if 'smiles' in syn_dict or 'formula' in syn_dict:
            # We only hardcode smiles and/or formula when rdkit/inchi isn't working
            mol_data.inchi_val = ''
            mol_data.inchikey = ''
            mol_data.smiles = ''
        smiles = syn_dict.get('smiles', mol_data.smiles)

        if pubchem_data:
            cid, iupac_name, names = pubchem_data
        else:
            cid, iupac_name, names = -1, iupac_standard_names.get(smiles, ''), []

        formula=serialize_formula(syn_dict.get('formula', mol_data.formula))
        name = iupac_name

        if CAS in self.pdf_data:
            pdf_dict = self.pdf_data[CAS]
            name = pdf_dict['Name']
            syns = pdf_dict.get('Other Names', [])
            MW = pdf_dict.get('MW')
            
            if not iupac_name:
                iupac_name = name
            else:
                syns.insert(0, name)
            
            other_CAS = []
            if 'Deleted CAS' in pdf_dict:
                other_CAS.extend(pdf_dict['Deleted CAS'])
            if 'Alternate CAS' in pdf_dict:
                other_CAS.extend(pdf_dict['Alternate CAS'])
            
            # Filter out duplicate names
            syns = [i for i in syns if i not in self.dup_names]
            names = syns + [i for i in names if i not in self.all_names] + other_CAS

        # Process names through user database
        synonyms = []
        for a_name in names:
            if a_name in self.all_user_names:
                if ('synonyms' in syn_dict and a_name in syn_dict['synonyms']):
                    synonyms.append(a_name)
            else:
                synonyms.append(a_name)

        if 'name' in syn_dict:
            synonyms.append(name)
            name = syn_dict['name']

        # Add user synonyms if available
        if 'synonyms' in syn_dict:
            for n in syn_dict['synonyms']:
                if n not in synonyms:
                    synonyms.append(n)

        synonyms = [i for i in synonyms if i]
        synonyms = deduplicate_names(synonyms)
        synonyms = remove_inchikeys(synonyms)

        if name in synonyms:
            synonyms.remove(name)
        if iupac_name in synonyms:
            synonyms.remove(iupac_name)
        if CAS in synonyms:
            synonyms.remove(CAS)
        if 'MW' in syn_dict:
            MW = syn_dict['MW']
        else:
            if formula:
                MW = molecular_weight(nested_formula_parser(formula))
            else:
                MW = (mol_data.molecular_weight if mol_data.molecular_weight is not None else MW)
        inchi_val = syn_dict.get('inchi', mol_data.inchi_val)
        inchikey = syn_dict.get('inchikey', mol_data.inchikey)

        iupac_name = syn_dict.get('iupac_name', iupac_name)

        if inchi_val in synonyms:
            synonyms.remove(inchi_val)
        if inchikey in synonyms:
            synonyms.remove(inchikey)
        synonyms = [v for v in synonyms if not 'inchi' in v.lower()]

        synonyms = [fix_synonym_case(v) for v in synonyms]

        # Good enough for now
        if name is not None and iupac_name is None:
            iupac_name = name
        return ChemicalData(
            cid=cid,
            iupac_name=iupac_name,
            name=name,
            names=synonyms,
            molecular_weight=MW,
            smiles=smiles,
            inchi_val=inchi_val,
            inchikey=inchikey,
            formula=formula
        )

    def _format_output(self, data: ChemicalData, CAS: str) -> str:
        """Format chemical data as tab-separated string"""
        inchi = data.inchi_val.replace('InChI=1S/', '').replace('InChI=1/', '')
        fields = [
            str(data.cid),
            CAS,
            data.formula,
            # f"{data.molecular_weight:.6f}",
            f"{round(data.molecular_weight,6)}",
            data.smiles,
            inchi,
            data.inchikey,
            data.iupac_name,
            data.name
        ]
        if not data.names:
            data.names = ['']
        return '\t'.join(fields + data.names)

    def process_mol_file(self, filepath: Path) -> Optional[str]:
        """Process a single molecular file and return formatted data"""
        CAS = filepath.stem
        if CAS in self.ignored_CASs:
            return None

        try:
            return self._process_compound(filepath, CAS)
        except Exception as e:
            logger.error(f"Failed to process {filepath}: {str(e)}")
            return None

    def remove_unwanted_compounds_after_processing(self, combined_data):
        if ((combined_data.name is not None and 'polymer' in combined_data.name)
            or (combined_data.iupac_name is not None and 'polymer' in combined_data.iupac_name)):
            return None
        return combined_data


    def _process_compound(self, filepath: Path, CAS: str) -> Optional[str]:
        """Internal method to process compound data"""
        try:
            mol_data = self._get_mol_data(filepath, CAS)
            if not mol_data:
                return None
                
            pubchem_data = self._get_pubchem_additional_data(mol_data)
            combined_data = self._combine_with_scifinder(CAS, mol_data, pubchem_data)
            combined_data = self.remove_unwanted_compounds_after_processing(combined_data)
            if combined_data is None:
                return None
            return self._format_output(combined_data, CAS)
            
        except Exception as e:
            import traceback
            logger.error(
                f"Error processing compound {CAS}:\n"
                f"Traceback:\n{traceback.format_exc()}"
            )            
            return None

    def process_files(self, input_files: List[Path], output_file: Path):
        """Process multiple input files and generate sorted output"""
        temp_output = output_file.with_suffix('.tmp')
        
        with open(temp_output, 'w', encoding='utf-8') as f:
            for filepath in input_files:
                if result := self.process_mol_file(filepath):
                    f.write(f"{result}\n")
        
        # Sort the output file
        sorted_lines = sorted(temp_output.read_text().splitlines(), 
                            key=lambda x: int(x.split('\t')[0]) if x.split('\t')[0].isdigit() else -1)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(sorted_lines) + '\n')
        
        temp_output.unlink()
