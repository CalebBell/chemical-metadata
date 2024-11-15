from __future__ import annotations
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
import logging
from dataclasses import dataclass

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from pubchempy import get_compounds, Compound
import json
from collections import Counter
from thermo import serialize_formula
from fcache.cache import FileCache
from chemical_metadata_tools import parse_CAS_data
from chemical_metadata_tools.parse_CAS_data import lower_case_first_letter_name
from chemical_metadata_tools import generate_db_tools
import appdirs
import os
main_cache_dir = appdirs.user_cache_dir(appname='chemical_metadata')
if not os.path.exists(main_cache_dir):
    os.mkdir(main_cache_dir)

cache_dir = os.path.join(main_cache_dir, 'inorganic')
if not os.path.exists(cache_dir):
    os.mkdir(cache_dir)

workdir_path = Path(__file__).parent


ignore_CASs = {'7782-42-5', # graphite/diamond/carbon issue
               '13550-53-3',
               '12024-10-1', 
               '12024-22-5',
               '12536-51-5', 
               '12536-52-6',
            '12536-51-5', '12536-52-6', '13966-04-6', '13940-21-1',
           '12033-49-7', '12061-70-0', '12122-00-8', '12143-17-8',
           '12184-84-8', '12184-90-6', '12505-77-0', '12596-60-0',
           '13494-92-3', '13770-40-6', '13774-94-2', '13940-21-1',
           '13966-04-6', '16904-65-7', '2074-87-5', '30664-12-1',
           '3170-80-7', '3352-57-6', '34518-80-4','3889-75-6',
           '51912-69-7', '19165-34-5', '30937-38-3', '12007-41-9', '12007-60-2'
        }

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# note 7558-80-7.mol should be H2NaO4P don't know why the mol from pubchem was wrong
def main():
    # First run the PDF parser
    logger.info("Running PDF parser...")
    parse_CAS_data.parse_pdfs_in_subfolder(workdir_path)
    
    # Parse command line arguments
    args = sys.argv[1:]
    # if not args:
    #     logger.error("No input files specified")
    #     sys.exit(1)
    
    output_file = Path(args.pop()) if args else workdir_path /'inorganic-tmp.tsv'
    input_files = [Path(f) for f in args]
    # input_files = [Path(workdir_path) / 'mol' / '13821-20-0.mol']

    # Include all files from syn_data if INCLUDE_EVERYTHING is True
    INCLUDE_EVERYTHING = True
    if INCLUDE_EVERYTHING:
        processor = generate_db_tools.ChemicalMetadataProcessor(cache_dir, workdir_path, require_charge=0, ignore_CASs=ignore_CASs)
        
        for CAS, d in processor.syn_data.items():
            if ('pubchem' in d or 'formula' in d) and not any(CAS in str(f) for f in input_files):
                mol_path =workdir_path / Path('mol') / f"{CAS}.mol"
                # even if the mol file doesn't exist we still append, will just result another rdkit failure and process OK
                if len(input_files) != 1:
                    input_files.append(mol_path)
    
    # Process files
    processor.process_files(input_files, output_file)
    logger.info(f"Processing complete. Output written to {output_file}")

if __name__ == '__main__':
    main()
