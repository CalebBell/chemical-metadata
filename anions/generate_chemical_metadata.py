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
from chemical_metadata_tools import parse_CAS_data
from chemical_metadata_tools import generate_db_tools
import appdirs
import os
main_cache_dir = appdirs.user_cache_dir(appname='chemical_metadata')
if not os.path.exists(main_cache_dir):
    os.mkdir(main_cache_dir)

cache_dir = os.path.join(main_cache_dir, 'anions')
if not os.path.exists(cache_dir):
    os.mkdir(cache_dir)


workdir_path = Path(__file__).parent


ignore_CASs = {
            '12339-27-4',  # H2P2+ can't draw it, not in any db
            '68111-10-4'   # AuBr2+ bromine charge problem, smiles is [Br-][Au+][Br-]
        }

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    # First run the PDF parser
    logger.info("Running PDF parser...")
    parse_CAS_data.parse_pdfs_in_subfolder(workdir_path)
    
    # Parse command line arguments
    args = sys.argv[1:]
    # if not args:
    #     logger.error("No input files specified")
    #     sys.exit(1)
    
    output_file = Path(args.pop()) if args else workdir_path /'cations-tmp.tsv'
    input_files = [Path(f) for f in args]
    # input_files = [Path(workdir_path) / 'mol' / '3198-32-1.mol']
    # input_files = [Path(workdir_path) / 'mol' / '18130-74-0.mol']
    # input_files = [Path(workdir_path) / 'mol' / '766-76-7.mol']
    # input_files = [Path(workdir_path) / 'mol' / '107596-48-5.mol']
    # input_files = [Path(workdir_path) / 'mol' / '11062-77-4.mol']

    # input_files = [Path(workdir_path) / 'mol' / '766-76-7.mol']
    # input_files = [Path(workdir_path) / 'mol' / '3198-32-1.mol']

    # Include all files from syn_data if INCLUDE_EVERYTHING is True
    INCLUDE_EVERYTHING = True
    if INCLUDE_EVERYTHING:
        processor = generate_db_tools.ChemicalMetadataProcessor(cache_dir, workdir_path, require_not_charge=0, ignore_CASs=ignore_CASs)
        
        for CAS, d in processor.syn_data.items():
            if ('pubchem' in d or 'formula' in d) and not any(CAS in str(f) for f in input_files):
                mol_path =workdir_path / Path('mol') / f"{CAS}.mol"
                # even if the mol file doesn't exist we still append, will just result another rdkit failure and process OK
                if len(input_files) != 1:
                    input_files.append(mol_path)
    
    # Process files
    processor.process_files(input_files, output_file)
    processor.write_preferences(str(output_file)[0:-4] + '_preferences.json')
    logger.info(f"Processing complete. Output written to {output_file}")

if __name__ == '__main__':
    main()
