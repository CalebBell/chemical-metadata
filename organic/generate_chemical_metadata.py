import sys
from pathlib import Path
from chemical_metadata_tools import parse_CAS_data
from chemical_metadata_tools import generate_db_tools
import appdirs
import os
main_cache_dir = appdirs.user_cache_dir(appname='chemical_metadata')
if not os.path.exists(main_cache_dir):
    os.mkdir(main_cache_dir)

cache_dir = os.path.join(main_cache_dir, 'organic')
if not os.path.exists(cache_dir):
    os.mkdir(cache_dir)

workdir_path = Path(__file__).parent

ignore_CASs = {'2025-56-1', '2122-48-7', '2143-68-2', '2226-96-2', '2370-18-5',
               '2597-44-6', '2896-70-0', '25085-53-4', '26063-22-9', 
                '38439-19-9', '29300-20-7',
               '25037-58-5', '9003-17-2',

               # polymers
               '25722-06-9', '25067-06-5', '33411-58-4', '25067-64-5', '24979-98-4', '30209-80-4',
               '28182-81-2', '26744-16-1', '52224-87-0', '2358-84-1', '24991-43-3', 
               '25068-01-3', '36486-90-5', '25722-18-3', '25267-51-0', '9002-83-9', '31213-03-3',  
               '25686-28-6', '25036-32-2', '26353-15-1', '27732-42-9',

                # di/tri etc mers which we don't have a good structure for
               '7756-94-7',
            }
def main():
    # First run the PDF parser
    parse_CAS_data.parse_pdfs_in_subfolder(workdir_path)
    
    # Parse command line arguments
    args = sys.argv[1:]
    
    output_file = Path(args.pop()) if args else workdir_path /'organic-tmp.tsv'
    input_files = [Path(f) for f in args]
    # input_files = [Path(workdir_path) / 'mol' / '73463-39-5.mol']

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

if __name__ == '__main__':
    main()
