from chemical_metadata import *
import chemicals
from chemicals import molecular_weight, nested_formula_parser

pdf_data = open('Parsed CAS metadata.json').read()
pdf_data = json.loads(pdf_data)

folder = os.path.dirname(__file__)

# Initialize with parsed scifinder data
fixed_data = {}
for k, d in pdf_data.items():
    
    # Collect the synonyms
    synonyms = []
    if 'Deleted CAS' in d:
        synonyms += d['Deleted CAS']
    if 'Other Names' in d:
        synonyms += d['Other Names']
    if 'Alternate CAS' in d:
        synonyms += d['Alternate CAS']
    
    new_d = {'CAS': k, 'synonyms': synonyms}
    
    # Handle the lowercase name
    if "Name" in d:
        new_d['name'] = d['Name']
    
    # Add 'mol' initialization file path
    mol_path = os.path.join(folder, 'mol', k + '.mol')
    if os.path.exists(mol_path):
        new_d['mol'] = mol_path
        
    fixed_data[k] = new_d

D2Se = {'formula': 'D2Se', 'MW': molecular_weight(nested_formula_parser('D2Se'))}
ammonium_hexafluorosilicate = {'pubchem': 28145}
CsBromate = {'pubchem': 23685550}
Br = {'pubchem': 5360770}
NaAlO4H4 = {'pubchem': 166673}
Na2HPO4 = {'pubchem': 24203}
NaH2PO4 = {'pubchem': 24204}
CuH2O4S = {'pubchem': 6536471}
NaSesquisulfate = {'pubchem': 19098072}
F6Na2Si = {'pubchem':  28127}
water = {'pubchem': 962}


custom_compounds = {'13536-95-3': D2Se, '16919-19-0': ammonium_hexafluorosilicate,
               '13454-75-6': CsBromate, '10097-32-2': Br, '2099990000-00-0': NaAlO4H4,
               '7558-79-4': Na2HPO4, '7558-80-7': NaH2PO4, '18939-61-2': CuH2O4S,
               '145226-95-5': NaSesquisulfate, '16893-85-9': F6Na2Si,
               '7732-18-5': water}
for k, v in custom_compounds.items():
    v['CAS'] = k 
fixed_data.update(custom_compounds)

# List of broken compounds
no_structure = ['12031-63-9', '12032-20-1', '12060-08-1', '12069-94-2', '12070-10-9', '1314-84-7',
                '1317-61-9', '1327-41-9', '1345-07-9', '13454-75-6', '13536-95-3', '16893-85-9',
                '16919-19-0', '7803-62-5', '93401-23-1']
for k in no_structure:
    if k in fixed_data:
        del fixed_data[k]


to_write = []
for d in list(fixed_data.values()):
    processed_json = process(d, use_cache=True)
    to_write.append(processed_json)

output = 'Inorganic db2.tsv'
write_database(to_write, output)

os.system('sort -n -o "%s" "%s"' %(output, output))
