from chemical_metadata import *
import chemicals
from chemicals import molecular_weight, nested_formula_parser

pdf_data = open('Parsed scifinder metadata.json').read()
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


custom_compounds = {}
for k, v in custom_compounds.items():
    v['CAS'] = k 
fixed_data.update(custom_compounds)

# List of broken compounds
no_structure = ['12135-76-1', '1309-37-1', '1327-53-3', # not even organic
    '75899-69-3',
    
    # no structure - improve?
    '9034-32-6',
    ]
for k in no_structure:
    if k in fixed_data:
        del fixed_data[k]


to_write = []
for d in list(fixed_data.values()):
    processed_json = process(d, use_cache=True)
    to_write.append(processed_json)

output = 'Organic_db2.tsv'
write_database(to_write, output)

os.system('sort -n -o "%s" "%s"' %(output, output))
