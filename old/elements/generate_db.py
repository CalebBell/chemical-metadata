from chemical_metadata import *
import chemicals

process_CASs = chemicals.elements.CAS_by_number_standard + chemicals.elements.homonuclear_elemental_singlets_CASs
fixed_data = {k: {'CAS': k} for k in process_CASs}

fixed_data['54083-77-1']['pubchem'] = 135476787 # Darmstadtium
fixed_data['54083-77-1']['inchi'] = 'InChI=1S/Ds'

fixed_data['54386-24-2']['pubchem'] = 135476786 # Roentgenium
fixed_data['54386-24-2']['inchi'] = 'InChI=1S/Rg'

fixed_data['54084-26-3']['pubchem'] = 135476785 # copernicium
fixed_data['54084-26-3']['inchi'] = 'InChI=1S/Cn'

#fixed_data['54084-70-7']['pubchem'] =    # Nihonium # Missing
fixed_data['54084-70-7']['inchi'] = 'InChI=1S/Nh'
fixed_data['54085-16-4']['inchi'] = 'InChI=1S/Fl'
fixed_data['54085-64-2']['inchi'] = 'InChI=1S/Mc'
fixed_data['54100-71-9']['inchi'] = 'InChI=1S/Lv'
fixed_data['54101-14-3']['inchi'] = 'InChI=1S/Ts'
fixed_data['54144-19-3']['inchi'] = 'InChI=1S/Og'


#pdf_data = open('Parsed CAS metadata.json').read()
#pdf_data = json.loads(pdf_data)



to_write = []
for d in fixed_data.values():
    processed_json = process(d)
    to_write.append(processed_json)

write_database(to_write, 'elements.tsv')
