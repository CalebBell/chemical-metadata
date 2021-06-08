from chemical_metadata import *
import chemicals

process_CASs = chemicals.elements.CAS_by_number_standard + chemicals.elements.homonuclear_elemental_singlets_CASs



pdf_data = open('Parsed scifinder metadata.json').read()
pdf_data = json.loads(pdf_data)

to_write = []
for CAS in process_CASs:
    processed_json = process_chemical(CAS)
    to_write.append(processed_json)
