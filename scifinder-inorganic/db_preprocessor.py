from chemicals import *
from chemicals.identifiers import ChemicalMetadataDB
from numpy.testing import assert_allclose
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import json

'''Generate a database of all dynonyms that are hardcoded to be included.
MW can be ard coded, as well as pubchem ID.

'''
db = ChemicalMetadataDB(elements=False, main_db='Inorganic db.tsv', user_dbs=[])
db.autoload_main_db()

good_syns = {}

for CAS, d in db.CAS_index.items():
    CAS = d.CASs
    if CAS in good_syns:
        pass
    else:
        good_syns[CAS] = {}
        good_syns[CAS]['synonyms'] = []

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



custom_ions = {'13536-95-3': D2Se, '16919-19-0': ammonium_hexafluorosilicate,
               '13454-75-6': CsBromate, '10097-32-2': Br, '2099990000-00-0': NaAlO4H4,
               '7558-79-4': Na2HPO4, '7558-80-7': NaH2PO4, '18939-61-2': CuH2O4S,
               '145226-95-5': NaSesquisulfate, '16893-85-9': F6Na2Si}

for CAS, d in custom_ions.items():
    if CAS in good_syns:
        good_syns[CAS].update(d)
    else:
        good_syns[CAS] = d


def write():
    f = open('Good synoynms by CAS.json', 'w')
    json.dump(good_syns, f, indent=2, separators=(',', ': '), sort_keys=True)
    f.close()
if __name__ == '__main__':
    write()
