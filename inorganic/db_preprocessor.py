from chemicals import *
from chemicals.identifiers import ChemicalMetadataDB
from numpy.testing import assert_allclose
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import json

'''Generate a database of all synonyms that are hardcoded to be included.
MW can be hard coded, as well as pubchem ID.

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

Cl2O2U = {'pubchem': 22033994}
H3N3NdO9 = {'pubchem': 159254}
AgPO3 = {'pubchem': 16172267} # can't find a CAS that doesn't have H
CsPO3 = {'pubchem': 23674445}
B2O3 = {'pubchem': 518682}
BeF2 = {'pubchem': 24589, 'name': 'beryllium fluoride'}
NaBH4 = {'pubchem': 4311764}
B2H6 = {'pubchem': 12544637}

TlNO2 = {'pubchem': 19753737, 'name': 'thallium(I) nitrite'}

SrP2O6 = {'pubchem': 16164193, 'name': 'strontium metaphosphate'}
RbPO3 = {'pubchem': 23674437, 'name': 'rubidium metaphosphate'} # https://i-systems.dechema.de/dbOverview.php?DSI=6479%7C8936
Na2S5 = {'pubchem': 129677176, 'name': 'sodium pentasulfide'} # https://www.wikidata.org/wiki/Q17419133
B4Na2O7 = {'pubchem': 10219853, 'name': 'sodium tetraborate'} # https://commonchemistry.cas.org/detail?cas_rn=1330-43-4

BH4Li = {'pubchem': 4148881, 'name': 'lithium borohydride'}
MnSiO3 = {'pubchem': 13932121, 'name': 'manganese silicate'}
MnF2 = {'pubchem': 24528, 'name': 'manganese difluoride'}
CsO3P = {'pubchem': 23674445, 'name': 'cesium metaphosphate'}
Br2Ga =  {'pubchem': 6394124, 'name': 'gallium(I,III) bromide'}


custom_compounds = {'13536-95-3': D2Se, '16919-19-0': ammonium_hexafluorosilicate,
               '13454-75-6': CsBromate, '10097-32-2': Br, '2099990000-00-0': NaAlO4H4,
               '7558-79-4': Na2HPO4, '7558-80-7': NaH2PO4, '18939-61-2': CuH2O4S,
               '145226-95-5': NaSesquisulfate, '16893-85-9': F6Na2Si,
               '7791-26-6': Cl2O2U, # https://commonchemistry.cas.org/detail?cas_rn=7791-26-6
               '10045-95-1': H3N3NdO9,

               '7787-49-7': BeF2, #https://commonchemistry.cas.org/detail?cas_rn=7787-49-7&search=7787-49-7
               '16949-15-8': BH4Li,
               '16940-66-2': NaBH4,
               '19287-45-7': B2H6,
               '13826-63-6': TlNO2,

               '16090-25-8': SrP2O6,
               '14640-60-9': RbPO3,
               '12034-40-1': Na2S5,
               '1330-43-4': B4Na2O7,
               '7759-00-4': MnSiO3,
               '7782-64-1': MnF2,
               '13783-11-4': CsO3P,
               '18933-31-8': Br2Ga,
               }

for CAS, d in custom_compounds.items():
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
