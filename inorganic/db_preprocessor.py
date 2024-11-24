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
K2S5 = {'name': 'K2S5', 'iupac_name': 'K2S5'}
K2S6 = {'name': 'K2S6', 'iupac_name': 'K2S6'}

H2 = {'pubchem': 783, 'synonyms': ['dihydrogen', 'molecular hydrogen', 'hydrogen (normal)', 'hydrogen gas', 'hydrogen molecule', 'hydrogen-molecule', 'hydrogen (H2)', 
                                    'mol. hydrogen', 'h-h', 'hydrogen [HSDB]', 'hydrogen [MI]', 'hydrogen [WHO-DD]', 'hydrogen, >=99.99%', 
                                    'hydrogen, >=99.999%', 'hydrogen, messer(r) cangas, 99.999%', 'UN 1049', 'UN 1966', 'E949', 'e 949', 'e-949', 
                                    'hidrogeno', 'hydrogene', 'wasserstoff', 'DTXSID9029643', 'CHEBI:18276', 'C00282', 'Q3027893',
                                    'equilibrium hydrogen', 'e-H2', 'eH2', 'e-hydrogen', 'equilibrium-H2', 'equilibrium form of hydrogen', 
                                    'equilibrium modification of hydrogen', 'equilibrium spin isomer of hydrogen', 
                                    'equilibrium nuclear spin isomer of molecular hydrogen', 'thermal equilibrium hydrogen',
                                    'temperature-equilibrated hydrogen', 'thermally equilibrated hydrogen', 'equilibrated hydrogen mixture', 
                                    'temperature-dependent hydrogen mixture', 'equilibrium ortho-para mixture', 'equilibrium o/p hydrogen'], 
                                    'hardcoded_synonyms': True}

H2S4 = {'pubchem': 82836, 'name': 'dihydrogen tetrasulfide'}
H2S5 = {'pubchem': 448216, 'name': 'dihydrogen pentasulfide'}


Magnesium_Hydroxide = {'pubchem': 73981, 'name': 'magnesium hydroxide'}

parahydrogen = {'pubchem': None, 'name': 'parahydrogen', 'iupac_name': 'parahydrogen', 'formula': 'H2','synonyms': ['pH2', 'p-H2', 'para-hydrogen', 'p-hydrogen', 
                        'para-H2','para form of hydrogen','para modification of hydrogen', 'para spin isomer of hydrogen',
                         'para nuclear spin isomer of molecular hydrogen', '800000-49-1', '1333-74-0p'], 'smiles': '', 'inchi': '', 'inchikey': '',
                         'hardcoded_synonyms': True}
orthohydrogen = {'pubchem': None, 'name': 'orthohydrogen', 'iupac_name': 'orthohydrogen', 'formula': 'H2','synonyms': ['oH2', 'o-H2', 'ortho-hydrogen', 'o-hydrogen', 
                        'ortho-H2','ortho form of hydrogen','ortho modification of hydrogen', 'ortho spin isomer of hydrogen',
                         'ortho nuclear spin isomer of molecular hydrogen', '800000-50-4', '1333-74-0o'], 'smiles': '', 'inchi': '', 'inchikey': '',
                         'hardcoded_synonyms': True}
normal_hydrogen =  {'pubchem': None, 'name': 'normal hydrogen', 'iupac_name': 'normal hydrogen', 'formula': 'H2','synonyms': ['nH2', 'n-H2', 'normal-hydrogen', 'n-hydrogen', 
                        'normal-H2','normal form of hydrogen', '800000-51-5',
                        'normal form of hydrogen', 'normal modification of hydrogen', 'room temperature hydrogen mixture',
                         'standard hydrogen mixture', '3:1 ortho:para hydrogen mixture', '75:25 ortho:para hydrogen mixture',
                         'hydrogen, normal'], 
                         'smiles': '', 'inchi': '', 'inchikey': '',
                         'hardcoded_synonyms': True}

H = {'pubchem': 5362549,
    'formula': 'H',
    'smiles': '[H]',
    'inchi': 'H',
    'inchikey': 'YZCKVEUIGOORGS-UHFFFAOYSA-N',
    'iupac_name': 'atomic hydrogen',
    'name': 'atomic hydrogen',
    'synonyms': ['hydrogen (atomic)',
    'hydrogen monatomic',
    'hydrogen atom',
    '1H',
    'H-1',
    'protium',
    'atomic hydrogen',
    'hydrogen radical',
    'hydrogen (H1)',
    'monatomic hydrogen',
    'atomic protium',
    'protium (atomic hydrogen)'],
    'hardcoded_synonyms': True
    }
paradeuterium = {
    'pubchem': None, 
    'name': 'paradeuterium', 
    'iupac_name': 'paradeuterium', 
    'formula': 'D2',
    'synonyms': [
        'pD2', 'p-D2', 'para-deuterium', 'p-deuterium',
        'para-D2', 'para form of deuterium', 'para modification of deuterium', 
        'para spin isomer of deuterium',
        'para nuclear spin isomer of molecular deuterium', 
    ],
    'smiles': '', 
    'inchi': '', 
    'inchikey': '',
    'hardcoded_synonyms': True
}

orthodeuterium = {
    'pubchem': None, 
    'name': 'orthodeuterium', 
    'iupac_name': 'orthodeuterium', 
    'formula': 'D2',
    'synonyms': [
        'oD2', 'o-D2', 'ortho-deuterium', 'o-deuterium',
        'ortho-D2', 'ortho form of deuterium', 'ortho modification of deuterium', 
        'ortho spin isomer of deuterium',
        'ortho nuclear spin isomer of molecular deuterium', 
    ],
    'smiles': '', 
    'inchi': '', 
    'inchikey': '',
    'hardcoded_synonyms': True
}

normal_deuterium = {
    'pubchem': None, 
    'name': 'normal deuterium', 
    'iupac_name': 'normal deuterium', 
    'formula': 'D2',
    'synonyms': [
        'nD2', 'n-D2', 'normal-deuterium', 'n-deuterium', 'deuterium, normal',
        'normal-D2', 'normal form of deuterium',
        'normal form of deuterium', 'normal modification of deuterium', 
        'room temperature deuterium mixture',
        'standard deuterium mixture', '2:1 ortho:para deuterium mixture', 
    ],
    'smiles': '', 
    'inchi': '', 
    'inchikey': '',
    'hardcoded_synonyms': True
}
custom_compounds = {    
    '13536-95-3': D2Se, '16919-19-0': ammonium_hexafluorosilicate,
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

               '12136-50-4': K2S5,
               '37188-07-1': K2S6,

                '1333-74-0': H2,
                '13845-25-5': H2S4,
                '13845-24-4': H2S5,

                '2099490000-00-0': parahydrogen,
                '2099479000-00-0': orthohydrogen,
                '2099474000-00-0': normal_hydrogen,
                '12385-13-6': H,

                '2099458000-00-0': paradeuterium,
                '2099453000-00-0': orthodeuterium,
                '2099437000-00-0': normal_deuterium,

                '1309-42-8': Magnesium_Hydroxide,
                '1317-43-7': {'pubchem': 14791},
                '14457-84-2': {'pubchem': -1}, # none
    
               }

custom_compounds['10022-50-1'] = {'synonyms': ['nitrogen dioxyfluoride']}
custom_compounds['10544-72-6'] = {'synonyms': ['dinitrogen tetroxide']}
custom_compounds['13536-59-9'] = {'synonyms': ['deuteriumbromide']}

custom_compounds['13537-15-0'] = {
    'pubchem': -1,
    'smiles': '[Eu+3].[Eu+3].[O-]S(=O)(=O)[O-].[O-]S(=O)(=O)[O-].[O-]S(=O)(=O)[O-]',
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
