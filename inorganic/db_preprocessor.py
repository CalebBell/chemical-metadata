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
B2H6 = {
    'pubchem': 12544637,
    'name': 'diborane',
    'iupac_name': 'diborane(6)',
    'formula': 'B2H6',
    'synonyms': [
        'diborane(6)',
        'diboron hexahydride',
        'diborane (B2H6)',
        'boron hydride (B2H6)',
        'hydrogen boride (H6B2)',
    ],
    'smiles': '[BH2]1[H][BH2][H]1',
    'inchi': 'InChI=1S/B2H6/c1-3-2-4-1/h1-2H2',
    'inchikey': 'KLDBIFITUCWVCC-UHFFFAOYSA-N',
    'hardcoded_synonyms': True,
    'preferred': True
}
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
                                    'hardcoded_synonyms': True, 'preferred': True}

H2S4 = {'pubchem': 82836, 'name': 'dihydrogen tetrasulfide'}
H2S5 = {'pubchem': 448216, 'name': 'dihydrogen pentasulfide'}


Magnesium_Hydroxide = {'pubchem': 73981, 'name': 'magnesium hydroxide'}

parahydrogen = {'pubchem': None, 'name': 'parahydrogen', 'iupac_name': 'parahydrogen', 'formula': 'H2','synonyms': ['pH2', 'p-H2', 'para-hydrogen', 'p-hydrogen', 
                        'para-H2','para form of hydrogen','para modification of hydrogen', 'para spin isomer of hydrogen',
                         'para nuclear spin isomer of molecular hydrogen', '800000-49-1', '1333-74-0p'], 'smiles': '', 'inchi': '', 'inchikey': '',
                         'hardcoded_synonyms': True, 'preferred': False}
orthohydrogen = {'pubchem': None, 'name': 'orthohydrogen', 'iupac_name': 'orthohydrogen', 'formula': 'H2','synonyms': ['oH2', 'o-H2', 'ortho-hydrogen', 'o-hydrogen', 
                        'ortho-H2','ortho form of hydrogen','ortho modification of hydrogen', 'ortho spin isomer of hydrogen',
                         'ortho nuclear spin isomer of molecular hydrogen', '800000-50-4', '1333-74-0o'], 'smiles': '', 'inchi': '', 'inchikey': '',
                         'hardcoded_synonyms': True, 'preferred': False}
normal_hydrogen =  {'pubchem': None, 'name': 'normal hydrogen', 'iupac_name': 'normal hydrogen', 'formula': 'H2','synonyms': ['nH2', 'n-H2', 'normal-hydrogen', 'n-hydrogen', 
                        'normal-H2','normal form of hydrogen', '800000-51-5',
                        'normal form of hydrogen', 'normal modification of hydrogen', 'room temperature hydrogen mixture',
                         'standard hydrogen mixture', '3:1 ortho:para hydrogen mixture', '75:25 ortho:para hydrogen mixture',
                         'hydrogen, normal'], 
                         'smiles': '', 'inchi': '', 'inchikey': '',
                         'hardcoded_synonyms': True, 'preferred': False}

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
    'hardcoded_synonyms': True,
    'preferred': False,
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
    'hardcoded_synonyms': True,
    'preferred': False
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
    'hardcoded_synonyms': True,
    'preferred': False
}

V = {
    'pubchem': 23990,
    'name': 'vanadium',
    'iupac_name': 'vanadium',
    'formula': 'V',
    'synonyms': [
        'atomic vanadium',
        'vanadium element',
        'vanadium-51',
        'elemental vanadium',
        'vanadium atom',
    ],
    'hardcoded_synonyms': True,
}
# synonyms conflicting with Nitrosyl https://commonchemistry.cas.org/detail?cas_rn=14452-93-8
NO = {
        'synonyms': [
        'nitric oxide',
        'nitrogen monoxide',
        'nitrogen oxide (NO)',
        'nitrogen(II) oxide',
        'mononitrogen monoxide',
        'nitrogen monooxide',
        'oxido nitrico',
        'oxyde nitrique',
        'Stickstoff(II)-oxid',
        'EDRF',
        'nitric oxide (NO)',
        'nitrogen oxide',
        'CCRIS 4319',
        'CHEBI:16480',
        'CHEMBL1200689',
        'C00533',
        'D00074',
        'DB00435',
        'DTXCID90938',
        'HSDB 1246',
        '31c4ky9esh'
    ],
    'hardcoded_synonyms': True
}
O = {
    'pubchem': 159832,
    'formula': 'O',
    'smiles': '[O]',
    'inchi': 'O',
    'inchikey': 'QVGXLLKOCUKJST-UHFFFAOYSA-N',
    'iupac_name': 'atomic oxygen',
    'name': 'atomic oxygen',
    'synonyms': [
        'oxygen atom',
        'monooxygen',
        'singlet oxygen',
        'atomic oxygen',
        'oxidanediyl',
        'O(3P)',
        'ground state oxygen atom',
        'CHEBI:29194',
        'monatomic oxygen',
    ],
    'hardcoded_synonyms': True
}

O2 = {
    'pubchem': 977,
    'formula': 'O2',
    'smiles': 'O=O',
    'inchi': 'O2/c1-2',
    'inchikey': 'MYMOFIZGZYHOMD-UHFFFAOYSA-N',
    'iupac_name': 'molecular oxygen',
    'name': 'oxygen',
    'synonyms': [
        'dioxygen',
        'molecular oxygen',
        'oxygen gas',
        'oxygen molecule',
        'pure oxygen',
        'diatomic oxygen',
        'triplet oxygen',
        'ground state oxygen',
        "oxygenium",
        "liquid oxygen",
        "oxygen-16",
        "sauerstoff",
        "oxigeno",
        "compressed oxygen",
        "CCRIS 1228",
        "HSDB 5054",
        "ins no.948",
        "UNII-S88TT14065",
        "CHEBI:15379",
        "E948",
        "oxygen 16",
        "singlet dioxygen",
        "triplet dioxygen",
        "singlet molecular oxygen",
        "triplet molecular oxygen",
        "CHEMBL1234886",
        "DTXCID0017681",
        "CHEBI:26689",
        "CHEBI:27140",
        "DB09140",
        "UN 1072",
        "UN 1073",
    ],
    'hardcoded_synonyms': True
}


# Atomic Nitrogen (N)
N = {
    'pubchem': 57370662,
    'formula': 'N',
    'smiles': '[N]',
    'inchi': 'N',
    'inchikey': 'QJGQUHMNIGDVPM-UHFFFAOYSA-N',
    'iupac_name': 'atomic nitrogen',
    'name': 'atomic nitrogen',
    'synonyms': [
        'nitrogen atom',
        'mononitrogen',
        'nitrogen, atomic',
        'nitrogen radical',
        'singlet nitrogen',
        'nitrogen(atom)',
        'monatomic nitrogen',
    ],
    'hardcoded_synonyms': True
}

# Molecular Nitrogen (N2)
N2 = {
    'pubchem': 947,
    'formula': 'N2',
    'smiles': 'N#N',
    'inchi': 'N2/c1-2',
    'inchikey': 'IJGRMHOSHXDMSA-UHFFFAOYSA-N',
    'iupac_name': 'molecular nitrogen',
    'name': 'nitrogen',
    'synonyms': [
        'diatomic nitrogen',
        'dinitrogen',
        'molecular nitrogen',
        'nitrogen gas',
        'nitrogen molecule',
    ],
    'hardcoded_synonyms': True
}

# Atomic Fluorine (F)
F = {
    'pubchem': 5360525,
    'formula': 'F',
    'smiles': '[F]',
    'inchi': 'F',
    'inchikey': 'YCKRFDGAMUMZLT-UHFFFAOYSA-N',
    'iupac_name': 'atomic fluorine',
    'name': 'atomic fluorine',
    'synonyms': [
        'fluorine atom',
        'monofluorine',
        'fluorine, atomic',
        'monatomic fluorine',
        'CHEBI:30239'
    ],
    'hardcoded_synonyms': True
}

# Molecular Fluorine (F2)
F2 = {
    'pubchem': 24524,
    'formula': 'F2',
    'smiles': 'FF',
    'inchi': 'F2/c1-2',
    'inchikey': 'PXGOKWXKJXAPGV-UHFFFAOYSA-N',
    'iupac_name': 'molecular fluorine',
    'name': 'fluorine',
    'synonyms': [
        'diatomic fluorine',
        'difluorine',
        'molecular fluorine',
        'fluorine gas',
        'fluorine-19',
        'mol. fluorine',
        'F2',
        'fluor',
        'fluoro',
        'fluorine-elemental',
        'fluorine 19',
        'fluorine, compressed',
        'CHEBI:30236',
        'CHEBI:36895',
        'UN1045',
    ],
    'hardcoded_synonyms': True
}
# Atomic Chlorine (Cl)
Cl = {
    'pubchem': 5360523,
    'formula': 'Cl',
    'smiles': '[Cl]',
    'inchi': 'Cl',
    'inchikey': 'ZAMOUSCENKQFHK-UHFFFAOYSA-N',
    'iupac_name': 'atomic chlorine',
    'name': 'atomic chlorine',
    'synonyms': [
        'chlorine atom',
        'monochlorine',
        'chlorine, atomic',
        'CHEBI:29311',
        'atomic Cl',
        'Cl atom',
        'monatomic chlorine',
    ],
    'hardcoded_synonyms': True
}

# Molecular Chlorine (Cl2)
Cl2 = {
    'pubchem': 24526,
    'formula': 'Cl2',
    'smiles': 'ClCl',
    'inchi': 'Cl2/c1-2',
    'inchikey': 'KZBUYRJDOAKODT-UHFFFAOYSA-N',
    'iupac_name': 'molecular chlorine',
    'name': 'chlorine',
    'synonyms': [
        'diatomic chlorine',
        'dichlorine',
        'molecular chlorine',
        'chlorine molecule',
        'chlorine mol.',
        'bertholite',
        'chlorinum',
        'chlore',
        'chlor',
        'cloro',
        'chloor',
        'liquid chlorine',
        'CHEBI:29310',
        'CHEBI:33432',
        'UN1017',
        'UNII-4R7X1O2820',
        'EINECS 231-959-5',
        'HSDB 206',
        'DTXSID1020273',
        'DB11109',
        'CCRIS 2280',
        'EC 231-959-5',
        'NS00076171'
    ],
    'hardcoded_synonyms': True
}
I = {
    'pubchem': 5360629,
    'formula': 'I',
    'smiles': '[I]',
    'inchi': 'I',
    'inchikey': 'ZCYVEMRRCGMTRW-UHFFFAOYSA-N',
    'iupac_name': 'atomic iodine',
    'name': 'atomic iodine',
    'synonyms': [
        'iodine atom',
        'iodine, atomic',
        'monoiodine',
        'CHEBI:33115',
        'AKOS024438077',
        'monatomic iodine',
    ],
    'hardcoded_synonyms': True
}

I2 = {
    'pubchem': 807,
    'formula': 'I2',
    'smiles': 'II',
    'inchi': 'I2/c1-2',
    'inchikey': 'PNDPGZBMCMUPRI-UHFFFAOYSA-N',
    'iupac_name': 'molecular iodine',
    'name': 'iodine',
    'synonyms': [
        'diatomic iodine',
        'diiodine',
        'molecular iodine',
        'iodine molecule',
        'iodine crystals',
        'iodine sublimed',
        'iodium',
        'iodum',
        'iode',
        'jood',
        'jod',
        'iodine (II)',
        'iodine [II]',
        'iodio',
        'iodine, crystalline',
        'iodine, elemental',
        'iodine, resublimed',
        'CHEBI:17606',
        'NSC 42355',
        'EC 231-442-4',
        'UNII-9679TC07X4',
        'CHEMBL1201225',
        'DB05382'
    ],
    'hardcoded_synonyms': True
}

# Atomic Bromine (Br)
Br = {
    'pubchem': 5360770,
    'formula': 'Br',
    'smiles': '[Br]',
    'inchi': 'Br',
    'inchikey': 'WKBOTKDWSSQWDR-UHFFFAOYSA-N',
    'iupac_name': 'atomic bromine',
    'name': 'atomic bromine',
    'synonyms': [
        'atomic bromine',
        'bromine atom',
        'bromine, atomic',
        'monobromine',
        'monatomic bromine',
    ],
    'hardcoded_synonyms': True
}

# Molecular Bromine (Br2)
Br2 = {
    'pubchem': 24408,
    'formula': 'Br2',
    'smiles': 'BrBr',
    'inchi': 'Br2/c1-2',
    'inchikey': 'GDTBXPJZTBHREO-UHFFFAOYSA-N',
    'iupac_name': 'molecular bromine',
    'name': 'bromine',
    'synonyms': [
        'diatomic bromine',
        'dibromine',
        'molecular bromine',
        'bromine molecule',
        'bromine element',
        'dibromane',
        'bromine liquid',
        # Database identifiers
        'CHEBI:29224',
        'HSDB 514',
        'EINECS 231-778-1',
        'MFCD00010896',
        'UNII-SBV4XY874G',
        'DTXSID1035238',
        'DTXCID9015238',
        'EC 231-778-1',
        'UN1744',
        'UN 1744',
        'AKOS015897100',
        'BCP26202',
        'NS00075687',
        # Multilingual names
        'brom',
        'brome',
        'bromo',
        'broom',
        'bromium',
        # Product identifiers
        'abws 20',
        'bromitize plus',
        'caswell no. 112',
        'B2414',
        'j-519936',
    ],
    'hardcoded_synonyms': True
}
O2Si = {
    'pubchem': 24261,
    'formula': 'O2Si',
    'smiles': 'O=[Si]=O',
    'inchi': 'O2Si/c1-3-2',
    'inchikey': 'VYPSYNLAJGMNEJ-UHFFFAOYSA-N',
    'iupac_name': 'dioxosilane',
    'name': 'silica',
    'synonyms': [
        'silicon dioxide',
        'silicon(IV) oxide',
        'silica gel',
        'dioxosilane',
        'silicon oxide',
        'silicic anhydride',
        'silicon dioxide, amorphous',
        'silica, amorphous',
        'silica, colloidal',
        'synthetic amorphous silica',
        'silicon(IV) oxide (SiO2)',
        'silicon dioxide (amorphous)'
    ],
    'hardcoded_synonyms': True,
    'preferred': True,
}

# CAS number mapping

custom_compounds = {    
    '13536-95-3': D2Se, '16919-19-0': ammonium_hexafluorosilicate,
               '13454-75-6': CsBromate,  '2099990000-00-0': NaAlO4H4,
               '7558-79-4': Na2HPO4, '7558-80-7': NaH2PO4, '18939-61-2': CuH2O4S,
               '145226-95-5': NaSesquisulfate, '16893-85-9': F6Na2Si,
               '7791-26-6': Cl2O2U, # https://commonchemistry.cas.org/detail?cas_rn=7791-26-6
               '10045-95-1': H3N3NdO9,

               '7787-49-7': BeF2, #https://commonchemistry.cas.org/detail?cas_rn=7787-49-7&search=7787-49-7
               '16949-15-8': BH4Li,
               '16940-66-2': NaBH4,
               '19287-45-7': B2H6,
               '13826-63-6': TlNO2,

               '7440-62-2': V,

               '7782-44-7': O2,
               '17778-80-2': O,

                '17778-88-0': N,
                '7727-37-9': N2,

                '14762-94-8': F,
                '7782-41-4': F2,

                '22537-15-1': Cl,
                '7782-50-5': Cl2,

                '14362-44-8': I,
                '7553-56-2': I2,

                '10097-32-2': Br,
                '7726-95-6': Br2,
                
                '1333-74-0': H2,
                '2099490000-00-0': parahydrogen,
                '2099479000-00-0': orthohydrogen,
                '2099474000-00-0': normal_hydrogen,
                '12385-13-6': H,

                '7782-39-0': {'preferred': True}, # Deuterium
                '2099458000-00-0': paradeuterium,
                '2099453000-00-0': orthodeuterium,
                '2099437000-00-0': normal_deuterium,

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
               
               '10102-43-9': NO,

                '13845-25-5': H2S4,
                '13845-24-4': H2S5,


                '1309-42-8': Magnesium_Hydroxide,
                '1317-43-7': {'pubchem': 14791},
                '14457-84-2': {'pubchem': -1}, # none


                '115-25-3': {'preferred': True}, # Octafluorocyclobutane
                '360-89-4': {'preferred': False}, # Perfluoro-2-butene

                '1309-42-8': {'preferred': True}, # Magnesium hydroxide
                '1317-43-7': {'preferred': False}, # Brucite (Mg(OH)2)

                '1318-23-6': {'preferred': True}, #Boehmite (Al(OH)O)
                '14457-84-2': {'preferred': False}, # Diaspore (AlHO2)

                '13776-62-0': {'preferred': False}, #trans Dinitrogen difluoride F2N2
                '13812-43-6': {'preferred': True}, #cis Dinitrogen difluoride F2N2 more stable
                '10578-16-2': {'preferred': False}, # generic F2N2

                '10294-56-1': {'preferred': True}, # Phosphorous acid
                '13598-36-2': {'preferred': False}, # Phosphonic acid


                '7631-86-9': O2Si,
                '14808-60-7': {'preferred': False}, # quartz
    
               }

custom_compounds['10022-50-1'] = {'synonyms': ['nitrogen dioxyfluoride']}
custom_compounds['10544-72-6'] = {'synonyms': ['dinitrogen tetroxide']}
custom_compounds['13536-59-9'] = {'synonyms': ['deuteriumbromide']}

custom_compounds['13537-15-0'] = {
    'pubchem': -1,
    'smiles': '[Eu+3].[Eu+3].[O-]S(=O)(=O)[O-].[O-]S(=O)(=O)[O-].[O-]S(=O)(=O)[O-]',
}

# Fe2O3 entries
# Prefer the general iron(III) oxide over specific mineral forms
custom_compounds['1309-37-1'] = {'preferred': True}  # Iron(III) oxide (Fe2O3)
custom_compounds['1317-60-8'] = {'preferred': False}  # Hematite (α-Fe2O3)
custom_compounds['12134-66-6'] = {'preferred': False}  # Maghemite (γ-Fe2O3)

# Fe3O4 entries
# Prefer the general iron(II,III) oxide over magnetite mineral form
custom_compounds['1317-61-9'] = {'preferred': True}  # Iron(II,III) oxide (Fe3O4)
custom_compounds['1309-38-2'] = {'preferred': False}  # Magnetite (Fe3O4)

# FeO entry
# Only one form so it's preferred by default
custom_compounds['1345-25-1'] = {'preferred': True}  # Iron(II) oxide/Wüstite (FeO)


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
