from thermo import *
from chemicals import *
from chemicals.identifiers import ChemicalMetadataDB, PUBCHEM_EXAMPLE_DB_NAME, folder
from numpy.testing import assert_allclose
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import json
import os

db = ChemicalMetadataDB(elements=False, main_db=os.path.join(folder, PUBCHEM_EXAMPLE_DB_NAME), user_dbs=[])
db.autoload_main_db()

good_syns = {}

for CAS, d in db.CAS_index.items():
    CAS = d.CASs
    if CAS in good_syns:
        pass
    else:
        good_syns[CAS] = {}
        good_syns[CAS]['synonyms'] = []


good_syns['502-52-3']['synonyms'].extend(['dipalmitin']) # There may be other chemicals with this name but biosteam was expecting this answer

good_syns['74-82-8']['synonyms'].extend(['C1', 'nC1', 'n-C1'])
good_syns['74-84-0']['synonyms'].extend(['C2', 'nC2', 'n-C2'])
good_syns['74-98-6']['synonyms'].extend(['C3', 'nC3', 'n-C3'])
good_syns['106-97-8']['synonyms'].extend(['C4', 'nC4', 'n-C4'])
good_syns['109-66-0']['synonyms'].extend(['C5', 'nC5', 'n-C5'])
good_syns['110-54-3']['synonyms'].extend(['C6', 'nC6', 'n-C6'])
good_syns['142-82-5']['synonyms'].extend(['C7', 'nC7', 'n-C7'])
good_syns['111-65-9']['synonyms'].extend(['C8', 'nC8', 'n-C8'])
good_syns['111-84-2']['synonyms'].extend(['C9', 'nC9', 'n-C9'])
good_syns['124-18-5']['synonyms'].extend(['C10', 'nC10', 'n-C10'])
good_syns['1120-21-4']['synonyms'].extend(['C11', 'nC11', 'n-C11'])
good_syns['112-40-3']['synonyms'].extend(['C12', 'nC12', 'n-C12'])
good_syns['629-50-5']['synonyms'].extend(['C13', 'nC13', 'n-C13'])
good_syns['629-59-4']['synonyms'].extend(['C14', 'nC14', 'n-C14'])
good_syns['629-62-9']['synonyms'].extend(['C15', 'nC15', 'n-C15'])
good_syns['544-76-3']['synonyms'].extend(['C16', 'nC16', 'n-C16'])
good_syns['629-78-7']['synonyms'].extend(['C17', 'nC17', 'n-C17'])
good_syns['593-45-3']['synonyms'].extend(['C18', 'nC18', 'n-C18'])
good_syns['629-92-5']['synonyms'].extend(['C19', 'nC19', 'n-C19'])
good_syns['112-95-8']['synonyms'].extend(['C20', 'nC20', 'n-C20'])
good_syns['629-94-7']['synonyms'].extend(['C21', 'nC21', 'n-C21'])
good_syns['629-97-0']['synonyms'].extend(['C22', 'nC22', 'n-C22'])
good_syns['638-67-5']['synonyms'].extend(['C23', 'nC23', 'n-C23'])
good_syns['646-31-1']['synonyms'].extend(['C24', 'nC24', 'n-C24'])
good_syns['629-99-2']['synonyms'].extend(['C25', 'nC25', 'n-C25'])
good_syns['630-01-3']['synonyms'].extend(['C26', 'nC26', 'n-C26'])
good_syns['593-49-7']['synonyms'].extend(['C27', 'nC27', 'n-C27'])
good_syns['630-02-4']['synonyms'].extend(['C28', 'nC28', 'n-C28'])
good_syns['630-03-5']['synonyms'].extend(['C29', 'nC29', 'n-C29'])
good_syns['638-68-6']['synonyms'].extend(['C30', 'nC30', 'n-C30'])
good_syns['630-04-6']['synonyms'].extend(['C31', 'nC31', 'n-C31'])
good_syns['544-85-4']['synonyms'].extend(['C32', 'nC32', 'n-C32'])
good_syns['630-05-7']['synonyms'].extend(['C33', 'nC33', 'n-C33'])
good_syns['14167-59-0']['synonyms'].extend(['C34', 'nC34', 'n-C34'])
good_syns['630-07-9']['synonyms'].extend(['C35', 'nC35', 'n-C35'])
good_syns['630-06-8']['synonyms'].extend(['C36', 'nC36', 'n-C36'])
good_syns['7194-84-5']['synonyms'].extend(['C37', 'nC37', 'n-C37'])
good_syns['7194-85-6']['synonyms'].extend(['C38', 'nC38', 'n-C38'])
good_syns['7194-86-7']['synonyms'].extend(['C39', 'nC39', 'n-C39'])
good_syns['4181-95-7']['synonyms'].extend(['C40', 'nC40', 'n-C40'])
good_syns['7194-87-8']['synonyms'].extend(['C41', 'nC41', 'n-C41'])
good_syns['7098-20-6']['synonyms'].extend(['C42', 'nC42', 'n-C42'])
good_syns['7098-21-7']['synonyms'].extend(['C43', 'nC43', 'n-C43'])
good_syns['7098-22-8']['synonyms'].extend(['C44', 'nC44', 'n-C44'])
good_syns['7098-23-9']['synonyms'].extend(['C45', 'nC45', 'n-C45'])
good_syns['7098-24-0']['synonyms'].extend(['C46', 'nC46', 'n-C46'])
good_syns['7098-25-1']['synonyms'].extend(['C47', 'nC47', 'n-C47'])
good_syns['7098-26-2']['synonyms'].extend(['C48', 'nC48', 'n-C48'])
good_syns['7098-27-3']['synonyms'].extend(['C49', 'nC49', 'n-C49'])
good_syns['6596-40-3']['synonyms'].extend(['C50', 'nC50', 'n-C50'])
good_syns['78-78-4']['synonyms'].extend(['isopentane', 'i-pentane', 'ipentane', 'iC5', 'i-C5'])
good_syns['75-28-5']['synonyms'].extend(['isobutane', 'i-butane', 'ibutane', 'iC4', 'i-C4'])
good_syns['107-83-5']['synonyms'].extend(['isohexane', 'i-hexane', 'ihexane', 'iC6', 'i-C6'])
good_syns['9005-53-2']['synonyms'].extend(['lignin'])

custom_chemicals = {}


custom_chemicals['15922-78-8'] = {'pubchem': 44782648}

custom_chemicals['12075-68-2'] = {'pubchem': 15977862}
custom_chemicals['65460-51-7'] = {'formula': "C64H138O10P2Ti", 'pubchem': 105197}


custom_chemicals['74-86-2'] = {'iupac_name': 'ethyne'}

# refprop names
custom_chemicals['112-63-0'] = {'synonyms': ['methyl (Z,Z)-9,12-octadecadienoate']}
custom_chemicals['301-00-8'] = {'synonyms': ['methyl (Z,Z,Z)-9,12,15-octadecatrienoate']}
custom_chemicals['616-38-6'] = {'synonyms': ['dimethyl ester carbonic acid']}
custom_chemicals['1885-48-9'] = {'synonyms': ['2,2,2-trifluoroethyl-difluoromethyl-ether']}
custom_chemicals['22410-44-2'] = {'synonyms': ['methyl-pentafluoroethyl-ether']}

custom_chemicals['75-69-4'] = {'synonyms': ['R-11']}
custom_chemicals['76-19-7'] = {'synonyms': ['R-218']}
custom_chemicals['359-35-3'] = {'synonyms': ['R-134']}
custom_chemicals['420-46-2'] = {'synonyms': ['R-143a']}
custom_chemicals['430-66-0'] = {'synonyms': ['R-143']}
custom_chemicals['431-89-0'] = {'synonyms': ['R-227ea']}


custom_chemicals['306-83-2'] = {'synonyms': ['R-123']}
custom_chemicals['354-23-4'] = {'synonyms': ['R-123a']}

# while adding Poling full table of names checks
custom_chemicals['92-52-4'] = {'synonyms': ['1,1`-biphenyl']}
custom_chemicals['355-25-9'] = {'synonyms': ['decafluoro-2-methylpropane']}
custom_chemicals['2207-04-7'] = {'synonyms': ['t-1,4-dimethylcyclohexane']}
custom_chemicals['75-31-0'] = {'synonyms': ['2-propanamine', 'methyl ethyl amine']}
custom_chemicals['95-65-8'] = {'synonyms': ['3,4-dimethylphenol', '3,4 xylenol']}
custom_chemicals['98-82-8'] = {'synonyms': ['1-methylethylbenzene', 'cumene']}
custom_chemicals['583-61-9'] = {'synonyms': ['2,3 lutidine']}
custom_chemicals['108-68-9'] = {'synonyms': ['3,5-dimethylphenol', '3,5 xylenol']}

# this is a trimer, its synonyms point to the actual chemical and replace it :(
custom_chemicals['26472-00-4'] = {'synonyms': ['3a,4,7,7a-tetrahydrodimethyl-4,7-methano-1h-indene',
'4,7-methano-1h-indene, 3a,4,7,7a-tetrahydrodimethyl-',
'4,7-methanoindene, 3a,4,7,7a-tetrahydrodimethyl-',
'methylcyclopentadienedimer'], 'hardcoded_synonyms': True}

for CAS, d in custom_chemicals.items():
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
