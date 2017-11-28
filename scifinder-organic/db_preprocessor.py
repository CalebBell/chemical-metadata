from thermo import *
from thermo.identifiers import ChemicalMetadataDB
from numpy.testing import assert_allclose
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import json

db = ChemicalMetadataDB(elements=False, main_db=('chemical identifiers example user db.tsv'), user_dbs=[])

good_syns = {}

for CAS, d in db.CAS_index.items():
    CAS = d.CASs
    if CAS in good_syns:
        pass
    else:
        good_syns[CAS] = {}
        good_syns[CAS]['synonyms'] = []


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


#D2Se = {'formula': 'D2Se', 'MW': molecular_weight(nested_formula_parser('D2Se'))}
#ammonium_hexafluorosilicate = {'pubchem': 28145}
#CsBromate = {'pubchem': 23685550}
#Br = {'pubchem': 5360770}
#NaAlO4H4 = {'pubchem': 166673}
#Na2HPO4 = {'pubchem': 24203}
#NaH2PO4 = {'pubchem': 24204}
#CuH2O4S = {'pubchem': 6536471}
#
#
#custom_ions = {'13536-95-3': D2Se, '16919-19-0': ammonium_hexafluorosilicate,
#               '13454-75-6': CsBromate, '10097-32-2': Br, '2099990000-00-0': NaAlO4H4,
#               '7558-79-4': Na2HPO4, '7558-80-7': NaH2PO4, '18939-61-2': CuH2O4S}

custom_chemicals = {}
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
