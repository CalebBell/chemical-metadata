#! /home/caleb/.anaconda3/bin/python
import sys

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from pubchempy import get_compounds, Compound
import json
from collections import Counter

pdf_data = open('Parsed scifinder metadata.json').read()
pdf_data = json.loads(pdf_data)
all_names = []
for CAS, d in pdf_data.items():
    if 'Other Names' in d:
        all_names.extend(d['Other Names'])
    if 'Name' in d:
        all_names.append(d['Name'])
    if 'Deleted CAS' in d:
        all_names.extend(d['Deleted CAS'])
    if 'Alternate CAS' in d:
        all_names.extend(d['Alternate CAS'])
    all_names.append(CAS)
dup_names =  [item for item, count in Counter(all_names).items() if count > 1]
all_names = set(all_names)

# TODO finish this
# Big database = hard to work on
# Small database = good to work on
failed_mol = set(['11062-77-4',# No charge
                  ])

arg = sys.argv[0:]
arg.pop(0)
for f in arg:
#    try:
    mol = Chem.MolFromMolFile(f)
    if mol is None:
        print('BAILING ON %s'%f)
        continue
    smi = Chem.MolToSmiles(mol, True)
    inchi_val = inchi.MolToInchi(mol)
    inchikey = inchi.InchiToInchiKey(inchi_val)
    mw = Descriptors.ExactMolWt(mol)
    formula = CalcMolFormula(mol)
    CAS = f.split('/')[1] if '/' in f else f
    CAS = CAS.split('.')[0]
    
    try:
        pc = get_compounds(inchikey, 'inchikey')[0]
        cid = pc.cid
        iupac_name = pc.iupac_name
        names = pc.synonyms
    except:
        cid = -1
        iupac_name = ''
        names = ['']
    
    other_CAS = []
#    print(CAS, CAS in pdf_data)
    if CAS in pdf_data:
#        print('hi')
        d = pdf_data[CAS]
        name = d['Name']
        if 'Other Names' in d:
            syns = d['Other Names']
        else:
            syns = []
        if not iupac_name:
            iupac_name = name
        else:
            syns.insert(0, name)
        if 'Deleted CAS' in d:
            other_CAS.extend(d['Deleted CAS'])
        if 'Alternate CAS' in d:
            other_CAS.extend(d['Alternate CAS'])
            
        syns = [i for i in syns if i not in dup_names]
        names =  syns + [i for i in names if i not in all_names] + other_CAS
    

    s = '%d\t%s\t%s\t%g\t%s\t%s\t%s\t%s\t' %(cid, CAS, formula, mw, smi, inchi_val, inchikey, iupac_name)
    s += '\t'.join(names)
    print(s)
#    except:
#        pass
    
    
#4	78-96-6	C3H9NO	75.10966	CC(CN)O	C3H9NO/c1-3(5)2-4/h3,5H,2,4H2,1H3	HXKKHQJGJAFBHI-UHFFFAOYSA-N	1-azanylpropan-2-ol	1-amino-2-propanol	
