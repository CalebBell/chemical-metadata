#! /home/caleb/.anaconda3/bin/python
import sys

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from pubchempy import get_compounds, Compound


arg = sys.argv[0:]
arg.pop(0)
for f in arg:
    try:
        mol = Chem.MolFromMolFile(f)
        smi = Chem.MolToSmiles(mol, True)
        inchi_val = inchi.MolToInchi(mol)
        inchikey = inchi.InchiToInchiKey(inchi_val)
        mw = Descriptors.ExactMolWt(mol)
        formula = CalcMolFormula(mol)
        CAS = f[0:-4]
        
        try:
            pc = get_compounds(inchikey, 'inchikey')[0]
            cid = pc.cid
            iupac_name = pc.iupac_name
            names = pc.synonyms
        except:
            cid = -1
            iupac_name = ''
            names = ['']

        s = '%d\t%s\t%s\t%g\t%s\t%s\t%s\t%s\t' %(cid, CAS, formula, mw, smi, inchi_val, inchikey, iupac_name)
        s += '\t'.join(names)
        print(s)
    except:
        pass
    
    
    
#4	78-96-6	C3H9NO	75.10966	CC(CN)O	C3H9NO/c1-3(5)2-4/h3,5H,2,4H2,1H3	HXKKHQJGJAFBHI-UHFFFAOYSA-N	1-azanylpropan-2-ol	1-amino-2-propanol	
