#! /home/caleb/.anaconda3/bin/python
import sys

from pubchempy import get_compounds, Compound


arg = sys.argv[0:]
arg.pop(0)
for cid in arg:
    cid = int(cid)
    pc = Compound.from_cid(cid)
    iupac_name = pc.iupac_name
    if iupac_name == None:
        iupac_name = ''
    names = pc.synonyms
    smi = pc.isomeric_smiles
    inchi_val = pc.inchi
    inchikey = pc.inchikey


    mw = pc.molecular_weight
    formula = pc.molecular_formula
    CAS = ''
    

    s = '%d\t%s\t%s\t%g\t%s\t%s\t%s\t%s\t' %(cid, CAS, formula, mw, smi, inchi_val, inchikey, iupac_name)
    s += '\t'.join(names)
    print(s)
    
    
    
#4	78-96-6	C3H9NO	75.10966	CC(CN)O	C3H9NO/c1-3(5)2-4/h3,5H,2,4H2,1H3	HXKKHQJGJAFBHI-UHFFFAOYSA-N	1-azanylpropan-2-ol	1-amino-2-propanol	
