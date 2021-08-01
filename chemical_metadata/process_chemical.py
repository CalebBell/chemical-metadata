# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2021 Caleb Bell
<Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

This module contains a complete periodic table, routines for working with
chemical formulas, computing molecular weight, computing mass fractions and
atom fractions, and assorted other tasks.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemical_metadata/>`_.

"""

from natsort import natsorted, realsorted

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.inchi import MolToInchi, MolFromInchi, InchiToInchiKey

from chemical_metadata import common_chemistry_data, find_pubchem_from_ids
from chemicals.identifiers import serialize_formula
from chemicals.elements import molecular_weight, nested_formula_parser

def write_database(data, output):
    '''Write a list of dictionaries of data to a file. The expected keys in the
    list are shown below.
    '''
    f = open(output, 'w')
    lines = generate_database_lines(data)
    f.writelines(lines)
    f.close()
    
def generate_database_lines(data):
    lines = []
    for dat in data:
        s = '%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' %(dat['cid'], dat['CAS'], dat['formula'], round(dat['MW'],9),
                                                 dat['smiles'], dat['inchi'], dat['inchikey'], dat['name'])
        s += '\t'.join(dat['synonyms'])
        s += '\n'
        lines.append(s)
    lines = realsorted(lines)
    return lines
    
def process(init_data, use_cache=True):
    '''
    
    Examples
    --------
    
    >>> res = process({'CAS': '10170-69-1', 'synonyms': ['14267-36-8', 'NSC 22319'], 'name': 'Manganese, decacarbonyldi-, (Mn-Mn)'})
    >>> res['inchi'], res['smiles'], res['cid'], res['CAS']
    ('InChI=1S/10CO.2Mn/c10*1-2;;', '[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[Mn].[Mn]', 517769, '10170-69-1')
    '''
    # print(locals())
    init_data = init_data.copy()
    cc = cc_CAS = cc_name = cc_inchi = cc_inchikey = cc_smiles = cc_synonyms = cc_deprecated_CASs = None
    if 'CAS' in init_data:
        try:
            cc = common_chemistry_data(init_data['CAS'])
            cc_CAS, cc_name, cc_inchi, cc_inchikey, cc_smiles, cc_synonyms, cc_deprecated_CASs = cc
        except ValueError:
            # Compund is not in common chemistry; this is OK
            pass
    
    cid = iupac_name = p_MW = p_inchi = p_inchikey = p_smiles = p_formula = p_synonyms = None


    if init_data.get('mol', None) is not None:
        # If not in common chemistry or no InChi there, but if we have a mol file, get the inchi and inchikey for the
        # pubchem lookup
        mol = Chem.MolFromMolFile(init_data['mol'])
        if mol is not None:
            init_data['inchi'] = MolToInchi(mol)
            init_data['inchikey'] = InchiToInchiKey(init_data['inchi'])

    can_search_pubchem = (init_data.get('pubchem') is not None
                          or init_data.get('CASRN', cc_CAS) is not None
                          or init_data.get('inchi', cc_inchi) is not None
                          or init_data.get('inchikey', cc_inchikey) is not None
                          or init_data.get('smiles', cc_smiles) is not None)
    
    if can_search_pubchem:
        try:
            p = find_pubchem_from_ids(pubchem=init_data.get('pubchem'), 
                                      CASRN=init_data.get('CASRN', cc_CAS),
                                      inchi=init_data.get('inchi', cc_inchi),
                                      inchikey=init_data.get('inchikey', cc_inchikey),
                                      smiles=init_data.get('smiles', cc_smiles),
                                      use_cache=use_cache)
        except Exception as e:
            p = None
            print(e, 'exception')
        if p is not None:
            cid, iupac_name, p_MW, p_inchi, p_inchikey, p_smiles, p_formula, p_synonyms = p
    # print(locals())
    mol = None
    # Be aware some smiles descriptions are wrong
    # Start with user overridding
    if 'mol' in init_data:
        mol = Chem.MolFromMolFile(init_data['mol'])
    if mol is None and 'smiles' in init_data:
        mol = Chem.MolFromSmiles(init_data['smiles'])
    if mol is None and 'inchi' in init_data:
        mol = MolFromInchi(init_data['inchi']) if init_data['inchi'].startswith("InChI=1S/") else MolFromInchi("InChI=1S/" + init_data['inchi'])
    # Trust common chemistry next
    if mol is None and cc_smiles is not None:
        mol = Chem.MolFromSmiles(cc_smiles)
    if mol is None and cc_inchi is not None:
        mol = MolFromInchi(cc_inchi) if cc_inchi.startswith("InChI=1S/") else MolFromInchi("InChI=1S/" + cc_inchi)
    # Did we pull up the structure from pubchem??
    if mol is None and p_smiles is not None:
        mol = Chem.MolFromSmiles(p_smiles)
    if mol is None and p_inchi is not None:
        mol = MolFromInchi(p_inchi) if p_inchi.startswith("InChI=1S/") else MolFromInchi("InChI=1S/" + p_inchi)
    if mol is None:
        raise ValueError("No structure found")

    smiles = Chem.MolToSmiles(mol, True)
    inchi = MolToInchi(mol)
    inchikey = InchiToInchiKey(inchi)
    #MW = Descriptors.ExactMolWt(mol)
    formula = CalcMolFormula(mol, True, True)
    formula = serialize_formula(formula)
    MW = molecular_weight(nested_formula_parser(formula))

    # print(inchi, cc_inchi, p_inchi)
    # print(inchikey, cc_inchikey, p_inchikey)
    # print(smiles, cc_smiles, p_smiles)
    
    # output values
    if 'pubchem' in init_data:
        cid = init_data['pubchem']
    elif cid is None:
        cid = -1
    
    if cc_CAS is not None:
        CAS = cc_CAS
    elif 'CAS' in init_data:
        CAS = init_data['CAS']
    else:
        raise ValueError("CAS culd not be found")
        
    if 'formula' in init_data:
        # Override rdkit
        formula = init_data['formula']

    if 'MW' in init_data:
        # Override rdkit
        MW = init_data['MW']
        
    if 'smiles' in init_data:
        smiles = init_data['smiles']
    if 'inchi' in init_data:
        inchi = init_data['inchi']
    if 'inchikey' in init_data:
        inchikey = init_data['inchikey']
        
    if inchikey == '*' or smiles == '*' or inchi == '*':
        raise ValueError("Failure in rdkit")
        
    # Do we have a name specified in the settings?
    if 'name' in init_data:
        name = init_data['name']
    elif cc_name is not None:
        name = cc_name
    elif iupac_name is not None:
        name = iupac_name
    else:
        raise ValueError("There is no name for this compound")
        
    synonyms = []
    if cc_synonyms is not None:
        synonyms += cc_synonyms
    if cc_deprecated_CASs is not None:
        synonyms += cc_deprecated_CASs
    if p_synonyms is not None:
        synonyms += p_synonyms
    if 'synonyms' in init_data:
        synonyms += init_data['synonyms']
    synonyms = list(set(synonyms))
    if name in synonyms:
        synonyms.remove(name)
    if synonyms:
        def key_sort_str(s):
            return len(s), s.lower()
        synonyms = sorted(synonyms, key=key_sort_str)
        # synonyms = natsorted(synonyms)
    # synonyms = []
        
    
    return {'cid': cid, 'CAS': CAS, 'formula': formula, 'MW': MW, 'smiles': smiles,
            'inchi': inchi, 'inchikey': inchikey, 'name': name, 'synonyms': synonyms}
