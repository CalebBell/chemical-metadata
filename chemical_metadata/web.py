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
import json
import appdirs
import requests
import os
import lxml.html
from pubchempy import get_compounds, Compound
import hashlib
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi, MolFromInchi, InchiToInchiKey

cache_dir = appdirs.user_cache_dir(appname='chemical_metadata')
if not os.path.exists(cache_dir):
    os.mkdir(cache_dir)

common_chemistry_cache_dir = os.path.join(cache_dir, 'common_chemistry')
if not os.path.exists(common_chemistry_cache_dir):
    os.mkdir(common_chemistry_cache_dir)

pubchem_cache_dir = os.path.join(cache_dir, 'pubchem')
if not os.path.exists(pubchem_cache_dir):
    os.mkdir(pubchem_cache_dir)

def deterministic_hash(s):
    return hashlib.sha224(s.encode('utf-8')).hexdigest()

    
base_url = '''https://rboq1qukh0.execute-api.us-east-2.amazonaws.com/default/detail?cas_rn='''

def remove_html(text):
    '''
    Examples
    --------
    
    >>> remove_html('Boric acid (H<sub>3</sub>BO<sub>3</sub>), trimethyl ester')
    'Boric acid (H3BO3), trimethyl ester'
    
    '''
    if '<' in text:
        return lxml.html.fromstring(text).text_content()
    return text

def common_chemistry_data(CASRN):
    '''Cached lookup of a chemical's metadata on Common Chemistry.
    
    Parameters
    ----------
    CASRN : str
        CAS number, [-]

    Returns
    -------
    CASRN : str
        CAS number from lookup (MAY BE DIFFERENT THAN INPUT IF INPUT IS
        DEPRECATED!), [-]
    name : str
        Name as given in Common Chemistry, [-]
    InChI : str
        InChI identification string as given in Common Chemistry (there can be multiple
        valid InChI strings for a compound), [-]
    InChI_key : str
        InChI key identification string (meant to be unique to a compound), [-]        
    smiles : str
        SMILES identification string, [-]
    synonyms : list[str]
        List of synonyms of the compound, [-]
    deprecated_CASs : list[str]
        No longer in use CAS numbers which should be resolved to `CASRN`, [-]
        
        
    Examples
    --------
    >>> common_chemistry_data("7732-18-5")[0]
    '7732-18-5'
    >>> common_chemistry_data("7732-18-5")[1]
    'Water'
    >>> common_chemistry_data("7732-18-5")[6]
    ['558440-22-5', '558440-53-2', '652133-48-7', '1202864-49-0', '1371582-34-1']
    >>> common_chemistry_data("7732-18-5")[4]
    'O'
    >>> common_chemistry_data("7732-18-5")[3]
    'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
    
    >>> common_chemistry_data("2471-08-1")
    Traceback (most recent call last):
    ValueError: Compound is not in common chemistry
    
    >>> common_chemistry_data("56-87-1")[5]
    ['L-Lysine', 'Lysine, L-', 'Hexanoic acid, 2,6-diamino-, (S)-', 'Lysine acid', 'Lysine', '2,6-Diaminohexanoic acid', '(S)-α,ε-Diaminocaproic acid', 'L-(+)-Lysine', 'Aminutrin', '(S)-2,6-Diaminohexanoic acid', 'α-Lysine', '(S)-Lysine', '(+)-S-Lysine', 'L-Norleucine, 6-amino-', 'L-2,6-Diaminocaproic acid', 'h-Lys-oh', 'L-Lys', 'Malandil', 'Biolys', 'LLB 50', '38: PN: WO2020027640 SEQID: 101 claimed protein', '38: PN: KR20200014684 SEQID: 101 claimed sequence']
    
    >>> common_chemistry_data("108-39-4")[1]
    'm-Cresol'
    
    >>> common_chemistry_data("53850-36-5")[1]
    'Rutherfordium'
    '''
    cache_loc = os.path.join(common_chemistry_cache_dir, CASRN)
    cached = False
    if os.path.exists(cache_loc):
        f = open(cache_loc, 'r')
        json_data = json.loads(f.read())
        f.close()
    else:
        r = requests.get(base_url + CASRN)
        json_data = r.json()
        f = open(cache_loc, 'w')
        json.dump(json_data, f)
        f.close()
        
    if 'message' in json_data and 'Detail not found' == json_data['message']:
        raise ValueError("Compound is not in common chemistry")
#     print(json_data)

    synonyms = [remove_html(k) for k in json_data['synonyms']]
    
    inchi = json_data['inchi'] if json_data['inchi']  else None
    inchiKey = json_data['inchiKey'] if json_data['inchiKey']  else None
    canonicalSmile = json_data['canonicalSmile'] if json_data['canonicalSmile']  else None
        
    return (json_data['rn'], remove_html(json_data['name']), inchi, inchiKey,
            canonicalSmile, synonyms, json_data['replacedRns'])

def find_pubchem_from_ids(pubchem=None, CASRN=None, inchi=None, inchikey=None,
                          smiles=None, use_cache=True):
    '''Cached query of pubchem database, based on one of many identifiers.
    
    Parameters
    ----------
    pubchem : int, optional
        PubChem ID; prefered lookup, [-]
    CASRN : str, optional
        CAS number, [-]
    inchi : str, optional
        InChI identification string as given in Common Chemistry (there can be multiple
        valid InChI strings for a compound), [-]
    inchikey : str, optional
        InChI key identification string (meant to be unique to a compound), [-]        
    smiles : str, optional
        SMILES identification string, [-]
    use_cache : bool, optional
        Whether or not to use the cache, [-]
    
    Returns
    -------
    cid : intoxidane
        PubChem ID, [-]
    iupac_name : str
        IUPAC name as given in pubchem, [-]
    MW : float
        Molecular weight, [g/mol]
    InChI : str
        InChI identification string as given in Common Chemistry (there can be multiple
        valid InChI strings for a compound), [-]
    InChI_key : str
        InChI key identification string (meant to be unique to a compound), [-]        
    smiles : str
        SMILES identification string, [-]
    formula : str
        Formula, [-]
    synonyms : list[str]
        List of synonyms of the compound, [-]
        
    Examples
    --------
    
    >>> find_pubchem_from_ids(pubchem=962)[0]
    962
    >>> find_pubchem_from_ids(pubchem=962)[1]
    'oxidane'
    >>> find_pubchem_from_ids(pubchem=962)[2]
    18.015
    >>> find_pubchem_from_ids(pubchem=962)[3]
    'InChI=1S/H2O/h1H2'
    >>> find_pubchem_from_ids(pubchem=962)[4]
    'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
    >>> find_pubchem_from_ids(pubchem=962)[5]
    'O'
    >>> find_pubchem_from_ids(pubchem=962)[6]
    'H2O'
    >>> len(find_pubchem_from_ids(pubchem=962)[7]) > 100
    True
    
    >>> find_pubchem_from_ids(CASRN="53850-36-5")[0]
    56951715
    
    >>> find_pubchem_from_ids(CASRN="54084-70-7") # Nihonium is missing
    [None, None, None, None, None, None, None, None]
    
    try to use rdkit here to check the correct inchikey is found.
    
    >>> find_pubchem_from_ids(inchi='InChI=1S/Cl', inchikey="ZAMOUSCENKQFHK-UHFFFAOYSA-N")[0]
    5360523
    >>> find_pubchem_from_ids(inchi='InChI=1S/H2O/h1H2', inchikey="XLYOFNOQVPJJNP-UHFFFAOYSA-N")[0]
    962
    
    >>> find_pubchem_from_ids(inchi='InChI=1S/I2/c1-2')[0:5]
    [807, 'molecular iodine', 253.8089, 'InChI=1S/I2/c1-2', 'PNDPGZBMCMUPRI-UHFFFAOYSA-N']
    >>> find_pubchem_from_ids(inchi='InChI=1S/H2/h1H')[0:5]
    [783, 'molecular hydrogen', 2.016, 'InChI=1S/H2/h1H', 'UFHFLCQGNIYNRP-UHFFFAOYSA-N']
    '''
    abort = False
    key = (pubchem, CASRN, inchi, inchikey, smiles)
    hash_key = deterministic_hash(str(key))
    key_file = os.path.join(pubchem_cache_dir, hash_key)
    if os.path.exists(key_file) and use_cache:
        f = open(key_file, 'r')
        json_data = json.loads(f.read())
        f.close()
        return json_data
    
    if pubchem is not None:
        compound = Compound.from_cid(pubchem)
        cid = compound.cid
    else:
        if inchikey is not None:
            # Dup for chlorine atomic here
            # find_pubchem_from_ids(inchikey='ZAMOUSCENKQFHK-UHFFFAOYSA-N')[0]
             compounds = get_compounds(inchikey, 'inchikey')
        elif inchi is not None:
            # chlorine search "InChI=1S/Cl" finds HCl
            compounds = get_compounds(inchi, 'inchi')
        elif smiles is not None:
             compounds = get_compounds(smiles, 'smiles')
        elif CASRN is not None:
            compounds = get_compounds(CASRN, 'name')
        # maybe sort by ID in the future
        if not compounds:
            abort = True
            cid = None
        if not abort:
            compound = compounds[0]
            cid = compound.cid

    if cid is None:
        abort = True
    if abort:
        cid, iupac_name, mw, inchi_val, inchikey, smi, formula, names = [None]*8
    else:
        iupac_name = compound.iupac_name
        mw = float(compound.molecular_weight)
        smi = compound.canonical_smiles
        inchi_val = compound.inchi
        inchikey = compound.inchikey
        formula = compound.molecular_formula
        names = compound.synonyms
    ans = (cid, iupac_name, mw, inchi_val, inchikey, smi, formula, names)
    
    f = open(key_file, 'w')
    json.dump(ans, f)
    f.close()
    return ans
