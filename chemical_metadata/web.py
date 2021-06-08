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

cache_dir = appdirs.user_cache_dir(appname='chemical_metadata')
if not os.path.exists(cache_dir):
    os.mkdir(cache_dir)

common_chemistry_cache_dir = os.path.join(cache_dir, 'common_chemistry')
if not os.path.exists(common_chemistry_cache_dir):
    os.mkdir(common_chemistry_cache_dir)

pubchem_cache_dir = os.path.join(cache_dir, 'pubchem')
if not os.path.exists(pubchem_cache_dir):
    os.mkdir(pubchem_cache_dir)
    
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
        
    return (json_data['rn'], remove_html(json_data['name']), json_data['inchi'], json_data['inchiKey'],
            json_data['canonicalSmile'], synonyms, json_data['replacedRns'])

