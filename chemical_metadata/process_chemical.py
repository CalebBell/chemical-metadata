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

from chemical_metadata import common_chemistry_data, find_pubchem_from_ids

def write_database(data, output):
    '''Write a list of dictionaries of data to a file. The expected keys in the
    list are shown below.
    '''
    f = open(output, 'w')
    for dat in data:
        s = '%d\t%s\t%s\t%g\t%s\t%s\t%s\t%s\t' %(dat['cid'], dat['CAS'], dat['formula'], dat['MW'],
                                                 dat['smiles'], dat['inchi'], dat['inchikey'], dat['name'])
        s += '\t'.join(dat['synonyms'])
        f.write(s)
    f.close()
