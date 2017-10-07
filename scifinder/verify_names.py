#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 18:35:28 2017

@author: caleb
"""

import json
from collections import Counter

f = open('Parsed scifinder metadata.json').read()
data = json.loads(f)
all_names = []
for CAS, d in data.items():
    if 'Other Names' in d:
        all_names.extend(d['Other Names'])
    if 'Name' in d:
        all_names.append(d['Name'])
dup_names =  [item for item, count in Counter(all_names).items() if count > 1]
#print(len(all_names), len(set(all_names)))
#print(dup_names)

