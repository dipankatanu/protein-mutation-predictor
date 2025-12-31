#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# src/constants.py

# Hydropathy index from Kyte and Doolittle (1982)
# Positive = Hydrophobic (hates water, usually inside the protein)
# Negative = Hydrophilic (loves water, usually on the surface)
HYDROPATHY_SCALE = {
    'ALA': 1.8,  'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5,
    'CYS': 2.5,  'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4,
    'HIS': -3.2, 'ILE': 4.5,  'LEU': 3.8,  'LYS': -3.9,
    'MET': 1.9,  'PHE': 2.8,  'PRO': -1.6, 'SER': -0.8,
    'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}

# Average Residue Volume in cubic Angstroms (Å³) from Zamyatin (1972)
# This helps us predict "Steric Clashes" (if a mutation is too big for its spot)
VOLUME_SCALE = {
    'ALA': 88.6,  'ARG': 173.4, 'ASN': 114.1, 'ASP': 111.1,
    'CYS': 108.5, 'GLN': 143.8, 'GLU': 138.4, 'GLY': 60.1,
    'HIS': 153.2, 'ILE': 166.7, 'LEU': 166.7, 'LYS': 168.6,
    'MET': 162.9, 'PHE': 189.9, 'PRO': 112.7, 'SER': 89.0,
    'THR': 116.1, 'TRP': 227.8, 'TYR': 193.6, 'VAL': 140.0
}

