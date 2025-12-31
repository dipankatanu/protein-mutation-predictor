#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from Bio.PDB import NeighborSearch

def calculate_clash_risk(structure, chain_id, residue_id, new_aa_type, volume_scale):
    """
    Calculates if a mutation will likely cause a steric clash based on volume.
    """
    # 1. Get the residue object
    residue = structure[0][chain_id][residue_id]
    original_aa = residue.get_resname()
    
    # 2. Get the "local neighborhood" (atoms within 5 Angstroms)
    # We use the CA (Alpha Carbon) of our residue as the center point
    center_coord = residue['CA'].get_coord()
    atoms_in_structure = list(structure.get_atoms())
    ns = NeighborSearch(atoms_in_structure)
    
    # Find all atoms near this residue (excluding the residue's own atoms)
    neighbors = ns.search(center_coord, 5.0)
    nearby_atoms = [a for a in neighbors if a.get_parent() != residue]
    
    # 3. Calculate "Crowdedness"
    # We count how many neighbor atoms are already very close (within 3.5A)
    close_neighbors = [a for a in neighbors if np.linalg.norm(a.get_coord() - center_coord) < 3.5]
    num_neighbors = len(close_neighbors)
    
    # 4. Logic: If the new AA is much bigger AND the area is already crowded
    vol_diff = volume_scale.get(new_aa_type, 0) - volume_scale.get(original_aa, 0)
    
    # Heuristic: If volume increases by >40% and neighbors > 10, it's a high clash risk
    clash_score = (vol_diff / 100) * num_neighbors 
    
    return {
        "original": original_aa,
        "new": new_aa_type,
        "vol_increase": vol_diff,
        "neighbor_count": num_neighbors,
        "risk_level": "High" if clash_score > 5 else "Low"
    }


# In[ ]:




