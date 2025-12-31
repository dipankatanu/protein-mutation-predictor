#!/usr/bin/env python
# coding: utf-8

# # Setup

# In[2]:


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import py3Dmol
from Bio.PDB import PDBList, PDBParser, NeighborSearch

# Define workspace
base_path = os.path.dirname(os.path.abspath(__file__))

folders = ['data', 'src', 'notebooks', 'output']
for folder in folders:
    os.makedirs(os.path.join(base_path, folder), exist_ok=True)

os.chdir(base_path)

def fetch_pdb(pdb_id):
    """Downloads and renames PDB files for local analysis."""
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir='data', file_format='pdb')
    
    raw_file = os.path.join('data', f"pdb{pdb_id.lower()}.ent")
    clean_file = os.path.join('data', f"{pdb_id}.pdb")
    
    if os.path.exists(raw_file):
        os.rename(raw_file, clean_file)
    return clean_file

# Loading Hemoglobin (1A3N) as our primary model
pdb_filename = fetch_pdb("1A3N")
parser = PDBParser(QUIET=True)
structure = parser.get_structure("1A3N", pdb_filename)


# # Reference Scales (The "Physics" constants)
# ## We rely on two fundamental biophysical scales:
# ## Zamyatin Volume (1972): To measure the $Å^3$ space occupied by amino acids.
# ## Kyte-Doolittle Hydropathy (1982): To measure water-affinity.

# In[3]:


HYDROPATHY_SCALE = {
    'ALA': 1.8,  'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,  
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,  
    'LEU': 3.8,  'LYS': -3.9, 'MET': 1.9,  'PHE': 2.8,  'PRO': -1.6, 
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}

VOLUME_SCALE = {
    'ALA': 88.6,  'ARG': 173.4, 'ASN': 114.1, 'ASP': 111.1, 'CYS': 108.5, 
    'GLN': 143.8, 'GLU': 138.4, 'GLY': 60.1,  'HIS': 153.2, 'ILE': 166.7, 
    'LEU': 166.7, 'LYS': 168.6, 'MET': 162.9, 'PHE': 189.9, 'PRO': 112.7, 
    'SER': 89.0,  'THR': 116.1, 'TRP': 227.8, 'TYR': 193.6, 'VAL': 140.0
}


# # The Analytical Engines
# ## A. The Heuristic Engine (Geometry & Chemistry)
# ### This engine looks for "Rule Violations." For example, trying to fit a bulky amino acid into a crowded core, or swapping an oily (hydrophobic) residue for a charged (hydrophilic) one.

# In[4]:


def calculate_full_impact(structure, chain_id, residue_id, new_aa_type):
    residue = structure[0][chain_id][residue_id]
    original_aa = residue.get_resname()
    
    # Calculate geometric/volume shift
    vol_diff = VOLUME_SCALE.get(new_aa_type, 0) - VOLUME_SCALE.get(original_aa, 0)
    
    # Calculate chemical polarity shift
    hydro_diff = HYDROPATHY_SCALE.get(new_aa_type, 0) - HYDROPATHY_SCALE.get(original_aa, 0)
    
    # Check neighbors (Atomic Packing Density)
    center_coord = residue['CA'].get_coord()
    ns = NeighborSearch(list(structure.get_atoms()))
    neighbors = ns.search(center_coord, 5.0)
    neighbor_count = len([a for a in neighbors if a.get_parent() != residue])
    
    chem_risk = "Low"
    if neighbor_count > 15 and abs(hydro_diff) > 3.0:
        chem_risk = "High (Polarity Switch in Core)"
        
    return {
        "original": original_aa,
        "vol_impact": vol_diff,
        "hydro_impact": hydro_diff,
        "chem_risk": chem_risk,
        "neighbors": neighbor_count
    }


# ## B. The Physical Engine (Lennard-Jones Energy)
# ## To simulate the actual forces of atomic repulsion, we implement a Lennard-Jones 12-6 potential. This treats atoms not as hard blocks, but as soft spheres with energy fields.

# In[7]:


def calculate_lj_energy(distance, sigma=3.5, epsilon=0.1):
    if distance < 1.0: return 100.0  # Avoid asymptotic energy at overlaps
    ratio = sigma / distance
    return 4 * epsilon * ( (ratio**12) - (ratio**6) )

def advanced_mutation_scan(structure, chain_id, res_num, new_aa):
    residue = structure[0][chain_id][res_num]
    center_atom = residue['CA']
    ns = NeighborSearch(list(structure.get_atoms()))
    neighbors = ns.search(center_atom.coord, 5.0)
    
    total_energy = sum([calculate_lj_energy(np.linalg.norm(a.coord - center_atom.coord)) 
                        for a in neighbors if a.get_parent() != residue])

    # Add a penalty for buried hydrophilic residues
    if len(neighbors) > 15 and HYDROPATHY_SCALE[new_aa] < -2.0:
        total_energy += 15.0 
            
    return total_energy


# # 4. Comprehensive Analysis & Case Study
# ## By merging these engines, we can perform a Deep Mutational Scan. We will validate this by looking at Sickle Cell Anemia (HbS): a Glu6Val mutation.

# In[10]:


def analyze_mutation_comprehensive(structure, chain_id, res_id, new_aa):
    """
    Combines Heuristic (High/Low) and Physical (LJ Energy) analysis.
    Note: It uses the global VOLUME_SCALE and HYDROPATHY_SCALE.
    """
    # 1. Run the Physics Calculation
    energy_score = advanced_mutation_scan(structure, chain_id, res_id, new_aa)
    
    # 2. Run the Heuristic Calculation (Now correctly passing 4 arguments)
    heuristic = calculate_full_impact(structure, chain_id, res_id, new_aa)
    
    # 3. Determine overall risk
    risk = "High" if (energy_score > 10.0 or heuristic['chem_risk'] != "Low") else "Low"
    
    return {
        'Residue': f"{heuristic['original']}{res_id}",
        'LJ_Energy': round(energy_score, 2),
        'Heuristic_Risk': risk,
        'Chem_Risk': heuristic['chem_risk']
    }
# Running Clinical Validation (HbS)
hb_result = analyze_mutation_comprehensive(structure, 'B', 6, 'VAL')
print(f"Validation on Sickle Cell Mutation: {hb_result}")


# In[14]:


# Generate Heatmap Data
heatmap_data = []
for residue in structure[0]['A']:
    if residue.get_resname() in VOLUME_SCALE:
        res_id = residue.get_id()[1]
        score = advanced_mutation_scan(structure, 'A', res_id, 'TRP')
        heatmap_data.append({'Pos': res_id, 'Energy': score})

df = pd.DataFrame(heatmap_data)

# Visualization
plt.figure(figsize=(15, 3))
plt.plot(df['Pos'], df['Energy'], color='darkred', lw=1)
plt.fill_between(df['Pos'], df['Energy'], color='salmon', alpha=0.3)
plt.title("Physical Sensitivity Scan (1A3N - Chain A)")
plt.ylabel("Repulsion Energy (kcal/mol)")
plt.xlabel("Residue Position")
plt.savefig("output/physical_sensitivity_scan_1A3N_chainA.png", dpi=300, bbox_inches='tight')

plt.show()

print("Saved as: physical_sensitivity_scan_1A3N_chainA.png")


# In[15]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming VOLUME_SCALE and HYDROPATHY_SCALE are already imported/defined
amino_acids = sorted(list(VOLUME_SCALE.keys()))

# --- 1. THE SIZE-CHANGE MATRIX ---
vol_matrix = np.zeros((20, 20))
for i, start_aa in enumerate(amino_acids):
    for j, end_aa in enumerate(amino_acids):
        vol_matrix[i, j] = VOLUME_SCALE[end_aa] - VOLUME_SCALE[start_aa]

plt.figure(figsize=(10, 8))
sns.heatmap(vol_matrix, annot=False, cmap='RdBu_r', center=0,
            xticklabels=amino_acids, yticklabels=amino_acids)
plt.title("The 'Size-Change' Matrix: $\Delta$ Volume ($Å^3$)")
plt.xlabel("Mutant Amino Acid")
plt.ylabel("Wild-Type Amino Acid")
plt.savefig("output/delta_volume_matrix.png", dpi=300)

# --- 2. THE CHEMISTRY VS PHYSICS SCATTER ---
delta_v, delta_h, labels = [], [], []
for start_aa in amino_acids:
    for end_aa in amino_acids:
        if start_aa != end_aa:
            delta_v.append(VOLUME_SCALE[end_aa] - VOLUME_SCALE[start_aa])
            delta_h.append(HYDROPATHY_SCALE[end_aa] - HYDROPATHY_SCALE[start_aa])
            labels.append(f"{start_aa}->{end_aa}")

plt.figure(figsize=(10, 7))
plt.scatter(delta_v, delta_h, alpha=0.5, c=delta_h, cmap='coolwarm')
plt.axhline(0, color='grey', linestyle='--', alpha=0.3)
plt.axvline(0, color='grey', linestyle='--', alpha=0.3)
plt.xlabel("Change in Volume ($\Delta V$, $Å^3$)")
plt.ylabel("Change in Hydropathy ($\Delta H$)")
plt.title("The Mutation Landscape: Chemistry vs. Physics")
plt.savefig("output/mutation_landscape_scatter.png", dpi=300)

# --- 3. THE PHYSICS REFERENCE (LJ CURVE) ---
def lj_curve(r, sigma=3.5, epsilon=0.1):
    return 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)

r_values = np.linspace(3.0, 8.0, 200)
e_values = lj_curve(r_values)

plt.figure(figsize=(8, 5))
plt.plot(r_values, e_values, 'k-', lw=2)
plt.axhline(0, color='black', lw=1)
plt.fill_between(r_values, e_values, where=(e_values > 0), color='red', alpha=0.2, label='Repulsion')
plt.fill_between(r_values, e_values, where=(e_values < 0), color='blue', alpha=0.2, label='Attraction')
plt.ylim(-0.2, 0.5)
plt.xlabel("Distance $r$ ($Å$)")
plt.ylabel("Energy $V(r)$ (kcal/mol)")
plt.legend()
plt.savefig("output/lennard_jones_reference.png", dpi=300)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Define the segment to scan (e.g., first 30 residues of Chain A)
chain_id = 'A'
amino_acids_list = sorted(list(VOLUME_SCALE.keys()))
segment_residues = [r for r in structure[0][chain_id] if r.get_resname() in VOLUME_SCALE][:30]

# 2. Build the data matrix
dms_matrix = np.zeros((len(amino_acids_list), len(segment_residues)))

for j, res in enumerate(segment_residues):
    res_num = res.get_id()[1]
    for i, aa in enumerate(amino_acids_list):
        # Calculate the energy score using your physics engine
        dms_matrix[i, j] = advanced_mutation_scan(structure, chain_id, res_num, aa)

# 3. Create DataFrame for plotting
dms_df = pd.DataFrame(
    dms_matrix, 
    index=amino_acids_list, 
    columns=[f"{r.get_resname()}{r.get_id()[1]}" for r in segment_residues]
)

# 4. Plot the Heatmap
plt.figure(figsize=(16, 8))
sns.heatmap(dms_df, cmap="YlOrRd", cbar_kws={'label': 'Stability Penalty (kcal/mol)'})
plt.title("Deep Mutational Scan (DMS): Physical Stability Landscape (Segment 1-30)")
plt.xlabel("Wild-Type Residue Position")
plt.ylabel("Mutant Amino Acid")

# Save to output
plt.savefig("output/dms_heatmap_full.png", dpi=300, bbox_inches='tight')
plt.show()


# In[ ]:




