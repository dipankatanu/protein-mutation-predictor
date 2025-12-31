# Protein Mutation Predictor

A computational structural biology pipeline designed to evaluate the impact of **single-point mutations** on protein stability.

This tool bridges the gap between sequence-based prediction and high-compute molecular dynamics using a **Dual-Engine approach**:

> **Geometric Heuristics + Physical Energetics**  
> to estimate mutation risk at atomic resolution.

---

## Key Features

### Dual-Engine Mutation Analysis
| Engine | Description |
|--------|--------------|
| **Physical Engine** | Computes **Lennard-Jones 12-6 potential** to quantify van der Waals repulsion and steric clashes (kcal/mol). |
| **Heuristic Engine** | Uses **Zamyatin side-chain volume** & **Kyte-Doolittle hydropathy** to detect cavity formation and polarity mismatches. |

### High-Throughput Prediction
- *In-silico saturation mutagenesis* (19 substitutions per position)
- Whole-chain structural fragility scoring
- Mutation sensitivity index per residue
- ΔVolume × ΔHydropathy **risk quadrant map**

###  Advanced Visualizations
- DMS Heatmaps (residue × mutation matrix)
- Mutation fragility profiles
- Local steric energy peaks
- 3D structural mapping (via **py3Dmol**)

###  Clinical Relevance Example
- Identifies pathogenic signature of the **HbS E6V mutation** in Hemoglobin
- Predicts increased packing strain + polarity inversion

---

##  Scientific Foundation

### **1 Steric Hindrance & van der Waals Physics**
Short-range repulsion prevents atomic overlap due to electron cloud exclusion.

The model uses the Lennard-Jones 12-6 potential:
V(r) = 4ε[(σ/r)^12 - (σ/r)^6]


- \( r \) = interatomic distance  
- \( (σ/r)^{12} \) → **Pauli Exclusion / Steric Clash**  
- \( (σ/r)^{6} \) → **van der Waals attraction**

When a mutation introduces volume into an already packed core:

> \( r \to 0 \Rightarrow V(r) \to +\infty \)  
> **→ Highly destabilizing mutation**

---

### **2️ Chemistry-Driven Unfolding**
Folding is stabilized by **the hydrophobic effect**:

- Buried → Hydrophobic preferred  
- Surface → Polar preferred

The predictor flags:

 **Hidden Killers**  
Mutations that *fit geometrically* but **fail chemically** - e.g., introducing Asp, Glu, Lys, or His inside a hydrophobic core → desolvation penalty & unfolding tendency.

---

## Outputs & Example Results

| Output Type | Example Insight |
|-------------|-----------------|
| **DMS Heatmap** | Identifies structurally fragile positions |
| **Energy Profile** | “Energy Peaks” correspond to steric bottlenecks |
| **ΔV × ΔH Scatter** | Classifies mutations by physical & chemical disruption |
| **3D Visualization** | Shows spatial clustering of sensitive sites |

These enable **rapid triage prior to MD, Rosetta, or FoldX**.

---

##  Disclaimer

This tool is a heuristic rapid-screening model, not a substitute for:

- Molecular Dynamics (GROMACS / AMBER)
- Rosetta or FoldX ΔΔG
-Experimental DMS or cryo-EM refinement

It is intended for:
- Prioritization
- Early triage
- Hypothesis generation


## Case Study: Clinical Validation - **Sickle Cell Anemia (HbS)**

The classic **Hemoglobin S (HbS)** mutation - **Glu6Val (E6V)** - causes Sickle Cell Disease by inducing hemoglobin polymerization under low-oxygen conditions.

Using our Dual-Engine framework:

| Engine | Prediction | Why it Matters |
|--------|-------------|-----------------|
| **Physical Engine** | **Low Steric Penalty** → tolerated | Valine is *smaller* than Glutamate → no steric clash |
| **Heuristic Engine** | **High-Risk Polarity Switch** 🚨 | Replacing a *charged, hydrophilic Glu* with a *hydrophobic Val* creates an **oily patch** |

### Interpretation

- The E6V substitution creates a **sticky hydrophobic surface**
- This patch drives **protein–protein aggregation**
- The aggregation triggers **polymer fiber formation**
- Result: **red blood cell deformation → Sickle Cell phenotype**

This matches clinical and biochemical understanding:

> The mutation does **not** destabilize the structure mechanically -  
> it destabilizes it **chemically**, enabling abnormal interactions.

###  Why This Matters

This case validates that the tool can correctly identify:

- Mutations that are **benign in terms of geometry**
- But **pathogenic due to biochemical context**

These mutations are often missed by distance-only or volume-only screens.  
Our tool **catches them.**

---

**Conclusion:**  
> The predictor reproduces the known pathological behavior of the **Sickle Cell mutation**, supporting its validity as a **mutation triage tool** before running full MD / Rosetta / FoldX pipelines.

---

