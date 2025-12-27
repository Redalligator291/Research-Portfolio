# Drug Finder - Filter and screening Pipeline

# AutoFilter

A low-cost biocomputational framework for high-throughput screening of chemical databases and identification of novel drug inhibitors.

## Overview

AutoFilter is an integrated computational drug discovery pipeline that combines machine learning, molecular docking, and molecular dynamics simulations to efficiently screen large chemical databases. The framework was developed to address the high costs and lengthy timelines of traditional drug discovery, reducing both by approximately 75%.

Originally applied to identify malaria inhibitors targeting *Plasmodium falciparum* apicoplast DNA polymerase (apPOL), AutoFilter successfully screened 2.4 million compounds from the ChEMBL database and identified 5 promising candidates currently undergoing *in vitro* validation.

## Key Features

- **Sequential Filtering Pipeline**: Applies Lipinski's Rule of Five, Veber's rules, and PAINS filters
- **High-Throughput Molecular Docking**: Automated docking using AutoDock Vina with custom Python implementation
- **ADME Property Analysis**: Pharmacokinetic profiling using SwissADME
- **ML-Based Toxicity Prediction**: Toxicity and synthetic accessibility assessment using eToxPred
- **Molecular Dynamics Simulations**: Stability analysis using GROMACS
- **Cost and Time Efficient**: Reduces drug discovery costs and time by ~75%

## Workflow

```
ChEMBL Database (2.4M compounds)
         ↓
Chemical Filters (Lipinski, Veber, PAINS)
         ↓
Molecular Docking (AutoDock Vina)
         ↓
ADME Filtration (SwissADME)
         ↓
ML Toxicity & SA Prediction (eToxPred)
         ↓
Molecular Dynamics (GROMACS)
         ↓
Final Candidates
```

## Installation

### Prerequisites

```bash
# Python 3.8+
pip install pandas rdkit-pypi meeko vina mysql-connector-python

# Required external tools
- AutoDock Vina
- GROMACS
- Meeko
```

### Dependencies

```python
pandas
rdkit
meeko
vina
mysql-connector-python
```

## Usage

### 1. Database Extraction (`dataset.py`)

Extract drug-like compounds from ChEMBL database:

```python
python dataset.py
```

**Filters applied:**
- Molecular weight ≤ 500 Da
- AlogP ≤ 5
- Hydrogen bond acceptors ≤ 10
- Hydrogen bond donors ≤ 5
- Polar surface area ≤ 140 Ų
- Rotatable bonds ≤ 10

### 2. Chemical Property Filtering (`filtering.py`)

Apply Lipinski's Rule of Five and Veber's rules:

```python
python filtering.py
```

**Criteria:**
- **Lipinski's Rule of Five**: MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10
- **Veber's Rule**: PSA ≤ 140, Rotatable bonds ≤ 10
- **Additional**: MW ≥ 200, LogP ≥ -2

### 3. Molecular Docking (`docking.py`)

Perform high-throughput molecular docking:

```python
python docking.py
```

**Configuration:**
- Grid box: 20×20×20 Å
- Exhaustiveness: 32
- Number of poses: 20
- Output: Top 5 poses per ligand

**Required files:**
- `smiles_only.csv`: Input SMILES strings
- `ap_pol_only.pdbqt`: Prepared receptor file

## Project Structure

```
autofilter/
├── dataset.py              # ChEMBL database extraction
├── filtering.py            # Chemical property filtering
├── docking.py              # Molecular docking pipeline
├── chembl_drug_like_compounds.csv
├── smiles_filtered.csv
├── pdbqt_output/          # Docking poses
└── energy_output/         # Binding energies
```

## Results

When applied to malaria drug discovery:
- **Initial compounds**: 2,400,000 (ChEMBL)
- **After preprocessing**: 2,200,000
- **After docking**: 3,194 high-affinity compounds
- **After ADME/ML filtering**: 10 candidates
- **After MD simulations**: 5 final candidates

**Best compound (L1):**
- Binding affinity: -10.830 kcal/mol
- Toxicity score: 0.26 (low)
- Synthetic accessibility: 0.09 (easy)
- *In vitro* inhibition: >50% (preliminary)

## Validation

AutoFilter was validated using DNA polymerase theta (Pol θ) as a proof of concept:
- Successfully identified all 7 known inhibitors from ChEMBL
- Demonstrates high accuracy and reliability for virtual screening

## Applications

While developed for malaria, AutoFilter is a **universal framework** applicable to:
- HIV drug discovery
- Tuberculosis inhibitors
- Neurodegenerative disorders
- Cancer therapeutics
- Any disease with a known target protein


## Computational Requirements

- **CPU**: Multi-core processor recommended for parallel docking
- **RAM**: 16+ GB for large database screening
- **Storage**: ~50 GB for ChEMBL database and output files
- **Time**: ~24-48 hours for complete screening of 2M compounds

## Limitations

- Requires receptor structure (PDB file) with known or predicted binding site
- PDBQT format required for receptor and ligands
- MySQL database required for ChEMBL extraction
- MD simulations computationally intensive for large numbers of compounds

## Future Directions

- Integration with additional ML models for affinity prediction
- Support for ensemble docking with multiple receptor conformations
- Automated binding site prediction
- GPU acceleration for molecular dynamics
- Web interface for easier accessibility

## Acknowledgments

Special thanks to Dr. Kamal Singh for guidance on the chemistry and validation of results.


---

**Note**: This framework significantly accelerates early-stage drug discovery but does not replace experimental validation. All computational predictions should be confirmed through *in vitro* and *in vivo* studies.
