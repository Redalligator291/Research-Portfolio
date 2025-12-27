import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    mw = Descriptors.ExactMolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    psa = Descriptors.TPSA(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    
    return pd.Series({
        'MW': mw,
        'LogP': logp,
        'HBD': hbd,
        'HBA': hba,
        'PSA': psa,
        'RotatableBonds': rotatable_bonds
    })

def apply_filters(row):
    # Lipinski's Rule of Five
    lipinski = (
        (row['MW'] <= 500) &
        (row['LogP'] <= 5) &
        (row['HBD'] <= 5) &
        (row['HBA'] <= 10)
    )
    
    # Veber's Rule
    veber = (
        (row['PSA'] <= 140) &
        (row['RotatableBonds'] <= 10)
    )
    
    # Additional drug-likeness criteria
    drug_like = (
        (row['MW'] >= 200) &
        (row['LogP'] >= -2)
    )
    
    return lipinski and veber and drug_like

# Load the data
df = pd.read_csv(r'C:\Users\vaithi\Chemistry 42\smiles.csv')

# Calculate properties
properties = df['smiles'].apply(calculate_properties)
df = pd.concat([df, properties], axis=1)

# Apply filters
df['Passes_Filters'] = df.apply(apply_filters, axis=1)

# Filter the dataframe
filtered_df = df[df['Passes_Filters']]

# Save the filtered compounds
filtered_df.to_csv('smiles_filtered.csv', index=False)

print(f"Original compounds: {len(df)}")
print(f"Filtered compounds: {len(filtered_df)}")
print("Filtered compounds saved to 'filtered_drug_like_compounds.csv'")