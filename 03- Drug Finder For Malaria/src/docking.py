import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import meeko
from vina import Vina
import os
import sys

 

# Read the CSV file containing SMILES

try:

    smiles_df = pd.read_csv('smiles_only.csv')  # Replace with your CSV file name

except FileNotFoundError:

    print("Error: 'smiles_only.csv' file not found.")

    sys.exit(1)

 
# Ensure output directories exist

os.makedirs('pdbqt_output', exist_ok=True)

os.makedirs('energy_output', exist_ok=True)

 
# Check if receptor file exists

if not os.path.exists('ap_pol_only.pdbqt'):

    print("Error: 'ap_pol_only.pdbqt' file not found.")

    sys.exit(1)

 

# Function to process a single SMILES

def process_smiles(smiles, index):

    try:

        # Convert SMILES to RDKit molecule

        lig = Chem.MolFromSmiles(smiles)

 

        protonated_lig = Chem.AddHs(lig)

        AllChem.EmbedMolecule(protonated_lig)

 

            # Prepare ligand using Meeko

        meeko_prep = meeko.MoleculePreparation()

        meeko_prep.prepare(protonated_lig)

        lig_pdbqt = meeko_prep.write_pdbqt_string()

 

            # Set up Vina

        v = Vina(sf_name='vina')

        v.set_receptor('ap_pol_only.pdbqt')

        v.set_ligand_from_string(lig_pdbqt)

        v.compute_vina_maps(center=[9.2, -32.5, 112.7], box_size=[20, 20, 20])

    except Exception as e:

        print("Couldn't do Meeko")

 

        # Try to score and minimize, but continue if it fails

    try:

        energy = v.score()

        print(f'Score before minimization for ligand {index}: {energy[0]:.3f} (kcal/mol)')

 

        energy_minimized = v.optimize()

        print(f'Score after minimization for ligand {index}: {energy_minimized[0]:.3f} (kcal/mol)')

    except Exception as e:

        print("Couldn't score")


            # Dock the ligand

    try:      # Dock the ligand

        v.dock(exhaustiveness=32, n_poses=20)

 

                # Write poses to PDBQT file

        pdbqt_filename = f'pdbqt_output/ligand_{index}_vina_out.pdbqt'

        v.write_poses(pdbqt_filename, n_poses=5, overwrite=True)

 

                # Create energy dataframe

        column_names = ["total", "inter", "intra", "torsions", "intra best pose"]

        df = pd.DataFrame(v.energies(), columns=column_names)

 

                # Save energies to CSV

        csv_filename = f'energy_output/energies_{index}.csv'

        df.to_csv(csv_filename, index=False)

 

        print(f"Successfully processed ligand {index}")

    except Exception as e:

        print("Couldn't dock")

 

# Process each SMILES in the dataframe

for index, row in smiles_df.iterrows():

    smiles = row['smiles']  # Assuming 'smiles' is the column name in your CSV

    print(f"\nProcessing ligand {index}")

    process_smiles(smiles, index)