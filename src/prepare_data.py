import pandas as pd
import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np

def calculate_descriptors(smiles):
    """Calculate molecular descriptors for a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    descriptors = {
        'MolWt': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol)
    }
    return descriptors

def prepare_data(input_file, output_file):
    """Prepare dataset by calculating descriptors and saving to a new CSV."""
    # Load data
    df = pd.read_csv(input_file)
    
    # Calculate descriptors for each molecule
    descriptor_list = []
    for smiles in df['smiles']:
        desc = calculate_descriptors(smiles)
        if desc:
            descriptor_list.append(desc)
        else:
            descriptor_list.append({key: np.nan for key in calculate_descriptors(df['smiles'][0]).keys()})
    
    # Create descriptor DataFrame
    desc_df = pd.DataFrame(descriptor_list)
    
    # Combine descriptors with original data
    result_df = pd.concat([df, desc_df], axis=1)
    
    # Drop rows with missing values
    result_df = result_df.dropna()
    
    # Save processed data
    result_df.to_csv(output_file, index=False)
    print(f"Processed data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare data for ADMET prediction.")
    parser.add_argument("--input", required=True, help="Path to input CSV file")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    args = parser.parse_args()
    
    prepare_data(args.input, args.output)