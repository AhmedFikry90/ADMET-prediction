import pandas as pd
import argparse
import pickle
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

def make_predictions(model_file, input_file, output_file):
    """Make ADMET predictions for new compounds using the trained model."""
    # Load model
    with open(model_file, 'rb') as f:
        model = pickle.load(f)
    
    # Load new compounds
    df = pd.read_csv(input_file)
    
    # Calculate descriptors
    descriptor_list = []
    for smiles in df['smiles']:
        desc = calculate_descriptors(smiles)
        if desc:
            descriptor_list.append(desc)
        else:
            descriptor_list.append({key: np.nan for key in calculate_descriptors(df['smiles'][0]).keys()})
    
    desc_df = pd.DataFrame(descriptor_list)
    desc_df = desc_df.dropna()
    
    # Make predictions
    features = ['MolWt', 'LogP', 'NumHDonors', 'NumHAcceptors', 'TPSA', 'NumRotatableBonds']
    X = desc_df[features]
    predictions = model.predict(X)
    
    # Create output DataFrame
    pred_df = pd.DataFrame(predictions, columns=['logP_pred', 'logS_pred', 'bbb_prob_pred', 'herg_prob_pred', 'bioavailability_pred'])
    result_df = pd.concat([df.iloc[desc_df.index].reset_index(drop=True), pred_df], axis=1)
    
    # Save predictions
    result_df.to_csv(output_file, index=False)
    print(f"Predictions saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict ADMET properties for new compounds.")
    parser.add_argument("--model", required=True, help="Path to trained model file")
    parser.add_argument("--input", required=True, help="Path to input CSV file with SMILES")
    parser.add_argument("--output", required=True, help="Path to save predictions")
    args = parser.parse_args()
    
    make_predictions(args.model, args.input, args.output)