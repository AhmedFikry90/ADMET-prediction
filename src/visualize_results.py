import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def visualize_results(predictions_file, output_file):
    """Generate a bar plot of predicted ADMET properties for each compound."""
    # Load predictions
    df = pd.read_csv(predictions_file)
    
    # ADMET properties to visualize
    admet_props = ['logP_pred', 'logS_pred', 'bbb_prob_pred', 'herg_prob_pred', 'bioavailability_pred']
    labels = ['LogP', 'LogS', 'BBB Prob', 'hERG Prob', 'Bioavailability']
    
    # Create bar plot
    fig, ax = plt.subplots(figsize=(10, 6))
    x = range(len(df))
    width = 0.15
    
    for i, prop in enumerate(admet_props):
        ax.bar([xi + width * i for xi in x], df[prop], width, label=labels[i])
    
    # Customize plot
    ax.set_xlabel('Compounds')
    ax.set_ylabel('Predicted Values')
    ax.set_title('Predicted ADMET Properties')
    ax.set_xticks([xi + width * 2 for xi in x])
    ax.set_xticklabels(df['smiles'], rotation=45, ha='right')
    ax.legend()
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save plot
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Bar plot saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize predicted ADMET properties.")
    parser.add_argument("--predictions", required=True, help="Path to predictions CSV file")
    parser.add_argument("--output", required=True, help="Path to save the bar plot")
    args = parser.parse_args()
    
    visualize_results(args.predictions, args.output)