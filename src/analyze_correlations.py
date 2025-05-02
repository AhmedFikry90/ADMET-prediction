import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def analyze_correlations(predictions_file, output_file):
    """Generate scatter plots to analyze correlations between predicted and actual ADMET values."""
    # Load predictions
    df = pd.read_csv(predictions_file)
    
    # Check if actual values are present
    actual_columns = ['logP', 'logS', 'bbb_prob', 'herg_prob', 'bioavailability']
    predicted_columns = ['logP_pred', 'logS_pred', 'bbb_prob_pred', 'herg_prob_pred', 'bioavailability_pred']
    labels = ['LogP', 'LogS', 'BBB Prob', 'hERG Prob', 'Bioavailability']
    
    # Filter columns that exist in the DataFrame
    available_pairs = [(act, pred, label) for act, pred, label in zip(actual_columns, predicted_columns, labels) if act in df.columns]
    
    if not available_pairs:
        print("No actual ADMET values found in the predictions file for correlation analysis.")
        return
    
    # Create scatter plots
    n_plots = len(available_pairs)
    fig, axes = plt.subplots(nrows=n_plots, ncols=1, figsize=(8, 4 * n_plots))
    
    if n_plots == 1:
        axes = [axes]  # Ensure axes is iterable
    
    for ax, (act_col, pred_col, label) in zip(axes, available_pairs):
        ax.scatter(df[act_col], df[pred_col], alpha=0.5)
        ax.set_xlabel(f'Actual {label}')
        ax.set_ylabel(f'Predicted {label}')
        ax.set_title(f'Correlation: Actual vs Predicted {label}')
        
        # Add a diagonal line (y=x) for reference
        min_val = min(df[act_col].min(), df[pred_col].min())
        max_val = max(df[act_col].max(), df[pred_col].max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='Ideal (y=x)')
        ax.legend()
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save plot
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Correlation plot saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze correlations between predicted and actual ADMET values.")
    parser.add_argument("--predictions", required=True, help="Path to predictions CSV file")
    parser.add_argument("--output", required=True, help="Path to save the correlation plot")
    args = parser.parse_args()
    
    analyze_correlations(args.predictions, args.output)