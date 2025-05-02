# ADMET Prediction Repository

This repository contains a Python-based pipeline for predicting ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) properties of small molecules using machine learning. The project uses RDKit for molecular descriptor calculation and a Random Forest model for prediction.

## Features
- Predicts five key ADMET properties: LogP, Solubility (LogS), Blood-Brain Barrier (BBB) penetration, hERG inhibition, and Oral Bioavailability.
- Includes scripts for data preparation, model training, prediction, and visualization.
- Provides a sample dataset with SMILES strings and ADMET properties.

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/AhmedFikry90/ADMET-prediction.git
   cd admet-prediction
   ```
2. Create a virtual environment and activate it:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
4. Ensure you have the sample data in the `data/` directory or provide your own dataset.

## Usage
1. **Prepare Data**: Use the `prepare_data.py` script to process your dataset and calculate molecular descriptors.
   ```bash
   python src/prepare_data.py --input data/sample_admet_data.csv --output data/processed_data.csv
   ```
2. **Train Model**: Train the Random Forest model using `train_model.py`.
   ```bash
   python src/train_model.py --data data/processed_data.csv --model_output models/admet_model.pkl
   ```
3. **Make Predictions**: Use `predict.py` to predict ADMET properties for new compounds.
   ```bash
   python src/predict.py --model models/admet_model.pkl --input data/new_compounds.csv --output predictions/predictions.csv
   ```
4. **Visualize Results**: Use `visualize_results.py` to generate plots of the predicted ADMET properties.
   ```bash
   python src/visualize_results.py --predictions predictions/predictions.csv --output plots/admet_bar_plot.png
   ```
5. **Analyze Correlations**: Use `analyze_correlations.py` to plot correlations between predicted and actual ADMET values (if actual values are available in the predictions file).
   ```bash
   python src/analyze_correlations.py --predictions predictions/predictions.csv --output plots/correlation_plot.png
   ```

## Directory Structure
- `data/`: Contains sample data and processed datasets.
- `src/`: Python scripts for data preparation, training, prediction, and visualization.
- `models/`: Directory to save trained models.
- `predictions/`: Directory to save prediction results.
- `plots/`: Directory to save visualization plots.
- `requirements.txt`: List of dependencies.
- `LICENSE`: MIT License for the project.

## Sample Data
The `data/sample_admet_data.csv` file contains SMILES strings and five ADMET properties for training. You can replace it with your own dataset in the same format.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue for suggestions.

## Contact
For questions, please open an issue on GitHub.