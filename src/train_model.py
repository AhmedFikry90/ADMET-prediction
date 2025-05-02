import pandas as pd
import argparse
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

def train_model(data_file, model_output):
    """Train a Random Forest model for ADMET prediction and save it."""
    # Load processed data
    df = pd.read_csv(data_file)
    
    # Features (descriptors) and targets (ADMET properties)
    features = ['MolWt', 'LogP', 'NumHDonors', 'NumHAcceptors', 'TPSA', 'NumRotatableBonds']
    targets = ['logP', 'logS', 'bbb_prob', 'herg_prob', 'bioavailability']
    
    X = df[features]
    y = df[targets]
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Train model
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # Evaluate model
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    print(f"Model Evaluation - MSE: {mse:.4f}, R2: {r2:.4f}")
    
    # Save model
    with open(model_output, 'wb') as f:
        pickle.dump(model, f)
    print(f"Model saved to {model_output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train ADMET prediction model.")
    parser.add_argument("--data", required=True, help="Path to processed data CSV file")
    parser.add_argument("--model_output", required=True, help="Path to save the trained model")
    args = parser.parse_args()
    
    train_model(args.data, args.model_output)