# Import libraries
from flask import Flask, jsonify, request
import pandas as pd
import joblib

# Set port
port = 8001

# Initialize API
app = Flask(__name__)

# Define API service
@app.route('/predict', methods=['POST']) # API service path and method to call it
def predict(): # API service functionality
   
    # Get input data in JSON format
    dat = request.json
    
    # Cast data to pandas
    dat = pd.DataFrame(dat)
    
    # Load model
    model = joblib.load('model.pkl')
    
    # Compute predictions
    pred = model.predict_proba(dat).tolist()
    
    # Cast predictions to JSON format
    pred = jsonify({'pred': pred})
    
    # Return predictions
    return pred



# Run API
if __name__ == '__main__':
    app.run(host = '0.0.0.0', port=port, debug=True)