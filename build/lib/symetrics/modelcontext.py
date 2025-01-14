import pickle
import lightgbm as lgb
import numpy as np

class ModelContext:
    def __init__(self, model_path, col_feats):
        """
        Initialize the predictor class.

        Parameters:
        - model_path (str): Path to the saved LGBM model (pkl file).
        - col_feats (list): List of feature names used by the model.
        """
        self.model_path = model_path
        self.col_feats = col_feats
        self.model = self.load_model()

    def load_model(self):
        """
        Load the LGBM model from the pickle file.

        Returns:
        - model: Loaded LGBM model.
        """
        with open(self.model_path, 'rb') as file:
            model = pickle.load(file)
        return model

    def predict(self, input_data):
        """
        Make a prediction using the LGBM model.

        Parameters:
        - input_data (dict): Dictionary where keys are feature names and values are feature values.

        Returns:
        - prediction: The predicted class or probability.
        """
        # Ensure input_data contains all required features
        data = [input_data[feat] for feat in self.col_feats]
        data_array = np.array([data])

        # Make a prediction
        prediction = self.model.predict_proba(data_array)
        return prediction
