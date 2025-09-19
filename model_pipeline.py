import os
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import SGDClassifier, LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Normalizer
from sklearn.svm import LinearSVC
from tpot.builtins import StackingEstimator


class ManualStackingSGD(BaseEstimator, ClassifierMixin):
    """
    Custom stacked classifier combining a RandomForest and SGDClassifier.
    A RandomForest is trained first, its probability outputs are concatenated
    with the original features, normalized, and then used to train an SGDClassifier.
    """

    def __init__(self, rf_params=None, sgd_params=None, norm='l1'):
        """
        Initialize the stacked classifier.

        Parameters
        ----------
        rf_params : dict, optional
            Parameters for RandomForestClassifier. Defaults to a fixed set of values.
        sgd_params : dict, optional
            Parameters for SGDClassifier. Defaults to a fixed set of values.
        norm : str, default="l1"
            Normalization type used by sklearn.preprocessing.Normalizer.
        """
        self.rf_params = rf_params or {
            'bootstrap': True,
            'criterion': "entropy",
            'max_features': 0.45,
            'min_samples_leaf': 3,
            'min_samples_split': 2,
            'n_estimators': 100,
            'random_state': 42
        }
        self.sgd_params = sgd_params or {
            'alpha': 0.01,
            'eta0': 0.01,
            'fit_intercept': False,
            'l1_ratio': 0.75,
            'learning_rate': "constant",
            'loss': "squared_hinge",
            'penalty': "elasticnet",
            'power_t': 0.5,
            'random_state': 42
        }
        self.norm = norm
        self.rf = RandomForestClassifier(**self.rf_params)
        self.sgd = SGDClassifier(**self.sgd_params)
        self.normalizer = Normalizer(norm=self.norm)

    def fit(self, X, y):
        """
        Fit the RandomForest and then the SGDClassifier on normalized stacked features.

        Parameters
        ----------
        X : array-like
            Feature matrix.
        y : array-like
            Labels.

        Returns
        -------
        self
        """
        self.rf.fit(X, y)
        rf_scores = self.rf.predict_proba(X)[:, 1].reshape(-1, 1)
        X_stacked = np.hstack((X.values if hasattr(X, "values") else X, rf_scores))
        X_norm = self.normalizer.fit_transform(X_stacked)
        self.sgd.fit(X_norm, y)
        return self

    def predict(self, X):
        """
        Predict class labels for input samples.

        Parameters
        ----------
        X : array-like
            Feature matrix.

        Returns
        -------
        np.ndarray
            Predicted labels.
        """
        rf_scores = self.rf.predict_proba(X)[:, 1].reshape(-1, 1)
        X_stacked = np.hstack((X.values if hasattr(X, "values") else X, rf_scores))
        X_norm = self.normalizer.transform(X_stacked)
        return self.sgd.predict(X_norm)

    def decision_function(self, X):
        """
        Return decision scores for samples.

        Parameters
        ----------
        X : array-like
            Feature matrix.

        Returns
        -------
        np.ndarray
            Decision scores from the SGDClassifier.
        """
        rf_scores = self.rf.predict_proba(X)[:, 1].reshape(-1, 1)
        X_stacked = np.hstack((X.values if hasattr(X, "values") else X, rf_scores))
        X_norm = self.normalizer.transform(X_stacked)
        return self.sgd.decision_function(X_norm)


class ModelRunner:
    """
    Orchestrates training and evaluation of multiple classifiers
    on a gene-pair similarity dataset.
    """

    def __init__(self, train_X, train_y, full_X):
        """
        Initialize with training and full dataset splits.

        Parameters
        ----------
        train_X : pd.DataFrame
            Training feature matrix.
        train_y : pd.Series
            Training labels.
        full_X : pd.DataFrame
            Full feature matrix used for generating predictions and probabilities.
        """
        self.train_X = train_X
        self.train_y = train_y
        self.full_X = full_X

    def get_model_outputs(self, model):
        """
        Fit a model and return predictions, probabilities, and decision scores.

        Parameters
        ----------
        model : sklearn estimator
            Classifier implementing fit/predict.

        Returns
        -------
        dict
            Contains trained model, predictions, probabilities, and decision scores.
        """
        outputs = {'model': model}

        print("Fitting model...")
        model.fit(self.train_X, self.train_y)

        print("Generating predictions...")
        outputs['predictions'] = model.predict(self.full_X)

        if hasattr(model, "predict_proba"):
            try:
                outputs['probabilities'] = model.predict_proba(self.full_X)[:, 1]
            except Exception:
                outputs['probabilities'] = None
        else:
            outputs['probabilities'] = None

        if hasattr(model, "decision_function"):
            try:
                outputs['decision_scores'] = model.decision_function(self.full_X)
            except Exception:
                outputs['decision_scores'] = None
        else:
            outputs['decision_scores'] = None

        return outputs

    def run_all_models(self, models):
        """
        Train and evaluate all provided models.

        Parameters
        ----------
        models : dict
            Dictionary mapping model names to sklearn estimators.

        Returns
        -------
        dict
            Mapping of model names to output dictionaries.
        """
        results = {}
        for name, model in tqdm(models.items(), desc="Fitting models"):
            print(f"\nRunning model: {name}")
            try:
                results[name] = self.get_model_outputs(model)
            except Exception as e:
                print(f"Error with model '{name}': {e}")
        return results


def build_models():
    """
    Construct and return a dictionary of models to be tested.

    Returns
    -------
    dict
        Model name -> sklearn estimator.
    """
    gb_pipeline = GradientBoostingClassifier(
        learning_rate=0.5,
        max_depth=10,
        max_features=0.6,
        min_samples_leaf=4,
        min_samples_split=19,
        n_estimators=100,
        subsample=0.8
    )

    stacked_pipeline = make_pipeline(
        StackingEstimator(
            estimator=RandomForestClassifier(
                bootstrap=True,
                criterion="entropy",
                max_features=0.45,
                min_samples_leaf=3,
                min_samples_split=2,
                n_estimators=100
            )
        ),
        Normalizer(norm="l1"),
        SGDClassifier(
            alpha=0.01, eta0=0.01, fit_intercept=False, l1_ratio=0.75,
            learning_rate="constant", loss="squared_hinge",
            penalty="elasticnet", power_t=0.5
        )
    )

    return {
        "Logistic Regression": LogisticRegression(max_iter=10000),
        "Random Forest": RandomForestClassifier(n_estimators=10000, n_jobs=-1),
        "SGD": SGDClassifier(loss='log_loss', max_iter=10000, random_state=42),
        "LinearSVC": LinearSVC(max_iter=10000, verbose=0, random_state=42),
        "TPOT Gradient Boosting": gb_pipeline,
        "TPOT Stacked RF SGD": ManualStackingSGD(),
        "Pipeline RF+SGD": stacked_pipeline
    }


def main():
    """
    Example entry point for running all models and saving results.
    """
    df = pd.read_pickle("/home/ubuntu/full_dataset.pkl")

    features = [
        'scGPT_bc_embeddings_Cosine_Similarity',
        'scGPT_pancancer_embeddings_Cosine_Similarity',
        'scGPT_lung_embeddings_Cosine_Similarity',
        'scGPT_heart_embeddings_Cosine_Similarity',
        'scGPT_brain_embeddings_Cosine_Similarity',
        'scGPT_kidney_embeddings_Cosine_Similarity',
        'scGPT_human_embeddings_Cosine_Similarity',
        'GF-6L30M_HUMANemb_Cosine_Similarity',
        'GF-20L95M_HUMANemb_Cosine_Similarity',
        'GF-12L95M_HUMANemb_Cosine_Similarity',
        'GF-12L95MCANCER_UNIPROT_HUMANemb_Cosine_Similarity',
        'GF-12L30M_HUMANemb_Cosine_Similarity',
        'Correlation'
    ]

    X = df[features]
    y = df['Same_Complex']

    X_train = X[~df['Test']]
    y_train = y[~df['Test']]
    X_full = X

    runner = ModelRunner(X_train, y_train, X_full)
    models = build_models()
    results = runner.run_all_models(models)

    with open("model_results.pkl", "wb") as f:
        pickle.dump(results, f)
