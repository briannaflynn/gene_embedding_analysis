import sys
import os
import numpy as np
import pandas as pd
import pytest
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

# Add parent dir (where models.py lives) to sys.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from models import ManualStackingSGD, ModelRunner, build_models


@pytest.fixture
def small_classification_data():
    """Generate a tiny synthetic dataset for testing."""
    X, y = make_classification(
        n_samples=30, n_features=5, n_informative=3,
        n_redundant=0, random_state=42
    )
    X = pd.DataFrame(X, columns=[f"f{i}" for i in range(X.shape[1])])
    y = pd.Series(y)
    return X, y


def test_manual_stacking_fit_predict(small_classification_data):
    X, y = small_classification_data
    model = ManualStackingSGD()
    model.fit(X, y)

    preds = model.predict(X)
    scores = model.decision_function(X)

    assert len(preds) == len(y)
    assert set(np.unique(preds)).issubset({0, 1})
    assert scores.shape[0] == len(y)


def test_modelrunner_get_model_outputs(small_classification_data):
    X, y = small_classification_data
    runner = ModelRunner(X, y, X)

    outputs = runner.get_model_outputs(LogisticRegression(max_iter=500))

    assert "predictions" in outputs
    assert "probabilities" in outputs
    assert "decision_scores" in outputs
    assert len(outputs["predictions"]) == len(X)


def test_modelrunner_run_all_models(small_classification_data):
    X, y = small_classification_data
    runner = ModelRunner(X, y, X)

    models = {
        "LR": LogisticRegression(max_iter=500),
        "RF": RandomForestClassifier(n_estimators=10, random_state=42),
    }

    results = runner.run_all_models(models)

    assert set(results.keys()) == {"LR", "RF"}
    for out in results.values():
        assert "predictions" in out


def test_build_models_returns_expected_keys():
    models = build_models()
    expected_keys = {
        "Logistic Regression",
        "Random Forest",
        "SGD",
        "LinearSVC",
        "TPOT Gradient Boosting",
        "TPOT Stacked RF SGD",
        "Pipeline RF+SGD",
    }
    assert expected_keys.issubset(models.keys())

