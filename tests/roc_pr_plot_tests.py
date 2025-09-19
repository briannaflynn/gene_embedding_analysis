import sys
import os
import numpy as np
import pandas as pd
import pytest
import matplotlib
matplotlib.use("Agg")  # use non-interactive backend for testing

# Add parent dir (where plotting code lives) to sys.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from plotting import (  # adjust filename if your code is in a different script
    get_model_prefixes,
    get_valid_scores,
    plot_roc_curves,
    plot_pr_curves,
)


@pytest.fixture
def dummy_df():
    """Create a small DataFrame with fake model outputs and labels."""
    return pd.DataFrame({
        "Same_Complex": [0, 1, 0, 1],
        "ModelA_probabilities": [0.1, 0.9, 0.2, 0.8],
        "ModelB_decision_scores": [-1, 2, -0.5, 3],
    })


def test_get_model_prefixes(dummy_df):
    prefixes = get_model_prefixes(dummy_df)
    assert "ModelA" in prefixes
    assert "ModelB" in prefixes
    assert isinstance(prefixes, list)


def test_get_valid_scores_prefers_probabilities(dummy_df):
    scores = get_valid_scores(dummy_df, "ModelA")
    assert np.allclose(scores, [0.1, 0.9, 0.2, 0.8])


def test_get_valid_scores_falls_back_to_decision_scores(dummy_df):
    scores = get_valid_scores(dummy_df, "ModelB")
    assert np.allclose(scores, [-1, 2, -0.5, 3])


def test_get_valid_scores_returns_none_if_missing(dummy_df):
    scores = get_valid_scores(dummy_df, "NonExistentModel")
    assert scores is None


def test_get_valid_scores_skips_none_values():
    df = pd.DataFrame({
        "Same_Complex": [0, 1],
        "ModelC_probabilities": [None, None]
    })
    scores = get_valid_scores(df, "ModelC")
    assert scores is None


def test_plot_roc_curves_creates_file(dummy_df, tmp_path):
    savepath = tmp_path / "roc.png"
    plot_roc_curves(dummy_df, "Same_Complex", savepath=savepath)
    assert savepath.exists()


def test_plot_pr_curves_creates_file(dummy_df, tmp_path):
    savepath = tmp_path / "pr.png"
    plot_pr_curves(dummy_df, "Same_Complex", savepath=savepath)
    assert savepath.exists()

