import os
import pandas as pd
import pytest
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from dataset_builder import ComplexDatasetBuilder

@pytest.fixture
def dummy_data(tmp_path):
    """Create dummy pickle files and directories for testing."""
    # Complex labels
    complex_df = pd.DataFrame({
        "GeneAB": ["A_B", "C_D"],
        "Same_Complex": [1, 0]
    })
    complex_path = tmp_path / "complex.pkl"
    complex_df.to_pickle(complex_path)

    # Pairwise correlation
    pw_df = pd.DataFrame({
        "Gene_AB": ["A_B", "C_D"],
        "Correlation": [0.9, -0.2]
    })
    pw_path = tmp_path / "pw.pkl"
    pw_df.to_pickle(pw_path)

    # scGPT directory with cosine similarity
    scgpt_dir = tmp_path / "scgpt"
    scgpt_dir.mkdir()
    scgpt_df = pd.DataFrame({"Cosine_Similarity": [0.8, 0.3]})
    scgpt_file = scgpt_dir / "scGPT_dummy_output_XXXXXXXXXXXX.pkl"
    scgpt_df.to_pickle(scgpt_file)

    # Geneformer directory with cosine similarity
    gf_dir = tmp_path / "gf"
    gf_dir.mkdir()
    gf_df = pd.DataFrame({"Cosine_Similarity": [0.5, 0.7]})
    gf_file = gf_dir / "GF_dummy_output_XXXXXXXXXXXX.pkl"
    gf_df.to_pickle(gf_file)

    output_path = tmp_path / "out.pkl"

    return {
        "complex_path": str(complex_path),
        "pw_path": str(pw_path),
        "scgpt_dir": str(scgpt_dir),
        "gf_dir": str(gf_dir),
        "output_path": str(output_path),
    }


def test_load_complex_df(dummy_data):
    builder = ComplexDatasetBuilder(**dummy_data)
    df = builder.load_complex_df()
    assert "GeneAB" in df.columns
    assert "Same_Complex" in df.columns
    assert len(df) == 2


def test_load_pw_df(dummy_data):
    builder = ComplexDatasetBuilder(**dummy_data)
    df = builder.load_pw_df()
    # Should rename Gene_AB -> GeneAB
    assert "GeneAB" in df.columns
    assert "Correlation" in df.columns


def test_load_similarity_data_scgpt(dummy_data):
    builder = ComplexDatasetBuilder(**dummy_data)
    df = builder.load_similarity_data(dummy_data["scgpt_dir"], "scGPT")
    assert not df.empty
    assert "scGPT_dummy_output" in df.columns


def test_build_dataset(dummy_data):
    builder = ComplexDatasetBuilder(**dummy_data)
    df = builder.build_dataset()
    # Should merge everything correctly
    assert "Same_Complex" in df.columns
    assert "Correlation" in df.columns
    assert any(col.startswith("scGPT") for col in df.columns)
    assert any(col.startswith("GF") for col in df.columns)
    assert len(df) == 2


def test_save_dataset(dummy_data):
    builder = ComplexDatasetBuilder(**dummy_data)
    df = builder.build_dataset()
    builder.save_dataset(df)
    # File should exist and be loadable
    loaded = pd.read_pickle(dummy_data["output_path"])
    assert loaded.equals(df)

