import json
from pathlib import Path
import pandas as pd

def extract_confidence_scores_single(json_path):
    """
    Extract ipTM, pTM, and ranking_confidence from a single scores_rank_001_*.json file.

    Args:
        json_path (str or Path): Full path to the JSON file.

    Returns:
        dict with keys: json_path, iptm, ptm, ranking_confidence, bin, pair_dir
    """
    json_path = Path(json_path)
    with open(json_path) as f:
        data = json.load(f)
    print(data.keys())
    iptm = data.get("iptm")
    ptm = data.get("ptm")

    if iptm is None or ptm is None:
        return None  # Skip incomplete results

    ranking_conf = 0.8 * iptm + 0.2 * ptm

    parent = json_path.parent.name
    bin_str = parent.split("_bin")[-1]
    try:
        bin_num = int(bin_str)
    except ValueError:
        bin_num = None

    return {
        "json_path": str(json_path),
        "iptm": iptm,
        "ptm": ptm,
        "ranking_confidence": ranking_conf,
        "bin": bin_num,
        "pair_dir": parent
    }

def extract_confidence_scores_single(json_path):
    """
    Extracts ptm and fallback ranking_confidence from a single JSON file
    when iptm is missing. Assumes ptm is still informative.

    Args:
        json_path (str or Path): Full path to scores_rank_001_*.json

    Returns:
        dict or None
    """
    json_path = Path(json_path)
    with open(json_path) as f:
        data = json.load(f)

    ptm = data.get("ptm")

    if ptm is None:
        print(f"[SKIPPED] Missing ptm in: {json_path}")
        print("Available keys:", list(data.keys()))
        return None

    # Use ptm as proxy for ranking_confidence when iptm is unavailable
    ranking_conf = ptm

    parent = json_path.parent.name
    bin_str = parent.split("_bin")[-1]
    try:
        bin_num = int(bin_str)
    except ValueError:
        bin_num = None

    return {
        "json_path": str(json_path),
        "ptm": ptm,
        "ranking_confidence": ranking_conf,
        "bin": bin_num,
        "pair_dir": parent
    }


from glob import glob

json_paths = glob("output/**/scores_rank_001_*.json", recursive=True)
print(json_paths)

from pathlib import Path

json_paths = list(Path("output").rglob("scores_rank_001_*.json"))
print(json_paths)

json_paths = list(Path("/home/ubuntu/output").rglob("*scores_rank_001_alphafold2_multimer_v3_model_*.json"))
print(json_paths)

records = []
for path in json_paths:
    result = extract_confidence_scores_single(path)
    if result:
        records.append(result)

df = pd.DataFrame(records)
print(df)
df.to_csv('RF_ptm_by_bin.csv')
