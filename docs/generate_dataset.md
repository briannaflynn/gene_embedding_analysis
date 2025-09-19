# Complex Dataset Builder

`dataset_builder.py`

This module provides the `ComplexDatasetBuilder` class, which constructs a unified dataset of protein/gene pairs by combining:

- **Complex membership labels** (binary indicator of whether a pair belongs to the same complex).  
- **Pairwise correlation values** (e.g., Spearman correlations).  
- **Cosine similarity values** from two embedding models: **scGPT** and **Geneformer (GF)**.  

The combined dataset can be saved as a pickle file for downstream analysis.

---

## `ComplexDatasetBuilder`

### Purpose
Aggregates and merges multiple data sources into a single dataset, indexed by gene pair identifier (`GeneAB`).

### Initialization Parameters
- `complex_path` (`str`): Path to pickle file containing complex membership info with at least:
  - `Same_Complex`: Binary label (`1` if same complex, else `0`).
  - `GeneAB`: Unique identifier for the gene pair.  
- `pw_path` (`str`): Path to pickle file with pairwise correlation values. Must contain `GeneAB`.  
- `scgpt_dir` (`str`): Directory containing scGPT pickle files with a `Cosine_Similarity` column.  
- `gf_dir` (`str`): Directory containing Geneformer pickle files with a `Cosine_Similarity` column.  
- `output_path` (`str`): Path to save the combined dataset pickle.

---

### Methods

#### `load_complex_df()`
**Purpose:** Load complex membership DataFrame.  
**Output:** `pd.DataFrame` with at least `Same_Complex` and `GeneAB`.

---

#### `load_pw_df()`
**Purpose:** Load pairwise correlation values and standardize column names.  
- Renames `Gene_AB` → `GeneAB` if necessary.  
**Output:** `pd.DataFrame` with correlation values and `GeneAB`.

---

#### `load_similarity_data(directory, prefix)`
**Purpose:** Load cosine similarity values from a directory of pickle files.  
- Iterates through files in the given directory that start with the `prefix`.  
- Each file must contain a column `Cosine_Similarity`.  
- Extracts values and assigns them to a column named by the file prefix.  

**Parameters:**  
- `directory` (`str`): Path to directory of pickle files.  
- `prefix` (`str`): Filename prefix (e.g., `"scGPT"`, `"GF"`).  

**Output:**  
- `pd.DataFrame` with one column per file, named by file prefix.  

---

#### `build_dataset()`
**Purpose:** Merge all data sources into a single dataset.  

**Steps:**  
1. Load complex membership (`load_complex_df`).  
2. Load pairwise correlation data (`load_pw_df`).  
3. Load cosine similarity features from scGPT and GF (`load_similarity_data`).  
4. Concatenate similarities and align with `Same_Complex` and `GeneAB`.  
5. Merge with pairwise correlation values.  
6. Drop rows with missing values.  

**Output:**  
- `pd.DataFrame` containing:  
  - `Same_Complex` (binary label).  
  - `GeneAB` (pair ID).  
  - Pairwise correlation values.  
  - Cosine similarity values from scGPT and GF.  

---

#### `save_dataset(df)`
**Purpose:** Save the final combined dataset.  
**Input:** `pd.DataFrame`.  
**Output:** Pickle file saved to `output_path`.

---

## Example Usage

```python
builder = ComplexDatasetBuilder(
    complex_path='/home/ubuntu/complex_label.pkl',
    pw_path='/home/ubuntu/allpairs_spearman_correlation.pkl',
    scgpt_dir='/home/ubuntu/scgpt_split_outputs',
    gf_dir='/home/ubuntu/gf_split_outputs',
    output_path='full_dataset.pkl'
)

dataset = builder.build_dataset()
builder.save_dataset(dataset)

