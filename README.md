## Installation

This repository includes a convenience script, `install.sh`, to set up a clean Python environment with all required dependencies.

### Steps

1. **Clone the repository** (if you haven’t already):

```bash
   git clone https://github.com/briannaflynn/gene_embedding_analysis.git
   cd gene_embedding_analysis
   git checkout GOBP-classifier
```

2. **Run the install script**:
   
   `bash install.sh`

   This will:
   - Create a new virtual environment in `./venv` (if one does not already exist).
   - Activate the environment.
   - Upgrade `pip` to the latest version.
   - Install all dependencies listed in `requirements.txt`.
   - Install this repository as a local package in *editable* mode (`pip install -e .`), so any code changes you make are reflected immediately without reinstalling.

3. **Activate the environment** whenever you want to use the package:
   `source venv/bin/activate`

4. **Deactivate the environment** when you’re done:
   `deactivate`

---

## Requirements

- Python 3.9+  
- pandas, numpy, scikit-learn  
- matplotlib, tqdm  
- tpot (for `StackingEstimator`)

Use `bash install.sh` or use pip to install dependencies (without automatic venv setup):

`pip install -r requirements.txt`

---

# Protein Complex (or GO BP term) Prediction Pipeline

This repository provides a full workflow for building datasets, creating balanced train/test splits, training multiple machine learning models, and visualizing evaluation metrics. The project is structured into four main components:

1. **Dataset Builder**  
2. **Stratified Group Split**  
3. **Model Pipeline**  
4. **ROC & PR Plotting**

---

## 1. Dataset Builder

Script: `dataset_builder.py`  
Class: `ComplexDatasetBuilder`

### Purpose
Combines multiple data sources into a unified dataset:
- Complex membership labels (`Same_Complex`, `GeneAB`)
- Pairwise correlation values
- Cosine similarity values from **scGPT** and **Geneformer**

### Usage

```python
builder = ComplexDatasetBuilder(
    complex_path='complex_label.pkl',
    pw_path='allpairs_spearman_correlation.pkl',
    scgpt_dir='scgpt_split_outputs',
    gf_dir='gf_split_outputs',
    output_path='full_dataset.pkl'
)

dataset = builder.build_dataset()
builder.save_dataset(dataset)
```

### Output
- A single pickle file (`full_dataset.pkl`) containing:
  - `GeneAB` identifiers
  - `Same_Complex` labels
  - Correlation values
  - Cosine similarities from embedding models

---

## 2. Stratified Group Split

Script: `stratified_split.py`  
Functions:  
- `assign_same_complex_multi`  
- `ProteinComplexGroupSplitter`  
- `stratified_group_split_balanced_per_fold`

### Purpose
Performs **balanced group-aware train/test splits**:
- Keeps all rows of the same gene pair (`pair_id`) together.
- Stratifies by positive/negative labels.
- Balances the number of positives and negatives per fold.

### Usage

```python
folds = stratified_group_split_balanced_per_fold(
    df=prepared_df,
    label_col='Same_Complex',
    group_col='group_id',
    n_splits=5,
    random_state=42
)

for i, (train_idx, test_idx) in enumerate(folds):
    y_train = prepared_df.loc[train_idx, 'Same_Complex']
    y_test = prepared_df.loc[test_idx, 'Same_Complex']
    print(f"Fold {i+1}: Train {len(train_idx)} / Test {len(test_idx)}")
```

### Output
- A list of `(train_idx, test_idx)` tuples for cross-validation.
- Optionally, CSVs for each fold: `balanced_group_stratified_fold_{i}.csv`.

---

## 3. Model Pipeline

Script: `model_pipeline.py`  
Classes & Functions:  
- `ManualStackingSGD` (custom stacked RF+SGD model)  
- `ModelRunner` (train/evaluate models)  
- `build_models()` (constructs a dictionary of candidate models)

### Purpose
Trains and evaluates multiple machine learning models on the prepared dataset, including:
- Logistic Regression  
- Random Forest  
- SGDClassifier  
- LinearSVC  
- Gradient Boosting (TPOT-inspired)  
- Custom stacked RF+SGD (manual and pipeline versions)

### Usage

```python
df = pd.read_pickle("full_dataset.pkl")

features = ["scGPT_bc_embeddings_Cosine_Similarity", ..., "Correlation"]
X = df[features]
y = df["Same_Complex"]

X_train = X[~df['Test']]
y_train = y[~df['Test']]
X_full = X

runner = ModelRunner(X_train, y_train, X_full)
models = build_models()
results = runner.run_all_models(models)

with open("model_results.pkl", "wb") as f:
    pickle.dump(results, f)
```

### Output
- `model_results.pkl`: dictionary of model predictions, probabilities, and decision scores for the full dataset.

---

## 4. ROC & PR Plotting

Script: `roc_pr_plotting.py`  
Functions:  
- `get_model_prefixes`  
- `get_valid_scores`  
- `plot_roc_curves`  
- `plot_pr_curves`

### Purpose
Visualizes model performance across all classifiers with:
- **ROC Curves** (AUC values)  
- **Precision–Recall Curves** (Average Precision values)

### Usage

```python
with open("model_results.pkl", "rb") as f:
    results = pickle.load(f)

# Add outputs to dataset
for mod in results:
    for p in ["predictions", "probabilities", "decision_scores"]:
        df[f"{mod}_{p}"] = results[mod][p]

# Select test set and plot where xf is test set
plot_roc_curves(xf, "Same_Complex")
plot_pr_curves(xf, "Same_Complex")
```

### Output
- `full_train_test_roc.png`: ROC curves for all models.  
- `full_train_test_auprc.png`: Precision–Recall curves for all models.  
- Figures are displayed and also saved to disk.

---

## End-to-End Workflow

1. **Build dataset**  
   Run `dataset_builder.py` to produce `full_dataset.pkl`.  

2. **Split data**  
   Use `stratified_split.py` to create train/test folds.  

3. **Train models**  
   Run `model_pipeline.py` to train models and save results to `model_results.pkl`.  

4. **Visualize results**  
   Run `roc_pr_plotting.py` to generate ROC and PR plots for all models.  


