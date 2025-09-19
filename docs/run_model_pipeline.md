# Model Training and Stacking Framework

`model_pipeline.py`

This module provides classes and utilities for training, evaluating, and stacking multiple classifiers on a gene-pair similarity dataset. It includes a custom stacked model (`ManualStackingSGD`), a runner for orchestrating training/evaluation (`ModelRunner`), and a function to build a suite of candidate models.

---

## `ManualStackingSGD`

### Purpose
Implements a two-stage stacked classifier:
1. A **RandomForestClassifier** is trained first.
2. Its probability outputs are concatenated with the original feature matrix.
3. The combined feature set is normalized and passed to an **SGDClassifier**.

This allows tree-based feature extraction followed by a linear classifier.

### Initialization Parameters
- `rf_params` (`dict`, optional): Parameters for the RandomForest. Defaults include `entropy` criterion, `n_estimators=100`, etc.  
- `sgd_params` (`dict`, optional): Parameters for the SGDClassifier. Defaults include `elasticnet` penalty, `squared_hinge` loss.  
- `norm` (`str`, default = `"l1"`): Normalization type for `sklearn.preprocessing.Normalizer`.

### Methods
- `fit(X, y)`: Trains RF, generates probability scores, stacks with input features, normalizes, and fits SGD.  
- `predict(X)`: Returns predicted class labels from the stacked pipeline.  
- `decision_function(X)`: Returns decision scores from the SGD stage.

---

## `ModelRunner`

### Purpose
Automates the process of fitting multiple classifiers, generating predictions, and collecting outputs (predictions, probabilities, decision scores).

### Initialization
- `train_X` (`pd.DataFrame`): Training features.  
- `train_y` (`pd.Series`): Training labels.  
- `full_X` (`pd.DataFrame`): Full feature matrix for downstream predictions.

### Methods
#### `get_model_outputs(model)`
- **Input:** Any sklearn-compatible estimator with `.fit()` and `.predict()`.  
- **Output:** `dict` containing:
  - `'model'`: trained model object.  
  - `'predictions'`: predicted labels on `full_X`.  
  - `'probabilities'`: predicted probabilities (if available).  
  - `'decision_scores'`: decision scores (if available).

#### `run_all_models(models)`
- **Input:** `dict` of `{model_name: sklearn_estimator}`.  
- **Output:** `dict` mapping model names to results from `get_model_outputs`.  
- Uses `tqdm` for progress tracking and prints diagnostics.

---

## `build_models()`

### Purpose
Constructs a dictionary of baseline and TPOT-inspired models to test.

### Returns
A dictionary mapping descriptive names to instantiated models, including:
- Logistic Regression  
- Random Forest  
- SGDClassifier (logistic loss)  
- LinearSVC  
- TPOT Gradient Boosting (custom GB pipeline)  
- TPOT Stacked RF+SGD (via `ManualStackingSGD`)  
- Pipeline RF+SGD (via sklearn `make_pipeline` and `StackingEstimator`)

---

## `main()`

### Purpose
Provides an example workflow for running all models on a pre-built dataset.

### Workflow
1. Load dataset: `/home/ubuntu/full_dataset.pkl`.  
2. Define **features**: scGPT similarity columns, Geneformer similarity columns, and correlation.  
3. Split into:
   - `X_train`, `y_train` (training subset where `Test == False`)  
   - `X_full` (entire feature set).  
4. Instantiate `ModelRunner` and build models via `build_models()`.  
5. Train and evaluate all models with `runner.run_all_models(models)`.  
6. Save results to `model_results.pkl`.

---

## Example Output Structure (`results`)

Each model’s entry in `results` contains:
```python
{
  'model': <trained sklearn estimator>,
  'predictions': np.ndarray of shape (n_samples,),
  'probabilities': np.ndarray of shape (n_samples,) or None,
  'decision_scores': np.ndarray of shape (n_samples,) or None
}

