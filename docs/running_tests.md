# Running the Test Suites

This project uses [**pytest**](https://docs.pytest.org/) for testing.  
All tests are located in the `tests/` directory and cover the main modules in this repository.

---

## 1. Install dependencies

Make sure you are in your project’s Python environment and install pytest:

```bash
pip install pytest
```
If your tests depend on additional libraries (like scikit-learn, pandas, numpy, etc.), install them too:

`pip install pandas numpy scikit-learn networkx tqdm tpot`

## 2. Project layout
A typical project structure looks like this:

```
project_root/
│
├── dataset_builder.py
├── protein_complex.py
├── models.py
├── ...
└── tests/
    ├── test_dataset_builder.py
    ├── test_protein_complex.py
    └── test_models.py
```

Each test_*.py file contains a pytest suite for one module.

## 3. Running all tests
From the project root (the same directory that contains dataset_builder.py):

`pytest`

This will:

* Discover all test_*.py files under tests/

* Run all test functions starting with test_

* Print a summary of passed and failed tests

## 4. Running a specific test file
To run only the tests for one module:

`pytest tests/test_dataset_builder.py`

or

`pytest tests/test_protein_complex.py`

## 5. Running a specific test function
To run a single test function, use the -k flag with part of the test name:

`pytest -k "test_manual_stacking_fit_predict"`

## 6. Getting more detailed output [RECOMMENDED]
Use the -v (verbose) flag for more detail:

`pytest -v`

Use the -s flag to show print statements (useful for debugging):

`pytest -v -s`

## 7. Typical workflow

* Edit your code (e.g., models.py).

* Run the relevant pytest file:

`pytest tests/test_models.py -v`

* Fix any failing tests.

* Rerun until everything passes :)
