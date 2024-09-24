# Gene to GO Term API Integration

This Python script allows you to retrieve UniProt IDs, GO Biological Process (BP) terms, and associated gene counts using the UniProt and QuickGO APIs. Additionally, it can check for shared GO terms between two genes and perform GO term data analysis for large datasets.

## Table of Contents

- [Requirements](#requirements)
- [Usage](#usage)
  - [Convert Gene Name to UniProt ID](#convert-gene-name-to-uniprot-id)
  - [Get GO Biological Process Terms](#get-go-biological-process-terms)
  - [Fetch GO Term Details](#fetch-go-term-details)
  - [Fetch All GO BP Terms and Gene Counts](#fetch-all-go-bp-terms-and-gene-counts)
  - [Check Shared GO Terms](#check-shared-go-terms)
  - [Check Shared Terms and Cosine Similarity](#check-shared-terms-and-cosine-similarity)

## Requirements

This script requires the following Python libraries:

- `requests`
- `bioservices`
- `numpy`
- `scipy`
- `time`

You can install them using `pip`:

```bash
pip install requests bioservices numpy scipy
```

## Usage

### Convert Gene Name to UniProt ID

This function converts a gene name to a UniProt ID.

```python
uniprot_id = gene_name_to_uniprot_id('BRCA1')
print(uniprot_id)
```
* Input: A gene name or symbol (e.g., BRCA1).
* Output: Corresponding UniProt ID (e.g., P38398 for BRCA1).

### Get GO Biological Process Terms

This function retrieves GO Biological Process terms associated with a given gene's UniProt ID using the QuickGO API.

```python
terms = get_gene_go_terms('P38398')
print(terms)
```
* Input: UniProt ID (e.g., P38398 for BRCA1).
* Output: List of GO terms associated with the gene.

### Fetch GO Term Details

This function fetches detailed information about a given GO term from QuickGO.

```python
term_data = fetch_go_term_details('GO:0008150')
print(term_data)
```
* Input: GO term ID (e.g., GO:0008150).
* Output: Dictionary of details for the GO term.

### Fetch All GO BP Terms and Gene Counts

This function fetches all GO Biological Process terms and their associated gene counts.

```python
term_gene_counts = fetch_all_go_bp_terms_and_counts()
print(term_gene_counts)
```
* Output: Dictionary where keys are GO term IDs and values are the counts of associated genes.

### Check Shared GO Terms

This function checks if two genes share a GO Biological Process term. It can also filter based on a threshold for the number of genes sharing the term.

```python
shared = check_shared_go_terms('P38398', 'P04637', n_shared_threshold=10)
print(shared)
```
* Input: Two UniProt IDs (e.g., P38398 for BRCA1 and P04637 for TP53).
* Output: 1 if the genes share at least one GO term, 0 otherwise. Can be filtered by the number of genes sharing the GO term.

### Check Shared Terms and Cosine Similarity

This function takes two dictionaries containing gene names and their embedding vectors and performs the following tasks:

1. Finds UniProt IDs for the gene names.
2. Checks if the genes share a GO Biological Process term.
3. Calculates the cosine similarity between the two embedding vectors.

```python
embedding_dict1 = {'BRCA1': np.random.rand(512)}
embedding_dict2 = {'TP53': np.random.rand(512)}

result = check_shared_terms_and_cosine_similarity(embedding_dict1, embedding_dict2)
print(result)
```

Takes two dictionaries as arguments, each with a gene name as the key and an embedding vector as the value.

Returns a dictionary where the key is a tuple of the two gene names and the value is another dictionary with two keys:
```cosine_similarity```: The cosine similarity between the embedding vectors (float).
```share_GOBP```: Whether the two genes share a GO Biological Process term (True/False).

```bash
{('BRCA1', 'TP53'): {'cosine_similarity': 0.7475048426524952, 'share_GOBP': True}}
```
