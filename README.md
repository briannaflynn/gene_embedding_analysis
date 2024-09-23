# gene_embedding_analysis

# Gene to GO Term API Integration

This Python script allows you to retrieve UniProt IDs, GO Biological Process (BP) terms, and associated gene counts using the UniProt and QuickGO APIs. Additionally, it can check for shared GO terms between two genes and perform GO term data analysis for large datasets.

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Convert Gene Name to UniProt ID](#convert-gene-name-to-uniprot-id)
  - [Get GO Biological Process Terms](#get-go-biological-process-terms)
  - [Fetch GO Term Details](#fetch-go-term-details)
  - [Fetch All GO BP Terms and Gene Counts](#fetch-all-go-bp-terms-and-gene-counts)
  - [Check Shared GO Terms](#check-shared-go-terms)
- [License](#license)

## Requirements

This script requires the following Python libraries:

- `requests`
- `bioservices`
- `time`

You can install them using `pip`:

```bash
pip install requests bioservices
