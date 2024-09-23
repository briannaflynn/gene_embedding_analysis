# -*- coding: utf-8 -*-

import requests
from bioservices import QuickGO
import time

def gene_name_to_uniprot_id(gene_name: str) -> str:
    """
    Converts a gene name to the corresponding UniProt ID using a direct UniProt API request.

    Args:
        gene_name (str): The gene name or symbol (e.g., "BRCA1").

    Returns:
        str: The UniProt ID associated with the gene name, or None if not found.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        'query': f'gene_exact:{gene_name} AND organism_id:9606',  # Filtering for human (organism ID 9606)
        'fields': 'accession',  # Fetch only the UniProt accession ID
        'format': 'tsv',        # Requesting tab-separated values format
        'size': 1               # Limit results to one entry
    }

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()  # Raise an error for HTTP error codes

        # Parse the response content
        lines = response.text.splitlines()
        if len(lines) > 1:  # First line is the header, next lines contain data
            uniprot_id = lines[1].split('\t')[0]  # Extract the UniProt ID from the first result
            return uniprot_id
        else:
            print(f"No UniProt ID found for gene: {gene_name}")
            return None

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching UniProt ID for {gene_name}: {e}")
        return None

def get_gene_go_terms(gene: str) -> list:
    """
    Returns a list of GO Biological Process terms associated with a gene using QuickGO.

    Args:
        gene (str): The gene symbol or identifier (e.g., UniProt ID).

    Returns:
        list: A list of GO Biological Process terms (GO IDs) associated with the gene.
    """
    quickgo = QuickGO()

    try:
        # Correct aspect value for Biological Process is 'P'
        print(f"Querying GO terms for gene: {gene}")
        results = quickgo.Annotation(geneProductId=gene, geneProductType="protein", aspect="P")

        # Check if the results are in the expected format
        if isinstance(results, dict) and 'results' in results:
            # Extract GO terms from the results
            go_terms = [annotation['goId'] for annotation in results['results']
                        if annotation.get('goAspect') == 'biological_process']
            return go_terms
        else:
            print("Unexpected response format or empty response:", results)
            return []

    except Exception as e:
        print(f"An error occurred: {e}")
        return []

# return a dictionary of all data associated with a given GO BP term

def fetch_go_term_details(term: str) -> dict:
    """
    Fetches and returns all available data for a single GO term using the QuickGO API.

    Args:
        term (str): A GO term ID (e.g., 'GO:0008150').

    Returns:
        dict: A dictionary containing the data for the GO term, or None if the request fails.
    """
    base_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"
    headers = {'Accept': 'application/json'}

    retries = 3
    for attempt in range(retries):
        try:
            # detailed information related to the GO term
            response = requests.get(f"{base_url}/{term}", headers=headers)
            response.raise_for_status()  # Raise an exception for HTTP errors

            # full response data as a dictionary
            data = response.json()
            print(f"Full data for {term}:\n", data)
            return data

        except requests.exceptions.RequestException as e:
            print(f"Attempt {attempt + 1} failed for {term}: {e}")
            time.sleep(2)  # wait two seconds

    # Return None if all attempts fail
    print(f"Failed to fetch data for {term} after {retries} attempts.")
    return None

# Go through all available GO BP terms in QuickGo, get gene counts

def fetch_all_go_bp_terms_and_counts() -> dict:
    """
    Fetches all GO Biological Process terms from QuickGO and retrieves their gene counts,
    printing the terms and counts as they are retrieved.

    Returns:
        dict: A dictionary with GO term IDs as keys and counts of associated genes as values.
    """
    base_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"
    headers = {'Accept': 'application/json'}
    go_terms = []
    params = {
        'aspect': 'biological_process',  # restrict to Biological Process terms
        'pageSize': 1000,
        'page': 1
    }
    term_gene_counts = {}

    try:
        while True:
            response = requests.get(base_url, headers=headers, params=params)
            response.raise_for_status()
            data = response.json()

            # collect GO term IDs
            terms = [term['id'] for term in data.get('results', [])]
            go_terms.extend(terms)

            # print each retrieved term and get counts
            for term in terms:
                print(f"Retrieved GO term: {term}")
                # fetch the gene count for this term
                counts = get_term_gene_counts([term])  # fetch counts for the current term
                term_gene_counts.update(counts)  # update the main dictionary with fetched counts

            # use total number of terms if available
            total_terms = data.get('pageInfo', {}).get('total', 0)
            current_page = data.get('pageInfo', {}).get('current', 1)
            results_per_page = data.get('pageInfo', {}).get('resultsPerPage', len(data.get('results', [])))

            # do more pages exist?
            if current_page * results_per_page >= total_terms:
                break  # exit if the last page is reached

            params['page'] += 1  # increment the page number

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching GO BP terms: {e}")

    return term_gene_counts

# count genes associated with GO BP term, if server error print out all avaiable data - most likely term is 'obsolete', in that case print 'Obsolete? True'

def get_term_gene_counts(go_terms: list) -> dict:
    """
    Fetches the count of gene products (genes) associated with each GO term using QuickGO's annotation endpoint,
    with retry logic and error handling.

    Args:
        go_terms (list): A list of GO term IDs (e.g., ['GO:0008150', 'GO:0009987']).

    Returns:
        dict: A dictionary with GO term IDs as keys and counts of associated genes as values.
    """
    base_url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
    headers = {'Accept': 'application/json'}
    term_gene_counts = {}

    for term in go_terms:
        retries = 3
        for attempt in range(retries):
            try:
                # query QuickGO for annotations related to the GO term
                params = {
                    'goId': term,
                    'limit': 0,  # Set limit to 0 because we only need the count
                    'facetField': 'geneProductId'
                }
                response = requests.get(base_url, headers=headers, params=params)
                response.raise_for_status()  # Raise an exception for HTTP errors

                data = response.json()

                # count of gene products (genes) from the response
                gene_count = data.get('numberOfHits', 0)

                # store in dictionary
                term_gene_counts[term] = gene_count
                print(f"Count for {term}: {gene_count}")
                break  # break if retry loop successful

            except requests.exceptions.RequestException as e:
                print(f"Attempt {attempt + 1} failed for {term}: {e}")
                time.sleep(2)  # wait a couple seconds

                # for last attempt, log the error and set the count to 0
                if attempt == retries - 1:
                    print(f"Failed to fetch data for {term} after {retries} attempts.")
                    term_gene_counts[term] = 0
                    details=fetch_go_term_details(term=term)
                    obs = details['results'][0]['isObsolete']
                    print(f'Obsolete term?: {obs}')
                    print(details)

    return term_gene_counts

# get counts for all GO BP terms found, if server returns error print out full data returned for that term from API query
# term_gene_counts = fetch_all_go_bp_terms_and_counts()

def check_shared_go_terms(gene_1: str, gene_2: str, n_shared_threshold: int = None) -> int:
    """
    Check if two genes share a GO Biological Process term. Optionally filter the terms
    by a threshold of how many genes share the term.

    Args:
        gene_1 (str): The first gene to compare. UNIPROT ID, not gene name.
        gene_2 (str): The second gene to compare. UNIPROT ID, not gene name.
        n_shared_threshold (int, optional): Minimum number of genes that must share a term
                                            for it to be considered in the comparison.

    Returns:
        int: 1 if the genes share at least one filtered GO term, 0 otherwise.
    """
    # get_gene_go_terms returns a list, get ordered set for comparison
    go_terms_1 = set(get_gene_go_terms(gene_1))
    go_terms_2 = set(get_gene_go_terms(gene_2))

    # if threshold is none, just compare - do they share any GO terms? Then 1, else - 0.
    if n_shared_threshold is None:
        shared_terms = go_terms_1 & go_terms_2
        return 1 if shared_terms else 0

    # Get counts of genes for each GO term
    terms_1 = get_term_gene_counts(go_terms_1)
    terms_2 = get_term_gene_counts(go_terms_2)

    def filter_terms(term_dict, thresh):
      filt = {key: value for key, value in term_dict.items() if value >= thresh}
      return set(list(filt.keys()))

    shared_filtered_terms = filter_terms(terms_1, n_shared_threshold) & filter_terms(terms_2, n_shared_threshold)

    return 1 if shared_filtered_terms else 0

"""
DEMO
"""

# Should print the UniProt ID for BRCA1, need Uniprot IDs for QuickGO API queries later
uniprot_id = gene_name_to_uniprot_id('BRCA1')
print(uniprot_id)

# Output will be a list of GO biological process terms for the given gene's Uniprot ID
terms = get_gene_go_terms('P38398')
print(terms)

print(f'There are {len(terms)} GO BP terms associated with P38398/BRCA1')

BRCA1 = "P38398"
TP53 = "P04637"

check_shared_go_terms(BRCA1, TP53)

