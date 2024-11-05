#!/usr/bin/python

def test_ensembl_to_gene_name():
    """
    Tests the ensembl_to_gene_name function with a known Ensembl ID and expected output.
    """
    # Known Ensembl ID for the BRCA1 gene in humans
    ensembl_id = "ENSG00000012048"
    expected_gene_name = "BRCA1"  # The expected gene name for this Ensembl ID
    
    # Call the function
    gene_name = ensembl_to_gene_name(ensembl_id)
    print(gene_name)
    
    # Check if the returned gene name matches the expected gene name
    assert gene_name == expected_gene_name, f"Expected {expected_gene_name}, but got {gene_name}"
    
    print("Test passed: ensembl_to_gene_name returns the expected gene name.")

# Run the test
test_ensembl_to_gene_name()
