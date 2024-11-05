def test1_convert_ensembl_ids_and_compute_similarity():
    # Create example input DataFrame
    data = {
        'ensembl_id_pair': ["('ENSG00000012048', 'ENSG00000012048')", 
                            "('ENSG00000012048', 'ENSG00000139618')"],
        'correlation': [0.85, 0.65]
    }
    df = pd.DataFrame(data)
    print(df)

    # Create example embeddings dictionary with dummy vectors
    embeddings_dict = {
        'BRCA1': np.array([1.0, 2.0, 3.0]),
        'BRCA2': np.array([4.0, 5.0, 6.0])
    }

    print(embeddings_dict)
    print('\nRunning func')
    convert_ensembl_ids_and_compute_similarity(df, embeddings_dict)

    # Run the function
    result_df = convert_ensembl_ids_and_compute_similarity(df, embeddings_dict)
    
    # Print results for inspection
    print(result_df)

    # Expected values for assertions
    expected_gene_names = [('BRCA1', 'BRCA1'), ('BRCA1', 'BRCA2')]
    expected_cosine_similarities = [
        1 - cosine(embeddings_dict['BRCA1'], embeddings_dict['BRCA1']),
        1 - cosine(embeddings_dict['BRCA1'], embeddings_dict['BRCA2'])
    ]
    
    # Assertions to validate the output
    assert result_df['gene_name_1'].tolist() == [pair[0] for pair in expected_gene_names], \
        f"Expected gene_name_1: {[pair[0] for pair in expected_gene_names]}"
    assert result_df['gene_name_2'].tolist() == [pair[1] for pair in expected_gene_names], \
        f"Expected gene_name_2: {[pair[1] for pair in expected_gene_names]}"
    assert np.allclose(result_df['cosine_similarity'], expected_cosine_similarities, equal_nan=True), \
        f"Expected cosine similarities: {expected_cosine_similarities}"

    print("Test passed: convert_ensembl_ids_and_compute_similarity function works as expected.")

def test2_convert_ensembl_ids_and_compute_similarity():
    # Create example input DataFrame
    data = {
        'ensembl_id_pair': ["('ENSG00000012048', 'ENSG00000012048')", 
                            "('ENSG00000012048', 'ENSG00000200000')"],  # ENSG00000200000 has no associated gene name
        'correlation': [0.85, 0.65]
    }
    df = pd.DataFrame(data)
    print(df)

    # Create example embeddings dictionary with dummy vectors
    embeddings_dict = {
        'BRCA1': np.array([1.0, 2.0, 3.0]),
        'BRCA2': np.array([4.0, 5.0, 6.0])
    }

    print(embeddings_dict)
    print('\nRunning func')
    convert_ensembl_ids_and_compute_similarity(df, embeddings_dict)

    # Run the function
    result_df = convert_ensembl_ids_and_compute_similarity(df, embeddings_dict)
    
    # Print results for inspection
    print(result_df)

    # Expected values for assertions
    expected_gene_names = [('BRCA1', 'BRCA1'), ('BRCA1', 'BRCA2')]
    expected_cosine_similarities = [
        1 - cosine(embeddings_dict['BRCA1'], embeddings_dict['BRCA1']),
        1 - cosine(embeddings_dict['BRCA1'], embeddings_dict['BRCA2'])
    ]
    
    # Assertions to validate the output
    # Print the resulting DataFrame for inspection
    result_df = result_df.replace({None: np.nan})

    # Assertions to verify behavior
    assert result_df['gene_name_1'].iloc[0] == "BRCA1", "Expected gene_name_1 to be 'GENE1'"
    #print(result_df['gene_name_2'].iloc[1])
    assert pd.isna(result_df['gene_name_2'].iloc[1]), "Expected gene_name_2 to be NaN"
    assert pd.isna(result_df['cosine_similarity'].iloc[1]), "Expected cosine_similarity to be NaN"


    print("Test passed: convert_ensembl_ids_and_compute_similarity handles missing gene names as expected.")

