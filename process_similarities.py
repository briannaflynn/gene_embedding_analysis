import pandas as pd
import os
import pickle
from pathlib import Path

def process_similarity_pkls_from_files(pkl_files, output_dir, export_gene_pairs=False):
    """
    Process individual .pkl files to extract and store gene pair information and similarity values separately.
    
    Parameters:
        pkl_files (list): List of file paths to .pkl files (as strings or Path objects).
        reference_filename (str): Filename of the reference .pkl file (must be one of pkl_files). Take the first one in the list as default.
        output_dir (str or Path): Output directory where the processed files will be saved.
        export_gene_pairs (bool): Whether to export Gene_A and Gene_B columns from the reference file.
    """
    pkl_files = [Path(f) for f in pkl_files]
    ref_path = pkl_files[0]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if export_gene_pairs:
        print('Exporting gene pairs...')
        ref_df = pd.read_pickle(ref_path)
        gene_pairs = ref_df[['Gene_A', 'Gene_B']]
        gene_pairs.to_pickle(output_dir / "gene_pairs.pkl")
        print(f"Saved gene pairs from {ref_path} to gene_pairs.pkl")

    # Store similarity columns for each other file
    for file in pkl_files:
        print(file)
        df = pd.read_pickle(file)
        print(df.columns)
        if 'Cosine_Similarity' not in df.columns:
            
            print(f"Skipping {file.name} - no Similarity column")
            continue

        sim_only = df[['Cosine_Similarity']]
        out_path = output_dir / f"{file.stem}_similarity_only.pkl"
        sim_only.to_pickle(out_path)
        print(f"Saved similarity-only file: {out_path.name}")


scgpt=["./scgpt_embeddings/cosine_similarities/scGPT_CP_embeddings_similarities.pkl",
"./scgpt_embeddings/cosine_similarities/scGPT_bc_embeddings_similarities.pkl",
"./scgpt_embeddings/cosine_similarities/scGPT_brain_embeddings_similarities.pkl",
"./scgpt_embeddings/cosine_similarities/scGPT_kidney_embeddings_similarities.pkl",
"./scgpt_embeddings/cosine_similarities/scGPT_heart_embeddings_similarities.pkl",
"./scgpt_embeddings/cosine_similarities/scGPT_lung_embeddings_similarities.pkl",
"./scgpt_embeddings/cosine_similarities/scGPT_pancancer_embeddings_similarities.pkl",
"./scgpt_embeddings/cosine_similarities/scGPT_human_embeddings_similarities.pkl"]
geneformer = ["GF-12L30M/GF-12L30M_HUMANemb_similarities.pkl",
"GF-12L95M/GF-12L95M_HUMANemb_similarities.pkl",
"GF-12L95MCANCER/GF-12L95MCANCER_UNIPROT_HUMANemb_similarities.pkl",
"GF-20L95M/GF-20L95M_HUMANemb_similarities.pkl",
"GF-6L30M/GF-6L30M_HUMANemb_similarities.pkl"]

print('Begin process scgpt')
process_similarity_pkls_from_files(
    pkl_files=scgpt,
    output_dir="scgpt_split_outputs"
)

print('Begin process geneformer')
process_similarity_pkls_from_files(
    pkl_files=geneformer,
    output_dir="gf_split_outputs",
    export_gene_pairs=True
)