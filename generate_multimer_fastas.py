import argparse
import os
import time
import pandas as pd
from Bio import Entrez

Entrez.email = "brianna.flynn@utexas.edu"  # Replace with your email

def fetch_protein_fasta_from_gene(gene_name, organism="Homo sapiens"):
    print('FETCH')
    try:
        print("begin search")
        search_handle = Entrez.esearch(
            db="protein",
            term=f"{gene_name}[Gene Name] AND {organism}[Organism]",
            retmode="xml",
            retmax=1
        )
        search_results = Entrez.read(search_handle)
        print(search_results)
        search_handle.close()

        if not search_results["IdList"]:
            return None
        
        protein_id = search_results["IdList"][0]
        print(protein_id)
        fetch_handle = Entrez.efetch(
            db="protein",
            id=protein_id,
            rettype="fasta",
            retmode="text"
        )
        fasta_data = fetch_handle.read()
        fetch_handle.close()
        time.sleep(0.4)  # respectful pause for NCBI servers
        return fasta_data.strip()
    except Exception as e:
        print(f"Error retrieving {gene_name}: {e}")
        return None

def generate_multimer_fastas_from_dataframe(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    gene_to_fasta_dict = {}

    print('START gene 2 fasta')
    for idx, row in df.iterrows():
        gene_a = row["candidate_gene_1"]
        gene_b = row["candidate_gene_2"]

        if gene_a not in gene_to_fasta_dict:
            gene_to_fasta_dict[gene_a] = fetch_protein_fasta_from_gene(gene_a)
        if gene_b not in gene_to_fasta_dict:
            gene_to_fasta_dict[gene_b] = fetch_protein_fasta_from_gene(gene_b)

        fasta_a = gene_to_fasta_dict.get(gene_a)
        fasta_b = gene_to_fasta_dict.get(gene_b)

        if fasta_a and fasta_b:
            seq_a = ''.join(fasta_a.splitlines()[1:])
            seq_b = ''.join(fasta_b.splitlines()[1:])
            combined_header = f">{gene_a}_{gene_b}"
            combined_sequence = f"{seq_a}:\n{seq_b}"
            print('Writing OUTPUT')
            output_file = os.path.join(output_dir, f"pair_{idx}_{gene_a}_{gene_b}.fasta")
            with open(output_file, "w") as f:
                f.write(f"{combined_header}\n{combined_sequence}\n")
        else:
            print(f"Skipping {gene_a} and {gene_b} — missing sequence")

def main():
    parser = argparse.ArgumentParser(description="Generate AlphaFold2 multimer FASTA files from a CSV of gene pairs.")
    parser.add_argument("csv_file", help="Path to the CSV file containing 'candidate_gene_1' and 'candidate_gene_2' columns.")
    parser.add_argument("--output_dir", default="colabfold_inputs_output", help="Directory to write FASTA files to.")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)
    print(df)
    required_cols = {"candidate_gene_1", "candidate_gene_2"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")


    print('Generating fastas, exporting to', args.output_dir)

    generate_multimer_fastas_from_dataframe(df, args.output_dir)

if __name__ == "__main__":
    main()

