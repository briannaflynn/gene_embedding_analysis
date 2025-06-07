import os
import re
from pathlib import Path
from Bio import SeqIO
import sys

"""
Process existing fasta files based on localcolabfold documentation
Multimer needs two fasta sequences joined by a colon with a SINGLE header

EXAMPLE
>gene-a_gene-b_complex
MDKFJSLDKFNDK:
MNSQRHGSKDNL
"""

def extract_gene_name(header):
    """
    Extracts the primary gene name from a FASTA header.
    Handles both UniProt-style and NCBI-style headers.
    """
    # Try UniProt style: RecName: Full=GeneName;
    match = re.search(r'Full=([^\;\n]+)', header)
    if match:
        name = match.group(1).split()[0].upper()
        return name

    # Try NCBI RefSeq style: e.g. 'complement C1r subcomponent-like protein isoform 4 precursor [Homo sapiens]'
    match = re.search(r'\[.*\]', header)
    if match:
        words = header.split()
        if len(words) >= 2:
            return words[1].upper()

    return "UNKNOWN"

def process_and_write_fasta_pair(record1, record2, output_path):
    name1 = extract_gene_name(record1.description)
    name2 = extract_gene_name(record2.description)
    combined_header = f">{name1}_{name2}"
    combined_seq = f"{str(record1.seq)}:\n{str(record2.seq)}"
    
    with open(output_path, "w") as f:
        f.write(combined_header + "\n" + combined_seq + "\n")

def convert_fasta_directory(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(input_dir.glob("*.fa*"))

    for fasta_file in fasta_files:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if len(records) != 2:
            print(f"Skipping {fasta_file.name}, expected 2 sequences, found {len(records)}.")
            continue

        output_path = output_dir / fasta_file.name
        process_and_write_fasta_pair(records[0], records[1], output_path)

    print(f"Processed {len(fasta_files)} files. Output saved to: {output_dir.resolve()}")

INPUT = sys.argv[1]
OUTPUT = sys.argv[2]
convert_fasta_directory(INPUT, OUTPUT)

# EXAMPLE
#af = "alphafold2multimer"
#convert_fasta_directory(af + "_inputs", af + "_combined")


