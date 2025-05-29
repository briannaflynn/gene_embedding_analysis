import os
from pathlib import Path

def generate_run_parallel_script(input_dir, output_dir, num_gpus=8, model_type="alphafold2_multimer_v3", output_script="run_parallel.sh"):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    fasta_files = sorted(input_dir.glob("*.fasta"))

    if not fasta_files:
        raise ValueError("No .fasta files found in the input directory.")

    output_lines = ["#!/bin/bash\n"]

    for i, fasta_path in enumerate(fasta_files):
        gpu_id = i % num_gpus
        input_file = fasta_path
        output_subdir = output_dir / fasta_path.stem
        output_subdir.mkdir(parents=True, exist_ok=True)

        # Build command with optional model_type argument
        model_arg = f"--model-type {model_type}" if model_type else ""
        cmd = f'CUDA_VISIBLE_DEVICES={gpu_id} colabfold_batch {model_arg} "{input_file}" "{output_subdir}" &'
        output_lines.append(cmd)

    output_lines.append("wait\n")

    # Write to script
    with open(output_script, "w") as f:
        f.write("\n".join(output_lines))

    os.chmod(output_script, 0o755)
    print(f"Script '{output_script}' generated with {len(fasta_files)} jobs using model_type={model_type}")

import sys
input = sys.argv[1]
output = sys.argv[2]

generate_run_parallel_script(input, output, num_gpus=8)

