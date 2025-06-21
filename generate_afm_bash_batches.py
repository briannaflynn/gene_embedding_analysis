import os
import math
import glob

# ========== USER CONFIGURATION ==========
input_dir = "/home/ubuntu/NEW/AFM_CANDIDATES_182/"  # Directory of input FASTA files
output_base = "/home/ubuntu/NEW/AFM_CANDIDATES_182_OUTPUT"                 # Where to write ColabFold outputs
bash_dir = "/home/ubuntu/NEW/"          # Path to write the bash launcher
gpu_ids = list(range(8))                             # GPU IDs to cycle through
jobs_per_batch = 64                                  # Max number of jobs to generate
colabfold_args = "--model-type alphafold2_multimer_v3 --num-recycle 0 --num-models 1 --no-use-probs-extra"

print(colabfold_args)
# ========================================
# Ensure output directories exist
os.makedirs(output_base, exist_ok=True)
os.makedirs(bash_dir, exist_ok=True)

# Gather all .fasta files
fasta_files = sorted(glob.glob(os.path.join(input_dir, "*.fasta")))
total_batches = math.ceil(len(fasta_files) / jobs_per_batch)

for batch_idx in range(total_batches):
    batch_fasta_files = fasta_files[batch_idx * jobs_per_batch : (batch_idx + 1) * jobs_per_batch]
    bash_lines = ["#!/bin/bash\n", f"echo 'Launching AF2-Multimer batch {batch_idx}...'\n"]

    for i, fasta in enumerate(batch_fasta_files):
        gpu_id = gpu_ids[i % len(gpu_ids)]
        fasta_name = os.path.basename(fasta)
        output_path = os.path.join(output_base, os.path.splitext(fasta_name)[0])

        print(colabfold_args)
        line = f"CUDA_VISIBLE_DEVICES={gpu_id} colabfold_batch {colabfold_args} \"{fasta}\" \"{output_path}\" &\n"
        bash_lines.append(line)

    bash_lines.append("wait\n")
    bash_lines.append(f"echo 'Batch {batch_idx} complete.'\n")

    # Write to launch_afm_batch_X.sh
    script_path = os.path.join(bash_dir, f"launch_afm_batch_{batch_idx}.sh")
    with open(script_path, "w") as f:
        f.writelines(bash_lines)

    print(f"✅ Wrote: {script_path}")

print(f"\n🚀 {total_batches} bash scripts generated in {bash_dir}")
