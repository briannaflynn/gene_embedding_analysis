import os
import glob
import math

# ========== USER CONFIGURATION ==========
input_dir = "/home/ubuntu/NEW/AFM_CANDIDATES_182/"
output_base = "/home/ubuntu/NEW/AFM_CANDIDATES_182_OUTPUT"
bash_dir = "/home/ubuntu/NEW/"
gpu_ids = list(range(8))
jobs_per_batch = 64
colabfold_args = "--model-type alphafold2_multimer_v3 --num-recycle 0 --num-models 1 --no-use-probs-extra"
# ========================================

os.makedirs(output_base, exist_ok=True)
os.makedirs(bash_dir, exist_ok=True)

# STEP 1: Get all input FASTA files
all_fastas = {
    os.path.splitext(os.path.basename(f))[0]: f
    for f in glob.glob(os.path.join(input_dir, "*.fasta"))
}

# ✅ STEP 2: Find completed jobs (any file ending with 'done.txt')
completed = set()
for dirpath, dirnames, filenames in os.walk(output_base):
    for fname in filenames:
        if fname.endswith("done.txt"):
            completed.add(os.path.basename(dirpath))
            break

print(f"✅ Found {len(completed)} completed jobs.")
print(completed)
# STEP 3: Keep only unfinished ones
unfinished_fastas = [path for name, path in all_fastas.items() if name not in completed]
# print(unfinished_fastas)

# print(len(unfinished_fastas))
# STEP 4: Split into batches
total_batches = math.ceil(len(unfinished_fastas) / jobs_per_batch)

for batch_idx in range(total_batches):
    batch_fastas = unfinished_fastas[batch_idx * jobs_per_batch: (batch_idx + 1) * jobs_per_batch]
    bash_lines = ["#!/bin/bash\n", f"echo 'Launching AF2-Multimer batch {batch_idx}...'\n"]

    for i, fasta in enumerate(batch_fastas):
        gpu_id = gpu_ids[i % len(gpu_ids)]
        fasta_name = os.path.basename(fasta)
        output_path = os.path.join(output_base, os.path.splitext(fasta_name)[0])
        line = f"CUDA_VISIBLE_DEVICES={gpu_id} colabfold_batch {colabfold_args} \"{fasta}\" \"{output_path}\" &\n"
        bash_lines.append(line)

    bash_lines.append("wait\n")
    bash_lines.append(f"echo 'Batch {batch_idx} complete.'\n")

    script_path = os.path.join(bash_dir, f"launch_afm_batch_{batch_idx}.sh")
    with open(script_path, "w") as f:
        f.writelines(bash_lines)

    print(f"✅ Wrote: {script_path}")

print(f"\n🚀 {total_batches} new bash script(s) generated for unfinished jobs.")

