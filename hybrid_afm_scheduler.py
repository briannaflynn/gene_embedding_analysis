import os
import glob
import subprocess
import time
from datetime import datetime
import pandas as pd
from tqdm import tqdm

# ========== USER CONFIGURATION ==========
input_dir = "/home/ubuntu/NEW/AFM_CANDIDATES_182/"  # Directory of input FASTA files
output_base = "/home/ubuntu/NEW/AFM_CANDIDATES_182_OUTPUT"               # Directory where ColabFold will write results
log_dir = "/home/ubuntu/NEW/LOGS"                         # <-- Replace with your preferred log directory
model_type = "alphafold2_multimer_v3"             # Or other model if desired
max_concurrent_jobs = 64                          # Number of parallel jobs allowed (adjust if needed)

# Ensure output directories exist
os.makedirs(output_base, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

# ========================================

# Gather all input FASTA files
fasta_files = sorted(glob.glob(os.path.join(input_dir, "*.fasta")))
active_processes = []
logs = []
pbar = tqdm(total=len(fasta_files), desc="AF2-Multimer Jobs", ncols=100)

def launch_job(fasta_path):
    output_dir = os.path.join(output_base, os.path.splitext(os.path.basename(fasta_path))[0])
    log_file = os.path.join(log_dir, os.path.basename(fasta_path).replace(".fasta", ".log"))
    cmd = [
        "colabfold_batch",
        "--model-type", model_type,# "--num-recycle 0 --num-models 1 --no-use-probs-extra",
        fasta_path,
        output_dir
    ]
    env = os.environ.copy()
    log_handle = open(log_file, "w")
    proc = subprocess.Popen(cmd, env=env, stdout=log_handle, stderr=subprocess.STDOUT)
    return proc, fasta_path, log_file, datetime.now(), log_handle

# Main loop
i = 0
print(f"[{datetime.now()}] Launching up to {max_concurrent_jobs} jobs at once...\n")
while i < len(fasta_files) or active_processes:
    # Launch new jobs if under limit
    while len(active_processes) < max_concurrent_jobs and i < len(fasta_files):
        fasta = fasta_files[i]
        proc, fasta, log_file, start_time, log_handle = launch_job(fasta)
        active_processes.append((proc, fasta, log_file, start_time, log_handle))
        print(f"[{datetime.now()}] ⏳ Launched {fasta}")
        i += 1

    # Check for completed jobs
    still_running = []
    for proc, fasta, log_file, start_time, log_handle in active_processes:
        if proc.poll() is None:
            still_running.append((proc, fasta, log_file, start_time, log_handle))
        else:
            end_time = datetime.now()
            runtime = (end_time - start_time).total_seconds()
            log_handle.close()
            logs.append({
                "fasta": os.path.basename(fasta),
                "log_file": log_file,
                "runtime_sec": runtime,
                "completed_at": end_time.strftime("%Y-%m-%d %H:%M:%S")
            })
            pbar.update(1)
            print(f"[{end_time}] ✅ Completed {fasta} in {runtime:.1f} seconds")
    active_processes = still_running
    time.sleep(2)

pbar.close()

# Write final summary log
log_df = pd.DataFrame(logs)
summary_path = os.path.join(log_dir, "afm_job_log_summary.csv")
log_df.to_csv(summary_path, index=False)
print(f"[{datetime.now()}] All jobs complete. Summary written to {summary_path}")

