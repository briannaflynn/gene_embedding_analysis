# scGPT Embeddings Download and Set-Up Documentation

This guide provides step-by-step instructions to resolve compilation issues when installing scGPT with flash-attn on Lambda.

---
### 1. Setting Up the Instance

1. Download and install Anaconda:
```
wget <https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh> 
bash Anaconda3-2024.06-1-Linux-x86_64.sh source ~/.bashrc
```

2. Create and activate the scGPT Conda environment:

```conda create -n "scgpt_conda" python=3.9 conda activate scgpt_conda```

3. Install Jupyter kernel:

```conda install ipykernel python -m ipykernel install --user --name scgpt_conda --display-name "scgpt_conda 3.9"```

4. Install necessary Python packages:

```pip install scgpt pip install gseapy pip install gdown pip install PyMuPdf pip install fitz```

Note: Installing scGPT without flash attention is fine when not training.

5. Download required models using gdown:

```gdown --folder <https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU>```

 `Additional model downloads:`

```
gdown --folder <https://drive.google.com/drive/folders/1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB> # CP gdown --folder <https://drive.google.com/drive/folders/1vf1ijfQSk7rGdDGpBntR5bi5g6gNt-Gx> # Brain gdown --folder <https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU> # BC gdown --folder <https://drive.google.com/drive/folders/1GcgXrd7apn6y4Ze_iSCncskX3UsWPY2r> # Heart gdown --folder <https://drive.google.com/drive/folders/16A1DJ30PT6bodt4bWLa4hpS7gbWZQFBG> # Lung gdown --folder <https://drive.google.com/drive/folders/1S-1AR65DF120kNFpEbWCvRHPhpkGK3kK> # Kidney gdown --folder <http://drive.google.com/drive/folders/13QzLHilYUd0v3HTwa_9n4G4yEF-hdkqa> # Pan Cancer gdown --folder <https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y> # Human (all 33 million)
```

 6. Clone the gene embedding analysis repository:

```git clone <https://github.com/briannaflynn/gene_embedding_analysis.git>```

### 2. Installing Dependencies

1\. Update package lists and install build tools:

```sudo apt-get update sudo apt-get install build-essential ninja-build```

2. Install Python development headers and libraries:

```sudo apt-get install python3-dev```

### 3. Verifying Compiler Setup

1. Create a test program file named `test_pybind11.cpp`:

#include <pybind11/pybind11.h> int main() { return 0; }

2. Compile the test program to check if the compiler can find headers:

```g++ test_pybind11.cpp -o test_pybind11 -I${CONDA_PREFIX}/include/python3.8 -L${CONDA_PREFIX}/lib -lpython3.8```

### 4. Setting Environment Variables

Set necessary environment variables to ensure correct library linking:

```bash
export CPLUS_INCLUDE_PATH=${CONDA_PREFIX}/include:${CONDA_PREFIX}/include/python3.8:$CPLUS_INCLUDE_PATH export C_INCLUDE_PATH=${CONDA_PREFIX}/include:${CONDA_PREFIX}/include/python3.8:$C_INCLUDE_PATH export LIBRARY_PATH=${CONDA_PREFIX}/lib:$LIBRARY_PATH export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH
```

### 5. Ensuring pip Uses the Correct Environment

#### Issue:
When using Lambda, pip may still point to the base user Python installation instead of the Conda environment.

#### Fix:
1. Add the following to `~/.bashrc`:

```export PATH="/home/ubuntu/anaconda3/envs/scgpt_conda/bin:$PATH"```

2. Apply the changes:

```source ~/.bashrc```

3. Verify the correct pip path:

```which pip```

`Expected output:`

```/home/ubuntu/anaconda3/envs/scgpt_conda/bin/pip```

### 6. Installing Flash Attention

The original instructions reference `flash-attn<1.0.5`, but this conflicts with CUDA 12. Instead, install version 1.0.6:

```pip install "flash-attn==1.0.6" --no-build-isolation```

 `- `--no-build-isolation`: Ensures that the package is built using the current environment's installed packages instead of an isolated build.
- If this causes issues in the future, downgrading to CUDA 11.7 may be necessary.

---
### 7. Installing Missing Dependencies

1. Install `wandb` for logging and visualization:`

```pip install wandb```

### 8. Issues with scvi-tools and optax (JAX)

The documentation states that scGPT supports Python 3.8+, but certain dependencies required for fine-tuning need **Python 3.9**.

To resolve this:
- Set up a clean Python 3.9 environment.
- Reinstall scGPT and dependencies.
- Verify the exact installation versions used in the scGPT Docker branch (see the relevant pull request for details).

---
This guide should help you get scGPT working with flash-attn on Lambda or other cloud service while avoiding common pitfalls.
