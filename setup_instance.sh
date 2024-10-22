#!/bin/bash

# Download Anaconda installer
echo "Downloading Anaconda..."
wget https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh

# Run the Anaconda installer
echo "Installing Anaconda..."
bash Anaconda3-2024.06-1-Linux-x86_64.sh -b  # -b flag for silent installation (accepting default settings)

# Initialize conda
echo "Initializing Anaconda..."
source ~/.bashrc

# Create the 'scgpt_conda' environment with Python 3.9
echo "Creating conda environment 'scgpt_conda'..."
conda create -n "scgpt_conda" python=3.9 -y

# Activate the environment
echo "Activating 'scgpt_conda' environment..."
source activate scgpt_conda

# Install ipykernel
echo "Installing ipykernel..."
conda install ipykernel -y

# Add kernel to Jupyter
echo "Setting up Jupyter kernel..."
python -m ipykernel install --user --name scgpt_conda --display-name "scgpt_conda 3.9"

# Install additional packages
echo "Installing scgpt, gseapy, and gdown..."
pip install scgpt gseapy gdown PyMuPdf fitz

# Download files from Google Drive
echo "Downloading data from Google Drive folder..."
gdown --folder https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU

echo "Setup complete!"
