#!/bin/bash

set -e

wget -q https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh

bash Anaconda3-2024.06-1-Linux-x86_64.sh -b -p $HOME/anaconda3

eval "$($HOME/anaconda3/bin/conda shell.bash hook)"

conda init

source ~/.bashrc

conda create -y -n scgpt_conda python=3.9

conda activate scgpt_conda

conda install -y ipykernel

python -m ipykernel install --user --name scgpt_conda --display-name "scgpt_conda 3.9"

pip install scgpt

pip install gseapy==1.1.5

pip install pandas==1.5.3

pip install gdown

pip install biopython

pip install pandas==1.5.3

gdown --folder https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU

pip install PyMuPdf

pip install fitz

gdown --fuzzy https://drive.google.com/file/d/1dlWcqKVVuD0FnfU7Ak-k-KltKQCv0AfK/view?usp=sharing

git clone https://github.com/briannaflynn/gene_embedding_analysis.git

gdown --fuzzy https://drive.google.com/file/d/1XTTo5RfvLP6r-SRb4u36j7JGFizR7rP6/view?usp=sharing

gdown --fuzzy https://drive.google.com/file/d/1Po_Zjz4FtIkzDAIKxGFUhBCy1Qh3bn1n/view?usp=sharing

pip install goatools obonet networkx

wget http://current.geneontology.org/annotations/goa_human.gaf.gz

wget http://purl.obolibrary.org/obo/go/go-basic.obo

gunzip goa_human.gaf.gz

sudo apt install r-base

pip install rpy2

pip install pymer4

conda install -c conda-forge scipy=1.13.0

git clone https://github.com/facebookresearch/esm.git

cd esm

pip install -e .

cd ../

mv esm esm_repo

pip install torch==2.0.0 torchtext==0.15.1 --force-reinstall

pip install "numpy<2.0"

gdown --folder https://drive.google.com/drive/folders/1hGM6aYw9_HNepZfYJgdV96gBsORtG0HM?usp=sharing

echo "Setup complete! Your environment 'scgpt_conda' is ready."