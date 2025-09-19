#!/usr/bin/env python3
"""
Setup script for gene_embedding_analysis package
"""

from setuptools import setup, find_packages

# Read dependencies from requirements.txt
try:
    with open('requirements.txt', 'r') as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]
except FileNotFoundError:
    requirements = []

setup(
    name="gene_embedding_analysis",
    version="0.1.0",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=requirements,
)
