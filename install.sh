#!/usr/bin/env bash
set -e

# Create a fresh virtual environment if none exists
if [ ! -d "venv" ]; then
  echo "Creating virtual environment..."
  python3 -m venv venv
fi

# Activate the virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install requirements
echo "Installing dependencies..."
pip install -r requirements.txt

# Install the package in editable mode
echo "Installing this package..."
pip install -e .

echo "Installation complete. Activate the environment with: source venv/bin/activate"

