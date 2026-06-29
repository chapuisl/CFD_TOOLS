#!/bin/bash

set -e

echo "===================================="
echo "  Setup project: ModifSolutionh5"
echo "===================================="

# Nom de l'environnement
ENV_NAME="envh5Modif"

# Chemin du projet
PROJECT_PATH="/Users/lilian_chapuis/Documents/zz01-GITLAB/PYTHON"

# 1. Vérification Python
echo "Checking Python..."
python3 --version

# 2. Aller dans le dossier projet
echo "Going to project directory..."
cd "$PROJECT_PATH"

# 3. Création venv (dans le dossier projet)
echo "Creating venv: $ENV_NAME"
python3 -m venv "$ENV_NAME"

# 4. Activation
echo "Activating..."
source "$ENV_NAME/bin/activate"

# 5. Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# 6. Installation des dépendances
echo "Installing packages..."
pip install numpy h5py

# 7. Fin
echo "===================================="
echo "Setup complete ✔"
echo "Now you have to source the new env:" 
echo "source $PROJECT_PATH/$ENV_NAME/bin/activate"
