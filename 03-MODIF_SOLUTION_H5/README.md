# ModifSolutionH5

## Overview

This project provides a set of Python tools to modify AVBP HDF5 solution files. It allows users to manipulate species fields, select geometric regions, and generate the corresponding XMF file for visualization.

---

## Project Structure

```text
.
├── RUN_main.py              # Main executable script
├── Makefile                 # Installation helper
├── setup.sh                 # Creates the Python virtual environment
├── envh5Modif/              # Python virtual environment (created after setup)
├── SCRIPTUSED/
│   ├── __init__.py
│   ├── constant.py
│   ├── geometry.py
│   ├── mesh.py
│   ├── patches_selection.py
│   ├── predefined_patches.py
│   ├── species.py
│   ├── utils.py
│   └── xmf_writer.py
├── *.mesh.h5                # Input mesh files
├── *.h5                     # Solution files
└── *.xmf                    # Generated XMF files
```

* **RUN_main.py** is the only script that should be executed by the user.
* **SCRIPTUSED/** contains all the project modules used internally by `RUN_main.py`.
* The `.mesh.h5` and `.h5` files are the simulation input/output files.
* The generated `.xmf` file allows visualization of the modified solution in ParaView.

---

# Installation

Before running the project, edit **setup.sh** and update the project path so that it matches the location of your repository.

Example:

```bash
PROJECT_PATH="/path/to/your/project"
```

Replace this path with the absolute path to your local project directory.

---

## Create the Python environment

Run:

```bash
make setup
```

This command will:

* create the Python virtual environment (`envh5Modif`)
* upgrade `pip`
* install the required Python packages (`numpy` and `h5py`)

---

## Activate the virtual environment

After the installation is complete, activate the environment manually:

```bash
source envh5Modif/bin/activate
```

You should now see something similar to:

```text
(envh5Modif)
```

at the beginning of your terminal prompt.

---

## Run the program

Once the environment is activated, execute:

```bash
python RUN_main.py
```

The program will then guide you through the different modification steps.

---

# Dependencies

The project requires:

* Python 3
* NumPy
* h5py

These dependencies are automatically installed by `make setup`.

---

# Notes

* The virtual environment only needs to be created once.
* Each time you open a new terminal, remember to activate it using:

```bash
source envh5Modif/bin/activate
```

before running:

```bash
python RUN_main.py
```

---

# Typical Workflow

```bash
# Clone the repository
git clone <repository_url>

# Enter the project directory
cd <repository_name>

# Edit setup.sh and set PROJECT_PATH

# Install dependencies
make setup

# Activate the virtual environment
source envh5Modif/bin/activate

# Run the program
python RUN_main.py
```

