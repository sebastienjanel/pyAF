# pyAF : An AFM Force Curve Analysis Software Dedicated to Biologists

pyAF is a Python software to analyze force curves & force maps acquired by atomic force microscopes (up to now Bruker and JPK).
Developed at CMPI lab (https://cmpi.cnrs.fr).

It computes:
- Elasticity with several models and bottom effect correction (BEC)
- Stiffness tomography
- Viscoelasticity (Ting, Garcia)
- Microrheology (Alcaraz)
- Work of detachment and rupture force
- Unbinding events and loading rates (SMFS)
- Zero-force topography, flattening, plane fitting
- Regions of interest for substrate definition and stats
- Curve tilt correction and smoothing
- 3D map rendering
- Text and JSON export
- Full statistical analysis: violin plots and superplots
- Multi files and multiprocessing

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info".

The plots use the Inter font (https://rsms.me/inter), shipped in src/fonts under the SIL Open Font License 1.1 (src/fonts/OFL.txt).

Main initial coding made by Michka Popoff (https://theses.fr/api/v1/document/2014LIL10220).
Modified by Antoine Dujardin, Simone Bovio, Javier Lopez-Alonso, Nuno Duarte, Sébastien Janel.

When using the program, cite the original article : Stiffness tomography of eukaryotic intracellular compartments by atomic force microscopy (Nanoscale, 2019,11, 10320-10328).
DOI	https://doi.org/10.1039/C8NR08955H

DISCLAIMER: pyAF is a lab-developed software, it is not a commercial software.
Check as much as possible results provided vs other software (e.g. JPK DP) to rule out any no supported edge cases! And please report to us any wrong results.

## Installation from the sources

pyAF runs on macOS, Windows and Linux with Python 3.11 or newer (3.13 recommended).

### 1. Install Python and Git

- **Python**: download it from https://www.python.org/downloads/ (on Windows, tick "Add python.exe to PATH" in the installer). On Linux, use the Python of your distribution (with the `venv` module, e.g. `sudo apt install python3-venv` on Ubuntu/Debian).
- **Git**: https://git-scm.com/downloads (already installed on most Linux systems; on macOS, `git` proposes to install it at its first use).

### 2. Get the code

In a terminal (Terminal on macOS, PowerShell on Windows):

    git clone https://github.com/sebastienjanel/pyAF.git
    cd pyAF

The folder of the code must be named `pyAF` (default name of the clone): pyAF imports itself under this name.

### 3. Create a virtual environment and install the packages

macOS and Linux:

    python3 -m venv .venv
    source .venv/bin/activate
    pip install --upgrade pip
    pip install -r requirements.txt

Windows (PowerShell):

    py -m venv .venv
    .venv\Scripts\Activate.ps1
    python -m pip install --upgrade pip
    pip install -r requirements.txt

If PowerShell refuses to run the activation script, run once `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`.

The packages take about 3 GB (Qt, VTK and ITK for the 3D rendering).

### 4. Run pyAF

From the `pyAF` folder, with the virtual environment activated:

    python main.py

Next times: open a terminal, go to the `pyAF` folder, activate the environment (step 3, second line) and run `python main.py`.

### Update

    git pull
    pip install -r requirements.txt

### With PyCharm (free for academics)

- File / Project from Version Control: https://github.com/sebastienjanel/pyAF.git
- Settings / Project / Python Interpreter: add a new virtual environment with Python 3.11 or newer
- in the terminal of PyCharm: `pip install -r requirements.txt`
- Run / Edit Configurations: script `main.py` of the project, with the interpreter of the virtual environment, then Run

### Known issues

- Linux: if Qt does not start ("Could not load the Qt platform plugin xcb"), install the library missing for Qt 6: `sudo apt install libxcb-cursor0` (Ubuntu/Debian).
- `ModuleNotFoundError: No module named 'pyAF'`: the folder of the code is not named `pyAF`: rename it.

Please report bugs (https://github.com/sebastienjanel/pyAF/issues).
