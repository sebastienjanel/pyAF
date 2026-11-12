# pyAF : An AFM Force Curve Analysis Software Dedicated to Biologists

pyAF is a Python software to analyze force curves & force maps acquired by atomic force microscopes (up to now Bruker and JPK, including JPK QI and Smart Mapping files).
Developed at CMPI lab (https://cmpi.cnrs.fr).

It computes:
- Elasticity with several models and bottom effect correction (BEC)
- Stiffness tomography
- Viscoelasticity (Ting, Garcia)
- Microrheology (Alcaraz)
- Work of detachment and rupture force
- Unbinding events and loading rates (SMFS)
- Plant cells: cell wall modulus and turgor pressure (Yamasaki, Beauzamy)
- Zero-force topography, flattening (with the z motor steps of Smart Mapping), plane fitting
- Regions of interest (ROI manager) for substrate definition, cell geometry and stats
- Curve tilt correction and smoothing
- 3D map rendering
- Text and JSON export
- Full statistical analysis: violin plots and superplots
- Multi files and multiprocessing

The formulas and references of all the models are in [docs/models.md](docs/models.md).

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info".

The plots use the Inter font (https://rsms.me/inter), shipped in src/fonts under the SIL Open Font License 1.1 (src/fonts/OFL.txt).

Main initial coding made by Michka Popoff (https://theses.fr/api/v1/document/2014LIL10220).
Modified by Antoine Dujardin, Simone Bovio, Javier Lopez-Alonso, Nuno Duarte, Sébastien Janel.

When using the program, cite the original article : Stiffness tomography of eukaryotic intracellular compartments by atomic force microscopy (Nanoscale, 2019,11, 10320-10328).
DOI	https://doi.org/10.1039/C8NR08955H

DISCLAIMER: pyAF is a lab-developed software, it is not a commercial software.
Check as much as possible results provided vs other software (e.g. JPK DP) to rule out any no supported edge cases! And please report to us any wrong results.

## Applications (no Python needed)

Ready-to-use applications for macOS (Apple Silicon and Intel), Windows and Linux are attached to the releases: https://github.com/sebastienjanel/pyAF/releases (download the file of your system).

Minimum systems (64 bits): macOS 14, Windows 10, Ubuntu 22.04 or another Linux with glibc 2.35. Older systems (e.g. Windows 7 or Ubuntu 18.04 on some AFM control computers) cannot run pyAF: analyze the data on another computer.

- **macOS** (`.dmg`, `AppleSilicon` or `Intel` according to your Mac: Apple menu / About This Mac): open it and drag pyAF to Applications. The application is not signed by Apple: at the first opening, macOS refuses to open it. Open System Settings / Privacy & Security, then click "Open Anyway" in front of pyAF (or, in a terminal: `xattr -dr com.apple.quarantine /Applications/pyAF.app`).
- **Windows** (`.zip`): unzip it (e.g. in Documents) and run `pyAF.exe` in the `pyAF` folder. The application is not signed: if Windows SmartScreen shows a warning, click "More info", then "Run anyway".
- **Linux** (`.tar.xz`): extract it where you want to keep it (`tar -xJf pyAF-*.tar.xz`), then add pyAF to the applications menu once with `pyAF/install.sh` (no administrator rights needed; run it again if you move the folder, `pyAF/install.sh --uninstall` removes the launcher). You can also run `pyAF/pyAF` directly. Qt may need `libxcb-cursor0` (see Known issues).

The applications take about 0.5 GB once installed (Qt and VTK for the 3D view) and a few seconds to start.

Projects (.pyaf) saved with versions before pyAF 5 cannot be opened: load the AFM files again and compute them.

To check an installation (applications or sources), run pyAF with the environment variable `PYAF_SMOKE_TEST=1`: it loads its libraries and resources, prints the results and quits.

## Installation from the sources

pyAF runs on macOS, Windows and Linux with Python 3.11 or newer (3.13 recommended), on the same minimum systems as the applications.

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

The packages take about 2 GB (Qt and VTK for the 3D view).

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
