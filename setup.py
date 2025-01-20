"""
Python application for fitting Hertz model to different force curves
"""

import sys
import subprocess
import traceback

from setuptools import setup

CURRENT_PYTHON = sys.version_info[:2]
REQUIRED_PYTHON = (3, 6)

# This check and everything above must remain compatible with Python 2.7.
if CURRENT_PYTHON < REQUIRED_PYTHON:
    sys.stderr.write("""
    =============================
    Unsupported Python version!!!
    =============================
    This version of PyAF requires Python {}.{}, but you're trying to
    install it on Python {}.{}.
    
    """.format(*(REQUIRED_PYTHON + CURRENT_PYTHON)))
    sys.exit(1)

try:

    sys.stderr.write("""
    ==========================
    SETUP PHASE
    ==========================
    The following elements are required:
    - create-dmg --> JavaScript command line program to generate .dmg files from .app files.
                     Documentation: https://github.com/sindresorhus/create-dmg
                     Requirements: Node.js 8 or later.
                     
                     Please run $npm install --global create-dmg to install create-dmg
                     
    - Python modules stated in the requirements.txt
    
    - PyInstaller development version --> Documentation: https://www.pyinstaller.org/index.html
    
    """)

    while True:
        user_input = input("Continue with python modules installation? (y/n)")

        if user_input.lower() == 'y':

            # Install packages dependencies.
            subprocess.run("pip  install -r requirements.txt", shell=True, check=True)

            # Install PyInstaller 4.0 development version. At this time there are issues with previous versions.
            subprocess.run("pip install https://github.com/pyinstaller/pyinstaller/archive/develop.tar.gz",
                           shell=True, check=True)

            sys.stderr.write("""
            ================================
                    Setup Complete!!        
            ================================
            
            """)

            break

        elif user_input.lower() == 'n':
            break


except:

    sys.stderr.write("""
    ============================================
    AN ERROR DURING THE INSTALLATION OCCURRED!!!
    ============================================
    """)

    traceback.print_exc()

    sys.exit(1)




