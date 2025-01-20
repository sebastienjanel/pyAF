# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
#
# This software is a computer program whose purpose is to analyze force curves
# recorded with an atomic force microscope.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use, modify
# and/ or redistribute the software under the terms of the CeCILL license as
# circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

"""
The constants for PYAF are defined here.

In python you don't have global constants which can be shared through the
different modules. This file lets you define constants and use them as
**consts.CONSTANT** in your code.

During developement the version number contains the git sha to have a way
to know which version of the file we have in case of a bug (this version number
is written at the beggining of the saved file and displayed in the mainwindow's
title. For release versions or the master branch, VERSION should not contain
these informations but only a normal version number, like x.x.x
Note that the short_sha is only defined if we are in a git directory.

The other constants are useful when testing or debuggin PYAF. Often you may
want to relaunch the software with the same file, so just define a LOAD_FILE
constant and the TEST_PATH, and set AUTOLOAD to True.

"""

import os
from .tools import git_tools

MAINVERSION = "2.0.0"
branch, short_sha = git_tools.get_git_info()
if short_sha is None:
    VERSION = MAINVERSION
else:
    VERSION = MAINVERSION + "." + branch + "." + short_sha


DEBUG = False
DISABLE_ERROR_SENDER = False
AUTOLOAD = False
AUTOCALC_STIFFNESS = False
AUTOCALC_WORK = False
AUTOCALC_EVENTS = False
AUTOCALC_LOADING_RATE = False
AUTO3D = False
AUTO_LOAD_LAYER_IN_3D = False
AUTO_LOAD_TIFF_IN_3D = False
AUTO_LOAD_STACK_IN_3D = False
OPENLOGSENDER = False
UNIT_TESTING = False  # Is set to True by the unit_test
AUTOTAB = 1  # Select a specific tab (DEBUG needs to be set to True)
ADVANCED = True  # Special options, True on develop branch

# Some packages wich are difficult to install an not mandatory
ALLOW_VTK = True
ALLOW_ITK = True

# Test files
TEST_PATH = os.path.expanduser("~") + "/"
LOAD_FILE = None
# LOAD_FILE = "test.pyaf"

# Pyaf file Encoding
ENCODING = "latin1"
