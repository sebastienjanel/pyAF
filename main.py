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

import multiprocessing
from platform import system, release
import sys
import os
import logging

if getattr(sys, 'frozen', False):
    # Running in a bundle
    base_path = sys._MEIPASS
else:
    # Running in normal Python
    base_path = os.path.abspath(os.path.dirname(__file__))

# We need to add the parent folder of the inner "pyAF" package.
parent_path = os.path.abspath(os.path.join(base_path, '..'))
if parent_path not in sys.path:
    sys.path.insert(0, parent_path)

# Now import using the absolute package path:
from pyAF.src.pyaf import main as pyaf_main

if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)
    logger = logging.getLogger(__name__)

    logger.debug("Starting main.py")

    # If running macOS catalina use forkserver instead of fork for creating new processes.
    # See https://bugs.python.org/issue33725
    if system() == "Darwin" and "19.0" <= release()[:-2] < "20.0":
        multiprocessing.set_start_method("forkserver")

    # Add support for multiprocessing in frozen app
    # See http://docs.python.org/3/library/multiprocessing.html
    multiprocessing.freeze_support()

    logger.debug("Calling pyaf_main")
    pyaf_main()
    logger.debug("Finished pyaf_main")
