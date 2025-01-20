# Copyright Antoine Dujardin (2016-2018) toine.dujardin@gmail.com
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

"""Base Class for pyAF-related unit tests.

See http://docs.python.org/2/library/unittest.html for more details
on unit testing.
"""

import faulthandler
import numpy
import sys
import unittest
from PyQt5 import QtWidgets

from pyAF.src import consts
from pyAF.src import pyaf

from base_test_tools import dlfile

# During unit testing, fetch all numpy warnings as actual errors.
# They should fixed, we don't want unplanned stuff to happen.
numpy.seterr(all="raise")
faulthandler.enable()


class pyAFTestCase(unittest.TestCase):
    """Class for pyAF-based test cases.

    Should be inherited instead of unittest.TestCase.
    Start pyAF with the requested files before performing the
    tests.
    """
    _requested_files_info = []  # File info of requested files

    @classmethod
    def base_setUpClass(cls):
        """Start pyAF with the requested files.

        Must be called in setUpClass, after requesting the files.
        """
        # Disable some qmessageboxes in the GUI
        consts.UNIT_TESTING = True

        # Start application
        cls.app = QtWidgets.QApplication(sys.argv)
        cls.form = pyaf.MainWindow(cls._requested_files_info)
        cls.form.show()

    @classmethod
    def use_files(cls, file_infos):
        for file_info in file_infos:
            dlfile(file_info["filename"])
        cls._requested_files_info = file_infos

    @classmethod
    def tearDownClass(cls):
        cls.app.quit()
