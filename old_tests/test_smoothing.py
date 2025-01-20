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

import sys
import unittest
from PyQt5.QtGui import QApplication, QtWidgets
from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt
from src.pyaf import MainWindow
from src import widgets_list
from src import shared
from src import consts
from src.tools import misc_tools
from .base_test_tools import dlfile




class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        """
        Load a file only once for this test.

        """

        reload(shared)
        reload(widgets_list)

        cls.folder = misc_tools.get_app_path() + "/../tests/files/"

        dlfile("https://bitbucket.org/cmip/pyaf-files/src/master/tests/single_curve.001")

        file1 = cls.folder + "single_curve.001"

        # Load the default files
        file1 = {"version": None, "checked": True, "error": "", "file_type": "Nanoscope (Single File)", "path": file1, "nbrcurves": 1, "enabled": True, "filename": "single_curve.001"}

        files_info = []
        files_info.append(file1)

        # Disable some qmessageboxes in the GUI
        consts.UNIT_TESTING = True

        cls.app = QApplication(sys.argv)

        cls.form = MainWindow(files_info)
        cls.form.show()


    @classmethod
    def tearDownClass(cls):

        cls.app.quit()


    def test_smoothing(self):

        """
        Compute stiffness with smoothing on and check the results.

        Enable smoothing. If the smoothing is not working correctly we will
        get different values for the stiffness.

        """

        # Enable smoothing
        shared.exp.list[0].sg_smoothing_enabled = True

        # Change fitting parameter for this file
        shared.exp.list[0].fitparam_poc_skip_start = 1000.0
        shared.exp.list[0].fitparam_poc_refit_option = 200.0

        # Compute the stiffness
        button = widgets_list.widget_compute.BT_compute_stiffness
        QTest.mouseClick(button, Qt.LeftButton)

        # Check the results
        mean, median, sd, mode = shared.single_values[0]

        self.assertAlmostEqual(mean, 0.24366447, places = 5)
        self.assertAlmostEqual(median, 0.24366447, places = 5)
        self.assertAlmostEqual(sd, 0.0, places = 5)
        self.assertAlmostEqual(mode, 0)
