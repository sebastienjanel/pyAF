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
from src.tools import apply_to_all
from .base_test_tools import dlfile




class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        """
        Load the same files twice and compute the stiffness.

        """

        reload(shared)
        reload(widgets_list)

        cls.folder = misc_tools.get_app_path() + "/../tests/files/"

        dlfile("https://bitbucket.org/cmip/pyaf-files/src/master/tests/hela.003")

        file1 = cls.folder + "hela.003"
        file2 = cls.folder + "hela.003"

        # Load the default files
        file1 = {"version": None, "checked": True, "error": "", "file_type": "Nanoscope (Force Volume)", "path": file1, "nbrcurves": 256, "enabled": True, "filename": "hela.003"}
        file2 = {"version": None, "checked": True, "error": "", "file_type": "Nanoscope (Force Volume)", "path": file2, "nbrcurves": 256, "enabled": True, "filename": "hela.003"}

        files_info = []
        files_info.append(file1)
        files_info.append(file2)

        # Disable some qmessageboxes in the GUI
        consts.UNIT_TESTING = True

        cls.app = QApplication(sys.argv)

        cls.form = MainWindow(files_info)
        cls.form.show()

        # Compute the first file with an indentation step of 0
        # and the second with a step of 100
        shared.exp.calc_all = True
        shared.exp.apply_to_all_compute = False
        shared.exp.list[0].indentation_step = 100
        shared.exp.list[1].indentation_step = 0
        # Compute the stiffness
        button = widgets_list.widget_compute.BT_compute_stiffness
        QTest.mouseClick(button, Qt.LeftButton)


    @classmethod
    def tearDownClass(cls):

        cls.app.quit()


    def test_slice_meshgrid(self):

        """
        Go to first file, select slice meshgrid, then change to second file.

        Check that the meshgrid type changes to "stiffness"

        """

        shared.exp.list[0].meshgrid_type = "stiffness_slice"
        apply_to_all.apply_to_all("all")
        widgets_list.widget_results.update_widget()

        # Go to the second file
        shared.exp.id_selected = 1
        self.form.file_changed("box")

        # Check that the type was set to stiffness, as stiffness_slice is not
        # possible for the second file as the indentation step is set to 0.
        self.assertEqual(shared.exp.list[1].meshgrid_type, "stiffness")
