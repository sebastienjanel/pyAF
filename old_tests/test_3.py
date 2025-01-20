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

import os
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
from src.utils import tablewidget
from src.utils import meshgrid_menu
from src.load_and_save.save import Save
from .base_test_tools import dlfile




class Test(unittest.TestCase):

    """
    See http://docs.python.org/2/library/unittest.html for more details.

    """

    @classmethod
    def setUpClass(cls):

        """
        Load a file only once for this test.

        """

        reload(shared)
        reload(widgets_list)

        cls.folder = misc_tools.get_app_path() + "/../tests/files/"

        dlfile("https://bitbucket.org/cmip/pyaf-files/src/master/tests/hela.003")
        dlfile("https://bitbucket.org/cmip/pyaf-files/src/master/tests/map16ums.jpk-force-map")

        file1 = cls.folder + "hela.003"
        file2 = cls.folder + "map16ums.jpk-force-map"

        # Load the default files
        file1 = {"version": None, "checked": True, "error": "", "file_type": "Nanoscope (Force Volume)", "path": file1, "nbrcurves": 256, "enabled": True, "filename": "hela.003"}
        file2 = {"version": None, "checked": True, "error": "", "file_type": "JPK (Force Map)", "path": file2, "nbrcurves": 1024, "enabled": True, "filename": "map16ums.jpk-force-map"}

        files_info = []
        files_info.append(file1)
        files_info.append(file2)

        # Disable some qmessageboxes in the GUI
        consts.UNIT_TESTING = True

        cls.app = QApplication(sys.argv)

        cls.form = MainWindow(files_info)
        cls.form.show()

        # Compute just the first one
        button = widgets_list.widget_compute.BT_compute_stiffness

        shared.exp.list[0].indentation_start = 0
        shared.exp.list[0].indentation_stop = 0
        shared.exp.list[0].indentation_step = 0

        QTest.mouseClick(button, Qt.LeftButton)


    @classmethod
    def tearDownClass(cls):

        cls.app.quit()



    def test_display_force_curve_with_options(self):

        """
        Display force curve with all the options

        """

        shared.exp.list[0].display_curve_type = "force"
        shared.exp.list[0].display_trace_retrace = 0
        shared.exp.list[0].display_poc = True
        shared.exp.list[0].display_fit_poc = True
        shared.exp.list[0].display_segments = True
        shared.exp.list[0].display_fits_stiffness = True

        widgets_list.widget_results.update_MPL("MPL_canvas2")


    def test_change_file_and_go_back_and_forward_buttons(self):

        """
        Test the go back and go forward buttons.

        Change to second file, go to the curve (31, 31), go back to first file
        and use the go back button. Nothing should happen because the button
        should be disabled.

        """

        shared.exp.id_selected = 1
        self.form.file_changed("box")
        self.form.change_curve(31, 31)
        shared.exp.id_selected = 0
        self.form.file_changed("box")

        # Nothing should happen here because the buttons are disabled
        BT_back = widgets_list.widget_curve_selector_data.BT_back
        BT_forward = widgets_list.widget_curve_selector_data.BT_forward
        # Check if the buttons are disabled and the parameters reset
        self.assertFalse(BT_back.isEnabled())
        self.assertFalse(BT_forward.isEnabled())
        self.assertEqual(shared.exp.pos_in_last_ten_curves, 0)
        self.assertEqual(shared.exp.last_ten_curves, [[0, 0]])

        # Change curve and test the buttons
        self.form.change_curve(2, 3)
        self.form.change_curve(5, 4)
        self.form.change_curve(3, 1)

        # This should change the curves
        QTest.mouseClick(BT_back, Qt.LeftButton)
        QTest.mouseClick(BT_back, Qt.LeftButton)
        QTest.mouseClick(BT_forward, Qt.LeftButton)
        QTest.mouseClick(BT_forward, Qt.LeftButton)


    def test_copy(self):

        widgets_list.widget_data.button_clicked("copy_file")


    def test_add_with_multimeshgrid_open(self):

        widgets_list.widget_results.button_clicked("button_multi_meshgrids")

        # Load more
        self.form.load_more(path = "../misc/tests/files/hela.003")


if __name__ == "__main__":
    unittest.main()
