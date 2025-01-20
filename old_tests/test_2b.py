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
import warnings
from PyQt5.QtGui import QApplication, QtWidgets
from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt
from src.pyaf import MainWindow
from src import widgets_list
from src import shared
from src import consts
from src.tools import misc_tools
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

        # Disable some qmessageboxes in the GUI
        consts.UNIT_TESTING = True

        # Reopend to test the loading
        # Note : load first a single curve before loading the pyaf file; will test
        # the load more procedure. (And test for bug #333)

        file1 = cls.folder + "43_curve1.jpk-force"
        file2 = cls.folder + "computed2.pyaf"

        # Load the default file
        file1 = {"version": None, "checked": True, "error": "", "file_type": "JPK (Single File)", "path": file1, "nbrcurves": 1, "enabled": True, "filename": "43_curve1.jpk-force"}
        file2 = {"version": None, "checked": True, "error": "", "file_type": "pyAF", "path": file2, "nbrcurves": 0, "enabled": True, "filename": "computed2.pyaf"}

        files_info = []
        files_info.append(file1)
        files_info.append(file2)

        # Disable some qmessageboxes in the GUI
        consts.UNIT_TESTING = True

        cls.app = QApplication(sys.argv)
        cls.form = MainWindow(files_info)
        cls.form.show()

    @classmethod
    def tearDownClass(cls):
        cls.app.quit()

    def test_indentation_widget(self):

        """
        Test the indentation widget (change file)

        """

        # Open indentation widget
        meshgrid_menu.open_indentation_widget(widgets_list.widget_results)

        # Change to the second file
        shared.exp.id_selected = 1
        self.form.file_changed(option = "box")

    def test_display_surface_in_defl_ext(self):

        """
        Display the surface first in defl ext, with no joc displayed.

        """

        shared.exp.id_selected = 1
        self.form.file_changed(option = "box")

        widgets_list.widget_results.RB_defl_ext.setChecked(True)
        widgets_list.widget_results.button_clicked("button_group_type_of_curve")

        widgets_list.widget_results.CB_surface.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_surface")

    def test_events_results(self):

        """
        Test the results of the events computation.

        """

        # Check if the computation is done
        self.assertEqual(shared.exp.list[1].events_calculated, True)

        widgets_list.widget_results_chooser_single.list_choser.setCurrentIndex(3)
        widgets_list.widget_results_chooser_single.list_updated(None)

        # Check the results
        mean, median, sd, mode = shared.single_values[0]

        self.assertAlmostEqual(mean, 59.933037, places=2)
        self.assertAlmostEqual(median, 53.086143, places=2)
        self.assertAlmostEqual(sd, 38.004951, places=2)
        self.assertEqual(mode, 0)

    def test_lr_calculated(self):

        """
        Check if the loading rates were calculated.

        """

        self.assertEqual(shared.exp.list[1].loading_rates_calculated, True)

    def test_lr_values(self):

        """
        Check if the values of the loading rates are right.

        """

        # Go to loading rates results
        widg = widgets_list.widget_results_chooser_single
        widg.list_choser.setCurrentIndex(6)
        widg.list_updated(None)

        # Check the results
        mean, median, sd, mode = shared.single_values[0]

        self.assertAlmostEqual(mean, 8035.5366, places=2)
        self.assertAlmostEqual(median, 4959.0479, places=2)
        self.assertAlmostEqual(sd, 12155.03, places=2)
        self.assertEqual(mode, 0)

    def test_display_lr_scatter(self):

        """
        Test the scatter plot for the loading rates

        """

        # Go to loading rates results
        widg = widgets_list.widget_results_chooser_single
        widg.list_choser.setCurrentIndex(6)
        with warnings.catch_warnings():
            # Catch matplotlib warning for scatter plots:
            # UserWarning: Unable to find pixel distance along axis for
            # interval padding; assuming no interval padding needed.
            warnings.simplefilter("ignore")
            widg.list_updated(None)

        cb = widgets_list.widget_results_single.tableWidget.cellWidget(0, 1)
        cb.checkbox.setChecked(True)

        # Refresh
        button = widgets_list.widget_results_single.BT_refresh_hist
        QTest.mouseClick(button, Qt.LeftButton)

    def test_load_more_on_force_curve(self):

        """
        Displays the force curve and the fits, then loads more.

        """

        # Change back to the second file
        shared.exp.id_selected = 1
        self.form.file_changed(option = "box")

        # Display the options on the force curve
        widgets_list.widget_results.RB_force.setChecked(True)
        widgets_list.widget_results.button_clicked("button_group_type_of_curve")
        widgets_list.widget_results.CB_poc.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_poc")
        widgets_list.widget_results.CB_fit_poc.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_fit_poc")
        widgets_list.widget_results.CB_segments.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_segments")
        widgets_list.widget_results.CB_display_fits.setChecked(True)
        widgets_list.widget_results.button_clicked("display_fits_stiffness")

        widgets_list.widget_results.CB_joc.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_joc")
        widgets_list.widget_results.CB_fit_joc.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_fit_joc")
        widgets_list.widget_results.CB_surface.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_surface")
        widgets_list.widget_results.CB_force.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_force")

        widgets_list.widget_results.CB_events_filter_dist.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_events_filter_dist")
        widgets_list.widget_results.CB_events_display_joc.setChecked(True)
        widgets_list.widget_results.button_clicked("checkbox_events_display_joc")

        # Load more
        self.form.load_more(path = "../tests/files/hela.003")

    def test_get_pdf_rupture_force2(self):

        """
        Check for the pdf computation on the events rupture force.

        Bug #378

        """

        # Go to events rupture force
        widg = widgets_list.widget_results_chooser_single
        widg.list_choser.setCurrentIndex(5)
        widg.list_updated(None)

        # Get pdf
        button = widgets_list.widget_results_single.BT_get_pdf
        QTest.mouseClick(button, Qt.LeftButton)

    def test_display_events_force_curve(self):

        """
        Displays the fuzzy force curve.

        """

        shared.exp.id_selected = 1
        self.form.file_changed(option = "box")

        widgets_list.widget_results.RB_force_events.setChecked(True)
        widgets_list.widget_results.button_clicked(
            "button_group_type_of_curve")

    def test_recompute_loading_rates(self):
        """
        Computes the loading rates two times.

        """

        shared.exp.id_selected = 0
        self.form.file_changed(option = "box")

        button = widgets_list.widget_compute.BT_compute_events
        QTest.mouseClick(button, Qt.LeftButton)

        button = widgets_list.widget_compute.BT_get_lr
        QTest.mouseClick(button, Qt.LeftButton)

        # Check the results
        mean, median, sd, mode = shared.single_values[1]

        self.assertAlmostEqual(mean, 35280.195, places=2)
        self.assertAlmostEqual(median, 8155.2559, places=2)
        self.assertAlmostEqual(sd, 49155.527, places=2)
        self.assertEqual(mode, 0)

        # Change the loading rate coefficient
        shared.exp.current_data.lr_coef = 0.5

        button = widgets_list.widget_compute.BT_get_lr
        QTest.mouseClick(button, Qt.LeftButton)

        # Check the results
        mean, median, sd, mode = shared.single_values[1]

        self.assertAlmostEqual(mean, 10356.055, places=2)
        self.assertAlmostEqual(median, 5975.0981, places=2)
        self.assertAlmostEqual(sd, 7004.5493, places=2)
        self.assertEqual(mode, 0)

if __name__ == "__main__":
    unittest.main()
