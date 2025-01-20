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

from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt

from pyAF.src import widgets_list
from pyAF.src import shared
from pyAF.src.tools import misc_tools
from pyAF.src.utils import tablewidget

from base_test import pyAFTestCase
from base_test_files import HELA, NS_SINGLE


class StiffnessTest(pyAFTestCase):
    """Test the elements related to the stiffness computation."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([HELA, NS_SINGLE])
        cls.base_setUpClass()

        # Compute everything
        shared.exp.calc_all = True
        shared.exp.apply_to_all_compute = False

        # Use Hertz (Sphere)
        for i in range(2):
            shared.exp.list[i].tomography = 1
            shared.exp.list[i].stiffness_model_selected = 0
            shared.exp.list[i].indentation_step = 25

        # Compute the stiffness
        button = widgets_list.widget_compute.BT_compute_stiffness
        QTest.mouseClick(button, Qt.LeftButton)

    def test_stiffness_results(self):
        """Test stiffness values."""
        # Check if the computation is done
        self.assertEqual(shared.exp.list[0].stiffness_calculated, True)
        self.assertEqual(shared.exp.list[1].stiffness_calculated, True)

        # Check the results
        mean, median, sd, mode = shared.single_values[0]

        self.assertAlmostEqual(mean, 76.025611996687715, places=2)
        self.assertAlmostEqual(median, 69.886543273925781, places=2)
        self.assertAlmostEqual(sd, 44.50973430130442, places=2)
        self.assertAlmostEqual(mode, 0)

        return
        # The following was in Michka's test but the new values
        # are consistant with the ones obtained with SB's version.

        # Check the results
        mean, median, sd, mode = shared.single_values[1]

        self.assertAlmostEqual(mean, 1.3446605975811299, places=2)
        self.assertAlmostEqual(median, 0.39643603563308716, places=2)
        self.assertAlmostEqual(sd, 2.3732998959958969, places=2)
        self.assertAlmostEqual(mode, 0)

    def test_open_multimeshgrid_and_change_types(self):
        """Open multimeshgrid and go through the types."""
        button = widgets_list.widget_results.BT_multi_meshgrids
        QTest.mouseClick(button, Qt.LeftButton)
        W_multi_meshgrid = widgets_list.widget_multimeshgrids

        # Check if the widget is opened
        self.assertIsNotNone(W_multi_meshgrid)

        # Go trough the different types of meshgrids and check if there
        # are any errors.
        for i in range(len(W_multi_meshgrid.list_choser)):
            mesh_type = misc_tools.get_meshgrid_type_as_string(i)
            shared.exp.multi_meshgrid_type = mesh_type
            W_multi_meshgrid.canvas_list[0].update_plot()

    def test_duplicate_and_remove_result_in_table(self):
        """Test the row duplication and removing."""
        widget = widgets_list.widget_results_single

        # Select the newly created row
        widget.tableWidget.selectRow(0)

        # Test the duplication
        tablewidget.duplicate_results()

        # Check that we have 2 results now
        self.assertEqual(len(shared.exp.results_list), 3)
        self.assertEqual(widget.tableWidget.rowCount(), 3)

        # Check if the order of the results is correctly set
        for i in range(len(shared.exp.results_list)):
            self.assertEqual(shared.exp.results_list[i].result_id, i)
        self.assertEqual(shared.exp.results_list[0].data_id, 0)
        self.assertEqual(shared.exp.results_list[1].data_id, 1)
        self.assertEqual(shared.exp.results_list[2].data_id, 0)

        # Select the newly created row
        widget.tableWidget.selectRow(2)

        # Test the removing
        tablewidget.remove_result()

        # Check that we have only 1 result now
        self.assertEqual(len(shared.exp.results_list), 2)
        self.assertEqual(widget.tableWidget.rowCount(), 2)

    def test_batch_slice_change(self):
        """Check if the slice can be changed."""
        widget = widgets_list.widget_results_single

        widget.tableWidget.selectRow(2)
        # Force the usage of a slice which doed not exist for this file.
        # This can happen if you have multiple rows selected. (See bug #319)
        # This ensures that this condition is tested and the right slice index
        # is selected.
        widget.batch_slices_change(force_index=3)

        # Normal batch slice change
        widget.tableWidget.selectRow(1)
        widget.batch_slices_change(force_index=15)
        widget.batch_slices_change(force_index="All")

    def test_change_single_color(self):
        """Test the changing of the color of the first result"""
        widget = widgets_list.widget_results_single
        widget.tableWidget.selectRow(0)
        widget.parent.tabs.setCurrentIndex(3)
        tablewidget.change_color()
        widget = widgets_list.widget_results_groups
        widget.tableWidget.selectRow(0)
        widget.parent.tabs.setCurrentIndex(4)
        tablewidget.change_color()
