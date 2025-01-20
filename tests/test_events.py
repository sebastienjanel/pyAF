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

from base_test import pyAFTestCase
from base_test_files import NS130226001, NS130226002, NS130226003


class EventTest(pyAFTestCase):
    """Test the elements related to the events computation."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([NS130226001, NS130226002, NS130226003])
        cls.base_setUpClass()

    def test_events_per_scan(self):
        """Test events per scan computing and widget."""
        # Compute all the events
        shared.exp.calc_all = True
        #shared.exp.apply_to_all_compute = False
        button = widgets_list.widget_compute.BT_compute_events
        QTest.mouseClick(button, Qt.LeftButton)

        # Display the first events on the plots
        widgets_list.widget_results_single.tableWidget.selectRow(0)
        widgets_list.widget_results_single.display_selected()
        widgets_list.widget_results_single.button_clicked("refresh_hist")

        button = widgets_list.widget_results_single.BT_events_per_scan
        QTest.mouseClick(button, Qt.LeftButton)
