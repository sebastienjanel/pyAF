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

from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest

from pyAF.src import shared
from pyAF.src import widgets_list

from base_test import pyAFTestCase
from base_test_files import NONSQR1, NONSQR2


class NonsquareTest(pyAFTestCase):
    """Test the elements related to non-square files."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([NONSQR1, NONSQR2])
        cls.base_setUpClass()

    def test_open_flatten_non_square(self):
        """Tests the flatten widget with non square map.

        See error #480
        """
        button = widgets_list.widget_compute.BT_flatten
        QTest.mouseClick(button, Qt.LeftButton)

    def test_compute_non_square(self):
        """Test the computation of the stiffness on non square maps."""
        shared.exp.calc_all = True
        button = widgets_list.widget_compute.BT_compute_stiffness
        QTest.mouseClick(button, Qt.LeftButton)

    def test_approach_velocity_QI(self):
        """Check if the approach velocity is correctly extracted."""
        val = shared.exp.current_data.approach_velocity
        self.assertAlmostEqual(val, 30000)

    def test_retraction_velocity_QI(self):
        """Check if the retraction velocity is correctly extracted."""
        val = shared.exp.current_data.retraction_velocity
        self.assertAlmostEqual(val, 30000)

    def test_scan_rate_QI(self):
        """Check if the scan rate is correctly extracted."""
        self.assertAlmostEqual(shared.exp.current_data.scan_rate, 7.5)
