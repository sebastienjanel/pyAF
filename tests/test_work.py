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
from base_test_files import NS130226001


class StiffnessTest(pyAFTestCase):
    """Test the elements related to the stiffness computation."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([NS130226001])
        cls.base_setUpClass()

    def test_work_compute_failure(self):
        """Check the work and rupture force computation and display."""
        # Change parameters
        data = shared.exp.current_data
        data.fitparam_joc_skip_start = 0
        data.fitparam_joc_fit_length = 100.0
        data.fitparam_joc_noise_multiplicator = 1.0

        # Go the work and rupture force tab to display it
        widgets_list.widget_compute.tabs.setCurrentIndex(2)

        # Compute the work and rupture force
        button = widgets_list.widget_compute.BT_compute_work_and_rupture_force
        QTest.mouseClick(button, Qt.LeftButton)
