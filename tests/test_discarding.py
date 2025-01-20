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
from base_test_files import HELA


class StiffnessTest(pyAFTestCase):
    """Test the elements related to the stiffness computation."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([HELA])
        cls.base_setUpClass()

    def test_discarding(self):
        """Check the discarding of curves.

        Discard some curves and check the end result.
        """
        # Discard some curves
        shared.exp.list[0].discarded_curves[0][0] = 1
        shared.exp.list[0].discarded_curves[4][5] = 1
        shared.exp.list[0].discarded_curves[2][2] = 1
        shared.exp.list[0].discarded_curves[13][1] = 1
        shared.exp.list[0].discarded_curves[2][8] = 1

        # Compute the stiffness
        button = widgets_list.widget_compute.BT_compute_stiffness
        QTest.mouseClick(button, Qt.LeftButton)

        # Check the results
        mean, median, sd, mode = shared.single_values[0]

        self.assertAlmostEqual(mean, 19.802361, places=2)
        self.assertAlmostEqual(median, 16.76255, places=2)
        self.assertAlmostEqual(sd, 12.73009, places=2)
        self.assertAlmostEqual(mode, 0)
