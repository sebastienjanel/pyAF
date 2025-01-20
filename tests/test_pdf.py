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
from base_test_files import JPK_SINGLE


class PdfTest(pyAFTestCase):
    """Test the elements related to the pdf computation."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([JPK_SINGLE])
        cls.base_setUpClass()

    def test_pdf(self):
        """Test pdf computation.

        Use a single curve with very few values, this did trigger the error
        message. This test allows to see if the pdf computation in stat tools
        runs without error. See bug #250 for more.
        """
        # Defined a step
        shared.exp.list[0].indentation_step = 10
        shared.exp.list[0].fitparam_poc_noise_multiplicator = 1.0

        # Compute
        button = widgets_list.widget_compute.BT_compute_stiffness
        QTest.mouseClick(button, Qt.LeftButton)

        # Get pdf
        button = widgets_list.widget_results_single.BT_get_pdf
        QTest.mouseClick(button, Qt.LeftButton)
