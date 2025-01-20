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

from nose.plugins.skip import SkipTest
import os
from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt

from pyAF.src import widgets_list
from pyAF.src import shared
from pyAF.src.tools import utils

from base_test import pyAFTestCase
from base_test_files import HELA


class VTKTest(pyAFTestCase):
    """Test the vtk-related elements."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([HELA])
        cls.base_setUpClass()

        # Use Hertz (Sphere)
        shared.exp.list[0].tomography = 1
        shared.exp.list[0].stiffness_model_selected = 0
        shared.exp.list[0].indentation_step = 25

        # Compute the stiffness
        button = widgets_list.widget_compute.BT_compute_stiffness
        QTest.mouseClick(button, Qt.LeftButton)

    def test_open_vtk(self):
        """Open VTK widget."""
        if not utils.module_exists("vtk") or os.uname()[4].startswith("arm"):
            # Skip this test if VTK is not installed
            # Skip this test on ARM, they have no drivers for this ...
            raise SkipTest

        button = widgets_list.widget_results.BT_display_in_3D
        QTest.mouseClick(button, Qt.LeftButton)

        # Check if the vtk widget is open
        self.assertIsNotNone(widgets_list.widget_vtk)

    def test_add_afm_file_to_vtk(self):
        """Test the add afm file function of the VTK Widget."""
        # raise SkipTest("Doesn't currently work.")
        if not utils.module_exists("vtk") or os.uname()[4].startswith("arm"):
            # Skip this test if VTK is not installed
            # Skip this test on ARM, they have no drivers for this ...
            raise SkipTest

        button = widgets_list.widget_results.BT_display_in_3D
        QTest.mouseClick(button, Qt.LeftButton)

        action = widgets_list.widget_main.list_actions["add_afm_layer"]
        action.trigger()
        self.fail()
