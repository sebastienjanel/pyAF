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

from pyAF.src import shared
from pyAF.src import widgets_list

from base_test import pyAFTestCase
from base_test_files import HELA, NS_SINGLE
from base_test_tools import with_pyplot_interaction


class DataTabTest(pyAFTestCase):
    """Test the elements related to the data tab."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([HELA, NS_SINGLE])
        cls.base_setUpClass()

    def test_change_curve(self):
        """Go trough some curves to see if the change_curve works."""
        # Change curve
        self.form.change_curve(12, 7)
        # Check if it changed as planned
        # Array indexes are 0-based, change_curve is 1-based
        self.assertEqual(shared.exp.meshgrid_click_xpos, 11)
        self.assertEqual(shared.exp.meshgrid_click_ypos, 6)
        # Change at random
        self.form.change_curve(0, 0, rand=True)
        # Check that the change happened
        self.assertNotEqual(shared.exp.meshgrid_click_xpos, 11)
        self.assertNotEqual(shared.exp.meshgrid_click_ypos, 6)

    def test_change_file(self):
        """Go back and forth between the two files."""
        shared.exp.id_selected = 1
        self.form.file_changed("box")
        self.assertEqual(
            widgets_list.widget_compute.list_exp2.currentIndex(), 1)
        self.assertEqual(
            widgets_list.widget_results.list_exp.currentIndex(), 1)

        shared.exp.id_selected = 0
        self.form.file_changed("box")
        self.assertEqual(
            widgets_list.widget_compute.list_exp2.currentIndex(), 0)
        self.assertEqual(
            widgets_list.widget_results.list_exp.currentIndex(), 0)

    @with_pyplot_interaction
    def test_open_meshgrid(self):
        """Open meshgrid on data tab."""
        action = widgets_list.widget_main.list_actions["open_meshgrid"]
        action.trigger()

    @with_pyplot_interaction
    def test_open_curve(self):
        """Open curve on data tab."""
        action = widgets_list.widget_main.list_actions["open_curve"]
        action.trigger()

    def test_rename(self):
        """Test the rename method."""
        # Update gui for second file
        widgets_list.widget_main.update_names(force_id=1)
        # Update all names
        widgets_list.widget_main.update_names()
