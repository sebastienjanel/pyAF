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

"""Widget which allows to change the type of result to display."""


import logging
from ... import shared
from ... import widgets_list
from PyQt5 import QtWidgets
from ...tools.gui_tools import PYAFComboBox
from ...tools.PYAFWidget import PYAFWidget
from ...tools import misc_tools
from ...tools import gui_tools


class ResultsChooserWidget(PYAFWidget):
    """List of buttons to change the type of result to display."""

    def __init__(self, parent, name):
        super().__init__(parent, name)

        # Load the logger
        self.logger = logging.getLogger()

        self.VL = QtWidgets.QVBoxLayout()

        self.list_choser = PYAFComboBox(self, "list_choser")
        self.list_choser.addItem("Elasticity")
        self.list_choser.addItem("Detachment work")
        self.list_choser.addItem("Detachment force")
        self.list_choser.addItem("Events force")
        self.list_choser.addItem("Events per curve")
        self.list_choser.addItem("Event max. force")
        self.list_choser.addItem("Loading rates")
        self.list_choser.addItem("Distance to Joc")

        self.VL.addWidget(self.list_choser)

        self.setLayout(self.VL)

        self.update_widget()

    def update_widget(self):
        """Update the widget."""
        self.update_GUI("all")

    def list_updated(self, _):
        """Called when a list is updated."""
        # Update value
        val = self.list_choser.currentIndex()
        shared.exp.results_type = misc_tools.get_results_type_from_id(val)
        # Reset the columns which are displayed or hidden
        self.parent.tw_labels_is_checked = None
        # Refresh histograms and tables completely (will also refresh groups)
        widgets_list.widget_results_single.reset_table()
        widgets_list.widget_results_single.button_clicked("refresh_hist")

    def update_GUI(self, what):
        """Update the GUI."""
        if what == "button_type_hist" or what == "all":
            found_stiffness = False
            found_work = False
            found_events = False
            found_lr = False

            for item in shared.exp.list:
                if item.stiffness_array is not None:
                    found_stiffness = True
                if item.jocs1_indices is not None:
                    found_work = True
                if item.events_calculated:
                    found_events = True
                if item.loading_rates_calculated:
                    found_lr = True

            if found_stiffness:
                gui_tools.enable_combobox_item(self.list_choser, 0)
            else:
                gui_tools.disable_combobox_item(self.list_choser, 0)

            if found_work:
                gui_tools.enable_combobox_item(self.list_choser, 1)
                gui_tools.enable_combobox_item(self.list_choser, 2)
            else:
                gui_tools.disable_combobox_item(self.list_choser, 1)
                gui_tools.disable_combobox_item(self.list_choser, 2)

            if found_events:
                gui_tools.enable_combobox_item(self.list_choser, 3)
                gui_tools.enable_combobox_item(self.list_choser, 4)
                gui_tools.enable_combobox_item(self.list_choser, 5)
                gui_tools.enable_combobox_item(self.list_choser, 7)
            else:
                gui_tools.disable_combobox_item(self.list_choser, 3)
                gui_tools.disable_combobox_item(self.list_choser, 4)
                gui_tools.disable_combobox_item(self.list_choser, 5)
                gui_tools.disable_combobox_item(self.list_choser, 7)

            if found_lr:
                gui_tools.enable_combobox_item(self.list_choser, 6)
            else:
                gui_tools.disable_combobox_item(self.list_choser, 6)

            if shared.exp.results_type == "stiffness":
                self.list_choser.setCurrentIndex(0)
            elif shared.exp.results_type == "work":
                self.list_choser.setCurrentIndex(1)
            elif shared.exp.results_type == "rupture_force":
                self.list_choser.setCurrentIndex(2)
            elif shared.exp.results_type == "events_forces":
                self.list_choser.setCurrentIndex(3)
            elif shared.exp.results_type == "events_per_curve":
                self.list_choser.setCurrentIndex(4)
            elif shared.exp.results_type == "events_rupture_force":
                self.list_choser.setCurrentIndex(5)
            elif shared.exp.results_type == "loading_rates":
                self.list_choser.setCurrentIndex(6)
            elif shared.exp.results_type == "events_distance":
                self.list_choser.setCurrentIndex(7)
