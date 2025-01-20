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

"""Widget used to go through the curves."""

from PyQt5 import QtGui, QtWidgets

from ... import widgets_list
from ...tools.gui_tools import PYAFButton
from ...tools.PYAFWidget import PYAFWidget


class PlotSelectorWidget(PYAFWidget):
    """The plot selector allows to go through different plot types."""

    def __init__(self, parent, name):
        super().__init__(parent, name)

        if self.name == "plot_selector_data":
            self.parent_widget = widgets_list.widget_data

        self.plot_types = ("Deflection vs Scanner Extension", "Deflection vs Time", "Scanner Extension vs Time")

        self.BT_back = PYAFButton(self, "back", "<", size=60)
        self.BT_forward = PYAFButton(self, "forward", ">", size=60)

        self.list_plot_types = QtWidgets.QComboBox()

        for plot_type in self.plot_types:
            self.list_plot_types.addItem(plot_type)

        self.list_plot_types.activated.connect(lambda: self.list_updated("plot_types"))

        self.HL_select_curve = QtWidgets.QHBoxLayout()
        self.HL_select_curve.addStretch(1)
        self.HL_select_curve.addWidget(self.BT_back)
        self.HL_select_curve.addWidget(self.list_plot_types)
        self.HL_select_curve.addWidget(self.BT_forward)
        self.HL_select_curve.addStretch(1)

        self.setLayout(self.HL_select_curve)

        self.update_widget()

    def update_widget(self):
        """Update the widget."""
        if self.list_plot_types.currentIndex() == 0:
            self.BT_back.setEnabled(False)

        else:
            self.BT_back.setEnabled(True)

        if self.list_plot_types.currentIndex() == len(self.plot_types) - 1:
            self.BT_forward.setEnabled(False)

        else:
            self.BT_forward.setEnabled(True)

    def button_clicked(self, button):
        """Called when a button is clicked."""
        index = self.list_plot_types.currentIndex()

        if button == "back":
            new_index = index - 1
            self.list_plot_types.setCurrentIndex(new_index)
            self.parent_widget.canvas_stack.setCurrentIndex(new_index)
            self.update_widget()

        if button == "forward":
            new_index = index + 1
            self.list_plot_types.setCurrentIndex(new_index)
            self.parent_widget.canvas_stack.setCurrentIndex(new_index)
            self.update_widget()

    def list_updated(self, name):
        """Called when a list is updated."""

        if name == "plot_types":
            self.parent_widget.canvas_stack.setCurrentIndex(self.list_plot_types.currentIndex())
            self.update_widget()
