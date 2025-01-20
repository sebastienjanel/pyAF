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

"""Widget allowing to change the scales of the plot in the slices widget."""

from PyQt5 import QtWidgets
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFInput
from .. import shared


class SlicesAxesWidget(QtWidgets.QWidget):
    """Options for the scales of the slices widget plot."""

    def __init__(self, parent):
        super().__init__()
        self.open = True
        self.parent = parent

        data = shared.exp.current_data

        self.gridlayout = QtWidgets.QGridLayout()
        self.list = PYAFComboBox(self, None)
        self.list.addItem("Auto")
        self.list.addItem("Manual")
        if data.slices_z_mode == "auto":
            self.list.setCurrentIndex(0)
        elif data.slices_z_mode == "manual":
            self.list.setCurrentIndex(1)
        self.list.activated.connect(self.list_updated)

        self.input_min = PYAFInput(self, "input_min", "Min", width=70)
        self.input_max = PYAFInput(self, "input_max", "Max", width=70)
        self.button_apply = PYAFButton(self, "button_apply", "Apply")

        self.gridlayout.addWidget(self.list, 0, 0)
        self.setLayout(self.gridlayout)

        self.update_GUI("input")

    def closeEvent(self, _):
        """Called on close."""
        self.open = False

    def _input_updated(self, val):
        """Called when an input is updated."""
        # Can be called only if the widget still exists
        if self.open:
            self.input_updated(val)

    def list_updated(self):
        """Method called when the list is updated.

        Displays or hides the two inputs and refreshes parent slices plot.
        """
        data = shared.exp.list[shared.exp.id_selected]

        if self.list.currentIndex() == 0:
            data.slices_z_mode = "auto"
        elif self.list.currentIndex() == 1:
            data.slices_z_mode = "manual"

        self.update_GUI("input")

        # Update only on auto mode, if manual use the apply button
        if self.list.currentIndex() == 0:
            self.parent.update_MPL()

    def update_GUI(self, element):
        """Method called to update the GUI layout of the widget.

        Input fields are displayed or hidden.
        """
        if element == "input":
            data = shared.exp.list[shared.exp.id_selected]
            if data.slices_z_mode == "manual":
                self.gridlayout.addWidget(self.input_min, 1, 0)
                self.gridlayout.addWidget(self.input_max, 2, 0)
                self.gridlayout.addWidget(self.button_apply, 3, 0)
                self.input_min.changeValue(data.slices_z_min)
                self.input_max.changeValue(data.slices_z_max)
            else:
                self.input_min.setParent(None)
                self.input_max.setParent(None)
                self.button_apply.setParent(None)

    def input_updated(self, field):
        """Method called when a value is updated."""
        data = shared.exp.list[shared.exp.id_selected]
        if field == "input_min":
            if self.input_min.input.text() != "":
                data.slices_z_min = self.input_min.get_float_value()
        elif field == "input_max":
            if self.input_max.input.text() != "":
                data.slices_z_max = self.input_max.get_float_value()

    def button_clicked(self, button):
        """Methode called when a button is clicked."""
        if button == "button_apply":
            self.input_updated("input_min")
            self.input_updated("input_max")
            self.parent.update_MPL()
