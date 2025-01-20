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

"""Widget used to change the height of the slice for the stiffness meshgrid."""

from PyQt5 import QtWidgets
from ..tools.gui_tools import PYAFInput
from .. import widgets_list
from .. import shared


class ChangeSliceDepthWidget(QtWidgets.QWidget):
    """Small widget with an input field to change the value of the slice."""

    def __init__(self, parent):
        super().__init__()

        self.open = True
        self.parent = parent
        self.data = shared.exp.current_data

        VL = QtWidgets.QVBoxLayout()
        self.input = PYAFInput(self, "input", "Value", 200)
        self.input.changeValue(self.data.stiffness_slice_depth)
        VL.addWidget(self.input)
        self.setLayout(VL)

    def closeEvent(self, _):
        """Called on close."""
        self.open = False

    def _input_updated(self, val):
        """Called when an input is updated."""
        # Can be called only if the widget still exists
        if self.open:
            self.input_updated(val)

    def input_updated(self, field):
        """Called when an input is updated."""
        if field == "input":
            if self.input.input.text() != "":
                val = self.input.get_float_value()
                self.data.stiffness_slice_depth = val
                self.parent.update_MPL("MPL_canvas1")

                if widgets_list.widget_slices is not None:
                    widgets_list.widget_slices.update_MPL()
