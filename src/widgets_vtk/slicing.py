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

"""Widget with options for the AFM tomographies."""

from PyQt5 import QtWidgets

from .. import shared
from .. import widgets_list
from ..tools.gui_tools import PYAFButton, PYAFInput
from ..tools.PYAFWidget import PYAFWidget


class SlicingWidget(PYAFWidget):
    """Widget with options to slice the tomographies."""

    def __init__(self, parent):
        super().__init__(parent, "widget_slicing")

        self.VL = QtWidgets.QVBoxLayout()

        self.BX_slicing = QtWidgets.QGroupBox("Slicing")
        self.VL_slicing = QtWidgets.QVBoxLayout()
        self.HL1_slicing = QtWidgets.QHBoxLayout()
        self.HL2_slicing = QtWidgets.QHBoxLayout()

        self.IN_slicing_top = PYAFInput(self, "input_slicing_top", "Top ")
        name = "input_slicing_bottom"
        self.IN_slicing_bottom = PYAFInput(self, name, "Bottom ")
        self.IN_slicing_left = PYAFInput(self, "input_slicing_left", "Left ")
        self.IN_slicing_right = PYAFInput(
            self, "input_slicing_right", "Right ")
        self.BT_slicing_apply = PYAFButton(self, "slicing_apply", "Apply", 250)

        self.HL1_slicing.addWidget(self.IN_slicing_top)
        self.HL1_slicing.addWidget(self.IN_slicing_bottom)
        self.HL1_slicing.addStretch(1)

        self.HL2_slicing.addWidget(self.IN_slicing_left)
        self.HL2_slicing.addWidget(self.IN_slicing_right)
        self.HL2_slicing.addStretch(1)

        self.VL_slicing.addLayout(self.HL1_slicing)
        self.VL_slicing.addLayout(self.HL2_slicing)
        self.VL_slicing.addWidget(self.BT_slicing_apply)

        self.BX_slicing.setLayout(self.VL_slicing)

        self.VL.addWidget(self.BX_slicing)

        self.setLayout(self.VL)

        self.update_widget()

    def button_clicked(self, button):
        """Called when a button is clicked."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]
        data = shared.exp.list[layer.afm_id]

        if button == "slicing_apply":
            val = self.IN_slicing_top.get_int_value()
            data.opengl_slice_top = val
            val = self.IN_slicing_bottom.get_int_value()
            data.opengl_slice_bottom = val
            val = self.IN_slicing_left.get_int_value()
            data.opengl_slice_left = val
            val = self.IN_slicing_right.get_int_value()
            data.opengl_slice_right = val
            layer.update_tomography()

    def update_widget(self):
        """Update the GUI."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]
        data = shared.exp.list[layer.afm_id]

        # Slicing inputs
        self.IN_slicing_top.changeValue(str(data.opengl_slice_top))
        self.IN_slicing_bottom.changeValue(str(data.opengl_slice_bottom))
        self.IN_slicing_left.changeValue(str(data.opengl_slice_left))
        self.IN_slicing_right.changeValue(str(data.opengl_slice_right))

    def input_updated(self, val):
        """Do nothing."""
        pass
