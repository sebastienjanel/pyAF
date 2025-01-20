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

"""Small widget allowing to set the options for the binary thresholding."""

from PyQt5 import QtWidgets

from ... import shared
from ... import widgets_list
from ...tools.gui_tools import PYAFInput, PYAFCheckBox
from ...tools.PYAFWidget import PYAFWidget


class ThresholdOptionsWidget(PYAFWidget):
    """Small widget with options for the thresholding."""

    def __init__(self, parent):
        name = "widget_vtk_binary_threshold"
        super().__init__(parent, name)

        self.parent = parent

        VL = QtWidgets.QVBoxLayout()

        self.CB_use_thresh = PYAFCheckBox(self, "use_th", "Use")
        self.IN_lower_th = PYAFInput(self, "lower_th", "Lower Threshold")
        self.IN_upper_th = PYAFInput(self, "upper_th", "Upper Threshold")
        self.IN_lower_val = PYAFInput(self, "lower_val", "Lower Value")
        self.IN_upper_val = PYAFInput(self, "upper_val", "Upper Value")

        VL.addWidget(self.IN_lower_th)
        VL.addWidget(self.IN_upper_th)
        VL.addWidget(self.IN_lower_val)
        VL.addWidget(self.IN_upper_val)

        VL.addWidget(self.CB_use_thresh)

        self.setLayout(VL)

        self.update_widget()

    def update_widget(self):
        """Updates the GUI."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        self.CB_use_thresh.setChecked(layer.use_binary_thresholding)

        self.IN_lower_th.changeValue(layer.lower_thresh)
        self.IN_upper_th.changeValue(layer.upper_thresh)
        self.IN_lower_val.changeValue(layer.lower_thresh_value)
        self.IN_upper_val.changeValue(layer.upper_thresh_value)

    def input_updated(self, name):
        """Called when the input is updated."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if name == "lower_th":
            val = self.IN_lower_th.get_int_value()
            layer.lower_thresh = val
        elif name == "upper_th":
            val = self.IN_upper_th.get_int_value()
            layer.upper_thresh = val
        elif name == "lower_val":
            val = self.IN_lower_val.get_int_value()
            layer.lower_thresh_value = val
        elif name == "upper_val":
            val = self.IN_upper_val.get_int_value()
            layer.upper_thresh_value = val

        layer.update_actor()

    def checkbox_clicked(self, name):
        """Called whenever a checkox is clicked."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if name == "use_th":
            layer.use_binary_thresholding = self.CB_use_thresh.isChecked()
            layer.update_actor()
