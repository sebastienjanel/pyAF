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

"""Widget used to flatten the piezo image."""

from PyQt5 import QtWidgets
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFCheckBox
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from .. import widgets_list
from .. import shared


class FlattenWidget(PYAFWidget):
    """Allows to flatten the piezo image with order 1, 2, 3."""

    def __init__(self, parent):
        super().__init__(parent, "widget_flatten")

        if shared.exp.current_data.applied_flatten_order is not None:
            flatten_order = shared.exp.current_data.applied_flatten_order
        else:
            flatten_order = 1
        shared.exp.selected_flatten_order = flatten_order

        self.box_flatten = QtWidgets.QGroupBox("Flatten")
        self.VL_box_flatten = QtWidgets.QVBoxLayout()

        self.HL_flatten_canvas = QtWidgets.QHBoxLayout()
        self.W_original = QtWidgets.QWidget()
        self.W_original.setFixedSize(216, 216)
        mesh_type = "small_meshgrid_for_correction_original"
        self.MPL_canvas_original = \
            PYAFPlot(self, mesh_type, self.W_original, [3, 3, 72])
        self.canvas_flatten_widget = QtWidgets.QWidget()
        self.canvas_flatten_widget.setFixedSize(216, 216)
        mesh_type = "small_meshgrid_for_correction"
        self.MPL_canvas_flatten = \
            PYAFPlot(self, mesh_type, self.canvas_flatten_widget, [3, 3, 72])
        if shared.exp.current_data.is_single is False:
            self.update_MPL("MPL_canvas_flatten")

        self.HL_flatten_canvas.addStretch(1)
        self.HL_flatten_canvas.addWidget(self.W_original)
        self.HL_flatten_canvas.addWidget(self.canvas_flatten_widget)
        self.HL_flatten_canvas.addStretch(1)

        self.flatten_buttons_HL = QtWidgets.QHBoxLayout()
        self.BT_apply_flatten = PYAFButton(self, "apply_flatten", "Apply")
        self.BT_remove_flatten = \
            PYAFButton(self, "remove_flatten", "Go back to original")
        self.update_GUI("buttons_flatten")

        self.CB_apply_flatten_on_all = \
            PYAFCheckBox(self, "apply_flatten_on_all", "Apply on all")

        if len(shared.exp.list) < 2:
            self.CB_apply_flatten_on_all.setEnabled(False)
        if shared.exp.flatten_correction_all:
            self.CB_apply_flatten_on_all.setChecked(True)

        self.label_radiobuttons_order = QtWidgets.QLabel()
        self.label_radiobuttons_order.setText("Order : ")
        self.RV_flatten_order1 = QtWidgets.QRadioButton()
        self.RV_flatten_order2 = QtWidgets.QRadioButton()
        self.RV_flatten_order3 = QtWidgets.QRadioButton()
        self.RV_flatten_order1.setText("1")
        self.RV_flatten_order2.setText("2")
        self.RV_flatten_order3.setText("3")
        self.update_GUI("radiobuttons_flatten_order")
        self.RV_flatten_order1.clicked.connect(
            lambda: self.button_clicked("radiobutton_flatten_order1"))
        self.RV_flatten_order2.clicked.connect(
            lambda: self.button_clicked("radiobutton_flatten_order2"))
        self.RV_flatten_order3.clicked.connect(
            lambda: self.button_clicked("radiobutton_flatten_order3"))

        self.flatten_buttons_HL.addWidget(self.BT_apply_flatten)
        self.flatten_buttons_HL.addWidget(self.BT_remove_flatten)
        self.flatten_buttons_HL.addWidget(self.label_radiobuttons_order)
        self.flatten_buttons_HL.addWidget(self.RV_flatten_order1)
        self.flatten_buttons_HL.addWidget(self.RV_flatten_order2)
        self.flatten_buttons_HL.addWidget(self.RV_flatten_order3)
        self.flatten_buttons_HL.addStretch(1)

        self.VL_box_flatten.addLayout(self.HL_flatten_canvas)
        self.VL_box_flatten.addLayout(self.flatten_buttons_HL)
        self.VL_box_flatten.addWidget(self.CB_apply_flatten_on_all)
        self.VL_box_flatten.addStretch(1)
        self.box_flatten.setLayout(self.VL_box_flatten)

        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.addWidget(self.box_flatten, 0, 0)
        self.setLayout(self.gridLayout)

    def button_clicked(self, name):
        """Called when a button is clicked."""
        data = shared.exp.current_data

        if name == "apply_flatten":
            if shared.exp.flatten_correction_all:
                for i in range(len(shared.exp.list)):
                    data = shared.exp.list[i]
                    data.applied_flatten_order = \
                        shared.exp.selected_flatten_order
                    data.flatten_applied = True

            else:
                data.applied_flatten_order = shared.exp.selected_flatten_order
                data.flatten_applied = True
            widgets_list.widget_compute.update_MPL("MPL_meshgrid")
            self.BT_remove_flatten.setEnabled(True)
            self.BT_apply_flatten.setEnabled(False)

        elif name == "remove_flatten":
            if shared.exp.flatten_correction_all:
                for i in range(len(shared.exp.list)):
                    shared.exp.list[i].applied_flatten_order = None
                    data.flatten_applied = False
            else:
                data.applied_flatten_order = None
                data.flatten_applied = False

            widgets_list.widget_compute.update_MPL("MPL_meshgrid")
            self.BT_remove_flatten.setEnabled(False)
            self.BT_apply_flatten.setEnabled(True)

        elif name == "radiobutton_flatten_order1":
            shared.exp.selected_flatten_order = 1
            self.update_MPL("MPL_canvas_flatten")

        elif name == "radiobutton_flatten_order2":
            shared.exp.selected_flatten_order = 2
            self.update_MPL("MPL_canvas_flatten")

        elif name == "radiobutton_flatten_order3":
            shared.exp.selected_flatten_order = 3
            self.update_MPL("MPL_canvas_flatten")

    def checkbox_clicked(self, name):
        """Called when a checkbox is clicked."""
        if name == "apply_flatten_on_all":
            state = self.CB_apply_flatten_on_all.isChecked()
            if state:
                shared.exp.flatten_correction_all = True
            else:
                shared.exp.flatten_correction_all = False

    def update_MPL(self, name):
        """Updates the plot."""
        if name == "MPL_canvas_flatten":
            # Original piezo_image
            self.MPL_canvas_original.update_plot()

            # Flattened image
            self.MPL_canvas_flatten.update_plot()

    def update_widget(self):
        """Update the widget (GUI and plots)."""
        self.update_MPL("MPL_canvas_flatten")
        self.update_GUI("radiobuttons_flatten_order")
        self.update_GUI("buttons_flatten")

    def update_GUI(self, what):
        """Update the GUI."""
        order = shared.exp.current_data.applied_flatten_order
        if order is None:
            order = shared.exp.selected_flatten_order

        if what == "radiobuttons_flatten_order":
            if order == 1 or order is None:
                self.RV_flatten_order1.setChecked(True)
            elif order == 2:
                self.RV_flatten_order2.setChecked(True)
            elif order == 3:
                self.RV_flatten_order3.setChecked(True)

        elif what == "buttons_flatten":
            if shared.exp.current_data.flatten_applied is False:
                self.BT_remove_flatten.setEnabled(False)
                self.BT_apply_flatten.setEnabled(True)
            else:
                self.BT_apply_flatten.setEnabled(False)
                self.BT_remove_flatten.setEnabled(True)
