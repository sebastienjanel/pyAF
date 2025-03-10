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

"""
Widget used to do some clean up on curves.

Curves can be replaced by another. This is useful when you have a missing
pixel. A list of corrupted curves is also displayed.

"""

from .. import shared
from .. import widgets_list
from PyQt5 import QtGui, QtWidgets
from ..tools.gui_tools import PYAFButton, PYAFInput
from ..tools.gui_tools import ClearWidgetsFromLayout
from ..tools.PYAFWidget import PYAFWidget


class CleanUpWidget(PYAFWidget):
    """The cleanup widget."""

    def __init__(self, parent):
        super().__init__(parent, "widget_clean_up")

        self.VL_empty = None

        self.data = shared.exp.current_data

        self.VL = QtWidgets.QVBoxLayout()

        # BOX corrupted curves ------------------------------------------------
        self.box_corrupted_curves = QtWidgets.QGroupBox("Corrupted curves")
        self.VL_labels = QtWidgets.QVBoxLayout()

        self.gridlayout_box_labels = QtWidgets.QGridLayout()
        self.update_widget()

        self.scrollArea = QtWidgets.QScrollArea()
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setMinimumSize(150, 150)
        self.widget_labels = QtWidgets.QWidget()
        self.widget_labels.setLayout(self.gridlayout_box_labels)
        self.scrollArea.setWidget(self.widget_labels)
        self.VL_labels.addWidget(self.scrollArea)
        self.box_corrupted_curves.setLayout(self.VL_labels)

        # BOX Replace single curve --------------------------------------------

        text = "Replace single curve (Trace and retrace)"
        self.box_replace_single = QtWidgets.QGroupBox(text)
        self.HL_box_replace_single = QtWidgets.QHBoxLayout()
        self.input_single_original_x = PYAFInput(
            self,
            "input_single_original_x",
            "Replace")
        self.input_single_original_y = PYAFInput(
            self,
            "input_single_original_y", ",")
        self.input_single_new_x = PYAFInput(self, "input_single_new_x", "with")
        self.input_single_new_y = PYAFInput(self, "input_single_new_y", ",")
        self.BT_apply_replace_single = PYAFButton(
            self,
            "button_apply_replace_single",
            "Apply")
        self.HL_box_replace_single.addWidget(self.input_single_original_x)
        self.HL_box_replace_single.addWidget(self.input_single_original_y)
        self.HL_box_replace_single.addWidget(self.input_single_new_x)
        self.HL_box_replace_single.addWidget(self.input_single_new_y)
        self.HL_box_replace_single.addStretch(1)
        self.HL_box_replace_single.addWidget(self.BT_apply_replace_single)
        self.box_replace_single.setLayout(self.HL_box_replace_single)

        # BOX Replace line ----------------------------------------------------

        text = "Replace line (Trace and retrace)"
        self.box_replace_line = QtWidgets.QGroupBox(text)
        self.input_line_original_x = PYAFInput(
            self,
            "input_line_original_x",
            "Replace vertical line")
        self.input_line_new_x = PYAFInput(self, "input_line_new_x", "with")
        self.BT_apply_rep_line_x = PYAFButton(
            self,
            "button_apply_replace_line_x",
            "Apply")
        self.input_line_original_y = PYAFInput(
            self,
            "input_line_original_y",
            "Replace horizontal line")
        self.input_line_new_y = PYAFInput(self, "input_line_new_y", "with")
        self.BT_apply_rep_line_y = PYAFButton(
            self,
            "button_apply_replace_line_y",
            "Apply")

        self.gridLayout_replace_line = QtWidgets.QGridLayout()
        self.gridLayout_replace_line.addWidget(
            self.input_line_original_x,
            0,
            0)
        self.gridLayout_replace_line.addWidget(self.input_line_new_x, 0, 1)
        self.gridLayout_replace_line.addWidget(self.BT_apply_rep_line_x, 0, 2)
        self.gridLayout_replace_line.addWidget(
            self.input_line_original_y,
            1,
            0)
        self.gridLayout_replace_line.addWidget(self.input_line_new_y, 1, 1)
        self.gridLayout_replace_line.addWidget(self.BT_apply_rep_line_y, 1, 2)
        self.box_replace_line.setLayout(self.gridLayout_replace_line)

        # Layout --------------------------------------------------------------

        self.VL.addWidget(self.box_corrupted_curves)
        self.VL.addWidget(self.box_replace_single)
        self.VL.addWidget(self.box_replace_line)
        self.setLayout(self.VL)

    def button_clicked(self, what):
        """Called when a button is clicked."""
        # Note : The positions need also to be replaced, to be up to date !
        if what == "button_apply_replace_single":
            # Curve to be replaced
            orig_x = self.input_single_original_x.get_int_value() - 1
            orig_y = self.input_single_original_y.get_int_value() - 1
            # The new curve
            new_x = self.input_single_new_x.get_int_value() - 1
            new_y = self.input_single_new_y.get_int_value() - 1
            # Trace
            self.data.curves_approach[orig_x, orig_y] = \
                self.data.curves_approach[new_x][new_y]
            # Retrace
            self.data.curves_retraction[orig_x, orig_y] = \
                self.data.curves_retraction[new_x][new_y]
            # Piezo_image needs also to be replaced
            self.data.piezo_image[orig_x, orig_y] = \
                self.data.piezo_image[new_x][new_y]
            # Update positions
            self.data.approach_positions[orig_x, orig_y] = \
                self.data.approach_positions[new_x][new_y]

        elif what == "button_apply_replace_line_x":
            # Line to be replaced (vertical)
            orig_x = self.input_line_original_x.get_int_value() - 1
            # New line (vertical)
            new_x = self.input_line_new_x.get_int_value() - 1
            for j in range(self.data.nbr_pixels_y):
                # Trace
                self.data.curves_approach[orig_x, j] = \
                    self.data.curves_approach[new_x][j]

                # Retrace
                self.data.curves_retraction[orig_x, j] = \
                    self.data.curves_retraction[new_x][j]
                # Piezo image needs also to be replaced
                self.data.piezo_image[orig_x, j] = \
                    self.data.piezo_image[new_x][j]
                # Update positions
                self.data.approach_positions[orig_x, j] = \
                    self.data.approach_positions[new_x][j]

        elif what == "button_apply_replace_line_y":
            # Line to be replaced (horizontal)
            orig_y = self.input_line_original_y.get_int_value() - 1
            # New line (hoizontal)
            new_y = self.input_line_new_y.get_int_value() - 1
            for i in range(self.data.nbr_pixels_x):
                # Trace
                self.data.curves_approach[i, orig_y] = \
                    self.data.curves_approach[i][new_y]
                # Retrace
                self.data.curves_retraction[i, orig_y] = \
                    self.data.curves_retraction[i][new_y]
                # Piezo image needs also to be replaced
                self.data.piezo_image[i, orig_y] = \
                    self.data.piezo_image[i][new_y]
                # Update positions
                self.data.approach_positions[i, orig_y] = \
                    self.data.approach_positions[i][new_y]

        # Upate the plots in the compute widget
        widgets_list.widget_compute.update_MPL("MPL_meshgrid")
        widgets_list.widget_compute.update_MPL("MPL_canvas")

    def input_updated(self, name):
        """Do nothing here."""
        pass

    def update_widget(self):
        """Update the GUI."""
        # Change dataset
        self.data = shared.exp.current_data
        # Get the list of corrupted curves
        corrupted_curves = self.data.corrupted_curves
        # Clear the layout with the labels
        ClearWidgetsFromLayout(self.gridlayout_box_labels)
        # Fill the labels
        u = v = 0
        if corrupted_curves is not None:
            for curve in corrupted_curves:
                label_value = QtWidgets.QLabel()
                label_value.setFont(QtGui.QFont("defaultFamily", 12))
                text = "[" + str(curve[0] + 1) + "," + str(curve[1] + 1) + "]"
                label_value.setText(text)
                self.gridlayout_box_labels.addWidget(label_value, u, v)
                v = v + 1
                if v == 10:
                    v = 0
                    u = u + 1
            self.VL_empty = QtWidgets.QVBoxLayout()
            self.VL_empty.addStretch(1)
            self.gridlayout_box_labels.addLayout(self.VL_empty, u + 1, 0)
        else:
            label_value = QtWidgets.QLabel()
            label_value.setFont(QtGui.QFont("defaultFamily", 12))
            label_value.setText("None")
            self.gridlayout_box_labels.addWidget(label_value, 0, 0)
            self.VL_empty = QtWidgets.QVBoxLayout()
            self.VL_empty.addStretch(1)
            self.gridlayout_box_labels.addLayout(self.VL_empty, 1, 0)
