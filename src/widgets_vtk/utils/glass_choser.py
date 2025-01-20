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

"""Small widget asking the user to chose the glass layer to be added."""

from PyQt5 import QtWidgets

from ... import shared
from ... import widgets_list
from ...tools.gui_tools import PYAFButton, PYAFCheckBox
from ...tools.PYAFWidget import PYAFWidget


class VTKGlassChoserWidget(PYAFWidget):
    """VTKGlassChoserWidget asks the user which glass layer to add."""

    def __init__(self, parent, mode=None):
        name = "widget_vtk_glass_choser"
        super().__init__(parent, name)

        self.parent = parent
        self.CB_list = []
        self.mode = mode

        VL = QtWidgets.QVBoxLayout()

        # Label
        label = QtWidgets.QLabel()
        label.setText("Choose the ROI used for the glass layer :")

        # Display a list of checkboxes with the names in a scroll area
        scroll_area = QtWidgets.QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setMinimumSize(150, 250)
        widget = QtWidgets.QWidget()
        VL_list = QtWidgets.QVBoxLayout()

        current = widgets_list.widget_vtk.get_current_layer()
        data = shared.exp.list[shared.layer_list[current].afm_id]
        for i in range(len(data.roi_list)):
            roi = data.roi_list[i]
            if roi.glass_coeffs is not None:
                CB = PYAFCheckBox(self, str(i), "Roi " + str(roi.roi_id))
                CB.setChecked(False)
                self.CB_list.append(CB)
                VL_list.addWidget(CB)

        VL_list.addStretch(1)
        widget.setLayout(VL_list)

        scroll_area.setWidget(widget)

        HL_buttons = QtWidgets.QHBoxLayout()
        BT_continue = PYAFButton(self, "add", "Add")
        BT_cancel = PYAFButton(self, "cancel", "Cancel")
        HL_buttons.addWidget(BT_cancel)
        HL_buttons.addStretch(1)
        HL_buttons.addWidget(BT_continue)

        if data.roi_list == []:
            # Disable the button if there is no ROI to add
            BT_continue.setEnabled(False)

        VL.addWidget(label)
        VL.addWidget(scroll_area)
        VL.addLayout(HL_buttons)

        self.setLayout(VL)

    def button_clicked(self, what):
        """Called when a button is clicked."""
        if what == "add":
            # Make a list of rois which will be loaded.
            roi_list = []
            for CB in self.CB_list:
                if CB.isChecked():
                    roi_list.append(int(CB.name))

            self.close()

            if roi_list == []:
                # Close and do nothing if no roi was selected.
                return False

            # Add the layers
            menu = widgets_list.widget_vtk.menu_layer_list
            current = widgets_list.widget_vtk.get_current_layer()
            afm_id = shared.layer_list[current].afm_id
            canvas = widgets_list.widget_vtk.canvas
            for roi_id in roi_list:
                new_id = len(shared.layer_list)
                # Create the new layer
                canvas.addGlassLayer(new_id, afm_id, roi_id)
                # Add line to tablewidget
                menu.tableWidget.insertRow(new_id)
                menu.add_line_in_tablewidget(new_id, "glass")

            # Select the added line and select glass tab
            menu.tableWidget.setCurrentCell(new_id, 0)
            menu.tableWidget_clicked(new_id, 0)

        elif what == "cancel":
            self.close()

    def checkbox_clicked(self, what):
        """Called when a checkbox is clicked.

        Is not needed here as the values are only saved when the Add button is
        pushed.
        """
        pass
