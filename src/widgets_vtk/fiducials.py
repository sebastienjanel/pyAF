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

"""Small widget allowing to manage fiducial markers."""

import logging
import numpy
from .. import shared
from .. import widgets_list
from PyQt5 import QtCore, QtWidgets
from ..tools import gui_tools
from ..tools import misc_tools
from ..tools.gui_tools import QLineEditForTableWidet
from ..tools.gui_tools import PYAFButton
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFInput


class FiducialsWidget(PYAFWidget):
    """In this widget you can manage the fiducials.

    You can load fiducials from a text file, or add them manually.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_fiducials")

        self.logger = logging.getLogger()
        self.logger.debug("Opening fiducials widget")

        # Load the settings
        self.settings = QtCore.QSettings()

        self.VL = QtWidgets.QVBoxLayout()

        HL = QtWidgets.QHBoxLayout()

        self.tableWidget = QtWidgets.QTableWidget(0, 5)
        policy = QtCore.Qt.ScrollBarAlwaysOff
        self.tableWidget.setHorizontalScrollBarPolicy(policy)
        behavior = QtWidgets.QAbstractItemView.SelectRows
        self.tableWidget.setSelectionBehavior(behavior)
        sel = QtWidgets.QAbstractItemView.ExtendedSelection
        self.tableWidget.setSelectionMode(sel)
        self.tableWidget.setHorizontalHeaderLabels(["", "", "X", "Y", "Z"])
        self.tableWidget.setFixedWidth(380)
        self.tableWidget.setColumnWidth(0, 40)
        self.tableWidget.setColumnWidth(1, 40)
        self.tableWidget.setColumnWidth(2, 100)
        self.tableWidget.setColumnWidth(3, 100)
        self.tableWidget.setColumnWidth(4, 100)
        self.tableWidget.verticalHeader().hide()
        self.tableWidget.cellClicked.connect(tableWidget_clicked)

        VL_options = QtWidgets.QVBoxLayout()

        self.BT_load_fiducials = PYAFButton(self, "load", "Load fiducials")
        self.BT_add_fiducial = PYAFButton(self, "add", "Add fiducial")
        self.BT_remove_fiducial = PYAFButton(
            self, "remove", "Remove fiducials")
        self.BT_colors = PYAFButton(self, "color", "Color")
        self.BT_change_heights = PYAFButton(self, "heights", "Change heights")
        self.BT_center = PYAFButton(self, "center", "Center")
        self.IN_radius = PYAFInput(self, "radius", "Radius (nm)")

        VL_options.addWidget(self.BT_load_fiducials)
        VL_options.addWidget(self.BT_add_fiducial)
        VL_options.addWidget(self.BT_remove_fiducial)
        VL_options.addWidget(self.BT_colors)
        VL_options.addWidget(self.BT_change_heights)
        VL_options.addWidget(self.BT_center)
        VL_options.addWidget(self.IN_radius)
        VL_options.addStretch(1)

        HL.addWidget(self.tableWidget)
        HL.addLayout(VL_options)

        self.VL.addLayout(HL)

        self.setLayout(self.VL)

        self.update_widget()

        self.tableWidget.selectRow(0)

    def button_clicked(self, button):
        """Called upon click on a button."""
        row = widgets_list.widget_vtk.get_current_layer()
        current_layer = shared.layer_list[row]

        if button == "load":
            # Get saved path
            settings_path = misc_tools.get_user_path(self.settings)

            filename = QtWidgets.QFileDialog.getOpenFileName(directory=settings_path)

            if not filename:
                return
            else:
                fiducials_file = open(filename, "rb")

                # Skip first line
                fiducials_file.readline()

                positions = []
                display = []

                for line in fiducials_file:
                    values = line.split()
                    # Add X, Y positions
                    # Last line is ['\x00']; so it has a length of 1, skip it
                    if len(values) != 1:
                        val1 = float(values[1]) * 1e3
                        val2 = float(values[2]) * 1e3
                        pos = [val1, val2, 0.0]
                        positions.append(pos)
                        display.append(True)

                current_layer.fiducials.add_to_fiducials(positions, display)

        elif button == "add":
            current_layer.fiducials.add_empty_fiducial()

            self.tableWidget.selectRow(len(current_layer.fiducials.positions))

        elif button == "remove":
            indexes = []
            for item in self.tableWidget.selectedIndexes():
                indexes.append(item.row())

            # It depends which way you have selected so I need to sort the
            # index array
            indexes = numpy.sort(numpy.unique(indexes), -1)[::-1].tolist()
            if indexes != []:
                for index in indexes:
                    current_layer.fiducials.remove_fiducial(index)

            self.tableWidget.selectRow(0)

        elif button == "color":
            current_color = current_layer.fiducials.color
            color = misc_tools.ask_user_for_color(current_color)
            if color is not False:
                row = widgets_list.widget_vtk.get_current_layer()
                current_layer = shared.layer_list[row]
                current_layer.fiducials.color = color
                current_layer.fiducials.update_fiducials()

        elif button == "heights":
            text, ok = QtWidgets.QInputDialog.getDouble(self, "Relative heights",
                                                    "Change z to ", value=0)

            if ok:
                # Change the heights for the selected rows
                for item in self.tableWidget.selectionModel().selectedRows():
                    index = item.row()
                    current_layer.fiducials.positions[index][2] = float(text)

            # Update the fiducials
            current_layer.fiducials.update_fiducials()

        elif button == "center":
            text = ("The fiducials will be centered respectively to the "
                    "layer's size. Check that the layer has the right size "
                    "before centering your fiducials.")

            # Create a message box
            msg = QtWidgets.QMessageBox()
            msg.setText("Centering")
            msg.setInformativeText(text)
            msg.setIcon(QtWidgets.QMessageBox.Information)
            msg.addButton(QtWidgets.QMessageBox.Ok)
            msg.addButton(QtWidgets.QMessageBox.Abort)
            ret = msg.exec_()

            if ret == QtWidgets.QMessageBox.Ok:
                # Center the fiducials
                current_layer.fiducials.center_fiducials()

        self.update_widget()

    def update_widget(self):
        """Updates the GUI of the widget."""
        # Remove the old rows
        for i in range(self.tableWidget.rowCount(), -1, -1):
            self.tableWidget.removeRow(i)

        row = widgets_list.widget_vtk.get_current_layer()
        current_layer = shared.layer_list[row]
        positions = current_layer.fiducials.positions
        display = current_layer.fiducials.display_list

        # Fill new rows
        for i in range(len(positions)):
            self.tableWidget.insertRow(i)
            self.tableWidget.setRowHeight(i, 25)

            # Index
            label = QtWidgets.QLabel()
            label.setFont(self.parent.parent.smallfont)
            label.setText(str(i + 1))
            cw = gui_tools.CenteredCellwidget(label)
            self.tableWidget.setCellWidget(i, 0, cw)

            # Checkbox
            cb = gui_tools.CenteredCellCheckbox(self, "display", i)
            self.tableWidget.setCellWidget(i, 1, cb)
            self.tableWidget.cellWidget(i, 1).checkbox.setChecked(display[i])

            # X
            item = QLineEditForTableWidet(self, "posx", i)
            item.setText(str(positions[i][0]))
            item.setFont(self.parent.parent.smallfont)
            self.tableWidget.setCellWidget(i, 2, item)

            # Y
            item = QLineEditForTableWidet(self, "posy", i)
            item.setText(str(positions[i][1]))
            item.setFont(self.parent.parent.smallfont)
            self.tableWidget.setCellWidget(i, 3, item)

            # Z
            item = QLineEditForTableWidet(self, "posz", i)
            item.setText(str(positions[i][2]))
            item.setFont(self.parent.parent.smallfont)
            self.tableWidget.setCellWidget(i, 4, item)

        # Update the input field for the sizes
        self.IN_radius.changeValue(str(current_layer.fiducials.radius))

    def input_updated(self, name, row=None):
        """Called when an input field is updated.

        For input fields (QLineEditForTableWidet), name == col
        """
        layerrow = widgets_list.widget_vtk.get_current_layer()
        current_layer = shared.layer_list[layerrow]

        if name == "radius":
            current_layer.fiducials.radius = self.IN_radius.get_int_value()

        elif name == "posx":
            value = float(str(self.tableWidget.cellWidget(row, 2).text()))
            current_layer.fiducials.positions[row][0] = value

        elif name == "posy":
            value = float(str(self.tableWidget.cellWidget(row, 3).text()))
            current_layer.fiducials.positions[row][1] = value

        elif name == "posz":
            value = float(str(self.tableWidget.cellWidget(row, 4).text()))
            current_layer.fiducials.positions[row][2] = value

        current_layer.fiducials.update_fiducials()

    def checkbox_clicked(self, name, data_id=None):
        """Called when a checkbox is clicked."""
        row = widgets_list.widget_vtk.get_current_layer()
        current_layer = shared.layer_list[row]

        if name == "display":
            data_id = int(data_id)

            value = self.tableWidget.cellWidget(
                data_id, 1).checkbox.isChecked()
            current_layer.fiducials.display_list[data_id] = value
            current_layer.fiducials.hide_or_display()


def tableWidget_clicked(row, _):
    """Called when the user clicks on a row in the tablewidget.

    Will change the color of the current fiducial.
    """
    layer_row = widgets_list.widget_vtk.get_current_layer()
    current_layer = shared.layer_list[layer_row]

    current_layer.fiducials.set_current_fiducial(row)
