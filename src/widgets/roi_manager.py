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

"""Widget used to manage the regions of interest (ROIs)."""

import numpy
from .. import shared
from .. import widgets_list
import logging
from qtpy import QtCore, QtWidgets
from ..tools import gui_tools
from ..tools import math_tools
from ..tools import planefit
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFButtonGroup
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFInput
from .. import experiment


class RoiManagerWidget(PYAFWidget):
    """Widget used to manage the regions of interest (ROIs).

    The Rois can be used to define the substrate's surface, or to sort data
    in the results by excluding data from a ROI.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_roi_manager")

        self.BG_selected = None

        self.logger = logging.getLogger()
        self.logger.debug("Opening ROI manager")

        self.VL = QtWidgets.QVBoxLayout()
        self.HL = QtWidgets.QHBoxLayout()

        self.TW = QtWidgets.QTableWidget(0, 8)
        self.TW.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.TW.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.TW.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        labels = ["#", "Draw", "Show", "Pixels", "Area [µm²]", "Color",
                  "Substrate", "Data"]
        self.TW.setHorizontalHeaderLabels(labels)
        tooltips = ["", "ROI drawn or erased with the right click menu of "
                    "the map", "Display the ROI on the map", "", "", "",
                    "ROI defined as substrate", "Curves of the ROI "
                    "discarded with Discard data"]
        for col, tooltip in enumerate(tooltips):
            self.TW.horizontalHeaderItem(col).setToolTip(tooltip)
        widths = [30, 45, 45, 60, 85, 100, 75, 80]
        for col, width in enumerate(widths):
            self.TW.setColumnWidth(col, width)
        self.TW.setFixedWidth(sum(widths) + 4)
        self.TW.verticalHeader().hide()
        self.TW.cellClicked.connect(self.tableWidget_clicked)

        self.VL_options = QtWidgets.QVBoxLayout()

        self.input_size = PYAFInput(self, "size", "Size")
        self.BG = PYAFButtonGroup(self, "bg")
        self.radiobutton_square = QtWidgets.QRadioButton("Square")
        self.radiobutton_circle = QtWidgets.QRadioButton("Circle")
        self.radiobutton_dot = QtWidgets.QRadioButton("Dot")
        self.BG.addButton(self.radiobutton_square, 0)
        self.BG.addButton(self.radiobutton_circle, 1)
        self.BG.addButton(self.radiobutton_dot, 2)

        self.button_add = PYAFButton(self, "add", "Add ROI")
        self.button_remove = PYAFButton(self, "remove", "Remove")
        self.button_invert = PYAFButton(self, "invert", "Invert to new")
        self.BT_glass = PYAFButton(self, "glass", "Define as substrate")
        self.BT_planefit = PYAFButton(self, "planefit", "Planefit")
        self.BT_planefit.setToolTip(
            "Fit a plane (offset and tilt) on the selected substrate ROI and\n"
            "subtract it from the piezo height and topography maps of this\n"
            "file: the substrate becomes the zero height. Only the maps are\n"
            "corrected, not the computations. Define the ROI as substrate\n"
            "first.")
        self.BT_discard = PYAFButton(self, "discard", "Discard data")
        self.BT_discard.setToolTip(
            "Discard data: the curves of the selected ROI are not used in\n"
            "the results (crossed on the maps). The ROI is kept.\n"
            "Keep data: use these curves again.")
        self.BT_volume = PYAFButton(self, "volume", "Compute volume")
        self.BT_volume.setToolTip(
            "Volume between the zero-force topography and the plane of the\n"
            "selected substrate ROI, over the whole map. Define the ROI as\n"
            "substrate first.")

        self.VL_options.addWidget(self.input_size)
        self.VL_options.addWidget(self.radiobutton_square)
        self.VL_options.addWidget(self.radiobutton_circle)
        self.VL_options.addWidget(self.radiobutton_dot)
        self.VL_options.addWidget(self.button_add)
        self.VL_options.addWidget(self.button_remove)
        self.VL_options.addWidget(self.BT_discard)
        self.VL_options.addWidget(self.button_invert)
        self.VL_options.addWidget(self.BT_glass)
        self.VL_options.addWidget(self.BT_planefit)
        self.VL_options.addWidget(self.BT_volume)
        self.VL_options.addStretch(1)

        self.HL.addWidget(self.TW)
        self.HL.addLayout(self.VL_options)
        self.VL.addLayout(self.HL)
        self.setLayout(self.VL)

        self.update_widget()

        self.TW.selectRow(0)

    def tableWidget_clicked(self, row, col):
        """Updates the glass button when changing ROI."""
        # Do nothing with these arguments :
        _ = row
        _ = col

        self.update_GUI("BT_glass")

    def button_clicked(self, button):
        """Called when a button is clicked."""
        data = shared.exp.current_data

        if button == "add":
            if data.roi_selected_row is None:
                data.roi_selected_row = 0

            newid = len(data.roi_list)

            self.logger.debug("Adding ROI %s", newid)

            data.roi_list.append(new_roi(data, newid))

            self.update_widget()
            self.update_results()

            self.TW.selectRow(data.roi_selected_row)

            widgets_list.widget_results_single.reset_table()

        elif button == "remove":
            indexes = []
            for item in self.TW.selectedIndexes():
                indexes.append(item.row())

            # It depends which way you have selected so I need to sort the
            # index array
            indexes = numpy.sort(numpy.unique(indexes), -1)[::-1].tolist()

            # Check if there are no glas layers defined in the VTK widget
            # The ROI has to be from the same file and have the same id
            u_id = data.unique_id
            if widgets_list.widget_vtk is not None:
                for layer in shared.layer_list:
                    if layer.type == "glass" and layer.afm_id == u_id:
                        for index in indexes:
                            if index == layer.roi_id:
                                # Display an error message and abort removing
                                ask_user_to_delete_layer_first()
                                self.logger.debug(
                                    "Abort removing ROI %s",
                                    index)
                                return False

            self.logger.debug("Removing ROI %s", indexes)

            if indexes:
                self.update_planefit_roi(data, indexes)
                for i in indexes:
                    del data.roi_list[i]
                    self.TW.removeRow(i)

                # ROI used to draw: the same one if it is kept, else the
                # first one
                active = data.roi_selected_row
                if not data.roi_list:
                    data.roi_selected_row = None
                elif active is None or active in indexes:
                    data.roi_selected_row = 0
                else:
                    data.roi_selected_row = active - sum(
                        1 for i in indexes if i < active)

                if data.roi_selected_row is not None:
                    self.TW.selectRow(0)

                self.update_widget()
                self.parent.MPL_canvas1.canvas.update_blit("all")
                self.update_results()
                widgets_list.widget_compute.update_GUI("roi_list_height_corr")

            widgets_list.widget_results_single.reset_table()

            update_roi_on_multimeshgrids_widget()

        elif button == "bg":
            value = self.BG.checkedId()
            if value == 0:
                shape = "square"
            elif value == 1:
                shape = "circle"
            elif value == 2:
                shape = "dot"
            shared.exp.roi_cursor_shape = shape

            self.logger.debug("Changing roi shape to %s", shape)

        elif button == "invert":
            for item in self.TW.selectionModel().selectedRows():
                index = item.row()

                newid = len(data.roi_list)

                txt = (
                    "Adding ROI (by inverting). New ROI id = %s, "
                    "old ROI = %s")
                self.logger.debug(txt, newid, index)

                data.roi_list.append(new_roi(data, newid))

                # Get new values (can be slow)
                inverted_roi = []
                for i in range(data.nbr_pixels_x):
                    for j in range(data.nbr_pixels_y):
                        values = data.roi_list[index].values
                        if not math_tools.in_list(values, [i, j]):
                            inverted_roi.append([i, j])

                # Add data to ROI
                data.roi_list[newid].values = inverted_roi

            self.update_widget()
            self.parent.MPL_canvas1.canvas.update_blit("all")
            self.update_results()
            widgets_list.widget_results_single.reset_table()

            update_roi_on_multimeshgrids_widget()

        elif button == "button_group_selected":
            data.roi_selected_row = self.BG_selected.checkedId()
            self.logger.debug("Selecting ROI %s", data.roi_selected_row)

        elif button == "glass":
            x_size = data.scan_size_x / data.nbr_pixels_x
            y_size = data.scan_size_x / data.nbr_pixels_x

            for item in self.TW.selectionModel().selectedRows():
                index = item.row()
                roi = data.roi_list[index]
                if roi.glass_coeffs is None:
                    self.logger.debug("Set as glass. ROI %s", index)

                    x = []
                    y = []
                    z = []

                    for i, j in roi.values:
                        x.append(i * x_size)
                        y.append(j * y_size)
                        z.append(data.topography[i][j])
                    x = numpy.array(x)
                    y = numpy.array(y)
                    z = numpy.array(z)

                    # Ignore the pixels without topography (no POC found)
                    finite = numpy.isfinite(z)
                    x, y, z = x[finite], y[finite], z[finite]

                    # Get the coefficients of the 2D fit
                    co = math_tools.polyfit2d(
                        x,
                        y,
                        z,
                        order=1,
                        linear=True)
                    roi.glass_coeffs = co

                else:
                    self.logger.debug("Remove glass. ROI %s", index)

                    roi.glass_coeffs = None

                    # The planefit used this substrate
                    if data.planefit_applied and data.planefit_roi == index:
                        self.set_planefit(data, None)

            self.update_widget()
            widgets_list.widget_compute.update_GUI("roi_list_height_corr")

        elif button == "volume":
            self.show_volume(data)

        elif button == "planefit":
            rows = [item.row()
                    for item in self.TW.selectionModel().selectedRows()]
            if planefit.get_planefit_roi(data) is not None:
                # Remove the planefit
                self.set_planefit(data, None)
            elif rows:
                self.set_planefit(data, rows[0])
            self.update_GUI("BT_glass")

        elif button == "discard":
            # Discard (or keep again) the curves of the selected ROI
            rows = [item.row()
                    for item in self.TW.selectionModel().selectedRows()]
            if not rows:
                return
            roi = data.roi_list[rows[0]]
            if getattr(roi, "discarded", []):
                # Keep data: only the curves discarded by this ROI (not the
                # ones discarded otherwise)
                for i, j in roi.discarded:
                    data.discarded_curves[i][j] = 0
                roi.discarded = []
            else:
                roi.discarded = []
                for i, j in roi.values:
                    if data.discarded_curves[i][j] == 0:
                        data.discarded_curves[i][j] = 1
                        roi.discarded.append([i, j])
            self.update_widget()

            # Refresh the meshgrids and the results
            widgets_list.widget_data.meshgrid.update_plot()
            widgets_list.widget_results.update_MPL("MPL_canvas1")
            self.update_results()

    def set_planefit(self, data, roi_index):
        """Apply the planefit with a substrate ROI (None: remove it), then
        rescale and redraw the height maps."""
        self.logger.debug("Planefit with ROI %s", roi_index)
        data.planefit_applied = roi_index is not None
        data.planefit_roi = roi_index
        planefit.reset_color_scales(data)
        if data.meshgrid_type in ("piezo", "topo"):
            widgets_list.widget_results.update_MPL("MPL_canvas1")
        if widgets_list.widget_multimeshgrids is not None:
            widgets_list.widget_multimeshgrids.update_widget()
        self.update_3D(data)

    def show_volume(self, data):
        """Computes the volume above the selected substrate ROI and displays
        it."""
        rows = [item.row() for item in self.TW.selectionModel().selectedRows()]
        if not rows or data.roi_list[rows[0]].glass_coeffs is None:
            return
        res = planefit.get_volume(data, data.roi_list[rows[0]])
        x_size, y_size = res["pixel_size"]
        title = ("Volume relative to the plane of the substrate ROI " +
                 str(rows[0] + 1))
        text = ("Zero-force topography, whole map\n\n"
                "Net volume: %.4g µm³\n"
                "Above the plane: %.4g µm³\n"
                "Below the plane: %.4g µm³\n"
                "Maximum height: %.4g nm\n\n"
                "%d pixels of %.4g × %.4g nm") % (
                    res["net"] * 1e-9, res["above"] * 1e-9,
                    res["below"] * 1e-9, res["max_height"], res["pixels"],
                    x_size, y_size)
        if res["no_topography"]:
            text += ("\n%d pixels without topography (no contact point) are "
                     "not counted") % res["no_topography"]
        box = QtWidgets.QMessageBox(self)
        box.setWindowTitle("Volume - " + data.filename)
        box.setText(title)
        box.setInformativeText(text)
        box.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        box.exec()

    @staticmethod
    def update_3D(data):
        """Redraw the 3D layers of the file (tomography and substrates)."""
        if widgets_list.widget_vtk is None:
            return
        for layer in shared.layer_list:
            if layer.type not in ("afm", "glass"):
                continue
            if shared.exp.list[layer.afm_id] is not data:
                continue
            if layer.type == "afm":
                layer.update_tomography()
            else:
                layer.update_actor()

    def update_planefit_roi(self, data, removed):
        """Keep the index of the planefit ROI when ROIs are removed (no
        planefit if its ROI is removed)."""
        if data.planefit_roi is None:
            return
        if data.planefit_roi in removed:
            data.planefit_applied = False
            data.planefit_roi = None
            planefit.reset_color_scales(data)
            self.update_3D(data)
        else:
            data.planefit_roi -= sum(1 for i in removed
                                     if i < data.planefit_roi)

    def update_results(self):
        """Update the results in the results tables."""
        if self.parent.parent.tabs.currentIndex() != 3:
            # Update later if we go on results tab
            self.parent.parent.update_results = True
        elif self.parent.parent.tabs.currentIndex() == 3:
            # We are on the results tab, we should update directly
            widgets_list.widget_results_single.update_widget()
        elif self.parent.parent.tabs.currentIndex() == 4:
            # We are on the results tab, we should update directly
            widgets_list.widget_results_groups.update_widget()

    def checkbox_clicked(self, checkbox, cb_id):
        """Called when a checkbox is clicked."""
        data = shared.exp.current_data

        if checkbox == "display_roi":
            value = self.TW.cellWidget(cb_id, 2).checkbox.isChecked()
            self.logger.debug("Display/hide ROI %s, %s", cb_id, value)
            data.roi_list[cb_id].display = value
            self.parent.MPL_canvas1.canvas.update_blit("all")

            update_roi_on_multimeshgrids_widget()

    def update_widget(self):
        """Update the content of the ROI manager."""
        # Get previously selected rows
        selected_rows = []
        for item in self.TW.selectionModel().selectedRows():
            selected_rows.append(item.row())

        data = shared.exp.current_data
        for i in range(self.TW.rowCount(), -1, -1):
            self.TW.removeRow(i)

        self.BG_selected = PYAFButtonGroup(self, "button_group_selected")

        # Add the rows to the tablewidget
        for i in range(len(data.roi_list)):
            self.TW.insertRow(i)
            self.TW.setRowHeight(i, 32)

            # Index
            label = QtWidgets.QLabel()
            label.setFont(self.parent.parent.smallfont)
            label.setText(str(i + 1))
            self.TW.setCellWidget(i, 0, gui_tools.CenteredCellwidget(label))

            # Radiobutton
            self.TW.setCellWidget(i, 1, gui_tools.CenteredCellRadioButton())
            button = self.TW.cellWidget(i, 1).radio_button
            self.BG_selected.addButton(button, i)
            if data.roi_selected_row == i or data.roi_selected_row is None:
                button.setChecked(True)
                data.roi_selected_row = i

            # Checkbox
            cb = gui_tools.CenteredCellCheckbox(self, "display_roi", i)
            self.TW.setCellWidget(i, 2, cb)
            disp = data.roi_list[i].display
            self.TW.cellWidget(i, 2).checkbox.setChecked(disp)

            # Pixels
            label = QtWidgets.QLabel()
            label.setFont(self.parent.parent.smallfont)
            label.setText(str(len(data.roi_list[i].values)))
            self.TW.setCellWidget(i, 3, gui_tools.CenteredCellwidget(label))

            # Area
            xsize = data.scan_size_x / data.nbr_pixels_x
            ysize = data.scan_size_y / data.nbr_pixels_y
            one_area = xsize * ysize
            label = QtWidgets.QLabel()
            label.setFont(self.parent.parent.smallfont)
            area = len(data.roi_list[i].values) * one_area * 1e-6
            label.setText("%.4g" % area)
            self.TW.setCellWidget(i, 4, gui_tools.CenteredCellwidget(label))

            # Color
            colorlist = PYAFComboBox(self, "colorlist")
            colorlist.setFont(self.parent.parent.smallfont)
            for j in range(len(data.roi_color_names)):
                colorlist.addItem(data.roi_color_names[j])
            colorlist.setCurrentIndex(data.roi_list[i].color)
            # In a layout, to keep the size of the cell
            self.TW.setCellWidget(i, 5, gui_tools.FullCellWidget(colorlist))

            # Glass
            label = QtWidgets.QLabel()
            label.setFont(self.parent.parent.smallfont)
            if data.roi_list[i].glass_coeffs is not None:
                label.setText("Yes")
            else:
                label.setText("")
            self.TW.setCellWidget(i, 6, gui_tools.CenteredCellwidget(label))

            # Discarded data
            label = QtWidgets.QLabel()
            label.setFont(self.parent.parent.smallfont)
            if getattr(data.roi_list[i], "discarded", []):
                label.setText("Discarded")
            self.TW.setCellWidget(i, 7, gui_tools.CenteredCellwidget(label))

        # Select rows
        for row in selected_rows:
            self.TW.selectRow(row)

        # Select the first row
        if selected_rows == []:
            self.TW.selectRow(0)

        # Shape
        if shared.exp.roi_cursor_shape == "square":
            self.radiobutton_square.setChecked(True)
        elif shared.exp.roi_cursor_shape == "circle":
            self.radiobutton_circle.setChecked(True)
        elif shared.exp.roi_cursor_shape == "dot":
            self.radiobutton_dot.setChecked(True)

        # Size
        self.input_size.changeValue(shared.exp.roi_cursor_size)

        self.update_GUI("BT_glass")

    def update_GUI(self, element):
        """Update some elements of the ROI manager."""
        data = shared.exp.current_data

        if element == "BT_glass":
            found_item = False
            for item in self.TW.selectionModel().selectedRows():
                index = item.row()
                roi = data.roi_list[index]
                found_item = True

                if not data.stiffness_calculated:
                    self.BT_glass.setEnabled(False)
                else:
                    if roi.values == []:
                        self.BT_glass.setEnabled(False)
                    else:
                        self.BT_glass.setEnabled(True)

                if roi.glass_coeffs is None:
                    self.BT_glass.setText("Define as substrate")
                else:
                    self.BT_glass.setText("No substrate")

        if not found_item:
            # No Roi defined, disable glass button
            self.BT_glass.setEnabled(False)

        if element == "BT_glass":
            # Planefit: with the selected ROI if it is a substrate, or to
            # remove the current planefit
            applied = planefit.get_planefit_roi(data) is not None
            rows = [item.row()
                    for item in self.TW.selectionModel().selectedRows()]
            substrate = bool(rows) and \
                data.roi_list[rows[0]].glass_coeffs is not None
            self.BT_planefit.setText(
                "Remove planefit" if applied else "Planefit")
            self.BT_planefit.setEnabled(applied or substrate)
            self.BT_volume.setEnabled(substrate and data.stiffness_calculated)
            discarded = bool(rows) and \
                bool(getattr(data.roi_list[rows[0]], "discarded", []))
            self.BT_discard.setText("Keep data" if discarded
                                    else "Discard data")
            self.BT_discard.setEnabled(discarded or (
                bool(rows) and bool(data.roi_list[rows[0]].values)))

    def input_updated(self, name):
        """Called when an input is updated."""
        if name == "size":
            value = self.input_size.get_int_value()
            shared.exp.roi_cursor_size = value
            self.logger.debug("Change roi pencil size to %s", value)

    def list_updated(self, name):
        """Called when a list is updated."""
        data = shared.exp.current_data

        if name == "colorlist":
            for i in range(len(data.roi_list)):
                cell = self.TW.cellWidget(i, 5)
                data.roi_list[i].color = cell.widget.currentIndex()
            self.parent.MPL_canvas1.canvas.update_blit("all")

            update_roi_on_multimeshgrids_widget()


def new_roi(data, roi_id):
    """New ROI, with the first color not used by the other ROIs of the file
    (the colors are used again when all of them are used)."""
    roi = experiment.ROI(roi_id)
    used = [other.color for other in data.roi_list]
    nbr_colors = len(data.roi_color_names)
    free = [i for i in range(nbr_colors) if i not in used]
    roi.color = free[0] if free else len(data.roi_list) % nbr_colors
    return roi


def update_roi_on_multimeshgrids_widget():
    """Update the ROIs on the multimeshgrids widget."""
    # Update the ROIs if displayed on multimeshgrids widget
    if widgets_list.widget_multimeshgrids is not None:
        if shared.exp.display_roi_in_multimeshgrid:
            widgets_list.widget_multimeshgrids.update_MPL()


def ask_user_to_delete_layer_first():
    """Ask the user to delete the glass layer bound to the ROI first."""
    text = "Delete the glass layer bound to the ROI before deleting it."
    msg = QtWidgets.QMessageBox()
    msg.setText("Error : can not delete ROI")
    msg.setInformativeText(text)
    msg.setIcon(QtWidgets.QMessageBox.Critical)
    msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
    msg.exec_()
