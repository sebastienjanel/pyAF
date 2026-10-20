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

"""Menu for the 3D widget."""

import os
from .. import consts
from qtpy import QtCore, QtWidgets
from ..tools import gui_tools
from ..tools import misc_tools
from .utils.create_layer import create_layer
from .utils.glass_choser import VTKGlassChoserWidget
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFToggle
from .fiducials import FiducialsWidget
from .. import widgets_list
from .. import shared

from .utils.resample import resample_images

# Space between the two buttons of a line
PAIR_SPACING = 6


class MenuLayerListWidget(QtWidgets.QWidget):
    """Widget with options for the 3D widget."""

    def __init__(self, parent, load_list):
        super().__init__()
        policy = QtWidgets.QSizePolicy.Minimum
        self.setSizePolicy(policy, policy)

        self.parent = parent

        # Load the settings
        self.settings = QtCore.QSettings()

        self.smallfont = QtWidgets.QApplication.font()
        self.smallfont.setPointSize(12)

        self.VL = QtWidgets.QVBoxLayout()
        self.VL.setContentsMargins(0, 0, 0, 0)

        # List
        self.tableWidget = QtWidgets.QTableWidget()
        behavior = QtWidgets.QAbstractItemView.SelectRows
        self.tableWidget.setSelectionBehavior(behavior)
        mode = QtWidgets.QAbstractItemView.ExtendedSelection
        self.tableWidget.setSelectionMode(mode)
        self.tableWidget.horizontalHeader().setHighlightSections(False)
        policy = QtCore.Qt.ScrollBarAlwaysOff
        self.tableWidget.setHorizontalScrollBarPolicy(policy)
        self.tableWidget.verticalHeader().hide()

        self.tableWidget.setColumnCount(3)
        self.tableWidget.setHorizontalHeaderLabels(["", "AX", "Name"])
        self.tableWidget.setColumnWidth(0, 30)
        self.tableWidget.setColumnWidth(1, 30)
        self.tableWidget.setColumnWidth(2, 200)

        # Create new rows (only afm layers)
        for i in range(len(load_list)):
            self.tableWidget.insertRow(i)
            self.add_line_in_tablewidget(i, "afm")

        self.tableWidget.selectRow(0)
        self.tableWidget.cellClicked.connect(self.tableWidget_clicked)

        self.VL.addWidget(self.tableWidget)

        # --- Switches, aligned (same width of the labels) -------------
        self.TG_colorscale = PYAFToggle(self, "colorscale", "Color scale")
        self.TG_colorscale.setToolTip("Displays the color scale of the map "
                                      "of the selected file")
        self.TG_colorscale.toggle.setChecked(
            shared.exp.list[0].opengl_display_colorbar)
        self.TG_autospin = PYAFToggle(self, "autospin", "Let it spin !")
        self.TG_autospin.setToolTip("Turns the camera around the scene")
        self.TG_profiles = PYAFToggle(self, "profiles", "Display Profiles")
        self.TG_profiles.setToolTip(
            "Draws the profiles of the maps on the AFM surfaces")
        for toggle in (self.TG_colorscale, self.TG_autospin,
                       self.TG_profiles):
            toggle.label.setFixedWidth(110)

        self.BT_screenshot = PYAFButton(self, "screenshot", "Screenshot")
        self.BT_screenshot.setToolTip("Saves the 3D view as a PNG image")

        # Scales and background: buttons switching between black and white,
        # their label is the next color
        self.BT_scales = PYAFButton(self, "scales", "")
        self.BT_scales.setToolTip("Scales and color scale in black or white")
        self.BT_background = PYAFButton(self, "background", "")
        self.BT_background.setToolTip("Black or white background")
        self.update_color_buttons()

        self.BT_load = PYAFButton(self, "load", "Load tiffs")
        self.BT_load.setToolTip(
            "Adds 2D images (tiff, e.g. optical microscopy) as layers, to "
            "be aligned with the AFM maps")
        self.BT_load_stack = PYAFButton(self, "load_stack", "Load Stack")
        self.BT_load_stack.setToolTip(
            "Adds a 3D stack of images (tiff, e.g. confocal microscopy) as "
            "a volume layer (needs ITK)")
        # Only allowed if itk is imported
        if consts.ALLOW_ITK is False:
            self.BT_load_stack.setEnabled(False)
        self.BT_add_glass = PYAFButton(self, "add_glass", "Add substrate")
        self.BT_add_glass.setToolTip(
            "Adds the plane of a substrate ROI (ROI manager) of an AFM map "
            "as a layer")
        self.BT_fiducials = PYAFButton(self, "fiducials", "Fiducials")
        self.BT_fiducials.setToolTip(
            "Aligns a layer on another with pairs of points picked on both "
            "(fiducial markers)")
        self.BT_lock_top = PYAFButton(self, "lock_top", "View from top")
        self.BT_lock_top.setToolTip("Camera above the scene, looking down")
        self.BT_resample = PYAFButton(self, "resample", "Resample to map")
        self.BT_resample.setToolTip(
            "Cuts the selected tiff image on the area of an AFM map (with "
            "the positions and angles of the layers) and saves it; it can "
            "then be loaded on the map")
        # Can only be enabled for 2D images
        self.BT_resample.setEnabled(False)
        self.BT_isosurface = PYAFButton(self, "isosurface", "New Isosurface")
        self.BT_isosurface.setToolTip(
            "Surface of the selected stack at a threshold of intensity "
            "(binary threshold), as a new layer")
        # Can only be enabled for stacks
        self.BT_isosurface.setEnabled(False)

        self.VL.addWidget(self.TG_colorscale)
        self.VL.addWidget(self.TG_autospin)
        self.VL.addWidget(self.TG_profiles)
        self.add_pair(self.BT_screenshot, self.BT_lock_top)
        self.add_pair(self.BT_scales, self.BT_background)
        self.add_pair(self.BT_load,
                      self.BT_load_stack if consts.ADVANCED else None)
        self.add_pair(self.BT_add_glass, self.BT_fiducials)
        if consts.ADVANCED:
            self.add_pair(self.BT_resample, self.BT_isosurface)
        self.VL.addStretch(1)

        # All the buttons with the width of the largest text (with the two
        # labels of the color buttons), the list as wide as two buttons
        texts = ["Black scales", "White scales", "Black background",
                 "White background"]
        buttons = [self.BT_screenshot, self.BT_lock_top, self.BT_scales,
                   self.BT_background, self.BT_load, self.BT_load_stack,
                   self.BT_add_glass, self.BT_fiducials, self.BT_resample,
                   self.BT_isosurface]
        metrics = self.BT_scales.fontMetrics()
        width = max(metrics.horizontalAdvance(text) + 3 * metrics.xHeight()
                    for text in texts + [b.text() for b in buttons]) + 16
        for button in buttons:
            button.setFixedWidth(width)
        self.tableWidget.setColumnWidth(2, 2 * width + PAIR_SPACING - 60)
        self.resize_table()

        self.setLayout(self.VL)

    def add_pair(self, left, right):
        """Two buttons on a line."""
        HL = QtWidgets.QHBoxLayout()
        HL.setContentsMargins(0, 0, 0, 0)
        HL.setSpacing(PAIR_SPACING)
        HL.addWidget(left)
        if right is not None:
            HL.addWidget(right)
        HL.addStretch(1)
        self.VL.addLayout(HL)

    def add_position_box(self, box):
        """Position box (of the options menu) below the list, as wide."""
        box.setFixedWidth(self.tableWidget.maximumWidth())
        self.VL.insertWidget(1, box)

    def resize_table(self):
        """Height of the table: adapted to the layers, from 4 to 12
        rows."""
        rows = min(max(self.tableWidget.rowCount(), 4), 12)
        table_width = self.tableWidget.horizontalHeader().length()
        self.tableWidget.setFixedSize(table_width, rows * 25 + 25)

    def update_color_buttons(self):
        """Labels of the scales and background buttons: the color they
        switch to."""
        widget = widgets_list.widget_vtk
        white_scales = widget is None or widget.scales_color[0] != 0.0
        self.BT_scales.setText(
            "Black scales" if white_scales else "White scales")
        white_background = shared.exp.opengl_background_color[0] == 1.0
        self.BT_background.setText(
            "Black background" if white_background else "White background")

    def reset_table(self):
        """Recreates a new tableWidget."""

        # Delete all the rows
        for i in range(self.tableWidget.rowCount(), -1, -1):
            self.tableWidget.removeRow(i)

        # Create new rows
        for i in range(len(shared.layer_list)):
            self.tableWidget.insertRow(i)
            self.add_line_in_tablewidget(i, shared.layer_list[i].type)

        self.tableWidget.selectRow(0)

    def add_line_in_tablewidget(self, i, layer_type):
        """Adds a new line to the tablewidget."""
        layer = shared.layer_list[i]

        # Default row height is 30, but it is too big
        self.tableWidget.setRowHeight(i, 25)

        # Checkbox (hide or display all)
        cb = gui_tools.CenteredCellCheckbox(self, "all", str(i))
        self.tableWidget.setCellWidget(i, 0, cb)
        self.tableWidget.cellWidget(i, 0).checkbox.setChecked(layer.visible)

        # Checkbox (hide or display axes)
        cb = gui_tools.CenteredCellCheckbox(self, "axes", str(i))
        self.tableWidget.setCellWidget(i, 1, cb)
        val = layer.axes_visible
        self.tableWidget.cellWidget(i, 1).checkbox.setChecked(val)

        # Name
        label_name = QtWidgets.QLabel()
        label_name.setFont(self.smallfont)
        if layer_type == "afm":
            data_id = layer.afm_id
            label_name.setText(shared.exp.list[data_id].filename)
        elif layer_type == "single_image":
            label_name.setText(layer.filename)
        elif layer_type == "glass":
            roi_id = layer.roi_id
            nm = shared.exp.list[layer.afm_id].filename
            label_name.setText("Glass (Roi " + str(roi_id) + " , " + nm + ")")
        elif layer_type == "stack":
            label_name.setText(shared.layer_list[i].filename)
        cw = gui_tools.CenteredCellwidget(label_name)
        self.tableWidget.setCellWidget(i, 2, cw)

        self.resize_table()

    def checkbox_clicked(self, name, row=None):
        """Called when a checkbox is clicked."""
        if name == "all":
            row = int(row)
            val = self.tableWidget.cellWidget(row, 0).checkbox.isChecked()
            shared.layer_list[row].visible = val
            val = self.tableWidget.cellWidget(row, 1).checkbox.isChecked()
            shared.layer_list[row].axes_visible = val
            shared.layer_list[row].hide_or_display_layer()

        elif name == "axes":
            row = int(row)
            val = self.tableWidget.cellWidget(row, 1).checkbox.isChecked()
            shared.layer_list[row].axes_visible = val
            shared.layer_list[row].hide_or_display_layer()

    def input_updated(self, name):
        """Called when an input is updated (no input)."""

    def tableWidget_clicked(self, row, _):
        """Update GUI and 3D when clicking on a row in the layer list."""
        # Update the tabs
        widgets_list.widget_vtk.menu_options.update_widget()

        # Highlight the picked actor with a red box
        layer = shared.layer_list[row]
        self.parent.canvas.style.HighlightProp(layer.actor)

        # Color scale of the selected AFM layer
        if layer.type == "afm":
            self.parent.update_scalar_bar()

        # Update the current widget
        self.update_widget()

    def update_widget(self, index=None):
        """Update the GUI."""

        # Get the current layer (at load, use index = 0)
        if index is None:
            index = widgets_list.widget_vtk.get_current_layer()

        layer = shared.layer_list[index]

        self.BT_resample.setEnabled(layer.type == "single_image")
        self.BT_isosurface.setEnabled(layer.type == "stack")

        if widgets_list.widget_vtk_binary_threshold is not None:
            widgets_list.widget_vtk_binary_threshold.update_widget()

        val = widgets_list.widget_vtk.display_profiles
        self.TG_profiles.toggle.setChecked(val)

    def button_clicked(self, button):
        """Called when a button is clicked."""
        index = widgets_list.widget_vtk.get_current_layer()

        if button == "screenshot":
            # Get saved path
            settings_path = misc_tools.get_user_path(self.settings)

            filename, _ = QtWidgets.QFileDialog.getSaveFileName(self,
                self.tr("Choose a file name"), settings_path,
                self.tr("PNG (*.png)"))
            if not filename:
                return
            else:
                name = str(filename)

                # Save path for next use
                path = os.path.dirname(name)
                misc_tools.set_user_path(self.settings, path)

                widgets_list.widget_vtk.canvas.TakeScreenshot(name)

        elif button == "scales":
            # Black <-> white scales
            wl = widgets_list.widget_vtk
            black = wl.scales_color[0] != 0.0
            wl.scales_color = [0.0, 0.0, 0.0, 1.0] if black else \
                [1.0, 1.0, 1.0, 1.0]
            for layer in shared.layer_list:
                layer.updateAxesColors()
            # Text of the color scale
            wl.update_scalar_bar()
            self.update_color_buttons()

        elif button == "background":
            # Black <-> white background
            wl = widgets_list.widget_vtk
            value = 0.0 if shared.exp.opengl_background_color[0] == 1.0 \
                else 1.0
            shared.exp.opengl_background_color = [value, value, value, 1.0]
            wl.canvas.renderer.SetBackground(value, value, value)
            if shared.VTK_first is False:
                wl.canvas.interactor.GetRenderWindow().Render()
            self.update_color_buttons()

        elif button == "autospin":
            if widgets_list.widget_vtk.canvas.anim.running is False:
                widgets_list.widget_vtk.canvas.anim.play()
            else:
                widgets_list.widget_vtk.canvas.anim.stop()

        elif button == "load":
            create_layer("single_image")

        elif button == "load_stack":
            create_layer("stack")

        elif button == "add_glass":
            if widgets_list.widget_vtk_glass_choser is None:
                # Ask the user which ROIs he wants to chose
                VTKGlassChoserWidget(self)
                widgets_list.widget_vtk_glass_choser.show()

            else:
                widgets_list.widget_vtk_glass_choser.activateWindow()
                widgets_list.widget_vtk_glass_choser.raise_()

        elif button == "lock_top":
            widgets_list.widget_vtk.canvas.lock_top()

        elif button == "fiducials":
            if widgets_list.widget_fiducials is None:
                # Create new widget
                FiducialsWidget(self)
                widgets_list.widget_fiducials.resize(500, 250)
                widgets_list.widget_fiducials.show()
            else:
                # Bring to front
                widgets_list.widget_fiducials.activateWindow()
                widgets_list.widget_fiducials.raise_()

        elif button == "resample":
            self.resample(index)

        elif button == "isosurface":
            create_layer("isosurface")

    def resample(self, index):
        """Cuts the selected image on the area of an AFM map, saves it
        where the user wants and proposes to load it on the map."""
        if index is None or shared.layer_list[index].type != "single_image":
            QtWidgets.QMessageBox.warning(
                self, "Resample to map",
                "Select the image (tiff) layer to resample.")
            return
        layer = shared.layer_list[index]

        # AFM map: asked if there are several
        afm_layers = [lay for lay in shared.layer_list if lay.type == "afm"]
        if not afm_layers:
            return
        afm = afm_layers[0]
        if len(afm_layers) > 1:
            names = [shared.exp.list[lay.afm_id].filename
                     for lay in afm_layers]
            name, ok = QtWidgets.QInputDialog.getItem(
                self, "Resample to map", "AFM map:", names, 0, False)
            if not ok:
                return
            afm = afm_layers[names.index(name)]

        settings_path = misc_tools.get_user_path(self.settings)
        name = os.path.splitext(layer.filename)[0] + "_resampled.tiff"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save the resampled image",
            os.path.join(settings_path, name), "TIFF (*.tiff *.tif)")
        if not path:
            return
        if not path.lower().endswith((".tif", ".tiff")):
            path += ".tiff"
        misc_tools.set_user_path(self.settings, os.path.dirname(path))

        x, y, angle = resample_images(layer, afm, path)

        answer = QtWidgets.QMessageBox.question(
            self, "Resample to map",
            "The resampled image was saved to:\n" + path +
            "\n\nLoad it on the map in the 3D view?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if answer == QtWidgets.QMessageBox.Yes:
            count = len(shared.layer_list)
            create_layer("single_image", path)
            if len(shared.layer_list) > count:
                # On the map: same center and angle
                new = shared.layer_list[-1]
                new.rotate(None, None, new_angle=angle)
                new.translate(x - new.pos_x, y - new.pos_y, 0)

    def toggle_clicked(self, toogle):
        """Called when a toogle is clicked."""
        wl = widgets_list.widget_vtk
        if toogle == "autospin":
            if widgets_list.widget_vtk.canvas.anim.running is False:
                widgets_list.widget_vtk.canvas.anim.play()
            else:
                widgets_list.widget_vtk.canvas.anim.stop()

        elif toogle == "colorscale":
            # Display or hide the color scale for all the files
            display = self.TG_colorscale.toggle.isChecked()
            for data in shared.exp.list:
                data.opengl_display_colorbar = display
            wl.update_scalar_bar()

        elif toogle == "profiles":
            wl.display_profiles = self.TG_profiles.toggle.isChecked()
            # Refresh the profiles of the AFM layers
            for layer in shared.layer_list:
                if layer.type == "afm":
                    layer.update_profiles()
