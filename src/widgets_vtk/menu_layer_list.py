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
from PyQt5 import QtCore, QtWidgets
from ..tools import gui_tools
from ..tools import misc_tools
from .utils.create_layer import create_layer
from .utils.glass_choser import VTKGlassChoserWidget
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFToggle
from .fiducials import FiducialsWidget
from .. import widgets_list
from .. import shared

if consts.ALLOW_ITK:
    from .utils.itk_resample import resample_images


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

        self.HL_checkboxes = QtWidgets.QHBoxLayout()
        self.CB_disp_all = PYAFCheckBox(self, "disp_all", "")
        self.CB_disp_scales_all = PYAFCheckBox(self, "disp_scales_all", "")
        self.HL_checkboxes.insertSpacing(0, 12)
        self.HL_checkboxes.addWidget(self.CB_disp_all)
        self.HL_checkboxes.insertSpacing(2, 20)
        self.HL_checkboxes.addWidget(self.CB_disp_scales_all)
        self.HL_checkboxes.addStretch(0)
        self.CB_disp_all.setChecked(True)
        self.CB_disp_scales_all.setChecked(True)

        self.VL.addWidget(self.tableWidget)
        self.VL.addLayout(self.HL_checkboxes)

        self.VL_screenshot = QtWidgets.QVBoxLayout()
        self.BT_screenshot = PYAFButton(self, "screenshot", "Screenshot")
        label = "Magnification"
        self.input_magn = PYAFInput(self, "magn", label, width=50)
        self.input_magn.changeValue(str(shared.exp.opengl_screenshot_magn))
        self.VL_screenshot.addWidget(self.BT_screenshot)
        self.VL_screenshot.addWidget(self.input_magn)

        # --- Colors ---------------------------------------------------
        self.VL_colors = QtWidgets.QVBoxLayout()
        self.HL_scales = QtWidgets.QHBoxLayout()
        label = "Background color"
        name = "background_color"
        self.BT_background_color = PYAFButton(self, name, label, 200)
        label = "White scales"
        self.BT_scales_white = PYAFButton(self, "scales_white", label, 110)
        label = "Black scales"
        self.BT_scales_black = PYAFButton(self, "scales_black", label, 110)
        self.HL_scales.addWidget(self.BT_scales_white)
        self.HL_scales.addWidget(self.BT_scales_black)
        self.VL_colors.addLayout(self.HL_scales)
        self.VL_colors.addWidget(self.BT_background_color)

        # --------------------------------------------------------------

        self.TG_autospin = PYAFToggle(self, "autospin", "Let it spin !")
        self.BT_load = PYAFButton(self, "load", "Load tiffs")
        self.BT_load_stack = PYAFButton(self, "load_stack", "Load Stack")
        # Only allowed if itk is imported
        if consts.ALLOW_ITK is False:
            self.BT_load_stack.setEnabled(False)
        self.BT_load_hr = PYAFButton(self, "load_hr", "Load HR")
        self.BT_add_glass = PYAFButton(self, "add_glass", "Add substrate")
        self.BT_fiducials = PYAFButton(self, "fiducials", "Fiducials")
        self.BT_lock_top = PYAFButton(self, "lock_top", "View from top")
        self.BT_draw_square = PYAFButton(self, "draw_square", "Draw square")
        self.BT_resample = PYAFButton(self, "resample", "ITK resample")
        self.BT_resample.setEnabled(False)
        # Can be only enabled when drawing the rectangle
        self.BT_isosurface = PYAFButton(self, "isosurface", "New Isosurface")
        # Can only be enabled for stacks
        self.BT_isosurface.setEnabled(False)
        label = "Display/Hide Color Scale"
        self.BT_display_colorscale = PYAFButton(self, "colorscale", label, 200)

        name = "Display Profiles"
        self.CB_display_profiles = PYAFCheckBox(self, "profiles", name)

        self.VL.addLayout(self.VL_screenshot)
        self.VL.addLayout(self.VL_colors)
        self.VL.addWidget(self.TG_autospin)
        self.VL.addWidget(self.BT_load)
        if consts.ADVANCED:
            self.VL.addWidget(self.BT_load_stack)
            self.VL.addWidget(self.BT_load_hr)
        self.VL.addWidget(self.BT_add_glass)
        self.VL.addWidget(self.BT_fiducials)
        self.VL.addWidget(self.BT_lock_top)
        self.VL.addWidget(self.BT_display_colorscale)
        if consts.ADVANCED:
            self.VL.addWidget(self.BT_draw_square)
            self.VL.addWidget(self.BT_resample)
            self.VL.addWidget(self.BT_isosurface)
        self.VL.addWidget(self.CB_display_profiles)
        self.VL.addStretch(1)

        self.setLayout(self.VL)

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
        elif layer_type == "opticsHR":
            label_name.setText("path")
        elif layer_type == "glass":
            roi_id = layer.roi_id
            nm = shared.exp.list[layer.afm_id].filename
            label_name.setText("Glass (Roi " + str(roi_id) + " , " + nm + ")")
        elif layer_type == "stack":
            label_name.setText(shared.layer_list[i].filename)
        cw = gui_tools.CenteredCellwidget(label_name)
        self.tableWidget.setCellWidget(i, 2, cw)

        # Define the size of the table to hide the white space at the bottom
        if self.tableWidget.verticalHeader().length() >= 200:
            table_height = 200
        else:
            table_height = self.tableWidget.verticalHeader().length()
        table_width = self.tableWidget.horizontalHeader().length()
        self.tableWidget.setFixedSize(table_width, table_height + 23)

    def checkbox_clicked(self, name, row=None):
        """Called when a checkbox is clicked."""
        if name == "all":
            row = int(row)
            val = self.tableWidget.cellWidget(row, 0).checkbox.isChecked()
            shared.layer_list[row].visible = val
            val = self.tableWidget.cellWidget(row, 1).checkbox.isChecked()
            shared.layer_list[row].axes_visible = val
            shared.layer_list[row].hide_or_display_layer()

            # Uncheck main if all are unchecked
            found_checked = False
            for i in range(len(shared.layer_list)):
                if self.tableWidget.cellWidget(i, 0).checkbox.isChecked():
                    found_checked = True
                    break

            self.CB_disp_all.setChecked(found_checked)

        elif name == "axes":
            row = int(row)
            val = self.tableWidget.cellWidget(row, 1).checkbox.isChecked()
            shared.layer_list[row].axes_visible = val
            shared.layer_list[row].hide_or_display_layer()

            # Uncheck main if all are unchecked
            found_checked = False
            for i in range(len(shared.layer_list)):
                if self.tableWidget.cellWidget(i, 1).checkbox.isChecked():
                    found_checked = True
                    break

            self.CB_disp_scales_all.setChecked(found_checked)

        elif name == "disp_all":
            value = self.CB_disp_all.isChecked()

            for i in range(len(shared.layer_list)):
                shared.layer_list[i].visible = value
                val = self.tableWidget.cellWidget(i, 1).checkbox.isChecked()
                shared.layer_list[i].axes_visible = val
                shared.layer_list[i].hide_or_display_layer()
                self.tableWidget.cellWidget(i, 0).checkbox.setChecked(value)

        elif name == "disp_scales_all":
            value = self.CB_disp_scales_all.isChecked()

            for i in range(len(shared.layer_list)):
                shared.layer_list[i].axes_visible = value
                shared.layer_list[i].hide_or_display_layer()
                self.tableWidget.cellWidget(i, 1).checkbox.setChecked(value)

        elif name == "profiles":
            val = self.CB_display_profiles.isChecked()
            widgets_list.widget_vtk.display_profiles = val

            # Check all displayed tomogrpahies for profiles and refresh them
            # if needed
            for layer in shared.layer_list:
                if layer.type == "afm":
                    layer.update_profiles()

    def input_updated(self, name):
        """Called when an input is updated."""

        if name == "magn":
            shared.exp.opengl_screenshot_magn = self.input_magn.get_int_value()

    def tableWidget_clicked(self, row, _):
        """Update GUI and 3D when clicking on a row in the layer list."""
        # Update the tabs
        widgets_list.widget_vtk.menu_options.update_widget()

        # Highlight the picked actor with a red box
        layer = shared.layer_list[row]
        self.parent.canvas.style.HighlightProp(layer.actor)

        # Update the current widget
        self.update_widget()

    def update_widget(self, index=None):
        """Update the GUI."""

        # Get the current layer (at load, use index = 0)
        if index is None:
            index = widgets_list.widget_vtk.get_current_layer()

        layer = shared.layer_list[index]

        if consts.ALLOW_ITK is False:
            self.BT_resample.setEnabled(False)

        if layer.type == "stack":
            self.BT_isosurface.setEnabled(True)
        else:
            self.BT_isosurface.setEnabled(False)

        if widgets_list.widget_vtk_binary_threshold is not None:
            widgets_list.widget_vtk_binary_threshold.update_widget()

        val = widgets_list.widget_vtk.display_profiles
        self.CB_display_profiles.setChecked(val)

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

        elif button == "background_color":
            widg = widgets_list.widget_vtk
            current_color = shared.exp.opengl_background_color
            val = misc_tools.ask_user_for_color(current_color)
            if val is not False:
                shared.exp.opengl_background_color = val
                widg.canvas.renderer.SetBackground(val[0], val[1], val[2])

                if shared.VTK_first is False:
                    widg.canvas.interactor.GetRenderWindow().Render()

        elif button == "color_glass":
            index = widgets_list.widget_vtk.get_current_layer()
            dt = shared.exp.list[index]
            val = misc_tools.ask_user_for_color(dt.opengl_glass_color)
            if val is not False:
                dt.opengl_glass_color = [val[0], val[1], val[2], 1.0]

                if shared.VTK_first is False:
                    widg.canvas.interactor.GetRenderWindow().Render()

        elif button == "scales_white":
            widgets_list.widget_vtk.scales_color = [1.0, 1.0, 1.0, 1.0]
            for layer in shared.layer_list:
                layer.updateAxesColors()

        elif button == "scales_black":
            widgets_list.widget_vtk.scales_color = [0.0, 0.0, 0.0, 1.0]
            for layer in shared.layer_list:
                layer.updateAxesColors()

        elif button == "autospin":
            if widgets_list.widget_vtk.canvas.anim.running is False:
                widgets_list.widget_vtk.canvas.anim.play()
            else:
                widgets_list.widget_vtk.canvas.anim.stop()

        elif button == "load":
            create_layer("single_image")

        elif button == "load_stack":
            create_layer("stack")

        elif button == "load_hr":
            widgets_list.widget_vtk.tab_HR.button_clicked("load")

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

        elif button == "colorscale":
            # For the moment, depending on the first layer only
            # (not clean ...)
            if shared.exp.list[0].opengl_display_colorbar:
                shared.exp.list[0].opengl_display_colorbar = False
            else:
                shared.exp.list[0].opengl_display_colorbar = True

            widgets_list.widget_vtk.update_scalar_bar()

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

        elif button == "draw_square":
            wl = widgets_list.widget_vtk

            if wl.rectangle_actor is None:
                wl.canvas.lock_top()
                wl.action_mode = "draw_square"
                self.BT_draw_square.setText("Hide square")
                # self.BT_resample will be enabled after the drawing
            else:
                wl.canvas.renderer.RemoveActor(wl.rectangle_actor)
                wl.rectangle_actor = None
                wl.rectangle_pos = [0, 0, 0, 0]
                self.BT_draw_square.setText("Draw square")
                wl.canvas.interactor.GetRenderWindow().Render()
                self.BT_resample.setEnabled(False)

        elif button == "resample":
            tw = self.tableWidget
            for item in tw.selectionModel().selectedRows():
                index = item.row()
                resample_images(index)

        elif button == "isosurface":
            create_layer("isosurface")

    def toggle_clicked(self, toogle):
        """Called when a toogle is clicked."""
        index = widgets_list.widget_vtk.get_current_layer()

        if toogle == "autospin":
            if widgets_list.widget_vtk.canvas.anim.running is False:
                widgets_list.widget_vtk.canvas.anim.play()
            else:
                widgets_list.widget_vtk.canvas.anim.stop()