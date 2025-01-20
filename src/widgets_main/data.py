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
First tab of the PYAF software. Is created in the MainWindow class in main.py.

Allows the user to display the list of the files he loaded and to have an
overview of the informations contained in the different files.

"""

from .. import widgets_list
from .. import shared
import os
from .. import experiment
from .. import consts
from PyQt5 import QtCore, QtWidgets
from ..tools import gui_tools
from ..tools import misc_tools
from ..tools import math_tools
from ..tools import events_refresher
from ..utils import meshgrid_menu
from ..tools.gui_tools import PYAFInput
from ..widgets_main.subwidgets.curve_selector import ResultsSelectorWidget
from ..widgets_main.subwidgets.plot_selector import PlotSelectorWidget
from ..tools.gui_tools import PYAFButton
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from ..widgets.progressbar import Progressbar


class DataWidget(PYAFWidget):
    """Init the DataWidget UI.

    The two main plots (meshgrid and curve) are defined here.
    The list of loaded files is created with the create_list method.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_data")

        data = shared.exp.list[shared.exp.id_selected]
        self.file_type_flag = data.file_type in ("JPK (Single File)", "JPK (Force Map)", "JPK (QI)")

        self.list_boxes = None
        self.container = None
        self.list_sel_styles = None
        self.list_unsel_styles = None

        self.smallfont = QtWidgets.QApplication.font()
        self.smallfont.setPointSize(12)
        self.smallfontbold = QtWidgets.QApplication.font()
        self.smallfontbold.setPointSize(12)
        self.smallfontbold.setBold(True)

        self.VL = QtWidgets.QVBoxLayout()

        # List with the files -------------------------------------------------
        self.VL_list = QtWidgets.QVBoxLayout()

        self.list_widget = QtWidgets.QScrollArea()
        self.list_widget.setWidgetResizable(True)
        self.list_widget.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                       QtWidgets.QSizePolicy.Expanding)

        pol = QtCore.Qt.ScrollBarAlwaysOff
        self.list_widget.setHorizontalScrollBarPolicy(pol)
        self.create_new_list()

        self.HL_files_options = QtWidgets.QHBoxLayout()
        label = "Rename"
        self.BT_add_files = PYAFButton(self, "add_files", "+", size=60)
        self.BT_remove_file = PYAFButton(self, "remove_file", "-", size=60)
        self.BT_rename_file = PYAFButton(self, "rename_file", label, size=80)
        self.BT_copy_file = PYAFButton(self, "copy_file", "Copy", size=70)
        self.HL_files_options.addWidget(self.BT_add_files)
        self.HL_files_options.addWidget(self.BT_remove_file)
        self.HL_files_options.addWidget(self.BT_rename_file)
        self.HL_files_options.addWidget(self.BT_copy_file)

        self.VL_list.addWidget(self.list_widget)
        self.VL_list.addLayout(self.HL_files_options)

        self.GL_data = QtWidgets.QGridLayout()

        # General information -------------------------------------------------
        self.box_info1 = QtWidgets.QGroupBox()
        self.box_info1.setFixedSize(500, 200)
        self.HL_box_info = QtWidgets.QHBoxLayout()
        self.VL_box_info1 = QtWidgets.QVBoxLayout()
        self.VL_box_info2 = QtWidgets.QVBoxLayout()
        self.LB_microscope_name = QtWidgets.QLabel()
        self.LB_scanner_file = QtWidgets.QLabel()
        self.LB_version = QtWidgets.QLabel()
        self.VL_box_info1.addWidget(self.LB_microscope_name)
        self.VL_box_info1.addWidget(self.LB_scanner_file)
        self.VL_box_info1.addWidget(self.LB_version)

        self.LB_scan_size = QtWidgets.QLabel()
        self.LB_nbr_points = QtWidgets.QLabel()
        self.LB_ramp_size = QtWidgets.QLabel()
        self.LB_nbr_points_approach = QtWidgets.QLabel()
        self.LB_nbr_points_retraction = QtWidgets.QLabel()
        self.LB_ramp_size = QtWidgets.QLabel()
        self.VL_box_info1.addWidget(self.LB_scan_size)
        self.VL_box_info1.addWidget(self.LB_nbr_points)
        self.VL_box_info1.addWidget(self.LB_ramp_size)
        self.VL_box_info1.addWidget(self.LB_nbr_points_approach)
        self.VL_box_info1.addWidget(self.LB_nbr_points_retraction)

        self.LB_scan_angle = QtWidgets.QLabel()
        self.LB_frame_direction = QtWidgets.QLabel()
        self.LB_xy_closed_loop = QtWidgets.QLabel()

        self.LB_z_closed_loop = QtWidgets.QLabel()
        self.LB_scan_rate = QtWidgets.QLabel()
        self.LB_fwd_scan_velocity = QtWidgets.QLabel()
        self.LB_rev_scan_velocity = QtWidgets.QLabel()
        self.LB_trig_threshold = QtWidgets.QLabel()
        self.LB_retracted_delay = QtWidgets.QLabel()
        self.LB_extended_delay = QtWidgets.QLabel()
        self.VL_box_info2.addWidget(self.LB_scan_angle)
        self.VL_box_info2.addWidget(self.LB_frame_direction)
        self.VL_box_info2.addWidget(self.LB_xy_closed_loop)
        self.VL_box_info2.addWidget(self.LB_z_closed_loop)
        self.VL_box_info2.addWidget(self.LB_scan_rate)
        self.VL_box_info2.addWidget(self.LB_fwd_scan_velocity)
        self.VL_box_info2.addWidget(self.LB_rev_scan_velocity)
        self.VL_box_info2.addWidget(self.LB_trig_threshold)
        self.LB_trig_threshold.mousePressEvent = lambda event: \
            self.button_clicked("thresh_units")

        self.VL_box_info2.addWidget(self.LB_retracted_delay)
        self.VL_box_info2.addWidget(self.LB_extended_delay)

        self.VL_box_info1.addStretch(1)
        self.VL_box_info2.addStretch(1)

        self.HL_box_info.addLayout(self.VL_box_info1)
        self.HL_box_info.addLayout(self.VL_box_info2)

        self.box_info1.setLayout(self.HL_box_info)
        self.box_info1.setFont(self.parent.smallfont)

        # Cantilever & Microscope ---------------------------------------------
        self.box_cantilever = QtWidgets.QGroupBox()
        self.box_cantilever.setFixedSize(330, 200)

        self.IN_spring_constant = PYAFInput(
            self, "spring_constant", "Spring constant (N/m)", 90)
        self.IN_defl_sens = PYAFInput(
            self, "defl_sens", "Deflection sensitivity (nm/V)", 90)
        self.IN_temp = PYAFInput(
            self, "input_temp", "Temperature (C)", 90)
        self.IN_spring_constant.input.setValidator(misc_tools.validator("UF"))
        self.IN_defl_sens.input.setValidator(misc_tools.validator("UF"))
        self.IN_temp.input.setValidator(misc_tools.validator("UF"))

        self.IN_spring_constant.label.mousePressEvent = \
            lambda event: self.change_info_type("spring")
        self.IN_defl_sens.label.mousePressEvent = \
            lambda event: self.change_info_type("defl")
        self.IN_temp.label.mousePressEvent = \
            lambda event: self.change_info_type("temp")

        # Add event filters to watch clicks on the inputs
        # See the eventFilter() method below
        self.IN_spring_constant.input.installEventFilter(self)
        self.IN_defl_sens.input.installEventFilter(self)
        self.IN_temp.input.installEventFilter(self)

        self.BT_reset = PYAFButton(self, "reset", "X", size=60)
        labl = "Apply on all"
        self.BT_apply_on_all = PYAFButton(self, "apply_on_all", labl, size=100)
        self.BT_reset_all = PYAFButton(self, "reset_all", "Reset all", size=80)

        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.addWidget(self.IN_spring_constant, 0, 0, 1, 0)
        self.gridLayout.addWidget(self.IN_defl_sens, 1, 0, 1, 0)
        self.gridLayout.addWidget(self.IN_temp, 2, 0, 1, 0)
        self.gridLayout.addWidget(self.BT_reset, 3, 0)
        self.gridLayout.addWidget(self.BT_apply_on_all, 3, 1)
        self.gridLayout.addWidget(self.BT_reset_all, 3, 2)

        self.VL_box_cantilever = QtWidgets.QVBoxLayout()
        self.VL_box_cantilever.addLayout(self.gridLayout)
        self.VL_box_cantilever.addStretch(1)

        self.box_cantilever.setLayout(self.VL_box_cantilever)

        name = "widget_curve_selector_data"
        self.curve_selector = ResultsSelectorWidget(self, name)
        # self.curve_selector is added to the layout in update_widget()

        # Meshgrid
        self.mesh_size = [350, 350]
        self.canvas_mesh = QtWidgets.QWidget()
        sizes = [self.mesh_size[0] / 72.0, self.mesh_size[1] / 72.0, 72.0]
        self.meshgrid = PYAFPlot(
            self, "meshgrid_data", self.canvas_mesh, sizes)
        self.meshgrid.mpl_connect("button_press_event", self.mesh_press)
        self.meshgrid.tab_id = 0

        # Curve
        self.curve_size = [490, 350]
        sizes = [self.curve_size[0] / 72.0, self.curve_size[1] / 72.0, 72.0]
        if self.file_type_flag:
            self.canvas_stack = QtWidgets.QStackedWidget()
            self.canvas_curve = QtWidgets.QWidget()
            self.canvas_curve_a = QtWidgets.QWidget()
            self.canvas_curve_b = QtWidgets.QWidget()
            self.curve = PYAFPlot(self, "curve_data", self.canvas_curve, sizes)
            self.curve_1 = PYAFPlot(self, "curve_dataa", self.canvas_curve_a, sizes)
            self.curve_2 = PYAFPlot(self, "curve_datab", self.canvas_curve_b, sizes)

            self.canvas_stack.addWidget(self.curve)
            self.canvas_stack.addWidget(self.curve_1)
            self.canvas_stack.addWidget(self.curve_2)

            pol = QtWidgets.QSizePolicy.Minimum
            sizePolicy = QtWidgets.QSizePolicy(pol, pol)
            self.curve.setSizePolicy(sizePolicy)
            self.curve_1.setSizePolicy(sizePolicy)
            self.curve_2.setSizePolicy(sizePolicy)
            self.canvas_stack.setSizePolicy(sizePolicy)
            self.canvas_mesh.setSizePolicy(sizePolicy)
            self.meshgrid.setSizePolicy(sizePolicy)

        else:
            self.canvas_curve = QtWidgets.QWidget()
            self.curve = PYAFPlot(self, "curve_data", self.canvas_curve, sizes)

        self.HL_box_infos = QtWidgets.QHBoxLayout()
        self.HL_box_infos.addWidget(self.box_info1)
        self.HL_box_infos.addWidget(self.box_cantilever)
        self.HL_box_infos.addStretch(1)

        if self.file_type_flag:
            name = "widget_plot_selector_data"
            self.plot_selector = PlotSelectorWidget(self, name)

        self.GL_data.addLayout(self.HL_box_infos, 0, 0, 1, 0)

        if self.file_type_flag:
            self.GL_data.addWidget(self.plot_selector, 1, 1)

        self.GL_data.addWidget(self.meshgrid, 2, 0)

        if self.file_type_flag:
            self.GL_data.addWidget(self.canvas_stack, 2, 1)

        else:
            self.GL_data.addWidget(self.curve, 2, 1)

        self.HL = QtWidgets.QHBoxLayout()
        self.HL.addLayout(self.VL_list)
        self.HL.addLayout(self.GL_data)

        self.VL.addLayout(self.HL)

        self.setLayout(self.VL)

        self.update_widget()

    def eventFilter(self, obj, event):
        """Listens to the events.

        Triggers an update of the info fields when clicking in the input field
        for the spring constant, deflection sensitivity or temperature.
        """
        # Check if its a click
        if event.type() == event.MouseButtonPress:
            # Check if it's a click on an input field
            if self.IN_spring_constant.input == obj:
                self.change_info_type("spring")

            elif self.IN_defl_sens.input == obj:
                self.change_info_type("defl")

            elif self.IN_temp.input == obj:
                self.change_info_type("temp")

        # Just do the default
        return super().eventFilter(obj, event)

    def create_new_list(self):
        """Create a new list with the files on the right of the data tab."""
        self.container = QtWidgets.QFrame()

        # Reset list of boxes
        self.list_boxes = []

        # Set a new layout
        self.layout = QtWidgets.QVBoxLayout(self.container)
        self.layout.setContentsMargins(0, 0, 0, 0)

        # Add the widgets to the QScrollArea
        for i in range(len(shared.exp.list)):
            box = gui_tools.BoxData(i)
            self.list_boxes.append(box)
            self.layout.addWidget(box)

        self.layout.addStretch(1)
        self.list_widget.setWidget(self.container)

    def update_widget(self):
        """Update the widget."""
        data = shared.exp.list[shared.exp.id_selected]

        # For single files don't display the curve selector
        if data.is_single:
            self.curve_selector.setParent(None)
        else:
            self.GL_data.addWidget(self.curve_selector, 1, 0)
            self.curve_selector.update_widget()

        # Update the labels
        text = "Microscope : " + str(data.microscope_name)
        self.LB_microscope_name.setText(text)
        text = "Config. file : " + str(data.scanner_file)
        self.LB_scanner_file.setText(text)
        text = "Software version : " + str(data.version)
        self.LB_version.setText(text)

        x = str(round(data.scan_size_x / 1000.0, 4))
        y = str(round(data.scan_size_y / 1000.0, 4))
        self.LB_scan_size.setText("Scan size : " + x + " x " + y + " \u03bcm")
        self.LB_nbr_points.setText("XY samples : "
                                   + str(data.nbr_pixels_x)
                                   + " x " + str(data.nbr_pixels_y))
        r_size = str(round(data.ramp_size / 1000.0, 4))
        self.LB_ramp_size.setText("Ramp size : " + r_size + " \u03bcm")
        self.LB_nbr_points_approach.setText(
            "Nbr. points approach : " +
            str(data.nbr_points_per_curve_approach) + " points")
        self.LB_nbr_points_retraction.setText(
            "Nbr. points retraction : " +
            str(data.nbr_points_per_curve_retraction) + " points")

        if data.is_single:
            self.LB_scan_angle.setText("Scan Angle : None")
        else:
            self.LB_scan_angle.setText(
                "Scan Angle : " + str(data.scan_angle) + "\u00b0")
        self.LB_frame_direction.setText(
            "Scan direction : " + str(data.frame_direction))
        self.LB_xy_closed_loop.setText(
            "XY Closed Loop: " + str(data.xy_closed_loop))
        self.LB_z_closed_loop.setText(
            "Z Closed Loop : " + str(data.z_closed_loop))
        self.LB_scan_rate.setText(
            "Scan rate : " + str(data.scan_rate) + " Hz")
        app_vel = str(round(data.approach_velocity / 1000.0, 4))
        ret_vel = str(round(data.retraction_velocity / 1000.0, 4))
        self.LB_fwd_scan_velocity.setText(
            "Forward velocity : " + app_vel + " \u03bcm/s")
        self.LB_rev_scan_velocity.setText(
            "Retrace velocity : " + ret_vel + " \u03bcm/s")

        if data.retracted_delay is not None:
            # JPK QI files don't have these values
            val = data.retracted_delay
            text = "Retracted delay : " + str(val * 1000.0) + " ms"
            self.LB_retracted_delay.setText(text)
            val = data.extended_delay
            text = "Extended delay : " + str(val * 1000.0) + " ms"
            self.LB_extended_delay.setText(text)
        else:
            self.LB_retracted_delay.setText("Retracted delay : None")
            self.LB_extended_delay.setText("Extended delay : None")

        self.IN_spring_constant.changeValue(
            round(data.spring_constant * 1e9, 4))
        self.IN_defl_sens.changeValue(round(data.deflection_sensitivity, 2))
        self.IN_temp.changeValue(data.temperature - 273.15)

        self.update_GUI("all")

        # If there is only one file in the list, disable the remove button
        if len(shared.exp.list) <= 1:
            self.BT_remove_file.setEnabled(False)

        # Always reupdate the fonts
        self.setFont(self.smallfont)

        # Update the meshgrid
        self.meshgrid.update_plot()
        self.curve.update_plot()
        if self.file_type_flag:
            self.curve_1.update_plot()
            self.curve_2.update_plot()

    def change_info_type(self, newtype):
        """Change the type of information which is selected.

        Can be : spring constant, deflection sensitivity or temperature.
        Values : newtype = "spring", "defl", "temp"
        """
        shared.exp.info_selected = newtype
        self.update_GUI("info_selected")

    def update_GUI(self, what):
        """Update the GUI."""
        data = shared.exp.current_data

        if what == "single_or_forcevolume" or what == "all":
            if data.is_single:
                self.meshgrid.setParent(None)
            else:
                self.GL_data.addWidget(self.meshgrid, 2, 0)

        if what == "info_selected" or what == "all":
            info = shared.exp.info_selected

            if info == "spring":
                self.IN_defl_sens.setFont(self.smallfont)
                self.IN_temp.setFont(self.smallfont)
                self.IN_spring_constant.label.setFont(self.smallfontbold)
            elif info == "defl":
                self.IN_spring_constant.setFont(self.smallfont)
                self.IN_temp.setFont(self.smallfont)
                self.IN_defl_sens.label.setFont(self.smallfontbold)
            elif info == "temp":
                self.IN_spring_constant.setFont(self.smallfont)
                self.IN_defl_sens.setFont(self.smallfont)
                self.IN_temp.label.setFont(self.smallfontbold)

            if data._spring_constant is None and info == "spring":
                self.BT_reset.setEnabled(False)
            if data._spring_constant is not None and info == "spring":
                self.BT_reset.setEnabled(True)

            if data._deflection_sensitivity is None and info == "defl":
                self.BT_reset.setEnabled(False)
            if data._deflection_sensitivity is not None and info == "defl":
                self.BT_reset.setEnabled(True)

            if data._temperature is None and info == "temp":
                self.BT_reset.setEnabled(False)
            elif data._temperature is not None and info == "temp":
                self.BT_reset.setEnabled(True)

        if what == "thresh_units" or what == "all":
            # Trigger threshold
            if data.trig_threshold:
                value = round(data.trig_threshold, 1)  # Is stored in V

                if shared.exp.infowidget_threshold_units == "V":
                    text = str(value) + " V"

                elif shared.exp.infowidget_threshold_units == "nm":
                    val = value * data.deflection_sensitivity
                    text = str(round(val, 2)) + " nm"

                elif shared.exp.infowidget_threshold_units == "nN":
                    val = value * data.deflection_sensitivity * \
                        data.spring_constant * 1e9
                    text = str(round(val, 4)) + " nN"

            else:
                text = str(None)

            self.LB_trig_threshold.setText("Trigger threshold : " + text)

        if what == "list_widget" or what == "all":
            # Redefine all the stylesheets
            self.list_sel_styles = []
            self.list_unsel_styles = []
            for i in range(len(shared.exp.list)):
                color = " { background-color: #3daee9; }"
                style = "QWidget#data_" + str(i) + color
                self.list_sel_styles.append(style)
                if shared.exp.current_theme == "Dark":
                    color = " { background-color: #31363b; }"
                elif shared.exp.current_theme == "Light":
                    color = " { background-color: #EFF0F1; }"
                style = "QWidget#data_" + str(i) + color
                self.list_unsel_styles.append(style)

            # Update all the style sheets
            self.update_style_sheets()

    def update_style_sheets(self):
        """Update the stylesheets.

        Used for the boxes with the files in the file's list.
        """
        current_id = shared.exp.id_selected

        style = self.list_sel_styles[current_id]
        self.list_boxes[current_id].setStyleSheet(style)

        for i in range(len(shared.exp.list)):
            if i != current_id:
                self.list_boxes[i].setStyleSheet(self.list_unsel_styles[i])

    def mesh_press(self, event):
        """Called when the meshgrid is clicked."""
        if event.button == 1 and event.inaxes == self.meshgrid.canvas.axes:
            xpos, ypos = misc_tools.get_position_on_meshgrid(event)

            widgets_list.widget_main.change_curve(xpos, ypos)

        elif event.button == 3 and event.inaxes == self.meshgrid.canvas.axes:
            # Get position in screen coordinates (top left position of canvas)
            globalpos = self.meshgrid.mapToGlobal(QtCore.QPoint(0, 0))
            # Recalculate positions with canvas position
            x = globalpos.x() + event.x
            y = globalpos.y() + (self.mesh_size[1] - event.y)
            pos = QtCore.QPoint(x, y)
            meshgrid_menu.create_meshgrid_menu(self, "data", pos)

    def button_clicked(self, button):
        """Called when a button is clicked."""
        data = shared.exp.current_data
        no_update = False

        if button == "reset":
            if shared.exp.info_selected == "spring":
                data.spring_constant = None
                self.IN_spring_constant.changeValue(data.spring_constant)
                self.BT_reset.setEnabled(False)

            elif shared.exp.info_selected == "defl":
                data.deflection_sensitivity = None
                self.IN_defl_sens.changeValue(data.deflection_sensitivity)
                self.BT_reset.setEnabled(False)
                self.curve.update_plot()

            elif shared.exp.info_selected == "temp":
                data.temperature = None

                self.IN_temp.changeValue(data.temperature)
                self.BT_reset.setEnabled(False)

        elif button == "thresh_units":
            if shared.exp.infowidget_threshold_units == "V":
                shared.exp.infowidget_threshold_units = "nm"

            elif shared.exp.infowidget_threshold_units == "nm":
                shared.exp.infowidget_threshold_units = "nN"

            elif shared.exp.infowidget_threshold_units == "nN":
                shared.exp.infowidget_threshold_units = "V"

            no_update = True
            self.update_GUI("thresh_units")

        elif button == "apply_on_all":
            if shared.exp.info_selected == "spring":
                value = self.IN_spring_constant.get_float_value()
                for data in shared.exp.list:
                    data.spring_constant = value

            elif shared.exp.info_selected == "defl":
                value = self.IN_defl_sens.get_float_value()
                for data in shared.exp.list:
                    data.deflection_sensitivity = value

            elif shared.exp.info_selected == "temp":
                value = self.IN_temp.get_float_value()
                for data in shared.exp.list:
                    data.temperature = value + 273.15

        elif button == "reset_all":
            if shared.exp.info_selected == "spring":
                for data in shared.exp.list:
                    data.spring_constant = None

            elif shared.exp.info_selected == "defl":
                for data in shared.exp.list:
                    data.deflection_sensitivity = None

            elif shared.exp.info_selected == "temp":
                for data in shared.exp.list:
                    data.temperature = None

        elif button == "add_files":
            self.parent.load_more()

        elif button == "rename_file":
            newname, ok = QtWidgets.QInputDialog.getText(self,
                                                     "Input",
                                                     "Rename as : ",
                                                     QtWidgets.QLineEdit.Normal,
                                                     data.filename)
            if ok:
                if newname == data.filename:
                    text = "You should give another name to your copy !"
                    QtWidgets.QMessageBox.warning(self, "Error", text)
                else:
                    # Update the value
                    data.filename = str(newname)

                    # Update the GUI
                    index = shared.exp.id_selected
                    widgets_list.widget_main.update_names(force_id=index)

        elif button == "remove_file":
            # Due to caching in pytables, you will get some errors doing the
            # remove step. Flushing the file is not enough. Adding a print
            # statement here will do it, but this is not a real solution.
            # The best way is to close and to reopen the file.

            Progressbar()
            widgets_list.widget_progressbar.set_label("Removing ...")
            widgets_list.widget_progressbar.set_range(0, 10)

            # Close file (cleans up the nodes)
            shared.exp.temp_file.close_file()
            # Reopen temp_file and update values in exp
            shared.exp.temp_file.open_file()
            for unique_id in range(len(shared.exp.list)):
                shared.exp.list[unique_id].update()

            widgets_list.widget_progressbar.update()

            tfile = shared.exp.temp_file

            # The index of the item to be deleted
            index = shared.exp.current_data.unique_id
            # Remove node from temp file
            tfile.file.remove_node("/data/_" + str(index), recursive=True)

            widgets_list.widget_progressbar.update()

            list_of_changed_ids = []
            for i in range(index + 1, len(shared.exp.list), 1):
                # Rename nodes and set new unique ids
                tfile.file.renameNode("/data/_" + str(i), "_" + str(i - 1))
                ind = shared.exp.list[i].unique_id
                list_of_changed_ids.append([ind, ind - 1])
                shared.exp.list[i].unique_id = shared.exp.list[i].unique_id - 1

            widgets_list.widget_progressbar.update()

            # Flush file
            tfile.flush_file()

            # Delete the object
            del shared.exp.list[index]

            # If there is only one dataset remaining, disable the remove button
            if len(shared.exp.list) == 1:
                self.BT_remove_file.setEnabled(False)

            widgets_list.widget_progressbar.update()

            # Remove lines from results data tables and shared storage
            for i in range(len(shared.exp.results_list) - 1, -1, -1):
                if shared.exp.results_list[i].data_id == index:
                    del shared.exp.results_list[i]
                    del shared.single_data[i]
                    del shared.single_values[i]
                    if shared.single_factors is not None:
                        del shared.single_factors[i]
                    if shared.single_factors is not None:
                        del shared.single_factors[i]
                    if shared.single_pdfs_x is not None:
                        del shared.single_pdfs_x[i]
                        del shared.single_pdfs_y[i]
                    widgets_list.widget_results_single.tableWidget.removeRow(i)

            widgets_list.widget_progressbar.update()

            # The indices have to be shifted in the results_list
            for ids in list_of_changed_ids:
                for result in shared.exp.results_list:
                    if result.data_id == ids[0]:
                        result.data_id = ids[1]

            widgets_list.widget_progressbar.update()

            # Reorganize all the ids
            for i in range(len(shared.exp.results_list)):
                shared.exp.results_list[i].result_id = i

            widgets_list.widget_progressbar.update()

            # Close file (cleans up the nodes)
            shared.exp.temp_file.close_file()
            # Reopen temp_file and update values in exp
            shared.exp.temp_file.open_file()
            for unique_id in range(len(shared.exp.list)):
                shared.exp.list[unique_id].update()

            widgets_list.widget_progressbar.update()

            # Change value in the data choser (just go to file 0)
            shared.exp.id_selected = 0
            widgets_list.widget_results.list_exp.removeItem(index)
            widgets_list.widget_compute.list_exp2.removeItem(index)
            self.create_new_list()
            self.parent.file_changed("box")

            widgets_list.widget_progressbar.update()

            # Update the results
            widgets_list.widget_results_single.reset_table()
            widgets_list.widget_results_single.button_clicked("refresh_hist")
            # Perhaps there is nothing to display in the results, so we have
            # to grey out the tab
            self.parent.update_GUI("tabs")

            widgets_list.widget_progressbar.update()

            # Set this value to True. Will display a message before saving
            # telling the user saving could be slow. Then, PYAF will repack
            # the hdf5 temp file to remove definitively the dead nodes.
            shared.files_removed = True

            widgets_list.widget_progressbar.close()

        elif button == "copy_file":
            old_id = shared.exp.id_selected
            ids = []
            for dt in shared.exp.list:
                ids.append(dt.unique_id)
            new_id = max(ids) + 1

            if consts.UNIT_TESTING is False:
                newname, ok = QtWidgets.QInputDialog.getText(
                    self, "Input", "Copy as : ",
                    QtWidgets.QLineEdit.Normal, data.filename)
            else:
                # In case of unit testing just create an empty dataset
                ok = True
                newname = "test1"

            if ok:
                if newname == os.path.basename(str(data.filename)):
                    QtWidgets.QMessageBox.warning(
                        self,
                        "Error",
                        "You should give another name to your file !")
                else:
                    # The complete new name of the dataset (with path)
                    path = os.path.dirname(str(data.filename))
                    newname = str(path + newname)

                    # --- Copying ---

                    # Fill empty dataset
                    shared.exp.addData(newname, new_id)
                    newdata = shared.exp.list[new_id]
                    t_file = shared.exp.temp_file

                    # Copy the variables
                    for key in list(data.__dict__.keys()):
                        inlist = math_tools.in_list(data.dontsave, key)
                        if inlist is False and key != "dontsave":
                            misc_tools.setattr_special(newdata, key,
                                                       data.__dict__[key])
                    newdata.unique_id = new_id
                    newdata.filename = newname

                    # Create tables in temp file
                    t_file.create_storage_group(new_id)

                    # Create tables
                    t_file.create_data_tables(new_id)
                    t_file.create_tables_for_results(new_id, "all")

                    # Copy data
                    t_file.copy_data(data, old_id, new_id)

                    # Flush file
                    t_file.flush_file()

                    # Update the events array (with a progressbar)
                    events_refresher.update_events(data_id=new_id)

                    # Update the values in exp
                    newdata.update()

                    # Add a new result in the result's list
                    # (only if something has been computed)
                    if newdata.calculated:
                        experiment.add_result(new_id, newdata.filename)

                    # --- Displaying ---

                    # Select the new dataset
                    shared.exp.id_selected = new_id

                    # Select the new dataset
                    widg = widgets_list
                    name = str(os.path.basename(newname))
                    index = shared.exp.id_selected
                    widg.widget_results.list_exp.addItem(name)
                    widg.widget_results.list_exp.setCurrentIndex(index)
                    widg.widget_compute.list_exp2.addItem(name)
                    widg.widget_compute.list_exp2.setCurrentIndex(index)
                    self.create_new_list()
                    self.parent.file_changed("box")

                    # Scroll down to last position
                    pos = self.list_widget.verticalScrollBar().maximum()
                    self.list_widget.verticalScrollBar().setValue(pos)

                    # If there is more than one dataset we can delete it,
                    # so enable the button
                    if len(shared.exp.list) > 1:
                        self.BT_remove_file.setEnabled(True)

                    # Update the table with the data in the result's tab
                    widg.widget_results_single.reset_table()
                    widg.widget_results_single.button_clicked("refresh_hist")

        if not no_update:
            self.update_widget()

    def input_updated(self, field):
        """Called when an input field is updated."""
        data = shared.exp.current_data

        if field == "spring_constant":
            value = self.IN_spring_constant.get_str_value(True)
            # Compare strings, comparing floats won't work well ...
            if value != str(round(data.spring_constant * 1e9, 4)):
                data.spring_constant = float(value)
                self.update_GUI("info_selected")
                QtWidgets.QMessageBox.warning(
                    self,
                    "Info",
                    "You should recalculate the stiffness.")
        elif field == "defl_sens":
            value = self.IN_defl_sens.get_str_value(True)
            if value != str(round(data.deflection_sensitivity, 2)):
                data.deflection_sensitivity = float(value)
                self.update_GUI("info_selected")
                self.curve.update_plot()
        elif field == "input_temp":
            value = self.IN_temp.get_str_value(True)
            if value != str(data.temperature):
                data.temperature = float(value) + 273.15
                self.update_GUI("info_selected")

        self.update_widget()
