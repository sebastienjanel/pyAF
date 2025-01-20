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

"""Tools to modify the curves."""

import os
import logging
import platform
from PyQt5 import QtWidgets
from ..tools import math_tools
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFButtonGroup
from ..tools.gui_tools import PYAFCheckBox
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from ..widgets.progressbar import Progressbar
from ..widgets.tilt_check import CheckTiltWidget
from ..tools import apply_to_all
import numpy
import copy
from .. import widgets_list
from .. import shared


class CurveModWidget(PYAFWidget):
    """Widget with options to clean the curves or to smooth them."""

    def __init__(self, parent):
        super().__init__(parent, "widget_curve_mod")

        self.logger = logging.getLogger()

        self.sel_st1 = "QWidget#widget_tilt { background-color: #3daee9; }"
        self.sel_st2 = "QWidget#widget_stretch { background-color: #3daee9; }"

        if shared.exp.current_theme == "Dark":
            self.un_st1 = "QWidget#widget_tilt { background-color: #31363b; }"
            self.un_st2 = "QWidget#widget_stretch { background-color: #31363b; }"

        elif shared.exp.current_theme == "Light":
            self.un_st1 = "QWidget#widget_tilt { background-color: #EFF0F1; }"
            self.un_st2 = "QWidget#widget_stretch { background-color: #EFF0F1; }"

        self.gridLayout = QtWidgets.QGridLayout()

        # Buttons to chose curve for tilt correction
        self.widget_tilt_options = QtWidgets.QWidget()
        self.IN_curve_x = PYAFInput(self, "curve_x", "")
        self.IN_curve_y = PYAFInput(self, "curve_y", "Selected curve")
        self.IN_curve_x.changeValue(shared.exp.meshgrid_click_xpos + 1)
        self.IN_curve_y.changeValue(shared.exp.meshgrid_click_ypos + 1)
        self.BT_random = PYAFButton(
            self,
            "button_chose_random_curve",
            "Random curve")
        name = "checkbox_apply_on_all"
        self.CB_apply_on_all = PYAFCheckBox(self, name, "Apply to all")
        self.HL_tilt_options = QtWidgets.QHBoxLayout()
        self.HL_tilt_options.addWidget(self.BT_random)
        self.HL_tilt_options.addWidget(self.IN_curve_x)
        self.HL_tilt_options.addWidget(self.IN_curve_y)
        self.HL_tilt_options.addWidget(self.CB_apply_on_all)
        self.HL_tilt_options.addStretch(1)
        self.widget_tilt_options.setLayout(self.HL_tilt_options)

        # Box for canvas
        self.BX_canvas = QtWidgets.QGroupBox()
        self.HL_canvas = QtWidgets.QHBoxLayout()
        self.canvas_tilt_widget = QtWidgets.QWidget()
        self.canvas_tilt_widget.setFixedSize(504, 288)
        self.MPL_tilt_canvas = PYAFPlot(
            self,
            "tilt_curve_preview",
            self.canvas_tilt_widget,
            [7, 4, 72])
        self.update_MPL("MPL_canvas")
        self.HL_canvas.addWidget(self.canvas_tilt_widget)
        self.HL_canvas.addStretch(1)
        self.BX_canvas.setLayout(self.HL_canvas)

        # Box for tilt options
        self.BX_tilt = QtWidgets.QGroupBox("Tilt")
        self.VL_widget_tilt = QtWidgets.QVBoxLayout()
        self.widget_tilt = QtWidgets.QWidget()
        self.widget_tilt.setObjectName("widget_tilt")
        self.widget_tilt.mousePressEvent = lambda event: self.update_type(0)

        # Limits
        self.HL_tilt_limits = QtWidgets.QHBoxLayout()
        self.input_limit_1 = PYAFInput(self, "input_limit_1", ",")
        self.input_limit_2 = PYAFInput(self, "input_limit_2", "Segment")
        self.HL_tilt_limits.addWidget(self.input_limit_1)
        self.HL_tilt_limits.addWidget(self.input_limit_2)
        self.HL_tilt_limits.addStretch(1)
        self.BT_check_tilt = PYAFButton(self, "check_tilt", "Check")
        self.HL_tilt_limits.addWidget(self.BT_check_tilt)

        # Apply tilt corrections
        self.VL_tilt_apply = QtWidgets.QVBoxLayout()
        self.HL_tilt_apply = QtWidgets.QHBoxLayout()
        self.BT_apply_tilt_corr_trace = \
            PYAFButton(self, "apply_tilt_corr_trace", "Fit from trace")
        self.BT_apply_tilt_corr_retrace = \
            PYAFButton(self, "apply_tilt_corr_retrace", "Fit from retrace")
        self.button_reset_tilt_correction = \
            PYAFButton(self, "button_reset_tilt_correction", "Reset")
        self.HL_tilt_apply.addWidget(self.BT_apply_tilt_corr_trace)
        self.HL_tilt_apply.addWidget(self.BT_apply_tilt_corr_retrace)
        self.HL_tilt_apply.addWidget(self.button_reset_tilt_correction)
        self.HL_tilt_apply.addStretch(1)
        self.VL_tilt_apply.addLayout(self.HL_tilt_limits)
        self.VL_tilt_apply.addLayout(self.HL_tilt_apply)
        self.widget_tilt.setLayout(self.VL_tilt_apply)
        self.VL_widget_tilt.addWidget(self.widget_tilt)
        self.BX_tilt.setLayout(self.VL_widget_tilt)

        # Box for curve stretch
        self.BX_stretch = QtWidgets.QGroupBox("Stretch")
        self.VL_widget_stretch = QtWidgets.QVBoxLayout()
        self.widget_stretch = QtWidgets.QWidget()
        self.widget_stretch.setObjectName("widget_stretch")
        self.widget_stretch.mousePressEvent = lambda event: self.update_type(1)

        self.VL_stretch = QtWidgets.QVBoxLayout()
        self.HL_stretch1 = QtWidgets.QHBoxLayout()
        self.input_stretch1 = PYAFInput(self, "in_stretch1", ",")
        self.input_stretch2 = PYAFInput(self, "in_stretch2", "Segment")
        self.input_stretch_len = PYAFInput(self, "in_stretch_len", "Length")
        self.BTGRP_stretch_mode = PYAFButtonGroup(self, "grp_stretch_mode")
        self.RB_stretch_app = QtWidgets.QRadioButton("Approach")
        self.RB_stretch_ret = QtWidgets.QRadioButton("Retraction")
        self.BTGRP_stretch_mode.addButton(self.RB_stretch_app, 0)
        self.BTGRP_stretch_mode.addButton(self.RB_stretch_ret, 1)
        self.HL_stretch1.addWidget(self.input_stretch1)
        self.HL_stretch1.addWidget(self.input_stretch2)
        self.HL_stretch1.addWidget(self.input_stretch_len)
        self.HL_stretch1.addWidget(self.RB_stretch_app)
        self.HL_stretch1.addWidget(self.RB_stretch_ret)
        self.HL_stretch1.addStretch(1)

        self.HL_stretch2 = QtWidgets.QHBoxLayout()
        self.button_apply_stretch = PYAFButton(self, "apply_stretch", "Apply")
        self.button_reset_stretch = PYAFButton(self, "reset_stretch", "Reset")
        self.HL_stretch2.addWidget(self.button_apply_stretch)
        self.HL_stretch2.addWidget(self.button_reset_stretch)
        self.HL_stretch2.addStretch(1)

        self.VL_stretch.addLayout(self.HL_stretch1)
        self.VL_stretch.addLayout(self.HL_stretch2)
        self.widget_stretch.setLayout(self.VL_stretch)
        self.VL_widget_stretch.addWidget(self.widget_stretch)
        self.BX_stretch.setLayout(self.VL_widget_stretch)

        # Box Savitzky-Golay smoothing
        self.BX_sg_smoothing = QtWidgets.QGroupBox("Savitzky-Golay smoothing")
        self.HL_sg_smoothing = QtWidgets.QHBoxLayout()

        self.input_sg_order = PYAFInput(self, "sg_order", "Order")
        self.input_sg_width = PYAFInput(self, "sg_width", "Width")
        self.CB_sg_uniform = PYAFCheckBox(self, "uniform_sg", "Uniform SG")
        name = "enable_smoothing"
        self.CB_sg_enable = PYAFCheckBox(self, name, "Use smoothing")

        self.HL_sg_smoothing.addWidget(self.input_sg_order)
        self.HL_sg_smoothing.addWidget(self.input_sg_width)
        self.HL_sg_smoothing.addWidget(self.CB_sg_uniform)
        self.HL_sg_smoothing.addWidget(self.CB_sg_enable)
        self.HL_sg_smoothing.addStretch(1)

        self.BX_sg_smoothing.setLayout(self.HL_sg_smoothing)

        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.addWidget(self.widget_tilt_options, 0, 0)
        self.gridLayout.addWidget(self.CB_apply_on_all, 1, 0)
        self.gridLayout.addWidget(self.BX_canvas, 2, 0)
        self.gridLayout.addWidget(self.BX_tilt, 3, 0)
        self.gridLayout.addWidget(self.BX_stretch, 4, 0)
        self.gridLayout.addWidget(self.BX_sg_smoothing, 5, 0)
        self.setLayout(self.gridLayout)

        self.update_widget()

    def button_clicked(self, button):
        """Called when a button is clicked."""
        # Set Focus to set values in inputs
        self.setFocus()

        if button == "button_chose_random_curve":
            widgets_list.widget_main.change_curve(0, 0, "rand")

        elif button == "apply_stretch":
            if shared.exp.display_stretch_mode == "approach":
                mode = shared.exp.display_stretch_mode

                if shared.exp.tilt_correction_all:
                    for i in range(len(shared.exp.list)):
                        data = shared.exp.list[i]
                        if self.check_for_stretch(data, mode):
                            data.stretch_applied_app = True
                else:
                    data = shared.exp.current_data
                    if self.check_for_stretch(data, mode):
                        data.stretch_applied_app = True

            elif shared.exp.display_stretch_mode == "retraction":
                mode = shared.exp.display_stretch_mode

                if shared.exp.tilt_correction_all:
                    for i in range(len(shared.exp.list)):
                        data = shared.exp.list[i]
                        if self.check_for_stretch(data, mode):
                            data.stretch_applied_ret = True
                else:
                    data = shared.exp.current_data
                    if self.check_for_stretch(data, mode):
                        data.stretch_applied_ret = True

            apply_to_all.apply_to_all("calc_options")
            self.parent.update_MPL("MPL_canvas2")
            self.update_widget()

        elif button == "reset_stretch":
            if shared.exp.display_stretch_mode == "approach":
                shared.exp.current_data.stretch_applied_app = False

            elif shared.exp.display_stretch_mode == "retraction":
                shared.exp.current_data.stretch_applied_ret = False

            apply_to_all.apply_to_all("calc_options")
            self.parent.update_MPL("MPL_canvas2")
            self.update_widget()

        elif button == "apply_tilt_corr_trace":
            if shared.exp.tilt_correction_all:
                for i in range(len(shared.exp.list)):
                    data = shared.exp.list[i]
                    data.tilt_applied = "trace"

            else:
                data = shared.exp.current_data
                data.tilt_applied = "trace"

            widgets_list.widget_compute.update_MPL("MPL_canvas")
            self.update_widget()

        elif button == "apply_tilt_corr_retrace":
            if shared.exp.tilt_correction_all:
                for i in range(len(shared.exp.list)):
                    data = shared.exp.list[i]
                    data.tilt_applied = "retrace"

            else:
                data = shared.exp.current_data
                data.tilt_applied = "retrace"

            widgets_list.widget_compute.update_MPL("MPL_canvas")
            self.update_widget()

        elif button == "button_reset_tilt_correction":
            # Get the list of files to reset
            id_list = []
            if shared.exp.tilt_correction_all:
                for i in range(len(shared.exp.list)):
                    id_list.append(i)
            else:
                id_list.append(shared.exp.current_data.unique_id)

            for data_id in id_list:
                data = shared.exp.list[data_id]

                data.tilt_applied = False

            if widgets_list.widget_tilt_check is not None:
                widgets_list.widget_tilt_check.close()

            widgets_list.widget_compute.update_MPL("MPL_canvas")
            self.update_widget()

        elif button == "grp_stretch_mode":
            id_button = self.BTGRP_stretch_mode.checkedId()
            if id_button == 0:
                shared.exp.display_stretch_mode = "approach"
            elif id_button == 1:
                shared.exp.display_stretch_mode = "retraction"

            self.update_MPL("MPL_canvas")
            self.update_widget()

        elif button == "check_tilt":
            if widgets_list.widget_tilt_check is None:
                # Create new widget
                CheckTiltWidget(self)
                widgets_list.widget_tilt_check.show()
            else:
                # Bring to front
                widgets_list.widget_tilt_check.activateWindow()
                widgets_list.widget_tilt_check.raise_()

    def check_for_stretch(self, data, mode):
        """Checks if the stretch can be applied to all the curves of the file.

        If the limits given by the user are in a part of the curve which is
        missing (corrupted curve), the user is told either to discard the curve
        or to chose different fitting limits.

        The stretch is not directly applied to the curves. The stretch_applied
        variables are set to true, and the curves will have the stretch applied
        on the fly when they are called to be displayed or computed.
        """
        # Check if the stretch limits are ok,
        # display error messages or return True

        if mode == "approach":
            limit_1 = data.stretch_app_lim1
            curves = data.curves_approach
        elif mode == "retraction":
            limit_1 = data.stretch_ret_lim1
            curves = data.curves_retraction

        # Create a progressbar
        Progressbar("Stretching ...")
        text = "Stretching " + os.path.basename(data.filename) + " ..."
        widgets_list.widget_progressbar.set_label(text)
        widgets_list.widget_progressbar.set_range(0, len(curves[0]))

        if mode == "approach":
            positions = data.approach_positions
        elif mode == "retraction":
            positions = data.retraction_positions

        # Get the indexes of the limits
        xlimits = []

        _, start_1 = math_tools.get_x_size_and_start_1(curves, positions,
                                                       data.discarded_curves,
                                                       data.nbr_pixels_x,
                                                       data.nbr_pixels_y,
                                                       limit_1)

        for i in range(len(curves)):
            widgets_list.widget_progressbar.update()

            xlimits.append([])
            for j in range(len(curves)):
                if data.discarded_curves[i][j] == 0:
                    found1 = False

                    for k in range(start_1, len(curves[i][j][0]), 1):
                        if curves[i][j][0][k] >= limit_1:
                            xlimits[i].append(k)
                            found1 = True
                            break

                    if not found1:
                        xlimits[i].append(len(curves[i][j][0]) - 1)
                else:
                    xlimits[i].append(len(curves[i][j][0]) - 1)

        # Get the position which is the fartest
        maxposition = 0
        for i in range(len(curves)):
            for j in range(len(curves)):
                if data.discarded_curves[i][j] == 0:
                    current_pos = copy.deepcopy(positions[i][j][0])
                    if current_pos > maxposition:
                        maxposition = current_pos

        amincurves = numpy.amin(xlimits)

        list_curves = []

        if amincurves < maxposition:
            # Make a list of curves
            for i in range(len(curves)):
                for j in range(len(curves)):
                    if data.discarded_curves[i][j] == 0:
                        first = positions[i][j][0]
                        if xlimits[i][j] < first:
                            list_curves.append([i + 1, j + 1])

        # Close the progressbar
        widgets_list.widget_progressbar.close()

        if mode == "approach":
            z_size = data.ramp_size / data.nbr_points_per_curve_approach
        if mode == "retraction":
            z_size = data.ramp_size / data.nbr_points_per_curve_retraction
        curves_txt = ""

        if list_curves != []:
            if maxposition != 0:
                maxposition = str(round(maxposition * z_size, 2))
                txt = "greater than ~ " + maxposition + " nm."
                if list_curves != []:
                    curves_txt = str(list_curves)

            if curves_txt != "" and mode == "approach":
                curves_txt = "Check curves (Approach) : " + curves_txt
            elif curves_txt != "" and mode == "retraction":
                curves_txt = "Check curves (Retract) : " + curves_txt
            filename = os.path.basename(data.filename)

            text = "Some curve(s) have data missing in " + filename + \
                "\r\nYou should correct this with the curve cleaning tool," + \
                " discard them or chose fitting boundaries " + txt + \
                "\r\n" + curves_txt

            QtWidgets.QMessageBox.warning(self, "Error", text)
            return False

        return True

    def checkbox_clicked(self, what):
        """Called when a checkbox is clicked."""
        data = shared.exp.current_data

        if what == "checkbox_apply_on_all":
            state = self.CB_apply_on_all.isChecked()
            if state:
                # Reset everything and tell the user if needed
                found = False
                for i in range(len(shared.exp.list)):
                    data = shared.exp.list[i]

                    if data.tilt_applied:
                        found = True
                        break
                    if data.stretch_applied_app:
                        found = True
                        break
                    if data.stretch_applied_ret:
                        found = True
                        break
                    if data.stretch_applied_ret:
                        found = True
                        break
                    if data.sg_smoothing_enabled:
                        found = True
                        break

                if found:
                    QtWidgets.QMessageBox.warning(
                        self,
                        "Warning",
                        "All the modification were reseted.\r\n" +
                        "Please reapply them.")

                for i in range(len(shared.exp.list)):
                    data = shared.exp.list[i]

                    data.tilt_applied = False
                    data.stretch_applied_app = False
                    data.stretch_applied_ret = False
                    data.stretch_applied_ret = False
                    data.sg_smoothing_enabled = False

                widgets_list.widget_compute.update_MPL("MPL_canvas")
                self.update_widget()

                shared.exp.tilt_correction_all = True
            else:
                shared.exp.tilt_correction_all = False
            self.logger.debug("State, %s", str(state))
            self.update_widget()

        elif what == "uniform_sg":
            value = self.CB_sg_uniform.isChecked()
            data.sg_smoothing_uniform = value
            self.logger.debug("Uniform SG, %s", str(value))

            self.update_MPL("MPL_canvas")
            widgets_list.widget_compute.update_MPL("MPL_canvas")

        elif what == "enable_smoothing":
            check = self.check_for_os_x_bug()
            if check is False:
                # Reset the value and return
                value = shared.exp.current_data.sg_smoothing_width
                self.input_sg_width.changeValue(value)
                return False

            value = self.CB_sg_enable.isChecked()
            data.sg_smoothing_enabled = value
            self.logger.debug("Enable smoothing, %s", str(value))

            self.update_MPL("MPL_canvas")
            widgets_list.widget_compute.update_MPL("MPL_canvas")

    def check_for_os_x_bug(self):
        """Check for bug in smoothing on OS X 10.7.5.

        Check for OS X 10.7.5 and the value of the smoothing. A too high
        value will make the computation crash
        """
        if platform.mac_ver()[0][:-2] == "10.7":
            value = self.input_sg_width.get_int_value()
            if value > 50:
                text = "A too high smoothing width is not allowed on OS " \
                    + "10.7.x. (See FAQ in the docs). Maximum value is 50."
                # Create a message box
                msg = QtWidgets.QMessageBox()
                msg.setText("High smoothing value")
                msg.setInformativeText(text)
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.addButton(QtWidgets.QMessageBox.Ok)
                msg.exec_()
                return False
            else:
                return True
        else:
            return True

    def input_updated(self, what):
        """Called when an input is updated."""
        if what == "curve_x" or what == "curve_y":
            pos_x = self.IN_curve_x.get_int_value()
            pos_y = self.IN_curve_y.get_int_value()
            widgets_list.widget_main.change_curve(pos_x, pos_y)

        elif what == "input_limit_1":
            value = self.input_limit_1.get_int_value()
            shared.exp.current_data.tilt_limit_1 = value
            apply_to_all.apply_to_all("tilt_all")
            self.update_MPL("MPL_canvas")

        elif what == "input_limit_2":
            value = self.input_limit_2.get_int_value()
            shared.exp.current_data.tilt_limit_2 = value
            apply_to_all.apply_to_all("tilt_all")
            self.update_MPL("MPL_canvas")

        elif what == "in_stretch1":
            value = self.input_stretch1.get_int_value()
            if shared.exp.display_stretch_mode == "approach":
                shared.exp.current_data.stretch_app_lim1 = value
            else:
                shared.exp.current_data.stretch_ret_lim1 = value
            self.update_MPL("MPL_canvas")

        elif what == "in_stretch2":
            value = self.input_stretch2.get_int_value()
            if shared.exp.display_stretch_mode == "approach":
                shared.exp.current_data.stretch_app_lim2 = value
            else:
                shared.exp.current_data.stretch_ret_lim2 = value
            self.update_MPL("MPL_canvas")

        elif what == "in_stretch_len":
            value = self.input_stretch_len.get_int_value()
            if shared.exp.display_stretch_mode == "approach":
                shared.exp.current_data.stretch_len_app = value
            else:
                shared.exp.current_data.stretch_len_ret = value
            self.update_MPL("MPL_canvas")

        elif what == "sg_order":
            value = self.input_sg_order.get_int_value()
            shared.exp.current_data.sg_smoothing_order = value
            widgets_list.widget_compute.update_MPL("MPL_canvas")
            self.update_MPL("MPL_canvas")

        elif what == "sg_width":
            check = self.check_for_os_x_bug()
            if check is False:
                # Reset the value and return
                value = shared.exp.current_data.sg_smoothing_width
                self.input_sg_width.changeValue(value)
                return False

            value = self.input_sg_width.get_int_value()
            shared.exp.current_data.sg_smoothing_width = value
            widgets_list.widget_compute.update_MPL("MPL_canvas")
            self.update_MPL("MPL_canvas")

    def update_MPL(self, what):
        """Update the plot."""
        if what == "MPL_canvas":
            self.MPL_tilt_canvas.update_plot()

    def update_widget(self):
        """Update the widget."""
        # Set the stylesheets
        if shared.exp.display_curve_modif_mode == "tilt":
            self.widget_tilt.setStyleSheet(self.sel_st1)
            self.widget_stretch.setStyleSheet(self.un_st2)
        elif shared.exp.display_curve_modif_mode == "stretch":
            self.widget_tilt.setStyleSheet(self.un_st1)
            self.widget_stretch.setStyleSheet(self.sel_st2)

        # Single file ?
        if shared.exp.current_data.is_single:
            self.widget_tilt_options.setParent(None)
        else:
            self.gridLayout.addWidget(self.widget_tilt_options, 0, 0)

        # Apply on all checkbox
        if len(shared.exp.list) < 2:
            self.CB_apply_on_all.setEnabled(False)
        if shared.exp.tilt_correction_all:
            self.CB_apply_on_all.setChecked(True)

        # Update the inputs
        self.input_limit_1.changeValue(shared.exp.current_data.tilt_limit_1)
        self.input_limit_2.changeValue(shared.exp.current_data.tilt_limit_2)

        if shared.exp.display_stretch_mode == "approach":
            val1 = shared.exp.current_data.stretch_app_lim1
            val2 = shared.exp.current_data.stretch_app_lim2
            length = shared.exp.current_data.stretch_len_app
            self.RB_stretch_app.setChecked(True)
        else:
            val1 = shared.exp.current_data.stretch_ret_lim1
            val2 = shared.exp.current_data.stretch_ret_lim2
            length = shared.exp.current_data.stretch_len_ret
            self.RB_stretch_ret.setChecked(True)

        self.input_stretch1.changeValue(val1)
        self.input_stretch2.changeValue(val2)
        self.input_stretch_len.changeValue(length)

        # Radiobuttons
        if shared.exp.display_stretch_mode == "approach":
            self.RB_stretch_app.setChecked(True)
        else:
            self.RB_stretch_ret.setChecked(True)

        # Buttons
        if shared.exp.display_stretch_mode == "approach":
            if shared.exp.current_data.stretch_applied_app:
                self.button_apply_stretch.setEnabled(False)
                self.button_reset_stretch.setEnabled(True)
            else:
                self.button_apply_stretch.setEnabled(True)
                self.button_reset_stretch.setEnabled(False)
        elif shared.exp.display_stretch_mode == "retraction":
            if shared.exp.current_data.stretch_applied_ret:
                self.button_apply_stretch.setEnabled(False)
                self.button_reset_stretch.setEnabled(True)
            else:
                self.button_apply_stretch.setEnabled(True)
                self.button_reset_stretch.setEnabled(False)

        if shared.exp.current_data.tilt_applied:
            self.input_limit_1.setEnabled(False)
            self.input_limit_2.setEnabled(False)
        else:
            self.input_limit_1.setEnabled(True)
            self.input_limit_2.setEnabled(True)

        dt = shared.exp.current_data

        self.input_sg_order.changeValue(dt.sg_smoothing_order)
        self.input_sg_width.changeValue(dt.sg_smoothing_width)
        self.CB_sg_uniform.setChecked(dt.sg_smoothing_uniform)
        self.CB_sg_enable.setChecked(dt.sg_smoothing_enabled)

        if shared.exp.tilt_correction_all:
            found = False
            for dataset in shared.exp.list:
                if dataset.tilt_applied:
                    found = True
                    break
            if found:
                self.BT_apply_tilt_corr_trace.setEnabled(False)
                self.BT_apply_tilt_corr_retrace.setEnabled(False)
                self.button_reset_tilt_correction.setEnabled(True)
                self.BT_check_tilt.setEnabled(True)
            else:
                self.BT_apply_tilt_corr_trace.setEnabled(True)
                self.BT_apply_tilt_corr_retrace.setEnabled(True)
                self.button_reset_tilt_correction.setEnabled(False)
                self.BT_check_tilt.setEnabled(False)
        else:
            if shared.exp.current_data.tilt_applied:
                self.BT_apply_tilt_corr_trace.setEnabled(False)
                self.BT_apply_tilt_corr_retrace.setEnabled(False)
                self.button_reset_tilt_correction.setEnabled(True)
                self.BT_check_tilt.setEnabled(True)
            else:
                self.BT_apply_tilt_corr_trace.setEnabled(True)
                self.BT_apply_tilt_corr_retrace.setEnabled(True)
                self.button_reset_tilt_correction.setEnabled(False)
                self.BT_check_tilt.setEnabled(False)

        self.update_MPL("MPL_canvas")

    def update_type(self, newtype):
        """Update the type of correction."""
        if newtype == 0:
            shared.exp.display_curve_modif_mode = "tilt"
        elif newtype == 1:
            shared.exp.display_curve_modif_mode = "stretch"

        self.input_updated("input_limit_1")
        self.input_updated("input_limit_2")
        self.input_updated("in_stretch1")
        self.input_updated("in_stretch2")
        self.input_updated("in_stretch_len")

        self.update_widget()
        self.update_MPL("MPL_canvas")
