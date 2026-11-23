# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
# Copyright Sébastien Janel (2024- ) sebastien.janel@cnrs.fr
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
from qtpy import QtCore, QtWidgets
from ..tools import tilt_tools
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFToggle
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from ..widgets.progressbar import Progressbar
from ..widgets.tilt_check import CheckTiltWidget
from ..tools import apply_to_all
from .. import widgets_list
from .. import shared


class CurveModWidget(PYAFWidget):
    """Widget with options to clean the curves or to smooth them."""

    def __init__(self, parent):
        super().__init__(parent, "widget_curve_mod")

        self.logger = logging.getLogger()

        self.gridLayout = QtWidgets.QGridLayout()

        # Curve displayed (maps only), and the corrections applied to all
        # the files (tilt and smoothing), on one row
        self.widget_tilt_options = QtWidgets.QWidget()
        self.IN_curve_x = PYAFInput(self, "curve_x", "", 40)
        self.IN_curve_y = PYAFInput(self, "curve_y", "", 40)
        self.IN_curve_x.changeValue(shared.exp.meshgrid_click_xpos + 1)
        self.IN_curve_y.changeValue(shared.exp.meshgrid_click_ypos + 1)
        self.BT_random = PYAFButton(
            self, "button_chose_random_curve", "Random", size=0)
        self.HL_tilt_options = QtWidgets.QHBoxLayout()
        self.HL_tilt_options.setContentsMargins(0, 0, 0, 0)
        self.HL_tilt_options.setSpacing(4)
        self.HL_tilt_options.addWidget(self.IN_curve_x)
        self.HL_tilt_options.addWidget(self.IN_curve_y)
        self.HL_tilt_options.addWidget(self.BT_random)
        self.widget_tilt_options.setLayout(self.HL_tilt_options)
        # No labels after the small fields (no gap)
        for field in (self.IN_curve_x, self.IN_curve_y):
            field.label.hide()
        self.TG_apply_on_all = PYAFToggle(
            self, "checkbox_apply_on_all", "Apply to all")
        self.TG_apply_on_all.setToolTip(
            "Tilt correction and smoothing of all the files")
        self.HL_top = QtWidgets.QHBoxLayout()
        self.HL_top.addWidget(self.widget_tilt_options)
        self.HL_top.addStretch(1)
        self.HL_top.addWidget(self.TG_apply_on_all)

        # Curve (no box: no space above it)
        self.canvas_tilt_widget = QtWidgets.QWidget()
        self.canvas_tilt_widget.setFixedSize(504, 288)
        self.MPL_tilt_canvas = PYAFPlot(
            self,
            "tilt_curve_preview",
            self.canvas_tilt_widget,
            [7, 4, 72])
        self.update_MPL("MPL_canvas")

        # Box for tilt options
        self.BX_tilt = QtWidgets.QGroupBox("Tilt")
        self.VL_widget_tilt = QtWidgets.QVBoxLayout()
        self.widget_tilt = QtWidgets.QWidget()

        # Limits
        self.HL_tilt_limits = QtWidgets.QHBoxLayout()
        self.input_limit_1 = PYAFInput(self, "input_limit_1", ",")
        self.input_limit_2 = PYAFInput(self, "input_limit_2", "Segment")
        # Baseline of each curve found, slope from all the curves of the file
        self.TG_tilt_auto = PYAFToggle(self, "tilt_auto", "Auto")
        self.TG_tilt_auto.setToolTip(
            "Automatic correction: the baseline of each curve is found\n"
            "(no limits), its length follows the curve; the slope is\n"
            "the same for the whole file (from all the baselines), each\n"
            "curve keeps its own offset.")
        self.HL_tilt_limits.addWidget(self.TG_tilt_auto)
        self.HL_tilt_limits.addSpacing(12)
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

        # Box Savitzky-Golay smoothing
        self.BX_sg_smoothing = QtWidgets.QGroupBox("Savitzky-Golay smoothing")
        self.HL_sg_smoothing = QtWidgets.QHBoxLayout()

        self.CB_sg_enable = PYAFToggle(self, "enable_smoothing", "Smoothing")
        # The labels before the small fields
        self.input_sg_order = PYAFInput(self, "sg_order", "", 40)
        self.input_sg_width = PYAFInput(self, "sg_width", "", 40)
        self.CB_sg_uniform = PYAFToggle(self, "uniform_sg", "Uniform SG")
        for field in (self.input_sg_order, self.input_sg_width):
            field.label.hide()

        self.HL_sg_smoothing.addWidget(self.CB_sg_enable)
        self.HL_sg_smoothing.addSpacing(16)
        self.HL_sg_smoothing.addWidget(QtWidgets.QLabel("Order"))
        self.HL_sg_smoothing.addWidget(self.input_sg_order)
        self.HL_sg_smoothing.addSpacing(8)
        self.HL_sg_smoothing.addWidget(QtWidgets.QLabel("Width"))
        self.HL_sg_smoothing.addWidget(self.input_sg_width)
        self.HL_sg_smoothing.addSpacing(16)
        self.HL_sg_smoothing.addWidget(self.CB_sg_uniform)
        self.HL_sg_smoothing.addStretch(1)

        self.BX_sg_smoothing.setLayout(self.HL_sg_smoothing)

        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.addLayout(self.HL_top, 0, 0)
        self.gridLayout.addWidget(self.canvas_tilt_widget, 1, 0,
                                  QtCore.Qt.AlignHCenter)
        self.gridLayout.addWidget(self.BX_tilt, 2, 0)
        self.gridLayout.addWidget(self.BX_sg_smoothing, 3, 0)
        self.setLayout(self.gridLayout)

        self.update_widget()

    def button_clicked(self, button):
        """Called when a button is clicked."""
        # Set Focus to set values in inputs
        self.setFocus()

        if button == "button_chose_random_curve":
            widgets_list.widget_main.change_curve(0, 0, "rand")

        elif button in ("apply_tilt_corr_trace", "apply_tilt_corr_retrace"):
            mode = "trace" if button == "apply_tilt_corr_trace" \
                else "retrace"
            if shared.exp.tilt_correction_all:
                files = list(shared.exp.list)
            else:
                files = [shared.exp.current_data]
            if not self.apply_tilt(files, mode):
                return

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
                data.tilt_slope = None

            if widgets_list.widget_tilt_check is not None:
                widgets_list.widget_tilt_check.close()

            widgets_list.widget_compute.update_MPL("MPL_canvas")
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

    def apply_tilt(self, files, mode):
        """Tilt correction of the files fitted on the approach (mode
        "trace") or the retraction ("retrace"): with the limits, or the
        baselines found on each curve (automatic, slope of each file from
        all its curves). Reports the curves with a short baseline or a slope
        far from the others. False if no file can be corrected."""
        notes = []
        corrected = 0
        for data in files:
            if not data.tilt_auto:
                data.tilt_slope = None
                data.tilt_applied = mode
                corrected += 1
                continue
            name = os.path.basename(data.filename)
            Progressbar("Tilt")
            bar = widgets_list.widget_progressbar
            bar.set_label("Baselines of " + name)
            bar.set_range(0, data.nbr_pixels_x)
            try:
                baselines = tilt_tools.data_baselines(data, mode, bar.update)
            finally:
                if widgets_list.widget_progressbar is not None:
                    widgets_list.widget_progressbar.close()
            slope, short, far = tilt_tools.file_slope(
                [b[2:] for b in baselines])
            if slope is None:
                notes.append(name + ": no baseline found, not corrected")
                continue
            data.tilt_slope = slope
            data.tilt_applied = mode
            corrected += 1
            self.logger.debug("Tilt of %s: slope %s from %s curves", name,
                              slope, len(baselines))
            if short or far:
                parts = []
                if short:
                    parts.append("%d curve%s with a short baseline (< %d "
                                 "points)" % (len(short),
                                              "s" if len(short) > 1 else "",
                                              tilt_tools.SHORT_BASELINE))
                if far:
                    parts.append("%d curve%s with a slope far from the "
                                 "others" % (len(far),
                                             "s" if len(far) > 1 else ""))
                # A few positions (map coordinates, from 1)
                positions = [baselines[k][:2]
                             for k in sorted(set(short) | set(far))]
                listed = ", ".join("%d,%d" % (i + 1, j + 1)
                                   for i, j in positions[:10])
                if len(positions) > 10:
                    listed += ", ..."
                notes.append(name + ": " + " and ".join(parts) +
                             " (" + listed + ")")
        if notes:
            QtWidgets.QMessageBox.information(
                self, "Tilt correction",
                "The slope of each file comes from all its baselines; check "
                "these curves:\n\n" + "\n".join(notes))
        return corrected > 0

    def checkbox_clicked(self, what):
        """Called when a checkbox is clicked."""
        data = shared.exp.current_data

        if what == "checkbox_apply_on_all":
            state = self.TG_apply_on_all.isChecked()
            if state:
                # The tilt corrections reset (to apply again to all the
                # files), the smoothing of this file to all the files
                if any(data.tilt_applied for data in shared.exp.list):
                    QtWidgets.QMessageBox.warning(
                        self,
                        "Warning",
                        "The tilt corrections were reset: apply the "
                        "correction again (to all the files).")

                for data in shared.exp.list:
                    data.tilt_applied = False
                    data.tilt_slope = None

                shared.exp.tilt_correction_all = True
                apply_to_all.apply_to_all("tilt_all")
                widgets_list.widget_compute.update_MPL("MPL_canvas")
            else:
                shared.exp.tilt_correction_all = False
            self.logger.debug("State, %s", str(state))
            self.update_widget()

        elif what == "tilt_auto":
            data.tilt_auto = self.TG_tilt_auto.isChecked()
            apply_to_all.apply_to_all("tilt_all")
            self.update_widget()

        elif what == "uniform_sg":
            value = self.CB_sg_uniform.isChecked()
            data.sg_smoothing_uniform = value
            self.logger.debug("Uniform SG, %s", str(value))
            apply_to_all.apply_to_all("tilt_all")

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
            apply_to_all.apply_to_all("tilt_all")

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

        elif what == "sg_order":
            value = self.input_sg_order.get_int_value()
            shared.exp.current_data.sg_smoothing_order = value
            apply_to_all.apply_to_all("tilt_all")
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
            apply_to_all.apply_to_all("tilt_all")
            widgets_list.widget_compute.update_MPL("MPL_canvas")
            self.update_MPL("MPL_canvas")

    def update_MPL(self, what):
        """Update the plot."""
        if what == "MPL_canvas":
            self.MPL_tilt_canvas.update_plot()

    def update_widget(self):
        """Update the widget."""
        # Single file ?
        self.widget_tilt_options.setVisible(
            not shared.exp.current_data.is_single)

        # Apply on all checkbox
        if len(shared.exp.list) < 2:
            self.TG_apply_on_all.setEnabled(False)
        if shared.exp.tilt_correction_all:
            self.TG_apply_on_all.setChecked(True)

        # Update the inputs
        self.input_limit_1.changeValue(shared.exp.current_data.tilt_limit_1)
        self.input_limit_2.changeValue(shared.exp.current_data.tilt_limit_2)

        # Limits: not in automatic mode, nor once the correction applied
        auto = shared.exp.current_data.tilt_auto
        applied = bool(shared.exp.current_data.tilt_applied)
        self.TG_tilt_auto.setChecked(auto)
        self.TG_tilt_auto.setEnabled(not applied)
        self.input_limit_1.setEnabled(not applied and not auto)
        self.input_limit_2.setEnabled(not applied and not auto)

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
