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

"""A widget with parameters for the events computation."""

from ... import shared
from ... import widgets_list
from PyQt5 import QtCore, QtWidgets
from ...tools.PYAFWidget import PYAFWidget
from ...tools.gui_tools import PYAFInput
from ...tools.gui_tools import PYAFCheckBox
from ...tools import misc_tools
from ...tools import events_tools


class EventsParamsWidget(PYAFWidget):
    """PYAFWidget containing the input fields for the fits."""

    def __init__(self, parent, name):
        super().__init__(parent, name)

        # Threshold
        self.IN_threshold = PYAFInput(
            self, "threshold", "Threshold ", width=80)
        self.IN_threshold.input.setValidator(misc_tools.validator("UF"))

        # Fit - Size
        self.IN_fit_size = PYAFInput(
            self, "fit_size", "Fit size [pts]", width=80)
        self.IN_fit_size.input.setValidator(misc_tools.validator("UI"))

        # Fit - Display events
        self.CB_fit_event = PYAFCheckBox(
            self,
            "fit_event",
            "Display events")

        # Fit - Display fit event
        self.CB_fit_event_seg = PYAFCheckBox(
            self,
            "fit_event_seg",
            "Display fit")

        # Kernel - Kernel size input
        reg = QtCore.QRegExp("[1-9]\\d{0,1}")
        self.IN_kernel_size = PYAFInput(
            self, "kernel_size", "Kernel half size [pts]")
        self.IN_kernel_size.input.setValidator(
            misc_tools.validator("R", reg=reg))

        # Kernel - Adaptive threshold Option
        self.CB_adaptive_threshold_option = PYAFCheckBox(
            self,
            "adaptive_threshold_option",
            "Adaptive threshold")

        # Kernel - Adaptive kernel window size
        reg = QtCore.QRegExp("[1-9]\\d{0,2}")
        self.IN_adaptive_smoothing_window = PYAFInput(
            self, "adaptive_ws", "Smoothing window size [pts]")

        # Kernel - Adaptive kernel order
        reg = QtCore.QRegExp("[1-9]\\d{0,1}")
        self.IN_adaptive_smoothing_order = PYAFInput(
            self, "adaptive_order", "Smoothing mode")

        # MSF - Window size
        self.IN_window_size = PYAFInput(
            self, "window_size", "Window size [pts]")
        self.IN_window_size.input.setValidator(misc_tools.validator("UI"))

        GB_conv_params = QtWidgets.QGroupBox("Convolution parameters")
        VL_conv_params = QtWidgets.QVBoxLayout()
        GB_conv_params.setLayout(VL_conv_params)

        GB_threshold_params = QtWidgets.QGroupBox("Threshold parameters")
        VL_threshold_params = QtWidgets.QVBoxLayout()
        GB_threshold_params.setLayout(VL_threshold_params)

        GB_fit_params = QtWidgets.QGroupBox("Fit parameters")
        VL_fit_params = QtWidgets.QVBoxLayout()
        GB_fit_params.setLayout(VL_fit_params)

        if self.name == "msf":
            VL_conv_params.addWidget(self.IN_window_size)
        else:
            VL_conv_params.addWidget(self.IN_kernel_size)
            VL_threshold_params.addWidget(self.CB_adaptive_threshold_option)
            VL_threshold_params.addWidget(self.IN_adaptive_smoothing_window)
            VL_threshold_params.addWidget(self.IN_adaptive_smoothing_order)
        VL_threshold_params.addWidget(self.IN_threshold)
        VL_fit_params.addWidget(self.IN_fit_size)
        VL_fit_params.addWidget(self.CB_fit_event)
        VL_fit_params.addWidget(self.CB_fit_event_seg)

        # Align all on Top
        VL_conv_params.addStretch(1)
        VL_threshold_params.addStretch(1)
        VL_fit_params.addStretch(1)

        HL = QtWidgets.QHBoxLayout()
        HL.addWidget(GB_conv_params)
        HL.addWidget(GB_threshold_params)
        HL.addWidget(GB_fit_params)

        self.setLayout(HL)

        self.update_widget()

    def update_widget(self):
        """Generic method to update the GUI of the widget.

        The input fields are filled with their respective values and if needed
        the checkbox for point of contact recomputation is checked/unchecked
        and enabled/disabled.
        """
        data = shared.exp.current_data

        self.CB_fit_event.setChecked(data.display_fit_event)
        self.CB_fit_event_seg.setChecked(data.display_fit_event_seg)

        if self.name == "kernel":
            val = data.adaptive_threshold_option
            self.CB_adaptive_threshold_option.setChecked(val)

            self.IN_adaptive_smoothing_window.setEnabled(True)
            val = data.adaptive_smoothing_window
            self.IN_adaptive_smoothing_window.changeValue(val)

            self.IN_adaptive_smoothing_order.setEnabled(True)
            val = data.adaptive_smoothing_order
            self.IN_adaptive_smoothing_order.changeValue(val)

            if data.adaptive_threshold_option:
                threshold = data.events_kernel_adaptive_threshold
                self.IN_threshold.changeValue(threshold)
            else:
                threshold = data.events_kernel_detection_threshold
                self.IN_threshold.changeValue(threshold)
                self.IN_adaptive_smoothing_window.setEnabled(False)
                self.IN_adaptive_smoothing_order.setEnabled(False)

            kernel = data.kernel_size
            self.IN_kernel_size.changeValue(kernel)

        if self.name == "msf":
            self.IN_fit_size.setEnabled(False)
            window_size = data.events_msf_window_size
            self.IN_window_size.changeValue(window_size)
            threshold = data.events_msf_detection_threshold
            self.IN_threshold.changeValue(threshold)

        # Diseable fit GorupBox if the kernel is displayed
        if data.events_curve_type == "curve":
            self.IN_fit_size.setEnabled(True)
            self.CB_fit_event.setEnabled(True)
            self.CB_fit_event_seg.setEnabled(True)
        else:
            self.IN_fit_size.setEnabled(False)
            self.CB_fit_event.setEnabled(False)
            self.CB_fit_event_seg.setEnabled(False)

        fit_size = data.events_fit_size
        self.IN_fit_size.changeValue(fit_size)

    def input_updated(self, input_field):
        """Method called upon change in an input field of the GUI."""
        data = shared.exp.current_data

        if input_field == "threshold":
            threshold = self.IN_threshold.get_float_value()
            if self.name == "kernel":
                if self.CB_adaptive_threshold_option.isChecked():
                    data.events_kernel_adaptive_threshold = threshold
                else:
                    data.events_kernel_detection_threshold = threshold

            elif self.name == "msf":
                threshold = self.IN_threshold.get_float_value()
                data.events_msf_detection_threshold = threshold

        if input_field == "window_size":
            window_size = self.IN_window_size.get_int_value()
            data.events_msf_window_size = window_size

        if input_field == "fit_size":
            fit_size = self.IN_fit_size.get_int_value()
            data.events_fit_size = fit_size

        if input_field == "kernel_size":
            kernel_size = self.IN_kernel_size.get_int_value()
            data.kernel_size = kernel_size
            data.kernel = events_tools.create_kernel(kernel_size)

        if input_field == 'adaptive_ws':
            window_ws = self.IN_adaptive_smoothing_window.get_int_value()
            data.adaptive_smoothing_window = window_ws

        if input_field == 'adaptive_order':
            window_mode = self.IN_adaptive_smoothing_order.get_int_value()
            data.adaptive_smoothing_order = window_mode

        widgets_list.widget_compute.update_MPL("MPL_canvas")

    def checkbox_clicked(self, checkbox):
        """Called when a checkbox is clicked."""
        data = shared.exp.current_data

        if checkbox == "adaptive_threshold_option":
            option = self.CB_adaptive_threshold_option.isChecked()
            data.adaptive_threshold_option = option

            if self.CB_adaptive_threshold_option.isChecked():
                self.IN_adaptive_smoothing_window.setEnabled(True)
                val = data.adaptive_smoothing_window
                self.IN_adaptive_smoothing_window.changeValue(val)

                self.IN_adaptive_smoothing_order.setEnabled(True)
                val = data.adaptive_smoothing_order
                self.IN_adaptive_smoothing_order.changeValue(val)

                threshold = data.events_kernel_adaptive_threshold
                self.IN_threshold.changeValue(threshold)

            else:
                threshold = data.events_kernel_detection_threshold
                self.IN_threshold.changeValue(threshold)
                self.IN_adaptive_smoothing_window.setEnabled(False)
                self.IN_adaptive_smoothing_order.setEnabled(False)

        elif checkbox == "fit_event":
            data.display_fit_event = self.CB_fit_event.isChecked()

        elif checkbox == "fit_event_seg":
            data.display_fit_event_seg = self.CB_fit_event_seg.isChecked()

        widgets_list.widget_compute.update_MPL("MPL_canvas")
