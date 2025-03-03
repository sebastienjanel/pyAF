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

"""Displays the used parameters for the computation."""

from PyQt5 import QtWidgets
from ..tools.PYAFWidget import PYAFWidget
from .. import shared
from .. import widgets_list


class UsedParametersWidget(PYAFWidget):
    """Displays the used parameters for the computation."""

    def __init__(self, parent):
        name = "widget_used_parameters"
        super().__init__(parent, name)

        HL = QtWidgets.QHBoxLayout()

        box_curve_modif = QtWidgets.QGroupBox("General informations")
        box_elasticity = QtWidgets.QGroupBox("Computation method")
        box_work_and_rupture_force = QtWidgets.QGroupBox("Work and rupture force")
        box_events = QtWidgets.QGroupBox("Events")

        self.fontbold = QtWidgets.QApplication.font()
        self.fontbold.setBold(True)

        # General information
        self.LB_smoothing_title = QtWidgets.QLabel("Smoothing")
        self.LB_smoothing_title.setFont(self.fontbold)
        self.LB_smoothing_on = QtWidgets.QLabel()
        self.LB_smoothing_uniform = QtWidgets.QLabel()
        self.LB_smoothing_width = QtWidgets.QLabel()
        self.LB_smoothing_coef = QtWidgets.QLabel()
        self.LB_other_title = QtWidgets.QLabel("Other")
        self.LB_other_title.setFont(self.fontbold)
        self.LB_defl_sens = QtWidgets.QLabel()
        self.LB_spring_constant = QtWidgets.QLabel()

        VL = QtWidgets.QVBoxLayout()
        VL.addWidget(self.LB_smoothing_title)
        VL.addWidget(self.LB_smoothing_on)
        VL.addWidget(self.LB_smoothing_uniform)
        VL.addWidget(self.LB_smoothing_width)
        VL.addWidget(self.LB_smoothing_coef)
        VL.addWidget(self.LB_other_title)
        VL.addWidget(self.LB_defl_sens)
        VL.addWidget(self.LB_spring_constant)
        VL.addStretch(1)

        box_curve_modif.setLayout(VL)

        # Elasticity
        self.LB_model_selected = QtWidgets.QLabel("Not computed")
        self.LB_poc_skip_start = QtWidgets.QLabel()
        self.LB_poc_fit_length = QtWidgets.QLabel()
        self.LB_poc_refit_option = QtWidgets.QLabel()
        self.LB_poc_noise_multiplicator = QtWidgets.QLabel()
        self.LB_poc_refit_times = QtWidgets.QLabel()
        self.LB_indentation_start = QtWidgets.QLabel()
        self.LB_indentation_step = QtWidgets.QLabel()
        self.LB_indentation_stop = QtWidgets.QLabel()

        VL = QtWidgets.QVBoxLayout()
        VL.addWidget(self.LB_model_selected)
        VL.addWidget(self.LB_poc_skip_start)
        VL.addWidget(self.LB_poc_fit_length)
        VL.addWidget(self.LB_poc_noise_multiplicator)
        VL.addWidget(self.LB_poc_refit_option)
        VL.addWidget(self.LB_poc_refit_times)
        VL.addWidget(self.LB_indentation_start)
        VL.addWidget(self.LB_indentation_step)
        VL.addWidget(self.LB_indentation_stop)
        VL.addStretch(1)

        box_elasticity.setLayout(VL)

        # Work
        self.LB_joc_skip_start = QtWidgets.QLabel("Not computed")
        self.LB_joc_fit_length = QtWidgets.QLabel()
        self.LB_joc_refit_option = QtWidgets.QLabel()
        self.LB_joc_refit_times = QtWidgets.QLabel()
        self.LB_joc_noise_multiplicator = QtWidgets.QLabel()

        VL = QtWidgets.QVBoxLayout()
        VL.addWidget(self.LB_joc_skip_start)
        VL.addWidget(self.LB_joc_fit_length)
        VL.addWidget(self.LB_joc_refit_option)
        VL.addWidget(self.LB_joc_refit_times)
        VL.addWidget(self.LB_joc_noise_multiplicator)
        VL.addStretch(1)

        box_work_and_rupture_force.setLayout(VL)

        # Events
        self.LB_events_detection = QtWidgets.QLabel("Events detection")
        self.LB_events_detection.setFont(self.fontbold)
        self.LB_model_detection_selected = QtWidgets.QLabel("Not computed")
        self.LB_kernel_size = QtWidgets.QLabel()
        self.LB_msf_window_size = QtWidgets.QLabel()
        self.LB_threshold = QtWidgets.QLabel()
        self.LB_event_fit_size = QtWidgets.QLabel()
        self.LB_events_joc = QtWidgets.QLabel("JOC detection")
        self.LB_events_joc.setFont(self.fontbold)
        self.LB_events_joc_skip_start = QtWidgets.QLabel()
        self.LB_events_joc_fit_length = QtWidgets.QLabel()
        self.LB_events_joc_refit_times = QtWidgets.QLabel()
        self.LB_events_joc_refit_option = QtWidgets.QLabel()
        self.LB_events_joc_noise_multiplicator = QtWidgets.QLabel()
        self.LB_adaptive_threshold_option = QtWidgets.QLabel()
        self.LB_adaptive_smoothing_window = QtWidgets.QLabel()
        self.LB_adaptive_smoothing_order = QtWidgets.QLabel()

        VL = QtWidgets.QVBoxLayout()
        VL.addWidget(self.LB_model_detection_selected)
        VL.addWidget(self.LB_kernel_size)
        VL.addWidget(self.LB_msf_window_size)
        VL.addWidget(self.LB_threshold)
        VL.addWidget(self.LB_event_fit_size)
        VL.addWidget(self.LB_events_joc_skip_start)
        VL.addWidget(self.LB_events_joc_fit_length)
        VL.addWidget(self.LB_events_joc_refit_times)
        VL.addWidget(self.LB_events_joc_refit_option)
        VL.addWidget(self.LB_events_joc_noise_multiplicator)
        VL.addWidget(self.LB_adaptive_threshold_option)
        VL.addWidget(self.LB_adaptive_smoothing_window)
        VL.addWidget(self.LB_adaptive_smoothing_order)

        VL.addStretch(1)

        box_events.setLayout(VL)

        # Final Layout
        HL.addWidget(box_curve_modif)
        HL.addWidget(box_elasticity)
        HL.addWidget(box_work_and_rupture_force)
        HL.addWidget(box_events)

        self.setLayout(HL)

        self.update_widget()

    def update_widget(self):
        """Update the GUI."""
        dt = shared.exp.current_data

        # General informations
        # - Smoothing
        if dt.used_sg_smoothing_enabled:
            value_on = "Smoothing : " + str(dt.used_sg_smoothing_enabled)
            value = "Uniform : " + str(dt.used_sg_smoothing_uniform)
            self.LB_smoothing_uniform.setText(value)
            value = "Width : " + str(dt.used_sg_smoothing_width)
            self.LB_smoothing_width.setText(value)
            value = "Coefficient : " + str(dt.used_sg_smoothing_order)
            self.LB_smoothing_coef.setText(value)
        else:
            value_on = "Smoothing : " + str(dt.used_sg_smoothing_enabled)
            self.LB_smoothing_uniform.setText("")
            self.LB_smoothing_width.setText("")
            self.LB_smoothing_coef.setText("")

        self.LB_smoothing_on.setText(value_on)

        # - Other
        self.LB_defl_sens.setText(
            "Deflection sensitivity (nm/V) : " +
            str(dt.used_deflection_sensitivity))
        self.LB_spring_constant.setText(
            "Spring constant (N/m) : " + str(dt.used_spring_constant * 1e9))

        # Elasticity
        if dt.stiffness_calculated:
            val = dt.used_stiffness_model_selected
            value = "Model : " + str(widgets_list.widget_compute.models[val])
            self.LB_model_selected.setText(value)
            value = "Skip (start) (nm) : " + \
                str(dt.used_fitparam_poc_skip_start)
            self.LB_poc_skip_start.setText(value)
            value = "Length (nm) : " + str(dt.used_fitparam_poc_fit_length)
            self.LB_poc_fit_length.setText(value)
            value = "Fit noise : " + \
                str(dt.used_fitparam_poc_noise_multiplicator)
            self.LB_poc_noise_multiplicator.setText(value)
            value = "Refit (nm) : " + str(dt.used_fitparam_poc_refit_option)
            self.LB_poc_refit_option.setText(value)
            value = "Refit (times) : " + str(dt.used_fitparam_poc_refit_times)
            self.LB_poc_refit_times.setText(value)
            value = "Start (nm) : " + str(dt.used_indentation_start)
            self.LB_indentation_start.setText(value)
            value = "Step (nm) : " + str(dt.used_indentation_step)
            self.LB_indentation_step.setText(value)
            value = "Stop (nm) : " + str(dt.used_indentation_stop)
            self.LB_indentation_stop.setText(value)

        # Events detection
        if dt.events_calculated:
            value = "Method : " + str(dt.used_events_algorithm)
            self.LB_model_detection_selected.setText(value)
            if dt.used_events_algorithm == "kernel":
                value = "Threshold : " + \
                    str(dt.used_events_kernel_detection_threshold)
                self.LB_threshold.setText(value)
                value = "Kernel Size : " + str(dt.used_kernel_size)
                self.LB_kernel_size.setText(value)
            elif dt.used_events_algorithm == "msf":
                value = "Threshold : " + \
                    str(dt.used_events_msf_detection_threshold)
                self.LB_threshold.setText(value)
                value = "Window size : " + str(dt.used_events_msf_window_size)
                self.LB_msf_window_size.setText(value)
            value = "Event fit size : " + str(dt.used_events_fit_size)
            self.LB_event_fit_size.setText(value)
            value = "Fit start (nm) : " + \
                str(dt.used_fitparam_events_joc_skip_start)
            self.LB_events_joc_skip_start.setText(value)
            value = "Fit length (nm) : " + \
                str(dt.used_fitparam_events_joc_fit_length)
            self.LB_events_joc_fit_length.setText(value)
            value = "Refit (times) : " + \
                str(dt.used_fitparam_events_joc_refit_times)
            self.LB_events_joc_refit_times.setText(value)
            value = "Refit (nm) : " + \
                str(dt.used_fitparam_events_joc_refit_option)
            self.LB_events_joc_refit_option.setText(value)
            value = "Noise : " + \
                str(dt.used_fitparam_events_joc_noise_multiplicator)
            self.LB_events_joc_noise_multiplicator.setText(value)
            value = "Adaptive threshold : " + \
                str(dt.used_adaptive_threshold_option)
            self.LB_adaptive_threshold_option.setText(value)
            if dt.used_adaptive_threshold_option:
                value = "Smoothing window size : " + \
                    str(dt.used_adaptive_smoothing_window)
                self.LB_adaptive_smoothing_window.setText(value)
                value = "Smoothing order : " + \
                    str(dt.used_adaptive_smoothing_order)
                self.LB_adaptive_smoothing_order.setText(value)
            else:
                self.LB_adaptive_smoothing_window.hide()
                self.LB_adaptive_smoothing_order.hide()

        # Work
        if dt.work_and_rupture_force1_calculated:
            value = "Fit start (nm) : " + str(dt.used_fitparam_joc_skip_start)
            self.LB_joc_skip_start.setText(value)
            value = "Fit lenght (nm) : " + str(dt.used_fitparam_joc_fit_length)
            self.LB_joc_fit_length.setText(value)
            value = "Refit (nm) : " + str(dt.used_fitparam_joc_refit_option)
            self.LB_joc_refit_times.setText(value)
            value = "Refit (times)  : " + str(dt.used_fitparam_joc_refit_times)
            self.LB_joc_noise_multiplicator.setText(value)
            value = "Noise : " + str(dt.used_fitparam_joc_noise_multiplicator)
            self.LB_joc_refit_option.setText(value)
