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
This widget displays input fields for the fitting parameters.

Depending of the type of fitting parameters to display, the values are filled
in the inputs. Changing a value will update the plot in the compute tab.

"""

from ... import shared
import logging
from ... import widgets_list
from PyQt5 import QtWidgets
from ...tools.PYAFWidget import PYAFWidget
from ...tools.gui_tools import PYAFInput
from ...tools.gui_tools import PYAFCheckBox
from ...tools import misc_tools


class FitParamsWidget(PYAFWidget):
    """PYAFWidget containing the input fields for the fits.

    Is called by the widget_compute widget for each type of computation.
    """

    def __init__(self, parent, name):
        super().__init__(parent, name)

        # Load the logger
        self.logger = logging.getLogger()
        self.logger.debug("Creating widget : %s", self.name)

        VL = QtWidgets.QVBoxLayout()

        if self.name == "fit_param_pocs":
            title = "Point of contact"
        elif self.name == "fit_param_jocs":
            title = "Jump of contact"
        elif self.name == "fit_param_jocs_events":
            title = "Jump of contact"

        BOX = QtWidgets.QGroupBox(title)
        self.GL = QtWidgets.QGridLayout()
        label = "Fit skip (start) [nm]"

        self.IN_skip_start = PYAFInput(self, "skip_start", label)
        self.IN_fit_length = PYAFInput(
            self, "fit_length", "Fit length [nm]")
        self.IN_noise_multi = PYAFInput(self, "noise_multi", "Fit noise")
        self.IN_refit = PYAFInput(self, "refit", "Fit refit [nm]")
        self.IN_refit_times = PYAFInput(
            self, "refit_times", "Refit [times]")
        self.IN_skip_start.input.setValidator(misc_tools.validator("UI"))
        self.IN_fit_length.input.setValidator(misc_tools.validator("UI"))
        self.IN_noise_multi.input.setValidator(misc_tools.validator("UF"))
        self.IN_refit.input.setValidator(misc_tools.validator("UI"))
        self.IN_refit_times.input.setValidator(misc_tools.validator("UI"))

        # Is added to GL in update_widget
        self.CB_display_fit_events_joc2_preview = PYAFCheckBox(
            self, "display_fit_events_joc2_preview", "Display fit")

        self.GL.addWidget(self.IN_skip_start, 0, 0)
        self.GL.addWidget(self.IN_fit_length, 1, 0)
        self.GL.addWidget(self.IN_noise_multi, 2, 0)
        self.GL.addWidget(self.IN_refit, 0, 1)
        self.GL.addWidget(self.IN_refit_times, 1, 1)

        VL_empty = QtWidgets.QVBoxLayout()
        VL_empty.addStretch(1)
        self.GL.addLayout(VL_empty, 3, 0)

        BOX.setLayout(self.GL)
        VL.addWidget(BOX)

        self.setLayout(VL)

        self.update_widget()

    def update_widget(self):
        """Generic method to update the GUI of the widget.

        The input fields are filled with their respective values and if needed
        the checkbox for point of contact recomputation is checked/unchecked
        and enabled/disabled.
        """
        data = shared.exp.current_data

        if self.name == "fit_param_pocs":
            self.IN_skip_start.input.setText(str(data.fitparam_poc_skip_start))
            self.IN_fit_length.input.setText(str(data.fitparam_poc_fit_length))
            self.IN_noise_multi.input.setText(
                str(data.fitparam_poc_noise_multiplicator))
            self.IN_refit.input.setText(str(data.fitparam_poc_refit_option))
            self.IN_refit_times.input.setText(
                str(data.fitparam_poc_refit_times))
            self.CB_display_fit_events_joc2_preview.setParent(None)

        elif self.name == "fit_param_jocs":
            self.IN_skip_start.input.setText(str(data.fitparam_joc_skip_start))
            self.IN_fit_length.input.setText(str(data.fitparam_joc_fit_length))
            self.IN_noise_multi.input.setText(
                str(data.fitparam_joc_noise_multiplicator))
            self.IN_refit.input.setText(str(data.fitparam_joc_refit_option))
            self.IN_refit_times.input.setText(
                str(data.fitparam_joc_refit_times))
            self.CB_display_fit_events_joc2_preview.setParent(None)

        elif self.name == "fit_param_jocs_events":
            self.IN_skip_start.input.setText(
                str(data.fitparam_events_joc_skip_start))
            self.IN_fit_length.input.setText(
                str(data.fitparam_events_joc_fit_length))
            self.IN_noise_multi.input.setText(
                str(data.fitparam_events_joc_noise_multiplicator))
            self.IN_refit.input.setText(
                str(data.fitparam_events_joc_refit_option))
            self.IN_refit_times.input.setText(
                str(data.fitparam_events_joc_refit_times))
            self.GL.addWidget(self.CB_display_fit_events_joc2_preview, 2, 1)
            value = data.display_fit_events_joc2_preview
            self.CB_display_fit_events_joc2_preview.setChecked(value)

    def checkbox_clicked(self, name):
        """Called when a checkbox is clicked."""
        if name == "display_fit_events_joc2_preview":
            data = shared.exp.current_data
            val = self.CB_display_fit_events_joc2_preview.isChecked()
            data.display_fit_events_joc2_preview = val

            widgets_list.widget_compute.update_MPL("MPL_canvas")

    def input_updated(self, input_field):
        """Method called upon change in an input field of the GUI."""
        text = "Input updated (%s). Input : %s"
        self.logger.debug(text, self.name, input_field)

        data = shared.exp.current_data

        if input_field == "skip_start":
            value = self.IN_skip_start.get_int_value()
            if self.name == "fit_param_pocs":
                data.fitparam_poc_skip_start = value
            elif self.name == "fit_param_jocs":
                data.fitparam_joc_skip_start = value
            elif self.name == "fit_param_jocs_events":
                data.fitparam_events_joc_skip_start = value

        elif input_field == "fit_length":
            value = self.IN_fit_length.get_int_value()
            if self.name == "fit_param_pocs":
                data.fitparam_poc_fit_length = value
            elif self.name == "fit_param_jocs":
                data.fitparam_joc_fit_length = value
            elif self.name == "fit_param_jocs_events":
                data.fitparam_events_joc_fit_length = value

        elif input_field == "noise_multi":
            value = self.IN_noise_multi.get_float_value()
            if self.name == "fit_param_pocs":
                data.fitparam_poc_noise_multiplicator = value
            elif self.name == "fit_param_jocs":
                data.fitparam_joc_noise_multiplicator = value
            elif self.name == "fit_param_jocs_events":
                data.fitparam_events_joc_noise_multiplicator = value

        elif input_field == "refit":
            value = self.IN_refit.get_int_value()
            if self.name == "fit_param_pocs":
                data.fitparam_poc_refit_option = value
            elif self.name == "fit_param_jocs":
                data.fitparam_joc_refit_option = value
            elif self.name == "fit_param_jocs_events":
                data.fitparam_events_joc_refit_option = value

        elif input_field == "refit_times":
            value = self.IN_refit_times.get_int_value()
            if self.name == "fit_param_pocs":
                data.fitparam_poc_refit_times = value
            elif self.name == "fit_param_jocs":
                data.fitparam_joc_refit_times = value
            elif self.name == "fit_param_jocs_events":
                data.fitparam_events_joc_refit_times = value

        widgets_list.widget_compute.update_MPL("MPL_canvas")
