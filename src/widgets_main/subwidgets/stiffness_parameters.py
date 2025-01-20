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
from ... import consts
from PyQt5 import QtCore, QtWidgets
from ...tools.PYAFWidget import PYAFWidget
from ...tools.gui_tools import PYAFInput, PYAFComboBox
from ...tools.gui_tools import PYAFCheckBox
from ...tools import misc_tools


class StiffnessParamsWidget(PYAFWidget):
    """PYAFWidget containing the input fields for the Elasticity calculus/tomography and fits."""

    def __init__(self, parent, name):
        super().__init__(parent, name)

        self.name = name

        # Button group
        self.list_range_mode = PYAFComboBox(self, "list_range_mode")
        self.list_range_mode.addItem("Force")
        self.list_range_mode.addItem("Indentation")

        # Indentation start
        self.IN_start = PYAFInput(self, "input_start_in", "Start (nm)")
        self.IN_start.input.setValidator(misc_tools.validator("UI"))
        # Indentation end
        self.IN_stop = PYAFInput(self, "input_stop_in", "Stop (nm)")
        self.IN_stop.input.setValidator(misc_tools.validator("UI"))
        # Indentation step
        self.IN_step = PYAFInput(self, "input_step_in", "Segment size (nm)")
        self.IN_step.input.setValidator(misc_tools.validator("UI"))
        # Force start
        self.F_start = PYAFInput(self, "input_start_f", "Start (%)")
        self.F_start.input.setValidator(misc_tools.validator("UI"))
        # Force end
        self.F_stop = PYAFInput(self, "input_stop_f", "Stop (%)")
        self.F_stop.input.setValidator(misc_tools.validator("UI"))

        data = shared.exp.list[shared.exp.id_selected]

        self.IN_start.changeValue(data.indentation_start)
        self.IN_stop.changeValue(data.indentation_stop)
        self.IN_step.changeValue(data.indentation_step)

        # Stacked widget
        self.STK_data_selection = QtWidgets.QStackedWidget(self)

        # Force range selection
        self.force_selection = QtWidgets.QWidget(None)
        self.VL_force_selection = QtWidgets.QVBoxLayout()
        self.VL_force_selection.addWidget(self.F_start)
        self.VL_force_selection.addWidget(self.F_stop)
        self.VL_force_selection.addStretch(1)
        self.force_selection.setLayout(self.VL_force_selection)

        # Indentation range selection
        self.indentation_selection = QtWidgets.QWidget(None)
        self.VL_indentation_selection = QtWidgets.QVBoxLayout()
        self.VL_indentation_selection.addWidget(self.IN_start)
        self.VL_indentation_selection.addWidget(self.IN_stop)
        self.VL_indentation_selection.addStretch(1)
        self.indentation_selection.setLayout(self.VL_indentation_selection)

        self.STK_data_selection.addWidget(self.force_selection)
        self.STK_data_selection.addWidget(self.indentation_selection)

        # Checkboxes
        label = "Tomography"
        self.CB_tomography = PYAFCheckBox(self, "tomography", label)
        # self.CB_strict_stop = PYAFCheckBox(self, "strict_stop", "Strict stop")
        # self.CB_perform_fit = PYAFCheckBox(self, "perform_fit", "Perform fit")

        widgets_list.widget_progressbar.update()

        self.VL_params_stiff = QtWidgets.QVBoxLayout()
        # Display the start parameter only in advanced mode

        if self.name == "widget_calculus":
            self.VL_list_calc = QtWidgets.QVBoxLayout()
            if consts.ADVANCED:
                self.VL_list_calc.addWidget(self.IN_start)
                self.VL_list_calc.addWidget(self.IN_step)
                # By default tomography is set to False. We display the Kasas
                # model only in advanced mode.
                self.VL_list_calc.addWidget(self.CB_tomography)
            else:
                self.VL_list_calc.addWidget(self.IN_step)

            self.VL_list_calc.addStretch(1)
            self.VL_params_stiff.addLayout(self.VL_list_calc)

        elif self.name == "widget_fitting":

            self.HL_list_range_mode = QtWidgets.QHBoxLayout()
            self.label_range_mode = QtWidgets.QLabel("Range Mode : ")
            self.HL_list_range_mode.addWidget(self.label_range_mode)
            self.HL_list_range_mode.addWidget(self.list_range_mode)
            self.HL_list_range_mode.addStretch(1)

            self.VL_params_stiff.addLayout(self.HL_list_range_mode)
            self.VL_params_stiff.addWidget(self.STK_data_selection)

        self.setLayout(self.VL_params_stiff)

        self.update_widget()

    def update_widget(self):
        """Generic method to update the GUI of the widget.

        The input fields are filled with their respective values and if needed
        the checkbox for point of contact recomputation is checked/unchecked
        and enabled/disabled.
        """
        data = shared.exp.current_data

        if self.name == "widget_calculus":
            self.IN_start.changeValue(data.indentation_start)
            self.IN_step.changeValue(data.indentation_step)
            self.CB_tomography.setChecked(data.tomography)

        elif self.name == "widget_fitting":
            val = data.fit_range_type

            if val == 0:
                self.F_start.changeValue(data.force_start)
                self.F_stop.changeValue(data.force_stop)

            elif val == 1:
                self.IN_start.changeValue(data.indentation_start)
                self.IN_stop.changeValue(data.indentation_stop)

            self.list_range_mode.setCurrentIndex(val)
            self.STK_data_selection.setCurrentIndex(val)

    def list_updated(self, name):
        """Called when a list is updated."""
        data = shared.exp.current_data

        if name == "list_range_mode":
            # This list lets you chose the model for the stiffness computation.
            value = str(self.list_range_mode.currentText())

            if value == "Force":
                index = 0
                data.fit_range_type = index
                self.STK_data_selection.setCurrentIndex(index)

            elif value == "Indentation":
                index = 1
                data.fit_range_type = index
                self.STK_data_selection.setCurrentIndex(index)

            # Update Widget
            self.update_widget()

    def input_updated(self, field):
        """Method called upon change in an input field of the GUI."""
        data = shared.exp.current_data

        if field == "input_start_in":
            if self.IN_start.input.text() != "":
                data.indentation_start = self.IN_start.get_int_value()

        if field == "input_stop_in":
            if self.IN_stop.input.text() != "":
                data.indentation_stop = self.IN_stop.get_int_value()
                # self.update_GUI("strict_stop")

        if field == "input_step_in":
            if self.IN_step.input.text() != "":
                data.indentation_step = self.IN_step.get_int_value()

        if field == "input_start_f":
            if self.F_start.input.text() != "":
                data.force_start = self.F_start.get_int_value()

        if field == "input_stop_f":
            if self.F_stop.input.text() != "":
                data.force_stop = self.F_stop.get_int_value()

        widgets_list.widget_compute.update_MPL("MPL_canvas")

    def checkbox_clicked(self, checkbox):
        """Called when a checkbox is clicked."""

        if checkbox == "tomography":
            state = self.CB_tomography.isChecked()
            shared.exp.list[shared.exp.id_selected].tomography = state

        # widgets_list.widget_compute.update_MPL("MPL_canvas")
