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

"""Widget with options for the experiment results and statistic analysis."""

import logging

from ... import shared
from ... import widgets_list
from PyQt5 import QtWidgets
from ...tools import misc_tools
from ...tools import stat_tools
from ...tools.gui_tools import PYAFInput
from ...tools.gui_tools import PYAFComboBox
from ...tools.PYAFWidget import PYAFWidget
from ...tools.gui_tools import PYAFCheckBox, PYAFButton
from ...widgets.stats_advanced_options import AdvancedOptions


class ResultsExperimentOptionsWidget(PYAFWidget):
    """Options for the statistical analysis and experiment results plots."""

    def __init__(self, parent, name):
        super().__init__(parent, name)

        self.name = name

        # Load the logger
        self.logger = logging.getLogger()

        self.GL = QtWidgets.QGridLayout()

        self.BOX_general = QtWidgets.QGroupBox("Options")

        self.VL_general = QtWidgets.QVBoxLayout()

        # List to chose between normed or non normed histograms
        self.list_norm = PYAFComboBox(self, "norm")
        self.list_norm.addItem("PDF")
        self.list_norm.addItem("Counts")

        # List to chose plot view
        self.list_plot_view = PYAFComboBox(self, "plot_views")
        self.list_plot_view.addItem("Histograms")
        self.list_plot_view.addItem("Statistic plots")

        # Display or hide legend
        self.CB_legend = PYAFCheckBox(self, "legend", "Display legend")

        # Box general layout
        self.VL_general.addWidget(self.list_norm)
        self.VL_general.addWidget(self.list_plot_view)
        self.VL_general.addWidget(self.CB_legend)
        self.VL_general.addStretch(1)
        self.BOX_general.setLayout(self.VL_general)

        # Bins and bandwidth
        self.BOX_bins = QtWidgets.QGroupBox("Bins and bandwidth")
        self.VL_bins_and_bandwidth = QtWidgets.QVBoxLayout()
        self.list_bins = PYAFComboBox(self, "bins")
        self.list_bins.addItem("Auto")
        self.list_bins.addItem("Manual")
        self.input_bins = PYAFInput(self, "bins", "Nbr of Bins", width=70)
        self.input_bins.input.setValidator(misc_tools.validator("UI"))
        self.list_bandwidth = PYAFComboBox(self, "bandwidth")
        self.list_bandwidth.addItem("Auto")
        self.list_bandwidth.addItem("Manual")
        label = "PDF BW"
        self.input_bandwidth = PYAFInput(self, "bandwidth", label, width=70)
        self.input_bandwidth.input.setValidator(misc_tools.validator("UF"))
        self.VL_bins_and_bandwidth.addWidget(self.list_bins)
        self.VL_bins_and_bandwidth.addWidget(self.input_bins)
        self.VL_bins_and_bandwidth.addWidget(self.list_bandwidth)
        self.VL_bins_and_bandwidth.addWidget(self.input_bandwidth)
        self.VL_bins_and_bandwidth.addStretch(1)
        self.BOX_bins.setLayout(self.VL_bins_and_bandwidth)

        # Statistics
        self.BOX_stat_tests = QtWidgets.QGroupBox("Statistics")
        self.VL_stat_tests = QtWidgets.QVBoxLayout()
        self.list_stat_tests = PYAFComboBox(self, "stat_tests")

        for test in shared.stat_tests:
            self.list_stat_tests.addItem(test)

        self.list_plot_type = PYAFComboBox(self, "plot_type")

        for plot_type in shared.stat_plots:
            self.list_plot_type.addItem(plot_type)

        self.advanced_options = PYAFButton(self, "advanced_options", text="Advanced options")
        self.perform_test = PYAFButton(self, "perform_test", text="Compute")
        self.VL_stat_tests.addWidget(self.list_stat_tests)
        self.VL_stat_tests.addWidget(self.list_plot_type)
        self.VL_stat_tests.addWidget(self.advanced_options)
        self.VL_stat_tests.addWidget(self.perform_test)
        self.VL_stat_tests.addStretch(1)
        self.BOX_stat_tests.setLayout(self.VL_stat_tests)

        # X limits
        self.BOX_limits_x = QtWidgets.QGroupBox("X scale")
        self.VL_limits_x = QtWidgets.QVBoxLayout()
        self.list_limits_x = PYAFComboBox(self, "limits_x")
        self.list_limits_x.addItem("Auto")
        self.list_limits_x.addItem("0 to Max")
        self.list_limits_x.addItem("Manual")
        self.input_limits_min_x = PYAFInput(self, "min_x", "X min", width=70)
        self.input_limits_max_x = PYAFInput(self, "max_x", "X max", width=70)
        self.input_limits_min_x.input.setValidator(misc_tools.validator("F"))
        self.input_limits_max_x.input.setValidator(misc_tools.validator("F"))
        self.VL_limits_x.addWidget(self.list_limits_x)
        self.VL_limits_x.addWidget(self.input_limits_min_x)
        self.VL_limits_x.addWidget(self.input_limits_max_x)
        self.VL_limits_x.addStretch(1)
        self.BOX_limits_x.setLayout(self.VL_limits_x)

        # Y limits
        self.BOX_limits_y = QtWidgets.QGroupBox("Y scale")
        self.VL_limits_y = QtWidgets.QVBoxLayout()
        self.list_limits_y = PYAFComboBox(self, "limits_y")
        self.list_limits_y.addItem("Auto")
        self.list_limits_y.addItem("Manual")
        self.input_limits_min_y = PYAFInput(self, "min_y", "Y min", width=70)
        self.input_limits_max_y = PYAFInput(self, "max_y", "Y max", width=70)
        self.input_limits_max_y.input.setValidator(misc_tools.validator("UF"))
        self.input_limits_min_y.input.setValidator(misc_tools.validator("UF"))
        self.VL_limits_y.addWidget(self.list_limits_y)
        self.VL_limits_y.addWidget(self.input_limits_min_y)
        self.VL_limits_y.addWidget(self.input_limits_max_y)
        self.VL_limits_y.addStretch(1)
        self.BOX_limits_y.setLayout(self.VL_limits_y)

        self.BOX_general.setFixedSize(200, 210)
        self.BOX_bins.setFixedSize(200, 210)
        self.BOX_stat_tests.setFixedSize(200, 210)
        self.BOX_limits_x.setFixedSize(200, 160)
        self.BOX_limits_y.setFixedSize(200, 160)

        self.box_stack = QtWidgets.QStackedWidget()
        self.box_stack.addWidget(self.BOX_bins)
        self.box_stack.addWidget(self.BOX_stat_tests)

        # Finish layout
        self.GL.addWidget(self.BOX_general, 0, 0)
        self.GL.addWidget(self.BOX_limits_x, 1, 0)
        self.GL.addWidget(self.box_stack, 0, 1)
        self.GL.addWidget(self.BOX_limits_y, 1, 1)


        VL_options = QtWidgets.QVBoxLayout()
        HL_options = QtWidgets.QHBoxLayout()
        VL_options.addLayout(self.GL)
        VL_options.addStretch(1)
        HL_options.addLayout(VL_options)
        HL_options.addStretch(1)

        self.setLayout(HL_options)

        self.update_widget()

    def update_widget(self):
        """Update the widget."""
        self.update_GUI("all")

    def update_GUI(self, what):
        """Update the GUI."""
        if what == "norm" or what == "all":
            if shared.exp.norm_hist_experiment:
                self.list_norm.setCurrentIndex(0)
            else:
                self.list_norm.setCurrentIndex(1)

        if what == "legend" or what == "all":
            value = shared.exp.display_legend_hist_experiment
            self.CB_legend.setChecked(value)

        if what == "bins" or what == "all":
            mode = shared.exp.hist_experiment_bins_mode
            self.input_bins.changeValue(stat_tools.get_bins("experiment"))

            if mode == "auto":
                self.list_bins.setCurrentIndex(0)
                self.input_bins.setEnabled(False)
            elif mode == "manual":
                self.list_bins.setCurrentIndex(1)
                self.input_bins.setEnabled(True)

        if what == "bandwidth" or what == "all":
            mode = shared.exp.hist_experiment_bw_mode
            bw = stat_tools.get_bandwidth("experiment")
            self.input_bandwidth.changeValue(bw)

            if mode == "auto":
                self.list_bandwidth.setCurrentIndex(0)
                self.input_bandwidth.setEnabled(False)
            elif mode == "manual":
                self.list_bandwidth.setCurrentIndex(1)
                self.input_bandwidth.setEnabled(True)

        if what == "limits_x" or what == "all":
            mode = shared.exp.hist_experiment_x_mode
            value_max = shared.exp.hist_experiment_max_x
            value_min = shared.exp.hist_experiment_min_x

            if mode == "auto":
                self.list_limits_x.setCurrentIndex(0)
                self.input_limits_min_x.setEnabled(False)
                self.input_limits_max_x.setEnabled(False)
            elif mode == "max":
                self.list_limits_x.setCurrentIndex(1)
                self.input_limits_min_x.setEnabled(False)
                self.input_limits_max_x.setEnabled(False)
            elif mode == "manual":
                self.list_limits_x.setCurrentIndex(2)
                self.input_limits_min_x.setEnabled(True)
                self.input_limits_max_x.setEnabled(True)

            self.input_limits_min_x.changeValue(value_min)
            self.input_limits_max_x.changeValue(value_max)

        if what == "limits_y" or what == "all":
            mode = shared.exp.hist_experiment_y_mode
            value_max = shared.exp.hist_experiment_max_y
            value_min = shared.exp.hist_experiment_min_y

            if mode == "auto":
                self.list_limits_y.setCurrentIndex(0)
                self.input_limits_max_y.setEnabled(False)
                self.input_limits_min_y.setEnabled(False)
            elif mode == "manual":
                self.list_limits_y.setCurrentIndex(1)
                self.input_limits_max_y.setEnabled(True)
                self.input_limits_min_y.setEnabled(True)

            self.input_limits_min_y.changeValue(value_min)
            self.input_limits_max_y.changeValue(value_max)

    def button_clicked(self, button):
        """Method called upon a click on a button."""
        if button == "advanced_options":
            self.open_advanced_options_widget()
        if button == "perform_test":
            stat_tools.fetch_conditions_data()
            self.perform_statistical_analysis()

    def checkbox_clicked(self, checkbox):
        """Called when a checkbox is clicked."""
        if checkbox == "legend":
            value = self.CB_legend.isChecked()
            shared.exp.display_legend_hist_experiment = value
            widgets_list.widget_results_experiment.update_MPL()

    def list_updated(self, name):
        """Called when a list is updated."""
        if name == "norm":
            if self.list_norm.currentIndex() == 0:
                value = True
            else:
                value = False

            shared.exp.norm_hist_experiment = value
            widgets_list.widget_results_experiment.update_GUI("table_labels")
            widgets_list.widget_results_experiment.update_MPL()

        elif name == "bins":
            if self.list_bins.currentIndex() == 0:
                value = "auto"
            elif self.list_bins.currentIndex() == 1:
                value = "manual"

            shared.exp.hist_experiment_bins_mode = value
            widgets_list.widget_results_experiment.update_MPL()

        elif name == "bandwidth":
            if self.list_bandwidth.currentIndex() == 0:
                value = "auto"
            elif self.list_bandwidth.currentIndex() == 1:
                value = "manual"

            shared.exp.hist_experiment_bw_mode = value
            widgets_list.widget_results_experiment.update_widget()

        elif name == "limits_x":
            if self.list_limits_x.currentIndex() == 0:
                mode = "auto"
            elif self.list_limits_x.currentIndex() == 1:
                mode = "max"
            elif self.list_limits_x.currentIndex() == 2:
                mode = "manual"

            shared.exp.hist_experiment_x_mode = mode
            widgets_list.widget_results_experiment.update_MPL()

            self.update_GUI("limits_x")

        elif name == "limits_y":
            if self.list_limits_y.currentIndex() == 0:
                mode = "auto"
            elif self.list_limits_y.currentIndex() == 1:
                mode = "manual"

            shared.exp.hist_experiment_y_mode = mode
            widgets_list.widget_results_experiment.update_MPL()

            self.update_GUI("limits_y")

        elif name == "plot_views":
            widgets_list.widget_results_experiment.canvas_stack.setCurrentIndex(self.list_plot_view.currentIndex())
            self.box_stack.setCurrentIndex(self.list_plot_view.currentIndex())

        elif name == "plot_type":
            shared.exp.stat_plot_type = self.list_plot_type.currentText()
            widgets_list.widget_results_experiment.update_MPL()

        elif name == "stat_tests":
            shared.exp.selected_statistical_test = self.list_stat_tests.currentText()

    def input_updated(self, input_field):
        """Called when an input is updated."""
        noupdate = False

        if input_field == "bins":
            if self.input_bins.input.text() != "":
                value = self.input_bins.get_float_value()
                shared.exp.hist_experiment_bins = value

        elif input_field == "bandwidth":
            if self.input_bandwidth.input.text() != "":
                value = self.input_bandwidth.get_float_value()
                shared.exp.hist_experiment_bw = value
                widgets_list.widget_results_experiment.update_widget()

                # Update is already done
                noupdate = True

        elif input_field == "min_x":
            value = self.input_limits_min_x.get_float_value()
            shared.exp.hist_experiment_min_x = value

        elif input_field == "max_x":
            value = self.input_limits_max_x.get_float_value()
            shared.exp.hist_experiment_max_x = value

        elif input_field == "min_y":
            value = self.input_limits_min_y.get_float_value()
            shared.exp.hist_experiment_min_y = value

        elif input_field == "max_y":
            value = self.input_limits_max_y.get_float_value()
            shared.exp.hist_experiment_max_y = value

        if not noupdate:
            widgets_list.widget_results_experiment.update_MPL()

    def open_advanced_options_widget(self):
        """Open the advanced options for the statistical analysis Window."""
        if widgets_list.widget_advanced_options is None:
            # Create a new widget
            widgets_list.widget_advanced_options = AdvancedOptions(self)
            widgets_list.widget_advanced_options.setWindowTitle("Advanced Options")
            widgets_list.widget_advanced_options.show()
        else:
            # Bring to front
            widgets_list.widget_advanced_options.activateWindow()
            widgets_list.widget_advanced_options.raise_()

    def perform_statistical_analysis(self):
        """"Function to perform the selected statistical test when the compute button is clicked"""
        stat_tools.do_statistical_test(self.list_stat_tests.currentText())
        widgets_list.widget_results_experiment.update_MPL()
