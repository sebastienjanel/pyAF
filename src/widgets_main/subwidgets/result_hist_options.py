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

"""Widget with options for the results."""

import logging
from ... import shared
from ... import widgets_list
from PyQt5 import QtWidgets
from ...tools import misc_tools
from ...tools import stat_tools
from ...tools.gui_tools import PYAFInput
from ...tools.gui_tools import PYAFComboBox
from ...tools.PYAFWidget import PYAFWidget
from ...tools.gui_tools import PYAFCheckBox


class ResultsHistOptionsWidget(PYAFWidget):
    """Options for the results plots.

    The widget is used for the single results and the grouped results.
    (Plots and grouped plots)
    """

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

        # Display or hide legend
        self.CB_legend = PYAFCheckBox(self, "legend", "Display legend")

        # Display or hide fit
        self.CB_fit = PYAFCheckBox(self, "fit", "Display PDF fit")

        if self.name == "widget_results_hist_options_single":
            # Remove or display zeros
            self.CB_zeros = PYAFCheckBox(self, "zeros", "Remove zeros")

            # Use log scale
            self.CB_log = PYAFCheckBox(self, "log", "Use log scale")

        # Scatter (for loading rates)
        self.CB_scatter = PYAFCheckBox(self, "scatter", "Scatter")

        # Box general layout
        self.VL_general.addWidget(self.list_norm)
        self.VL_general.addWidget(self.CB_legend)
        self.VL_general.addWidget(self.CB_fit)
        if self.name == "widget_results_hist_options_single":
            self.VL_general.addWidget(self.CB_zeros)
            self.VL_general.addWidget(self.CB_log)
        self.VL_general.addWidget(self.CB_scatter)
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

        # Plots
        # self.BOX_plots = QtWidgets.QGroupBox("Plots")
        # self.VL_plots = QtWidgets.QVBoxLayout()
        # self.plot_list = PYAFComboBox(self, "plot_selector")
        # self.plot_list.addItem("Hist")
        # self.plot_list.addItem("Boxplot")
        # self.VL_plots.addWidget(self.plot_list)
        # self.BOX_plots.setLayout(self.VL_plots)

        self.BOX_general.setFixedSize(200, 210)
        self.BOX_bins.setFixedSize(200, 210)
        self.BOX_limits_x.setFixedSize(200, 160)
        self.BOX_limits_y.setFixedSize(200, 160)
        # self.BOX_plots.setFixedSize(100, 100)

        # Finish layout
        self.GL.addWidget(self.BOX_general, 0, 0)
        self.GL.addWidget(self.BOX_bins, 0, 1)
        self.GL.addWidget(self.BOX_limits_x, 1, 0)
        self.GL.addWidget(self.BOX_limits_y, 1, 1)
        # self.GL.addWidget(self.BOX_plots, 0, 2)

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
            if self.name == "widget_results_hist_options_single":
                if shared.exp.norm_hist_single:
                    self.list_norm.setCurrentIndex(0)
                    # self.plot_list.setCurrentIndex(0)
                else:
                    self.list_norm.setCurrentIndex(1)
                    # self.plot_list.setCurrentIndex(1)
            elif self.name == "widget_results_hist_options_groups":
                if shared.exp.norm_hist_groups:
                    self.list_norm.setCurrentIndex(0)
                else:
                    self.list_norm.setCurrentIndex(1)

        if what == "legend" or what == "all":
            if self.name == "widget_results_hist_options_single":
                value = shared.exp.display_legend_hist_single
            elif self.name == "widget_results_hist_options_groups":
                value = shared.exp.display_legend_hist_groups

            self.CB_legend.setChecked(value)

        if what == "fit" or what == "all":
            if self.name == "widget_results_hist_options_single":
                value = shared.exp.hist_single_display_pdf_fit
                if shared.single_pdfs_x is None:
                    self.CB_fit.setEnabled(False)
                else:
                    self.CB_fit.setEnabled(True)
            elif self.name == "widget_results_hist_options_groups":
                value = shared.exp.hist_groups_display_pdf_fit
                if shared.groups_pdfs_x is None:
                    self.CB_fit.setEnabled(False)
                else:
                    self.CB_fit.setEnabled(True)

            self.CB_fit.setChecked(value)

        if what == "zeros" or what == "all":
            if self.name == "widget_results_hist_options_single":
                self.CB_zeros.setChecked(shared.exp.hist_remove_zeros)

        if what == "log" or what == "all":
            if self.name == "widget_results_hist_options_single":
                self.CB_log.setChecked(shared.exp.hist_log)

        if what == "scatter" or what == "all":
            if self.name == "widget_results_hist_options_single":
                value = shared.exp.hist_lr_single_display_mode
            elif self.name == "widget_results_hist_options_groups":
                value = shared.exp.hist_lr_groups_display_mode

            if value == "scatter":
                value = True
            else:
                value = False

            self.CB_scatter.setChecked(value)

            # Enable this checkbox only for loading rates
            if shared.exp.results_type != "loading_rates":
                self.CB_scatter.setEnabled(False)
            else:
                self.CB_scatter.setEnabled(True)

        if what == "bins" or what == "all":
            if self.name == "widget_results_hist_options_single":
                mode = shared.exp.hist_single_bins_mode
                self.input_bins.changeValue(stat_tools.get_bins("single"))
            elif self.name == "widget_results_hist_options_groups":
                mode = shared.exp.hist_groups_bins_mode
                self.input_bins.changeValue(stat_tools.get_bins("groups"))

            if mode == "auto":
                self.list_bins.setCurrentIndex(0)
                self.input_bins.setEnabled(False)
            elif mode == "manual":
                self.list_bins.setCurrentIndex(1)
                self.input_bins.setEnabled(True)

        if what == "bandwidth" or what == "all":
            if self.name == "widget_results_hist_options_single":
                mode = shared.exp.hist_single_bw_mode
                bw = stat_tools.get_bandwidth("single")
                self.input_bandwidth.changeValue(bw)
            elif self.name == "widget_results_hist_options_groups":
                mode = shared.exp.hist_groups_bw_mode
                bw = stat_tools.get_bandwidth("groups")
                self.input_bandwidth.changeValue(bw)

            if mode == "auto":
                self.list_bandwidth.setCurrentIndex(0)
                self.input_bandwidth.setEnabled(False)
            elif mode == "manual":
                self.list_bandwidth.setCurrentIndex(1)
                self.input_bandwidth.setEnabled(True)

        if what == "limits_x" or what == "all":
            if self.name == "widget_results_hist_options_single":
                mode = shared.exp.hist_single_x_mode
                value_max = shared.exp.hist_single_max_x
                value_min = shared.exp.hist_single_min_x
            elif self.name == "widget_results_hist_options_groups":
                mode = shared.exp.hist_groups_x_mode
                value_max = shared.exp.hist_groups_max_x
                value_min = shared.exp.hist_groups_min_x

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
            if self.name == "widget_results_hist_options_single":
                mode = shared.exp.hist_single_y_mode
                value_max = shared.exp.hist_single_max_y
                value_min = shared.exp.hist_single_min_y
            elif self.name == "widget_results_hist_options_groups":
                mode = shared.exp.hist_groups_y_mode
                value_max = shared.exp.hist_groups_max_y
                value_min = shared.exp.hist_groups_min_y

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

        # Keep the plots combo box index value before the refresh
        if what == "plots" or what == "all":
            if self.name == "widget_results_hist_options_single":
                plot = shared.exp.hist_single_plot_selected
            elif self.name == "widget_results_hist_options_groups":
                plot = shared.exp.hist_groups_plot_selected

            # if plot == "hist":
            #    self.plot_list.setCurrentIndex(0)
            # elif plot == "box":
            #    self.plot_list.setCurrentIndex(1)


    def button_clicked(self, button):
        """Method called upon a click on a button."""
        if button == "button_type_hist":
            # Update value
            val = self.button_type_hist.checkedId()
            shared.exp.results_type = misc_tools.get_results_type_by_id(val)
            # Refresh histograms and tables completely
            widgets_list.widget_results_single.update_widget()
            widgets_list.widget_results_groups.update_widget()

    def checkbox_clicked(self, checkbox):
        """Called when a checkbox is clicked."""
        if checkbox == "legend":
            value = self.CB_legend.isChecked()

            if self.name == "widget_results_hist_options_single":
                shared.exp.display_legend_hist_single = value
                widgets_list.widget_results_single.update_MPL()
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.display_legend_hist_groups = value
                widgets_list.widget_results_groups.update_MPL()

        elif checkbox == "fit":
            value = self.CB_fit.isChecked()

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_display_pdf_fit = value
                widgets_list.widget_results_single.update_MPL()
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_display_pdf_fit = value
                widgets_list.widget_results_groups.update_MPL()

        elif checkbox == "zeros":
            shared.exp.hist_remove_zeros = self.CB_zeros.isChecked()

            # Refresh values
            widgets_list.widget_results_single.button_clicked("refresh_hist")
            widgets_list.widget_results_single.update_MPL()
            widgets_list.widget_results_groups.update_MPL()

        elif checkbox == "log":
            shared.exp.hist_log = self.CB_log.isChecked()

            # Refresh values
            widgets_list.widget_results_single.button_clicked("refresh_hist")
            widgets_list.widget_results_single.update_MPL()
            widgets_list.widget_results_groups.update_MPL()

        elif checkbox == "scatter":
            value = self.CB_scatter.isChecked()
            if value:
                value = "scatter"
            else:
                value = "points"

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_lr_single_display_mode = value
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_lr_groups_display_mode = value

            widgets_list.widget_results_single.update_MPL()
            widgets_list.widget_results_groups.update_MPL()

    def list_updated(self, name):
        """Called when a list is updated."""
        if name == "norm":
            if self.list_norm.currentIndex() == 0:
                value = True
            else:
                value = False

            if self.name == "widget_results_hist_options_single":
                shared.exp.norm_hist_single = value
                widgets_list.widget_results_single.update_GUI("table_labels")
                widgets_list.widget_results_single.update_MPL()
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.norm_hist_groups = value
                widgets_list.widget_results_groups.update_GUI("table_labels")
                widgets_list.widget_results_groups.update_MPL()

        elif name == "bins":
            if self.list_bins.currentIndex() == 0:
                value = "auto"
            elif self.list_bins.currentIndex() == 1:
                value = "manual"

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_bins_mode = value
                widgets_list.widget_results_single.update_MPL()
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_bins_mode = value
                widgets_list.widget_results_groups.update_MPL()

            self.update_GUI("bins")

        elif name == "bandwidth":
            if self.list_bandwidth.currentIndex() == 0:
                value = "auto"
            elif self.list_bandwidth.currentIndex() == 1:
                value = "manual"

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_bw_mode = value
                widgets_list.widget_results_single.update_widget()
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_bw_mode = value
                widgets_list.widget_results_groups.update_widget()

            self.update_GUI("bandwidth")

        elif name == "limits_x":
            if self.list_limits_x.currentIndex() == 0:
                mode = "auto"
            elif self.list_limits_x.currentIndex() == 1:
                mode = "max"
            elif self.list_limits_x.currentIndex() == 2:
                mode = "manual"

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_x_mode = mode
                widgets_list.widget_results_single.update_MPL()
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_x_mode = mode
                widgets_list.widget_results_groups.update_MPL()

            self.update_GUI("limits_x")

        elif name == "limits_y":
            if self.list_limits_y.currentIndex() == 0:
                mode = "auto"
            elif self.list_limits_y.currentIndex() == 1:
                mode = "manual"

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_y_mode = mode
                widgets_list.widget_results_single.update_MPL()
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_y_mode = mode
                widgets_list.widget_results_groups.update_MPL()

            self.update_GUI("limits_y")

        # Define the plot selected
        elif name == "plot_selector":
            if self.plot_list.currentIndex() == 0:
                plot_selected = "hist"
                self.BOX_limits_x.setEnabled(True)
            elif self.plot_list.currentIndex() == 1:
                plot_selected = "box"
                self.BOX_limits_x.setEnabled(False)
            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_plot_selected = plot_selected
                widgets_list.widget_results_single.update_MPL()
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_plot_selected = plot_selected
                widgets_list.widget_results_groups.update_MPL()


    def input_updated(self, input_field):
        """Called when an input is updated."""
        noupdate = False

        if input_field == "bins":
            if self.input_bins.input.text() != "":
                value = self.input_bins.get_float_value()

                if self.name == "widget_results_hist_options_single":
                    shared.exp.hist_single_bins = value
                elif self.name == "widget_results_hist_options_groups":
                    shared.exp.hist_groups_bins = value

        elif input_field == "bandwidth":
            if self.input_bandwidth.input.text() != "":
                value = self.input_bandwidth.get_float_value()

                if self.name == "widget_results_hist_options_single":
                    shared.exp.hist_single_bw = value
                    widgets_list.widget_results_single.update_widget()
                elif self.name == "widget_results_hist_options_groups":
                    shared.exp.hist_groups_bw = value
                    widgets_list.widget_results_groups.update_widget()

                # Update is already done
                noupdate = True

        elif input_field == "min_x":
            value = self.input_limits_min_x.get_float_value()

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_min_x = value
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_min_x = value

        elif input_field == "max_x":
            value = self.input_limits_max_x.get_float_value()

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_max_x = value
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_max_x = value

        elif input_field == "min_y":
            value = self.input_limits_min_y.get_float_value()

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_min_y = value
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_min_y = value

        elif input_field == "max_y":
            value = self.input_limits_max_y.get_float_value()

            if self.name == "widget_results_hist_options_single":
                shared.exp.hist_single_max_y = value
            elif self.name == "widget_results_hist_options_groups":
                shared.exp.hist_groups_max_y = value

        if not noupdate:
            if self.name == "widget_results_hist_options_single":
                widgets_list.widget_results_single.update_MPL()
            elif self.name == "widget_results_hist_options_groups":
                widgets_list.widget_results_groups.update_MPL()
