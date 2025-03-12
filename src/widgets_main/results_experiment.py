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

"""Widget displaying histograms and results (Grouped)."""

from itertools import combinations
import logging
from .. import shared
from .. import widgets_list
from PyQt5 import QtCore, QtWidgets
from ..tools import gui_tools
from ..tools import stat_tools
from ..tools.gui_tools import PYAFComboBox
from ..widgets_main.subwidgets.result_chooser import ResultsChooserWidget
from ..widgets_main.subwidgets.result_experiment_options import \
    ResultsExperimentOptionsWidget
from ..widgets_small.rename_conditions import RenameConditionsWidget
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFTableWidget, PYAFButton
from ..experiment import ConditionGroup, SampleGroup
from .subwidgets.result_header import QCheckableHeader
from ..widgets.sample_grouping_dialog import SampleGroupingDialog


class ResultsExperimentWidget(PYAFWidget):
    """Widget displaying histograms and results (Grouped).

    This widget is in the fifth tab of PYAF. It is called from the QMainWindow
    in main.py.
    """

    def __init__(self, parent):
        super().__init__(parent,
                         "widget_results_experiment")

        # Load the settings
        self.settings = QtCore.QSettings()

        # Load the logger
        self.logger = logging.getLogger()

        self.canvas_resolution = parent.canvas_resolution
        self.canvas_size = [780, 350]
        self.hist_checkboxes = []
        self.stats_checkboxes = []
        self.groups_sample_size_labels = []
        self.groups_means_labels = []
        self.groups_medians_labels = []
        self.groups_sds_labels = []
        self.groups_modes_labels = []
        self.groups_pvalue_labels = []
        self.groups_significance_labels = []

        sizes = [self.canvas_size[0] / self.canvas_resolution,
                 self.canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]

        # Canvas for the histograms (groups)
        self.canvas_stack = QtWidgets.QStackedWidget()

        self.canvas_hist = QtWidgets.QWidget()
        self.canvas_hist.setFixedSize(self.canvas_size[0], self.canvas_size[1])
        self.MPL_canvas_hist = PYAFPlot(self, "results_experiment", self.canvas_hist, sizes)
        self.MPL_canvas_hist.mpl_connect("button_press_event", self.canvas_press)

        self.canvas_stats = QtWidgets.QWidget()
        self.canvas_stats.setFixedSize(self.canvas_size[0], self.canvas_size[1])
        self.MPL_canvas_stats = PYAFPlot(self, "results_stats", self.canvas_stats, sizes)
        self.MPL_canvas_stats.mpl_connect("button_press_event", self.canvas_press)

        self.canvas_stack.addWidget(self.MPL_canvas_hist)
        self.canvas_stack.addWidget(self.MPL_canvas_stats)

        name = "widget_results_hist_options_experiment"
        self.widget_options = ResultsExperimentOptionsWidget(self, name)

        # The histogram's data and statistical results
        self.GL_tableWidget_groups = QtWidgets.QGridLayout()

        # Histogram data ###############################################################################################

        # Number of rows
        rows = len(shared.exp.conditions_list) - 1
        # Define small sizes, will for the tablewidget to shrink so that the
        # window fits for 13 inch screens.
        hist_w = 256
        hist_h = 50
        self.tableWidget = PYAFTableWidget(hist_w, hist_h, rows, 7)
        behavior = QtWidgets.QAbstractItemView.SelectRows
        self.tableWidget.setSelectionBehavior(behavior)
        self.tableWidget.horizontalHeader().setHighlightSections(False)
        policy = QtCore.Qt.ScrollBarAlwaysOff
        self.tableWidget.setHorizontalScrollBarPolicy(policy)
        labels = ["", "Condition", "Sample size", "Mean", "Median", "Sd", "Mode"]
        custom_header = QCheckableHeader(QtCore.Qt.Horizontal, self.tableWidget)
        self.tableWidget.setHorizontalHeader(custom_header)
        self.tableWidget.setHorizontalHeaderLabels(labels)
        hist_headers = self.tableWidget.horizontalHeader()
        hist_headers.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        hist_headers.customContextMenuRequested. \
            connect(self.popUpMenuTableHeaderGroups)

        self.tableWidget.verticalHeader().hide()
        for i in range(len(shared.exp.conditions_list) - 1):
            self.add_row_in_tablewidget_conditions(i)

        # Let the tablewidget expand
        pol = QtWidgets.QSizePolicy.MinimumExpanding
        sizePolicy = QtWidgets.QSizePolicy(pol, pol)
        self.tableWidget.setSizePolicy(sizePolicy)

        # Let only the name column be stretched
        self.tableWidget.setColumnWidth(0, 30)  # Minimal size for checkboxes
        hist_headers.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)

        # Stats table ##################################################################################################

        # Number of rows
        rows = len(shared.exp.sample_groups)
        # Define small sizes, will for the tablewidget to shrink so that the
        # window fits for 13 inch screens.
        stats_w = 60
        stats_h = 50
        self.stats_tableWidget = PYAFTableWidget(stats_w, stats_h, rows, 4)
        behavior = QtWidgets.QAbstractItemView.SelectRows
        self.stats_tableWidget.setSelectionBehavior(behavior)
        self.stats_tableWidget.horizontalHeader().setHighlightSections(False)
        policy = QtCore.Qt.ScrollBarAlwaysOff
        self.stats_tableWidget.setHorizontalScrollBarPolicy(policy)
        labels = ["", "Pair", "P-Value", "Significance"]
        custom_header = QCheckableHeader(QtCore.Qt.Horizontal, self.stats_tableWidget)
        self.stats_tableWidget.setHorizontalHeader(custom_header)
        self.stats_tableWidget.setHorizontalHeaderLabels(labels)
        stats_headers = self.stats_tableWidget.horizontalHeader()
        stats_headers.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        # headers.customContextMenuRequested. \
        # connect(self.popUpMenuTableHeaderGroups)

        self.stats_tableWidget.verticalHeader().hide()

        # Set maximum size for stats widget
        max_w = 400
        self.stats_tableWidget.setMaximumWidth(max_w)

        # Let the tablewidget expand
        pol = QtWidgets.QSizePolicy.MinimumExpanding
        sizePolicy = QtWidgets.QSizePolicy(pol, pol)
        self.stats_tableWidget.setSizePolicy(sizePolicy)

        # Let only the name column be stretched
        self.stats_tableWidget.setColumnWidth(0, 30)  # Minimal size for checkboxes
        stats_headers.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)

        widgets_list.widget_progressbar.update()

        # Buttons for calculations/actions on histograms
        self.HL_hist_buttons = QtWidgets.QHBoxLayout()

        name = "widget_results_chooser_groups"
        self.widget_result_chooser = ResultsChooserWidget(self, name)

        # self.list_values = PYAFComboBox(self, "values")
        # self.list_values.addItem("Values")
        # self.list_values.addItem("Frequencies")

        self.HL_hist_buttons.addWidget(self.widget_result_chooser)
        # self.HL_hist_buttons.addWidget(self.list_values)
        self.HL_hist_buttons.addStretch(1)

        # Buttons for calculations/actions on stats
        self.HL_stats_buttons = QtWidgets.QHBoxLayout()

        self.BT_add_rule = PYAFButton(self, "add_rule", "Add rule")
        self.BT_auto_group = PYAFButton(self, "auto_group", "Auto group")

        self.HL_stats_buttons.addWidget(self.BT_add_rule)
        self.HL_stats_buttons.addWidget(self.BT_auto_group)
        self.HL_stats_buttons.addStretch(1)

        # Leave some blank space on the right
        spacer = QtWidgets.QSpacerItem(30, 10)

        self.GL_tableWidget_groups.addWidget(self.tableWidget, 0, 0)
        self.GL_tableWidget_groups.addLayout(self.HL_hist_buttons, 1, 0)
        self.GL_tableWidget_groups.addWidget(self.stats_tableWidget, 0, 1)
        self.GL_tableWidget_groups.addLayout(self.HL_stats_buttons, 1, 1)
        self.GL_tableWidget_groups.addItem(spacer, 0, 2, 1, 0)

        # Connect header checkboxes
        # hist_headers.all_selected.connect(lambda: self.check_all(mode="hist"))
        # stats_headers.all_selected.connect(lambda: self.check_all(mode="stats"))

        self.GL_results_data = QtWidgets.QGridLayout()
        self.GL_results_data.addWidget(self.widget_options, 0, 0)
        self.GL_results_data.addWidget(self.canvas_stack, 0, 1)
        self.GL_results_data.addLayout(self.GL_tableWidget_groups, 1, 0, 1, 0)

        self.setLayout(self.GL_results_data)
        widgets_list.widget_progressbar.update()

        stat_tools.fetch_conditions_data()

        self.update_widget()

    def update_MPL(self):
        """Update the plot."""
        self.MPL_canvas_hist.update_plot()
        self.MPL_canvas_stats.update_plot()

    def update_widget(self):
        """Updates the widget.

        No fetch_group_data() here, it would overwrite the mode value when
        doing get_pdf in single results.
        """
        widgets_list.widget_results_chooser_groups.update_widget()
        widgets_list.widget_results_hist_options_groups.update_widget()
        self.update_GUI("all")
        self.update_MPL()

    def canvas_press(self, event):
        """Called when clicking on the plot."""
        if event.button == 3 and event.inaxes == self.MPL_canvas_hist.canvas.axes:
            # Get position in screen coordinates (top left position of canvas)
            globalpos = self.MPL_canvas_hist.mapToGlobal(QtCore.QPoint(0, 0))
            # Recalculate positions with canvas position
            x = globalpos.x() + event.x
            y = globalpos.y() + (self.canvas_size[1] - event.y)
            pos = QtCore.QPoint(x, y)
            self.popUpMenu(pos, "histogram_groups")

    def list_updated(self, name):
        """Called when a list is updated."""
        #if name == "values":
        #    if self.list_values.currentIndex() == 0:
        #        shared.exp.hist_values_or_frequencies = "values"
        #    else:
        #        shared.exp.hist_values_or_frequencies = "frequencies"

        self.update_GUI("table_labels")
        widgets_list.widget_results_groups.update_GUI("list_values")
        widgets_list.widget_results_groups.update_GUI("table_labels")

    def add_row_in_tablewidget_conditions(self, i):
        """Adds a row in the groups results table."""
        # Default row height is 30, but it is too big
        self.tableWidget.setRowHeight(i, 40)

        center = QtCore.Qt.AlignHCenter + QtCore.Qt.AlignVCenter
        condition = shared.exp.conditions_list[i + 1]

        # Checkbox (hide or display result)
        cb = gui_tools.CenteredCellCheckbox(self, "display_result_hist", i)
        self.tableWidget.setCellWidget(i, 0, cb)
        self.tableWidget.cellWidget(i, 0).checkbox.setChecked(condition.display)
        self.hist_checkboxes.append(cb)
        self.check_checkboxes_state(mode="hist")

        item = QtWidgets.QTableWidgetItem(str(condition.name))
        # Do not add Qt.IsEditable here, will disable editing capability
        item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        item.setFont(self.parent.smallfont)
        item.setTextAlignment(center)
        self.tableWidget.setItem(i, 1, item)

        label_sample_size = QtWidgets.QLabel()
        label_sample_size.setFont(self.parent.smallfont)
        self.groups_sample_size_labels.append(label_sample_size)
        label = gui_tools.CenteredCellwidget(label_sample_size)
        self.tableWidget.setCellWidget(i, 2, label)

        label_mean = QtWidgets.QLabel()
        label_mean.setFont(self.parent.smallfont)
        self.groups_means_labels.append(label_mean)
        label = gui_tools.CenteredCellwidget(label_mean)
        self.tableWidget.setCellWidget(i, 3, label)

        label_median = QtWidgets.QLabel()
        label_median.setFont(self.parent.smallfont)
        self.groups_medians_labels.append(label_median)
        label = gui_tools.CenteredCellwidget(label_median)
        self.tableWidget.setCellWidget(i, 4, label)

        label_sd = QtWidgets.QLabel()
        label_sd.setFont(self.parent.smallfont)
        self.groups_sds_labels.append(label_sd)
        label = gui_tools.CenteredCellwidget(label_sd)
        self.tableWidget.setCellWidget(i, 5, label)

        label_mode = QtWidgets.QLabel()
        label_mode.setFont(self.parent.smallfont)
        self.groups_modes_labels.append(label_mode)
        label = gui_tools.CenteredCellwidget(label_mode)
        self.tableWidget.setCellWidget(i, 6, label)

    def add_row_in_stats_tablewidget_conditions(self, i):
        """Adds a row in the experiment results statistics table."""
        # Default row height is 30, but it is too big
        self.stats_tableWidget.setRowHeight(i, 25)

        center = QtCore.Qt.AlignHCenter + QtCore.Qt.AlignVCenter
        pair = shared.exp.sample_groups[i]
        condition1 = shared.exp.conditions_list[pair[0]]
        condition2 = shared.exp.conditions_list[pair[1]]

        # Checkbox (hide or display result)
        cb = gui_tools.CenteredCellCheckbox(self, "display_result_stats", i)
        self.stats_tableWidget.setCellWidget(i, 0, cb)
        # self.tableWidget.cellWidget(i, 0).checkbox.setChecked(condition.display) # Implement stat result class
        self.stats_checkboxes.append(cb)
        self.check_checkboxes_state(mode="stats")

        item = QtWidgets.QTableWidgetItem(f"{condition1.name} - {condition2.name}")
        # Do not add Qt.IsEditable here, will disable editing capability
        item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        item.setFont(self.parent.smallfont)
        item.setTextAlignment(center)
        self.stats_tableWidget.setItem(i, 1, item)

        label_pvalue = QtWidgets.QLabel()
        label_pvalue.setFont(self.parent.smallfont)
        self.groups_pvalue_labels.append(label_pvalue)
        label = gui_tools.CenteredCellwidget(label_pvalue)
        self.stats_tableWidget.setCellWidget(i, 2, label)

        label_significance = QtWidgets.QLabel()
        label_significance.setFont(self.parent.smallfont)
        self.groups_significance_labels.append(label_significance)
        label = gui_tools.CenteredCellwidget(label_significance)
        self.stats_tableWidget.setCellWidget(i, 3, label)

    def checkbox_clicked(self, name, row):
        """Called when a checkbox is clicked."""
        if name == "display_result_hist":
            # +1 = empty group (None)
            condition_id = row + 1

            # Save the value
            val = self.tableWidget.cellWidget(row, 0).checkbox.isChecked()
            shared.exp.conditions_list[condition_id].display = val

            # Check hist checkboxes state
            self.check_checkboxes_state(mode="hist")

            # Refresh plot
            self.update_MPL()

        if name == "display_result_stats":
            # +1 = empty group (None)
            condition_id = row + 1

            # Save the value
            val = self.tableWidget.cellWidget(row, 0).checkbox.isChecked()
            shared.exp.conditions_list[condition_id].display = val

            # Check stats checkboxes state
            self.check_checkboxes_state(mode="stats")

            # Refresh plot
            self.update_MPL()

    def check_checkboxes_state(self, mode):
        """Checks the state of all the checkboxes in the experiment results widget"""
        if mode == "hist":
            states = [c.checkbox.isChecked() for c in self.hist_checkboxes]
            hist_headers = self.tableWidget.horizontalHeader()
            all_checked_flag = False not in states
            hist_headers.changeState(all_checked_flag)

        elif mode == "stats":
            states = [c.checkbox.isChecked() for c in self.stats_checkboxes]
            stats_headers = self.stats_tableWidget.horizontalHeader()
            all_checked_flag = False not in states
            stats_headers.changeState(all_checked_flag)

    def check_all(self, mode):
        """Checks or unchecks all checkboxes in the table based on the state of the header checkbox"""
        if mode == "hist":
            headers = self.tableWidget.horizontalHeader()
            checkboxes_list = self.hist_checkboxes

        elif mode == "stats":
            headers = self.stats_tableWidget.horizontalHeader()
            checkboxes_list = self.stats_checkboxes

        else:
            checkboxes_list = []

        for c in checkboxes_list:
            c.checkbox.setChecked(headers.select_all)
            shared.exp.groups_list[int(c.cell_id)].display = headers.select_all

        self.update_widget()

    def button_clicked(self, button):
        """Method called upon a click on a button."""
        if button == "add_rule":
            self.open_sample_grouping_dialog()

        if button == "auto_group":
            self.auto_group_conditions()

    def update_GUI(self, what):
        """Update the GUI."""
        if what == "table_labels" or what == "all":
            if shared.conditions_values is not None:
                for i in range(len(shared.conditions_values)):
                    norm = shared.exp.norm_hist_experiment
                    if shared.conditions_factors is not None and norm is False:
                        factor = shared.conditions_factors[i]
                        print(factor)
                    elif shared.conditions_factors is None or norm:
                        factor = 1.0

                    sample_size = str(0)
                    mean = str(0)
                    median = str(0)
                    std = str(0)
                    mode = str(0)

                    values = [0, 0, 0, 0]

                    if shared.exp.hist_values_or_frequencies == "values":
                        values = shared.conditions_values[i]
                        factor = 1.0
                    else:
                        if shared.conditions_frequencies is not None:
                            if shared.conditions_frequencies[i]:
                                values = shared.conditions_frequencies[i]

                    # Write the labels
                    if values[0] != 0:
                        sample_size = str(values[0])
                    if values[1] != 0:
                        mean = str(round(values[1] * factor, 6))
                    if values[2] != 0:
                        median = str(round(values[2] * factor, 6))
                    if values[3] != 0 and values[3] is not None:
                        std = str(round(values[3] * factor, 6))
                    if values[4] != 0:
                        mode = str(round(values[4] * factor, 6))

                    self.groups_sample_size_labels[i - 1].setText(sample_size)
                    self.groups_means_labels[i - 1].setText(mean)
                    self.groups_medians_labels[i - 1].setText(median)
                    self.groups_sds_labels[i - 1].setText(std)
                    self.groups_modes_labels[i - 1].setText(mode)

        if what == "stat_table_labels" or what == "all":
            pass

        # if what == "list_values" or what == "all":
        #    if shared.exp.hist_values_or_frequencies == "values":
        #        self.list_values.setCurrentIndex(0)
        #    elif shared.exp.hist_values_or_frequencies == "frequencies":
        #        self.list_values.setCurrentIndex(1)

    def add_condition(self):
        """Add a new condition."""
        val = len(shared.exp.conditions_list)
        shared.exp.conditions_list.append(ConditionGroup(val))
        self.tableWidget.insertRow(val - 1)
        self.add_row_in_tablewidget_conditions(val - 1)
        # Set the color. Go through the list and determine a color
        # periodically.
        # -2 because we don't take into account the first "None" group
        pos = len(shared.exp.groups_list) - 2
        while pos >= len(shared.colors_list):
            pos = pos - len(shared.colors_list)
        shared.exp.conditions_list[-1].color = shared.colors_list[pos]

        # Update the tablewidget (will also refresh this tab)
        widgets_list.widget_results_single.reset_table()  # Implement reset group table, but how? Do not eliminate all choices by user.

    def open_sample_grouping_dialog(self):
        """Open the advanced options for the statistical analysis Window."""
        if widgets_list.widget_sample_grouping_dialog is None:
            # Create new widget
            SampleGroupingDialog(self)
            widgets_list.widget_sample_grouping_dialog.activateWindow()
            widgets_list.widget_sample_grouping_dialog.show()
        else:
            # Bring to front
            widgets_list.widget_sample_grouping_dialog.activateWindow()
            widgets_list.widget_sample_grouping_dialog.raise_()

    def auto_group_conditions(self):
        group_size = 2    # Allow user to change it depending on the statistical test
        conditions_indices = [c.condition_id for c in shared.exp.conditions_list[1:]]
        groups = combinations(conditions_indices, group_size)   # Ignore index 0 corresponding to None.
        for grp in groups:
            self.add_sample_group(grp)

    def filter_sample_conditions(self):
        sample_pairs = []
        rows = len(shared.exp.conditions_list) - 1
        for i in range(rows):
            cond1 = i + 1
            for cond in shared.exp.conditions_list:
                if cond.name == self.tableWidget.cellWidget(i, 2).currentText():
                    cond2 = cond.condition_id
                    sample_pair = self.order_sample_pair((cond1, cond2))
                    sample_pairs.append(sample_pair)

        samp_groups_copy = shared.exp.sample_groups.copy()

        for i in range(len(shared.exp.sample_groups)):
            pair = shared.exp.sample_groups[i]
            if pair not in sample_pairs:
                samp_groups_copy.remove(pair)
                self.stats_tableWidget.removeRow(i)

        shared.exp.sample_groups = samp_groups_copy

    def add_sample_group(self, sample_pair):
        """Add a new sample group."""

        if sample_pair not in shared.exp.sample_groups:
            val = len(shared.exp.sample_groups)
            shared.exp.sample_groups_list.append(SampleGroup(val))
            shared.exp.sample_groups.append(sample_pair)
            self.stats_tableWidget.insertRow(val)
            self.add_row_in_stats_tablewidget_conditions(val)
            # Set the color. Go through the list and determine a color
            # periodically.
            pos = len(shared.exp.sample_groups_list) - 1
            while pos >= len(shared.colors_list):
                pos = pos - len(shared.colors_list)
            shared.exp.sample_groups_list[-1].color = shared.colors_list[pos]
            # Update the tablewidget (will also refresh this tab)
            widgets_list.widget_results_single.reset_table()  # Implement reset group table, but how? Do not eliminate all choices by user.

        else:
            return


    def change_lr_display_mode(self, ptype):
        """Change the type of plot which is displayed for the loading rates.

        Can be either a scatter plot or a single point per loading rate.
        """
        if ptype == "histogram_groups":
            if shared.exp.hist_lr_groups_display_mode == "scatter":
                shared.exp.hist_lr_groups_display_mode = "points"
            elif shared.exp.hist_lr_groups_display_mode == "points":
                shared.exp.hist_lr_groups_display_mode = "scatter"
            self.update_MPL()

    def rename_conditions(self):
        """Method to display a widget with inputs to rename the groups."""
        if widgets_list.widget_rename_conditions is None:
            # Create new widget
            RenameConditionsWidget(self)
            widgets_list.widget_rename_conditions.setWindowTitle("Rename Conditions")
            widgets_list.widget_rename_conditions.activateWindow()
            widgets_list.widget_rename_conditions.show()
        else:
            # Bring to front
            widgets_list.widget_rename_conditions.activateWindow()
            widgets_list.widget_rename_conditions.raise_()

    def popUpMenuTableHeaderGroups(self, pos):
        """Pop up menu displayed on the headers of the groups table."""
        # Get the global pos for the menu
        pos = self.tableWidget.mapToGlobal(pos)

        # Create menu
        menu = QtWidgets.QMenu()

        # Add group
        action = QtWidgets.QAction("Add Condition", self)
        action.triggered.connect(self.add_condition)
        menu.addAction(action)

        # Rename groups
        action = QtWidgets.QAction("Rename Conditions", self)
        action.triggered.connect(self.rename_conditions)
        menu.addAction(action)

        # Display menu
        menu.popup(pos, menu.menuAction())
        menu.exec_()

    def popUpMenu(self, pos, ptype):
        """Right click menu for the results plots."""
        # Create menu
        menu = QtWidgets.QMenu()

        if shared.exp.results_type == "loading_rates":
            if ptype == "histogram_single":
                if shared.exp.hist_lr_single_display_mode == "scatter":
                    txt = "Single point"
                elif shared.exp.hist_lr_single_display_mode == "points":
                    txt = "Scatter"
            elif ptype == "histogram_groups":
                if shared.exp.hist_lr_groups_display_mode == "scatter":
                    txt = "Single point"
                elif shared.exp.hist_lr_groups_display_mode == "points":
                    txt = "Scatter"

            # Display as scatter or single point with error bars
            action = QtWidgets.QAction(txt, self)
            action.triggered.connect(
                lambda: self.change_lr_display_mode(ptype))
            menu.addAction(action)

        # Display menu
        menu.popup(pos, menu.menuAction())
        menu.exec_()
