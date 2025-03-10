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

import logging
from .. import shared
from .. import widgets_list
from PyQt5 import QtCore, QtWidgets
from ..tools import gui_tools
from ..tools import stat_tools
from ..widgets_main.subwidgets.result_chooser import ResultsChooserWidget
from ..widgets_main.subwidgets.result_hist_options import \
    ResultsHistOptionsWidget
from ..widgets_small.rename_groups import RenameGroupsWidget
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFTableWidget, PYAFButton, PYAFComboBox
from ..experiment import ResultGroup
from .subwidgets.result_header import QCheckableHeader


class ResultsGroupsWidget(PYAFWidget):
    """Widget displaying histograms and results (Grouped).

    This widget is in the fifth tab of PYAF. It is called from the QMainWindow
    in main.py.
    """

    def __init__(self, parent):
        super().__init__(parent,
                                                  "widget_results_groups")

        # Load the settings
        self.settings = QtCore.QSettings()

        # Load the logger
        self.logger = logging.getLogger()

        self.canvas_resolution = parent.canvas_resolution
        self.canvas_size = [780, 350]
        self.checkboxes = []
        self.groups_means_labels = []
        self.groups_medians_labels = []
        self.groups_sds_labels = []
        self.groups_modes_labels = []

        sizes = [self.canvas_size[0] / self.canvas_resolution,
                 self.canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]

        # Canvas for the histograms (groups)
        self.canvas = QtWidgets.QWidget()
        self.canvas.setFixedSize(self.canvas_size[0], self.canvas_size[1])
        self.MPL_canvas = PYAFPlot(self, "results_groups", self.canvas, sizes)
        self.MPL_canvas.mpl_connect("button_press_event", self.canvas_press)

        name = "widget_results_hist_options_groups"
        self.widget_options = ResultsHistOptionsWidget(self, name)

        # The histogram's data
        self.HL_tableWidget_groups = QtWidgets.QHBoxLayout()

        # Number of rows
        rows = len(shared.exp.groups_list) - 1
        # Define small sizes, will for the tablewidget to shrink so that the
        # window fits for 13 inch screens.
        w = 256
        h = 50
        self.tableWidget = PYAFTableWidget(w, h, rows, 7)
        behavior = QtWidgets.QAbstractItemView.SelectRows
        self.tableWidget.setSelectionBehavior(behavior)
        self.tableWidget.horizontalHeader().setHighlightSections(False)
        policy = QtCore.Qt.ScrollBarAlwaysOff
        self.tableWidget.setHorizontalScrollBarPolicy(policy)
        labels = ["", "Group", "Condition", "Mean", "Median", "Sd", "Mode"]
        custom_header = QCheckableHeader(QtCore.Qt.Horizontal, self.tableWidget)
        self.tableWidget.setHorizontalHeader(custom_header)
        self.tableWidget.setHorizontalHeaderLabels(labels)
        headers = self.tableWidget.horizontalHeader()
        headers.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        headers.customContextMenuRequested.\
            connect(self.popUpMenuTableHeaderGroups)
        headers.all_selected.connect(self.check_all)

        self.tableWidget.verticalHeader().hide()
        for i in range(len(shared.exp.groups_list) - 1):
            self.add_row_in_tablewidget_groups(i)

        # Let the tablewidget expand
        pol = QtWidgets.QSizePolicy.MinimumExpanding
        sizePolicy = QtWidgets.QSizePolicy(pol, pol)
        self.tableWidget.setSizePolicy(sizePolicy)

        # Let only the name column be stretched
        self.tableWidget.setColumnWidth(0, 30)  # Minimal size for checkboxes
        headers.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)

        # Leave some blank space on the right
        spacer = QtWidgets.QSpacerItem(400, 10)

        self.HL_tableWidget_groups.addWidget(self.tableWidget)
        self.HL_tableWidget_groups.addItem(spacer)

        widgets_list.widget_progressbar.update()

        # Buttons for calculations/actions on histograms
        self.HL_hist_buttons = QtWidgets.QHBoxLayout()

        name = "widget_results_chooser_groups"
        self.widget_result_chooser = ResultsChooserWidget(self, name)

        self.list_values = PYAFComboBox(self, "values")
        self.list_values.addItem("Values")
        self.list_values.addItem("Frequencies")

        self.BT_refresh_hist = PYAFButton(self, "refresh_hist", "Refresh")

        self.HL_hist_buttons.addWidget(self.widget_result_chooser)
        self.HL_hist_buttons.addWidget(self.list_values)
        self.HL_hist_buttons.addWidget(self.BT_refresh_hist)
        self.HL_hist_buttons.addStretch(1)

        self.GL_results_data = QtWidgets.QGridLayout()
        self.GL_results_data.addWidget(self.widget_options, 0, 0)
        self.GL_results_data.addWidget(self.MPL_canvas, 0, 1)
        self.GL_results_data.addLayout(self.HL_tableWidget_groups, 1, 0, 1, 0)
        self.GL_results_data.addLayout(self.HL_hist_buttons, 2, 0, 1, 0)

        self.setLayout(self.GL_results_data)
        widgets_list.widget_progressbar.update()

        stat_tools.fetch_group_data()

        self.update_widget()

    def update_MPL(self):
        """Update the plot."""
        self.MPL_canvas.update_plot()

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
        if event.button == 3 and event.inaxes == self.MPL_canvas.canvas.axes:
            # Get position in screen coordinates (top left position of canvas)
            globalpos = self.MPL_canvas.mapToGlobal(QtCore.QPoint(0, 0))
            # Recalculate positions with canvas position
            x = globalpos.x() + event.x
            y = globalpos.y() + (self.canvas_size[1] - event.y)
            pos = QtCore.QPoint(x, y)
            self.popUpMenu(pos, "histogram_groups")

    def list_updated(self, name):
        """Called when a list is updated."""
        if name == "values":
            if self.list_values.currentIndex() == 0:
                shared.exp.hist_values_or_frequencies = "values"
            else:
                shared.exp.hist_values_or_frequencies = "frequencies"

            self.update_GUI("table_labels")
            widgets_list.widget_results_single.update_GUI("list_values")
            widgets_list.widget_results_single.update_GUI("table_labels")

    def add_row_in_tablewidget_groups(self, i):
        """Adds a row in the groups results table."""
        # Default row height is 30, but it is too big
        self.tableWidget.setRowHeight(i, 25)

        center = QtCore.Qt.AlignHCenter + QtCore.Qt.AlignVCenter
        group = shared.exp.groups_list[i + 1]

        # Checkbox (hide or display result)
        cb = gui_tools.CenteredCellCheckbox(self, "display_result", i)
        self.tableWidget.setCellWidget(i, 0, cb)
        self.tableWidget.cellWidget(i, 0).checkbox.setChecked(group.display)
        self.checkboxes.append(cb)
        self.check_checkboxes_state()

        item = QtWidgets.QTableWidgetItem(str(group.name))
        # Do not add Qt.IsEditable here, will disable editing capability
        item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        item.setFont(self.parent.smallfont)
        item.setTextAlignment(center)
        self.tableWidget.setItem(i, 1, item)

        # Condition combobox selector
        sumlist = PYAFComboBox(self, None)
        sumlist.setFont(self.parent.smallfont)
        for condition in shared.exp.conditions_list:
            sumlist.addItem(str(condition.name))
        sumlist.setCurrentIndex(group.condition)
        self.tableWidget.setCellWidget(i, 2, sumlist)

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

    def checkbox_clicked(self, name, row):
        """Called when a checkbox is clicked."""
        if name == "display_result":
            # +1 = empty group (None)
            group_id = row + 1

            # Save the value
            val = self.tableWidget.cellWidget(row, 0).checkbox.isChecked()
            shared.exp.groups_list[group_id].display = val

            # Check checkboxes state
            self.check_checkboxes_state()

            # Refresh plot
            self.update_MPL()

    def check_checkboxes_state(self):
        """Checks the state of all the checkboxes in the single results widget"""
        states = [c.checkbox.isChecked() for c in self.checkboxes]
        headers = self.tableWidget.horizontalHeader()
        all_checked_flag = False not in states
        headers.changeState(all_checked_flag)

    def check_all(self):
        """Checks or unchecks all checkboxes in the table based on the state of the header checkbox"""
        headers = self.tableWidget.horizontalHeader()
        for c in self.checkboxes:
            c.checkbox.setChecked(headers.select_all)
            shared.exp.groups_list[int(c.cell_id)].display = headers.select_all
        self.update_widget()

    def button_clicked(self, button):
        """Method called upon a click on a button."""
        if button == "refresh_hist":
            self.save_hist_prefs()

            stat_tools.get_values()

            self.update_widget()

            stat_tools.fetch_group_data()
            stat_tools.fetch_conditions_data()
            widgets_list.widget_results_experiment.update_widget()

    def update_GUI(self, what):
        """Update the GUI."""
        if what == "table_labels" or what == "all":
            if shared.groups_values is not None:
                for i in range(len(shared.groups_values)):
                    norm = shared.exp.norm_hist_groups
                    if shared.groups_factors is not None and norm is False:
                        factor = shared.groups_factors[i]
                    elif shared.groups_factors is None or norm:
                        factor = 1.0

                    mean = str(0)
                    median = str(0)
                    std = str(0)
                    mode = str(0)

                    values = [0, 0, 0, 0]

                    if shared.exp.hist_values_or_frequencies == "values":
                        values = shared.groups_values[i]
                        factor = 1.0
                    else:
                        if shared.groups_frequencies is not None:
                            if shared.groups_frequencies[i] != []:
                                values = shared.groups_frequencies[i]

                    # Write the labels
                    if values[0] != 0:
                        mean = str(round(values[0] * factor, 6))
                    if values[1] != 0:
                        median = str(round(values[1] * factor, 6))
                    if values[2] != 0 and values[2] is not None:
                        std = str(round(values[2] * factor, 6))
                    if values[3] != 0:
                        mode = str(round(values[3] * factor, 6))

                    self.groups_means_labels[i - 1].setText(mean)
                    self.groups_medians_labels[i - 1].setText(median)
                    self.groups_sds_labels[i - 1].setText(std)
                    self.groups_modes_labels[i - 1].setText(mode)

        if what == "list_values" or what == "all":
            if shared.exp.hist_values_or_frequencies == "values":
                self.list_values.setCurrentIndex(0)
            elif shared.exp.hist_values_or_frequencies == "frequencies":
                self.list_values.setCurrentIndex(1)

    def add_group(self):
        """Add a new group."""
        val = len(shared.exp.groups_list)
        shared.exp.groups_list.append(ResultGroup(val))
        self.tableWidget.insertRow(val - 1)
        self.add_row_in_tablewidget_groups(val - 1)
        # Set the color. Go through the list and determine a color
        # periodically.
        # -2 because we don't take into account the first "None" group
        pos = len(shared.exp.groups_list) - 2
        while pos >= len(shared.colors_list):
            pos = pos - len(shared.colors_list)
        shared.exp.groups_list[-1].color = shared.colors_list[pos]
        # Update the tablewidget (will also refresh this tab)
        widgets_list.widget_results_single.reset_table()


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

    def rename_groups(self):
        """Method to display a widget with inputs to rename the groups."""
        if widgets_list.widget_rename_groups is None:
            # Create new widget
            RenameGroupsWidget(self)
            widgets_list.widget_rename_groups.setWindowTitle("Rename groups")
            widgets_list.widget_rename_groups.activateWindow()
            widgets_list.widget_rename_groups.show()
        else:
            # Bring to front
            widgets_list.widget_rename_groups.activateWindow()
            widgets_list.widget_rename_groups.raise_()

    def save_hist_prefs(self):
        """Save the currently chosen conditions options in the groups results table."""

        for i in range(1, len(shared.exp.groups_list)):
            group = shared.exp.groups_list[i]

            # Conditions
            cell = self.tableWidget.cellWidget(i - 1, 2)
            group.condition = cell.currentIndex()

            for i in range(len(shared.exp.results_list)):
                result = shared.exp.results_list[i]
                if result.group == group.group_id:
                    result.condition = cell.currentIndex()


    def popUpMenuTableHeaderGroups(self, pos):
        """Pop up menu displayed on the headers of the groups table."""
        # Get the global pos for the menu
        pos = self.tableWidget.mapToGlobal(pos)

        # Create menu
        menu = QtWidgets.QMenu()

        # Add group
        action = QtWidgets.QAction("Add group", self)
        action.triggered.connect(self.add_group)
        menu.addAction(action)

        # Rename groups
        action = QtWidgets.QAction("Rename groups", self)
        action.triggered.connect(self.rename_groups)
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
