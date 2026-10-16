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

"""Widget displaying the violin plot and the results of the groups."""

import logging
import numpy
from .. import shared
from .. import widgets_list
from .. import experiment
from qtpy import QtCore, QtWidgets
from ..tools import stat_tools
from ..tools import gui_tools
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFComboBox
from ..widgets_main.subwidgets.result_chooser import ResultsChooserWidget
from ..widgets_main.subwidgets.violin_options import ViolinOptionsWidget
from ..widgets_main.results_single import VIOLIN_WIDTH, PLOT_MARGINS_WIDTH, \
    PLOT_MIN_HEIGHT
from ..widgets_small.rename_groups import RenameGroupsWidget
from ..widgets.dfs import DFSWidget
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFTableWidget
from ..experiment import ResultGroup
from .subwidgets.result_header import QCheckableHeader

# Columns of the table (the values share the remaining space)
COLUMNS = ("CB", "Group", "Files", "Condition", "Mean", "Median", "S.D.",
           "Mode")
WIDTHS = {"CB": 35, "Group": 350, "Files": 70, "Condition": 350}


def col(name):
    """Position of a column of the table."""
    return COLUMNS.index(name)


class ResultsGroupsWidget(PYAFWidget):
    """Widget displaying the violin plot and the results of the groups.

    This widget is in the fifth tab of PYAF (Grouped Plots). It looks like
    the plots tab, with the groups instead of the files: a row of the table
    is a group (the files are put in groups in the plots tab). The first
    group of shared.exp.groups_list (None, the files without group) is not
    displayed.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_results_groups")

        # Load the settings
        self.settings = QtCore.QSettings()

        # Load the logger
        self.logger = logging.getLogger()

        self.checkboxes = []
        self.means_labels = []
        self.medians_labels = []
        self.sds_labels = []
        self.mode_labels = []
        self.files_labels = []
        self.submenu_conditions = None

        self.canvas_resolution = parent.canvas_resolution
        self.canvas_size = [780, 350]

        sizes = [self.canvas_size[0] / self.canvas_resolution,
                 self.canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]

        # Violin plot of the groups, in a scroll area: the plot gets wider
        # with the number of groups (see update_MPL)
        self.canva = QtWidgets.QWidget()
        self.MPL_canva = PYAFPlot(self, "results_groups", self.canva, sizes)
        self.SA_plot = QtWidgets.QScrollArea()
        self.SA_plot.setWidget(self.MPL_canva)
        self.SA_plot.setWidgetResizable(False)
        self.SA_plot.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarAlwaysOff)
        self.SA_plot.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.SA_plot.setStyleSheet("QScrollArea { border: none; }")
        self.SA_plot.setMinimumHeight(PLOT_MIN_HEIGHT)
        self.SA_plot.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        # Redraw the plot at its new size when the tab is resized
        self.plot_resize_timer = QtCore.QTimer(self)
        self.plot_resize_timer.setSingleShot(True)
        self.plot_resize_timer.timeout.connect(self.update_MPL)
        self.SA_plot.viewport().installEventFilter(self)

        # The colors of the groups follow the palette of their violin plot
        experiment.apply_palette(shared.exp, "groups")

        name = "widget_violin_options_groups"
        self.widget_options = ViolinOptionsWidget(self, name, "groups")

        # The groups
        self.VL_tableWidget = QtWidgets.QVBoxLayout()
        # Define small sizes, will for the tablewidget to shrink so that the
        # window fits for 13 inch screens.
        w = 256
        h = 50
        self.tableWidget = PYAFTableWidget(w, h, 0, len(COLUMNS))
        behaviour = QtWidgets.QAbstractItemView.SelectRows
        self.tableWidget.setSelectionBehavior(behaviour)
        custom_header = QCheckableHeader(QtCore.Qt.Horizontal,
                                         self.tableWidget)
        self.tableWidget.setHorizontalHeader(custom_header)
        self.tableWidget.horizontalHeader().setHighlightSections(False)
        policy = QtCore.Qt.ScrollBarAlwaysOff
        self.tableWidget.setHorizontalScrollBarPolicy(policy)
        self.tableWidget.setHorizontalHeaderLabels(
            ["" if name == "CB" else name for name in COLUMNS])
        headers = self.tableWidget.horizontalHeader()
        headers.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        headers.customContextMenuRequested.connect(self.popUpMenuTableHeader)
        headers.all_selected.connect(self.check_all)

        for name in COLUMNS:
            if name in WIDTHS:
                headers.setSectionResizeMode(
                    col(name), QtWidgets.QHeaderView.Interactive)
                self.tableWidget.setColumnWidth(col(name), WIDTHS[name])
            else:
                headers.setSectionResizeMode(
                    col(name), QtWidgets.QHeaderView.Stretch)

        self.tableWidget.verticalHeader().hide()
        # Let the tablewidget expand
        pol = QtWidgets.QSizePolicy.MinimumExpanding
        sizePolicy = QtWidgets.QSizePolicy(pol, pol)
        self.tableWidget.setSizePolicy(sizePolicy)

        self.VL_tableWidget.addWidget(self.tableWidget)

        widgets_list.widget_progressbar.update()

        # Type of results and buttons
        self.HL_hist_buttons = QtWidgets.QHBoxLayout()
        self.BT_dfs = PYAFButton(self, "dfs", "DFS")

        name = "widget_results_chooser_groups"
        self.widget_result_chooser = ResultsChooserWidget(self, name)

        self.HL_hist_buttons.addWidget(self.widget_result_chooser)
        self.HL_hist_buttons.addStretch(1)
        self.HL_hist_buttons.addWidget(self.BT_dfs)

        # Plot and its options, then the results type and the buttons, then
        # the groups
        self.GL_results_data = QtWidgets.QGridLayout()
        self.GL_results_data.addWidget(self.SA_plot, 0, 0)
        self.GL_results_data.addWidget(
            self.widget_options, 0, 1, QtCore.Qt.AlignTop)
        self.GL_results_data.addLayout(self.HL_hist_buttons, 1, 0, 1, 2)
        self.GL_results_data.addLayout(self.VL_tableWidget, 2, 0, 1, 2)
        self.GL_results_data.setRowStretch(0, 3)
        self.GL_results_data.setRowStretch(2, 2)
        self.GL_results_data.setColumnStretch(0, 1)

        self.setLayout(self.GL_results_data)
        widgets_list.widget_progressbar.update()

        self.reset_table()

        stat_tools.fetch_group_data()

        self.update_widget()

    def update_widget(self):
        """Update only GUI, do not recalulate nothing."""
        widgets_list.widget_results_chooser_groups.update_widget()
        self.widget_options.update_widget()
        self.update_GUI("all")
        self.update_MPL()

    def update_MPL(self):
        """Update the plot.

        The plot fills the visible area, and gets wider (with a scroll bar)
        when there are too many groups for it.
        """
        viewport = self.SA_plot.viewport().size()
        nbr_groups = sum(1 for group in shared.exp.groups_list[1:]
                         if group.display)
        width = max(viewport.width(),
                    PLOT_MARGINS_WIDTH + nbr_groups * VIOLIN_WIDTH)
        height = max(viewport.height(), PLOT_MIN_HEIGHT)
        if self.MPL_canva.size() != QtCore.QSize(width, height):
            self.MPL_canva.resize(width, height)
        self.MPL_canva.update_plot()

    def eventFilter(self, obj, event):
        """Redraw the plot when the visible area is resized."""
        if obj is self.SA_plot.viewport() and \
                event.type() == QtCore.QEvent.Resize:
            self.plot_resize_timer.start(150)
        return super().eventFilter(obj, event)

    def reset_table(self):
        """Fill the table with the groups (without the None group)."""
        self.checkboxes = []
        self.means_labels = []
        self.medians_labels = []
        self.sds_labels = []
        self.mode_labels = []
        self.files_labels = []

        self.tableWidget.setRowCount(0)
        for row in range(len(shared.exp.groups_list) - 1):
            self.tableWidget.insertRow(row)
            self.add_row_in_tablewidget(row)

        self.check_checkboxes_state()
        self.update_GUI("table_labels")

    def add_row_in_tablewidget(self, row):
        """Adds a row in the table (row 0 is the group 1)."""
        # Default row height is 30, but it is too big
        self.tableWidget.setRowHeight(row, 40)
        group = shared.exp.groups_list[row + 1]

        # Checkbox (hide or display the group)
        cb = gui_tools.CenteredCellCheckbox(self, "display_group", row)
        cb.checkbox.setChecked(group.display)
        self.tableWidget.setCellWidget(row, col("CB"), cb)
        self.checkboxes.append(cb)

        # Name in the plots, editable
        line_edit = QtWidgets.QLineEdit(group.name)
        line_edit.setFont(self.parent.smallfont)
        line_edit.setAlignment(QtCore.Qt.AlignCenter)
        line_edit.setFrame(False)
        line_edit.setToolTip("Name of the group (edit it)")
        line_edit.editingFinished.connect(
            lambda r=row, e=line_edit: self.name_edited(r, e))
        self.tableWidget.setCellWidget(
            row, col("Group"), gui_tools.FullCellWidget(line_edit))

        # Condition (Experiment Plots tab)
        conditions = PYAFComboBox(self, None)
        conditions.setFont(self.parent.smallfont)
        for condition in shared.exp.conditions_list:
            conditions.addItem(str(condition.name))
        conditions.setCurrentIndex(group.condition)
        conditions.activated.connect(self.conditions_changed)
        self.tableWidget.setCellWidget(
            row, col("Condition"), gui_tools.FullCellWidget(conditions))

        for name, labels in (("Files", self.files_labels),
                             ("Mean", self.means_labels),
                             ("Median", self.medians_labels),
                             ("S.D.", self.sds_labels),
                             ("Mode", self.mode_labels)):
            label = QtWidgets.QLabel()
            label.setFont(self.parent.smallfont)
            labels.append(label)
            self.tableWidget.setCellWidget(
                row, col(name), gui_tools.CenteredCellwidget(label))
        self.files_labels[-1].setToolTip(
            "Number of files in the group (column Group of the Plots tab)")

    def name_edited(self, row, line_edit):
        """A group name was edited."""
        group = shared.exp.groups_list[row + 1]
        text = line_edit.text().strip()
        if text == group.name:
            return
        if not text:
            # A group needs a name
            line_edit.setText(group.name)
            return
        group.name = text
        # Names in the lists of the plots tab (keep the choices of its table)
        widgets_list.widget_results_single.save_hist_prefs()
        widgets_list.widget_results_single.reset_table()
        widgets_list.widget_results_single.update_GUI("table_labels")
        self.update_MPL()

    def conditions_changed(self, *args):
        """A condition was changed: update the experiment plots tab."""
        self.save_hist_prefs()
        stat_tools.fetch_conditions_data()
        widgets_list.widget_results_experiment.update_widget()

    def checkbox_clicked(self, name, row):
        """Called when a checkbox is clicked."""
        if name == "display_group":
            value = self.checkboxes[row].checkbox.isChecked()
            shared.exp.groups_list[row + 1].display = value
            self.check_checkboxes_state()
            self.update_MPL()

    def check_checkboxes_state(self):
        """Checks the state of all the checkboxes of the table."""
        states = [c.checkbox.isChecked() for c in self.checkboxes]
        headers = self.tableWidget.horizontalHeader()
        headers.changeState(False not in states)

    def check_all(self):
        """Checks or unchecks all checkboxes in the table based on the
        state of the header checkbox."""
        headers = self.tableWidget.horizontalHeader()
        for row, cb in enumerate(self.checkboxes):
            cb.checkbox.setChecked(headers.select_all)
            shared.exp.groups_list[row + 1].display = headers.select_all
        self.update_MPL()

    def button_clicked(self, button):
        """Method called upon a click on a button."""
        if button == "refresh_hist":
            # The values of the groups come from the ones of the files
            widgets_list.widget_results_single.button_clicked("refresh_hist")

        elif button == "dfs":
            shared.exp.dfs_mode = "groups"
            if widgets_list.widget_dfs is None:
                # Create new widget
                DFSWidget(self)
                label = "Dynamic force spectroscopy"
                widgets_list.widget_dfs.setWindowTitle(label)
                widgets_list.widget_dfs.resize(600, 420)
                widgets_list.widget_dfs.show()
            else:
                widgets_list.widget_dfs.update_widget()
                # Bring to front
                widgets_list.widget_dfs.activateWindow()
                widgets_list.widget_dfs.raise_()

    def update_GUI(self, what):
        """Update the GUI."""
        if what == "table_labels" or what == "all":
            groups = shared.exp.groups_list
            for row in range(min(len(groups) - 1, len(self.means_labels))):
                group_id = row + 1
                nbr = sum(1 for result in shared.exp.results_list
                          if result.group == group_id)
                self.files_labels[row].setText(str(nbr))

                values = [numpy.nan] * 4
                if shared.groups_values is not None and \
                        group_id < len(shared.groups_values):
                    values = shared.groups_values[group_id][:4]
                # NaN (no values): "-"
                texts = ["–" if v is None or not numpy.isfinite(v)
                         else str(round(v, 6)) for v in values]
                self.means_labels[row].setText(texts[0])
                self.medians_labels[row].setText(texts[1])
                self.sds_labels[row].setText(texts[2])
                self.mode_labels[row].setText(texts[3])

        if what == "BT_dfs" or what == "all":
            found = any(data.loading_rates_calculated
                        for data in shared.exp.list)
            self.BT_dfs.setEnabled(found)

    def add_group(self):
        """Add a new group."""
        group = ResultGroup(len(shared.exp.groups_list))
        shared.exp.groups_list.append(group)
        # Color: the palette of the groups, else a color of the list
        pos = (len(shared.exp.groups_list) - 2) % len(shared.colors_list)
        group.color = shared.colors_list[pos]
        experiment.apply_palette(shared.exp, "groups")
        self.reset_table()
        # The groups of the plots tab (keep the choices of its table)
        widgets_list.widget_results_single.save_hist_prefs()
        widgets_list.widget_results_single.reset_table()
        widgets_list.widget_results_single.update_GUI("table_labels")
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
        """Save the conditions chosen in the table (also in the results of
        the files of the groups)."""
        for row in range(self.tableWidget.rowCount()):
            group = shared.exp.groups_list[row + 1]
            cell = self.tableWidget.cellWidget(row, col("Condition"))
            if cell is None:
                continue
            group.condition = cell.widget.currentIndex()

            for result in shared.exp.results_list:
                if result.group == group.group_id:
                    result.condition = group.condition

    def display_selected(self, value):
        """Tick (value True) or untick the selected groups."""
        for item in self.tableWidget.selectionModel().selectedRows():
            row = item.row()
            self.checkboxes[row].checkbox.setChecked(value)
            shared.exp.groups_list[row + 1].display = value
        self.check_checkboxes_state()
        self.update_MPL()

    def batch_conditions_change(self, action):
        """Change the condition of the selected groups."""
        index = action.data()
        for item in self.tableWidget.selectionModel().selectedRows():
            cell = self.tableWidget.cellWidget(item.row(), col("Condition"))
            cell.widget.setCurrentIndex(index)
        self.conditions_changed()

    def popUpMenuTableHeader(self, pos):
        """Pop up menu displayed on the headers of the groups table."""
        # Get the global pos for the menu
        pos = self.tableWidget.mapToGlobal(pos)

        # Create menu
        menu = QtWidgets.QMenu()

        action = QtWidgets.QAction("Check selected", self)
        action.triggered.connect(lambda: self.display_selected(True))
        menu.addAction(action)
        action = QtWidgets.QAction("Uncheck selected", self)
        action.triggered.connect(lambda: self.display_selected(False))
        menu.addAction(action)
        menu.addSeparator()

        # Add group
        action = QtWidgets.QAction("Add group", self)
        action.triggered.connect(self.add_group)
        menu.addAction(action)

        # Rename groups
        action = QtWidgets.QAction("Rename groups", self)
        action.triggered.connect(self.rename_groups)
        menu.addAction(action)

        # Rename conditions (as in the experiment plots tab)
        action = QtWidgets.QAction("Rename conditions", self)
        action.triggered.connect(
            widgets_list.widget_results_experiment.rename_conditions)
        menu.addAction(action)

        # Batch change conditions
        action = QtWidgets.QAction("Change selected conditions", self)
        menu.addAction(action)
        self.submenu_conditions = QtWidgets.QMenu()
        for index, condition in enumerate(shared.exp.conditions_list):
            subaction = QtWidgets.QAction(str(condition.name), self)
            subaction.setData(index)
            self.submenu_conditions.addAction(subaction)
        self.submenu_conditions.triggered.connect(
            self.batch_conditions_change)
        action.setMenu(self.submenu_conditions)

        # Display menu
        menu.popup(pos, menu.menuAction())
        menu.exec_()
