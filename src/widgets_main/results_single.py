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

"""Widget displaying histograms and results (Single)."""

import logging
import numpy
from .. import shared
from .. import widgets_list
from .. import consts
from .. import experiment
from qtpy import QtCore, QtWidgets
from ..tools import stat_tools
from ..tools import gui_tools
from ..tools import superplot_tools
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFComboBox
from ..widgets_main.subwidgets.result_chooser import ResultsChooserWidget
from ..widgets_main.subwidgets.violin_options import ViolinOptionsWidget
from ..widgets.events_per_scan import EventsPerScanWidget
from ..widgets.events_filters import ChangeFilterWidget
from ..widgets.events_filters import DIST_TOOLTIP, SLOPE_TOOLTIP
from ..tools import event_filter_tools
from ..widgets.dfs import DFSWidget
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFTableWidget
from .subwidgets.result_header import QCheckableHeader


# Size of the violin plot: width per file and for the axis and its labels,
# minimum height
VIOLIN_WIDTH = 60
PLOT_MARGINS_WIDTH = 120
PLOT_MIN_HEIGHT = 300


class Columns:
    """Columns object

    Contains a list of columns.
    """

    def __init__(self):
        self.results_type = shared.exp.results_type

        # Create the single columns
        self.list = []
        names = ["CB", "File", "Short name", "Group", "Mean", "Median",
                 "Sd", "Mode"]

        if self.results_type == "stiffness":
            names += ["Slice", "ROI", "B.E.C."]

        elif self.results_type == "events_forces" or \
                self.results_type == "events_per_curve" or \
                self.results_type == "events_rupture_force" or \
                self.results_type == "loading_rates":
            names += ["Slope filter", "Distance filter", "ROI"]

        elif self.results_type == "events_distance":
            names += ["Slope filter", "Distance filter", "Event(s)", "ROI"]

        else:
            names += ["ROI"]

        for index, name in enumerate(names):
            self.list.append(Column(index, name))

    @property
    def len_all(self):
        """Number of columns."""
        return len(self.list)

    @property
    def len_visible(self):
        """Number of visible columns."""
        count = 0
        for column in self.list:
            if column.visibility:
                count += 1
        return count

    @property
    def names(self):
        """Column names."""
        name_list = []
        for column in self.list:
            if column.visibility:
                if column.name == "CB":
                    # Remove layer for the first column so it will not overlap with checkbox.
                    name_list.append("")
                elif column.name == "Sd":
                    name_list.append("S.D.")
                elif column.name == "Group":
                    # The groups are the replicates of the SuperPlots
                    name_list.append("Replicate")
                else:
                    name_list.append(column.name)
        return name_list

    def get_column(self, name):
        """Return a column object by name."""
        for column in self.list:
            if column.name == name:
                return column


class Column:
    """Column object"""

    def __init__(self, column_index, name):
        self.index = column_index
        self.position = column_index
        self.name = name
        self.visibility = True


class ResultsSingleWidget(PYAFWidget):
    """Widget displaying histograms and results (Single).

    The single results widget is in the fourth tab of PYAF. It is called from
    the QMainWindow in main.py.
    """

    def __init__(self, parent):
        super().__init__(
            parent, "widget_results_single")

        # Load the settings
        self.settings = QtCore.QSettings()

        # Set submenus to None
        self.submenu_slices = None
        self.submenu_rois = None
        self.submenu_events = None
        self.submenu_groups = None

        # Checkboxes
        self.checkboxes = None

        # Labels
        self.means_labels = None
        self.medians_labels = None
        self.sds_labels = None
        self.mode_labels = None

        # A menu for the pop up submenu
        self.submenu_columns = None

        # Columns object, containing references to which column is displayed.
        self.columns = Columns()

        # Load the logger
        self.logger = logging.getLogger()

        self.canvas_resolution = parent.canvas_resolution
        self.canvas_size = [780, 350]

        sizes = [self.canvas_size[0] / self.canvas_resolution,
                 self.canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]

        # Violin plot of the files, in a scroll area: the plot gets wider
        # with the number of files (see update_MPL)
        self.canva = QtWidgets.QWidget()
        self.MPL_canva = PYAFPlot(self, "results_single", self.canva, sizes)
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

        name = "widget_violin_options_single"
        self.widget_options = ViolinOptionsWidget(self, name)

        # The histogram's data and statistical results
        self.GL_tableWidget_groups = QtWidgets.QGridLayout()

        # The histogram's data
        self.VL_tableWidget = QtWidgets.QVBoxLayout()
        rows = len(shared.exp.results_list)
        # Define small sizes, will for the tablewidget to shrink so that the
        # window fits for 13 inch screens.
        w = 256
        h = 50
        self.tableWidget = PYAFTableWidget(w, h, rows, 1)
        behaviour = QtWidgets.QAbstractItemView.SelectRows
        self.tableWidget.setSelectionBehavior(behaviour)
        custom_header = QCheckableHeader(QtCore.Qt.Horizontal, self.tableWidget)
        self.tableWidget.setHorizontalHeader(custom_header)
        self.tableWidget.horizontalHeader().setHighlightSections(False)
        policy = QtCore.Qt.ScrollBarAlwaysOff
        self.tableWidget.setHorizontalScrollBarPolicy(policy)
        headers = self.tableWidget.horizontalHeader()
        headers.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        headers.customContextMenuRequested.connect(self.popUpMenuTableHeader)
        headers.all_selected.connect(self.check_all)

        self.tableWidget.verticalHeader().hide()
        # Let the tablewidget expand
        pol = QtWidgets.QSizePolicy.MinimumExpanding
        sizePolicy = QtWidgets.QSizePolicy(pol, pol)
        self.tableWidget.setSizePolicy(sizePolicy)

        self.VL_tableWidget.addWidget(self.tableWidget)

        widgets_list.widget_progressbar.set_label("Loading GUI ...")
        widgets_list.widget_progressbar.update()

        # The values are recomputed when an option of the table is changed
        # (once for several quick changes)
        self.refreshing = False
        self.refresh_timer = QtCore.QTimer(self)
        self.refresh_timer.setSingleShot(True)
        self.refresh_timer.setInterval(300)
        self.refresh_timer.timeout.connect(
            lambda: self.button_clicked("refresh_hist"))

        # Buttons for calculations/actions on histograms
        self.HL_hist_buttons = QtWidgets.QHBoxLayout()
        label = "Events per scan"
        self.BT_events_per_scan = PYAFButton(self, "events_per_scan", label)
        self.BT_dfs = PYAFButton(self, "dfs", "DFS")

        # self.list_values = PYAFComboBox(self, "values")
        # self.list_values.addItem("Values")
        # self.list_values.addItem("Frequencies")

        name = "widget_results_chooser_single"
        self.widget_result_chooser = ResultsChooserWidget(self, name)

        self.HL_hist_buttons.addWidget(self.widget_result_chooser)
        # self.HL_hist_buttons.addWidget(self.list_values)
        self.HL_hist_buttons.addStretch(1)
        self.HL_hist_buttons.addWidget(self.BT_events_per_scan)
        self.HL_hist_buttons.addWidget(self.BT_dfs)

        # Plot and its options, then the results type and the buttons, then
        # the files
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

        # The values are computed when the plots tabs are opened
        self.reset_table()

        self.update_widget()

    def update_widget(self):
        """Update only GUI, do not recalulate nothing."""
        widgets_list.widget_results_chooser_single.update_widget()
        self.widget_options.update_widget()
        self.update_GUI("all")
        self.update_MPL()

    def update_MPL(self):
        """Update the plot.

        The plot fills the visible area, and the violin plots get wider
        (with a scroll bar) when there are too many files for them (not the
        force vs loading rate plot).
        """
        viewport = self.SA_plot.viewport().size()
        nbr_files = sum(1 for result in shared.exp.results_list
                        if result.display)
        width = max(viewport.width(),
                    PLOT_MARGINS_WIDTH + nbr_files * VIOLIN_WIDTH)
        if shared.exp.results_type == "loading_rates":
            width = viewport.width()
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

    def list_updated(self, name):
        """Called when a list is updated."""
        if name == "values":
            if self.list_values.currentIndex() == 0:
                shared.exp.hist_values_or_frequencies = "values"
            else:
                shared.exp.hist_values_or_frequencies = "frequencies"

            self.update_GUI("table_labels")

            widgets_list.widget_results_groups.update_GUI("list_values")
            widgets_list.widget_results_groups.update_GUI("table_labels")

    def reset_table(self):
        """Updates the result's table widget.

        Will change the number of columns, depending of what type of results
        have to be displayed. It will then repopulate the tablewidget with
        the computed results and aviable options.
        """
        if shared.exp.results_type != self.columns.results_type:
            # The type has changed, create a new column object
            self.columns = Columns()

        # Create the headers
        self.tableWidget.setColumnCount(self.columns.len_visible)
        self.tableWidget.setHorizontalHeaderLabels(self.columns.names)
        superplot_tools.set_header_tooltips(self.tableWidget, "plots")

        # Set specific column widths
        headers = self.tableWidget.horizontalHeader()

        # Widths of the columns (the file names take the remaining space)
        widths = {"CB": 35, "Short name": 170, "Group": 140, "Mean": 90,
                  "Median": 90, "Sd": 90, "Mode": 90, "Slope filter": 110,
                  "Distance filter": 160}
        for column in self.columns.list:
            if not column.visibility:
                continue
            if column.name == "File":
                headers.setSectionResizeMode(
                    column.position, QtWidgets.QHeaderView.Stretch)
            else:
                headers.setSectionResizeMode(
                    column.position, QtWidgets.QHeaderView.Interactive)
                self.tableWidget.setColumnWidth(
                    column.position, widths.get(column.name, 120))

        self.means_labels = []
        self.medians_labels = []
        self.sds_labels = []
        self.mode_labels = []
        self.checkboxes = []
        self.means_labels = []
        self.medians_labels = []
        self.sds_labels = []
        self.mode_labels = []
        self.checkboxes = []

        if shared.exp.results_list:
            # Delete all the rows
            for i in range(self.tableWidget.rowCount(), -1, -1):
                self.tableWidget.removeRow(i)

        self.tableWidget.setCurrentCell(0, 0)  # Select first row

        # Create new rows
        for i in range(len(shared.exp.results_list)):
            self.tableWidget.insertRow(i)
            self.add_row_in_tablewidget(i)

    def add_row_in_tablewidget(self, row_index):
        """Adds a row in the result's tablewidget."""
        # Default row height is 30, but it is too big
        self.tableWidget.setRowHeight(row_index, 40)

        if self.columns.get_column("CB").visibility:
            self.add_cell_at_pos("CB", row_index)

        if self.columns.get_column("File").visibility:
            self.add_cell_at_pos("File", row_index)

        if self.columns.get_column("Short name").visibility:
            self.add_cell_at_pos("Short name", row_index)

        if self.columns.get_column("Group").visibility:
            self.add_cell_at_pos("Group", row_index)

        if self.columns.get_column("Mean").visibility:
            self.add_cell_at_pos("Mean", row_index)

        if self.columns.get_column("Median").visibility:
            self.add_cell_at_pos("Median", row_index)

        if self.columns.get_column("Sd").visibility:
            self.add_cell_at_pos("Sd", row_index)

        if self.columns.get_column("Mode").visibility:
            self.add_cell_at_pos("Mode", row_index)

        rtype = shared.exp.results_type

        if rtype == "stiffness":
            self.add_cell_at_pos("Slice", row_index)

            self.add_cell_at_pos("B.E.C.", row_index)

            self.add_cell_at_pos("ROI", row_index)

        elif rtype == "events_forces" or rtype == "events_per_curve" or \
                rtype == "events_rupture_force" or rtype == "loading_rates" \
                or rtype == "events_distance":
            self.add_cell_at_pos("Slope filter", row_index)

            self.add_cell_at_pos("Distance filter", row_index)

            if rtype == "events_distance":
                self.add_cell_at_pos("Event(s)", row_index)

                self.add_cell_at_pos("ROI", row_index)

            else:
                self.add_cell_at_pos("ROI", row_index)

        else:
            self.add_cell_at_pos("ROI", row_index)

    def add_cell_at_pos(self, name, row_index):
        """Add a cell to a colum, by name and row index."""
        # Get the result to be displayed
        result = shared.exp.results_list[row_index]
        dt = shared.exp.list[result.data_id]

        column_pos = self.columns.get_column(name).position

        if name == "File":
            label = QtWidgets.QLabel()
            label.setFont(self.parent.smallfont)
            label.setText(result.name)
            label = gui_tools.CenteredCellwidget(label)
            self.tableWidget.setCellWidget(row_index, column_pos, label)

            # Text for the tooltip
            name = shared.exp.list[result.data_id].filename
            cell = self.tableWidget.cellWidget(row_index, column_pos)
            cell.setToolTip("Copy from file " + name)

        elif name == "Short name":
            # Name in the plots, editable
            line_edit = QtWidgets.QLineEdit(experiment.get_short_name(result))
            line_edit.setFont(self.parent.smallfont)
            line_edit.setAlignment(QtCore.Qt.AlignCenter)
            line_edit.setFrame(False)
            line_edit.setToolTip(
                "Name in the plots (edit it, empty: automatic name)")
            line_edit.editingFinished.connect(
                lambda r=row_index, e=line_edit: self.short_name_edited(r, e))
            # In a layout: keeps the width of the cell (see FullCellWidget)
            self.tableWidget.setCellWidget(
                row_index, column_pos, gui_tools.FullCellWidget(line_edit))

        elif name == "CB":
            # Checkbox (hide or display histogram)
            cb = gui_tools.CenteredCellCheckbox(
                self, "display_hist", str(row_index))
            self.tableWidget.setCellWidget(row_index, column_pos, cb)
            self.tableWidget.cellWidget(
                row_index, column_pos).checkbox.setChecked(result.display)
            self.checkboxes.append(cb)
            self.check_checkboxes_state()

        elif name == "Group":
            sumlist = PYAFComboBox(self, None)
            sumlist.setFont(self.parent.smallfont)
            for group in shared.exp.groups_list:
                sumlist.addItem(str(group.name))
            sumlist.setCurrentIndex(result.group)
            sumlist.activated.connect(self.schedule_refresh)
            self.tableWidget.setCellWidget(
                row_index, column_pos, gui_tools.FullCellWidget(sumlist))

        elif name == "Mean":
            label_mean = QtWidgets.QLabel()
            label_mean.setFont(self.parent.smallfont)
            self.means_labels.append(label_mean)
            label = gui_tools.CenteredCellwidget(label_mean)
            self.tableWidget.setCellWidget(row_index, column_pos, label)

        elif name == "Median":
            label_median = QtWidgets.QLabel()
            label_median.setFont(self.parent.smallfont)
            self.medians_labels.append(label_median)
            label = gui_tools.CenteredCellwidget(label_median)
            self.tableWidget.setCellWidget(row_index, column_pos, label)

        elif name == "Sd":
            label_sd = QtWidgets.QLabel()
            label_sd.setFont(self.parent.smallfont)
            self.sds_labels.append(label_sd)
            label = gui_tools.CenteredCellwidget(label_sd)
            self.tableWidget.setCellWidget(row_index, column_pos, label)

        elif name == "Mode":
            label_mode = QtWidgets.QLabel()
            label_mode.setFont(self.parent.smallfont)
            self.mode_labels.append(label_mode)
            label = gui_tools.CenteredCellwidget(label_mode)
            self.tableWidget.setCellWidget(row_index, column_pos, label)

        elif name == "ROI":
            roi_list = PYAFComboBox(self, None)
            roi_list.setFont(self.parent.smallfont)
            roi_list.addItem("None")
            for j in range(len(shared.exp.list[result.data_id].roi_list)):
                roi_list.addItem(str(j + 1))
            roi_list.setCurrentIndex(result.roi)
            roi_list.activated.connect(self.schedule_refresh)
            rois = gui_tools.CenteredCellComboBox(roi_list)
            self.tableWidget.setCellWidget(row_index, column_pos, rois)

        elif name == "Slice":
            # Add a list of indentation depths for each row (first choice is
            # "All", which gives an histogram with all the slices)
            list_indent_hist = PYAFComboBox(self, None)
            list_indent_hist.setFont(self.parent.smallfont)
            if dt.stiffness_calculated:
                if dt.indentation_step != 0:
                    list_indent_hist.addItem("All")
                    list_indent_hist.setEnabled(True)
                    for item in dt.indentation_list:
                        list_indent_hist.addItem(item)
                else:
                    # Start to max
                    list_indent_hist.addItem(dt.indentation_list[0])
                    list_indent_hist.setEnabled(False)
                list_indent_hist.setCurrentIndex(result.slice)
            else:
                list_indent_hist.setEnabled(False)

            list_indent_hist.activated.connect(self.schedule_refresh)
            ind = gui_tools.CenteredCellComboBox(list_indent_hist)
            self.tableWidget.setCellWidget(row_index, column_pos, ind)

        elif name == "B.E.C.":
            # Add list for the stiffness correction
            corr_list = PYAFComboBox(self, None)
            corr_list.setFont(self.parent.smallfont)
            corr_list.addItem("False")
            corr_list.addItem("True")
            corr_list.setCurrentIndex(result.sub_type)
            if not shared.exp.list[result.data_id].stiffness_corrected:
                corr_list.setEnabled(False)

            corr_list.activated.connect(self.schedule_refresh)
            corr = gui_tools.CenteredCellComboBox(corr_list)
            self.tableWidget.setCellWidget(row_index, column_pos, corr)

        elif name == "Slope filter" or name == "Distance filter":
            # The filters are changed in a small window (click on the cell)
            button = gui_tools.ClickableCellLabel()
            button.setFont(self.parent.smallfont)
            filt = "slope" if name == "Slope filter" else "dist_filt"
            tooltip = SLOPE_TOOLTIP if filt == "slope" else DIST_TOOLTIP
            button.setToolTip(tooltip + "\nClick to change it.")
            button.clicked.connect(
                lambda b=button, f=filt: self.filter_cell_clicked(b, f))
            widget = gui_tools.CenteredCellwidget2([button])
            self.tableWidget.setCellWidget(row_index, column_pos, widget)
            self.update_filter_cells(row_index)

        elif name == "Event(s)":
            # Event to select (Distance to Joc)
            input_dist_to_joc = PYAFComboBox(self, None)
            input_dist_to_joc.setFont(self.parent.smallfont)
            input_dist_to_joc.addItem("All")
            input_dist_to_joc.addItem("Last")
            input_dist_to_joc.setFont(self.parent.smallfont)
            input_dist_to_joc.setCurrentIndex(result.dist_to_joc)
            input_dist_to_joc.activated.connect(self.schedule_refresh)
            dist = gui_tools.CenteredCellComboBox(input_dist_to_joc)
            self.tableWidget.setCellWidget(row_index, column_pos, dist)

    def update_filter_cells(self, row_index):
        """Display the filters of a row of the table."""
        result = shared.exp.results_list[row_index]
        col = self.columns.get_column("Slope filter")
        if col is not None and col.visibility:
            cell = self.tableWidget.cellWidget(row_index, col.position)
            if cell is not None:
                cell.items[0].setText(event_filter_tools.slope_filter_text(
                    result.filter_slope_min, result.filter_slope_max))
        col = self.columns.get_column("Distance filter")
        if col is not None and col.visibility:
            cell = self.tableWidget.cellWidget(row_index, col.position)
            if cell is not None:
                cell.items[0].setText(event_filter_tools.dist_filter_text(
                    result.filter_dist_min, result.filter_dist_max,
                    result.filter_dist_keep))

    def filter_cell_clicked(self, button, filt):
        """Change the filter of the row of the clicked cell (and of the
        other selected rows if it is one of them)."""
        position = button.parentWidget().pos()
        row = self.tableWidget.indexAt(position).row()
        selected = [item.row() for item in
                    self.tableWidget.selectionModel().selectedRows()]
        if row not in selected:
            self.tableWidget.selectRow(row)
        self.change_filters(filt)

    def save_hist_prefs(self):
        """Save the currently chosen options in the results table."""
        for i in range(len(shared.exp.results_list)):
            result = shared.exp.results_list[i]

            # Groups
            col = self.columns.get_column("Group")
            if col is not None and col.visibility:
                cell = self.tableWidget.cellWidget(i, col.position)
                result.group = cell.widget.currentIndex()

            # Checkboxes
            col = self.columns.get_column("CB")
            if col is not None and col.visibility:
                cell = self.tableWidget.cellWidget(i, col.position)
                result.display = cell.checkbox.isChecked()

            # Rois
            col = self.columns.get_column("ROI")
            if col is not None and col.visibility:
                cell = self.tableWidget.cellWidget(i, col.position)
                roi_id = cell.list.currentIndex()
                result.roi = roi_id

            col = self.columns.get_column("Slice")
            if col is not None and col.visibility:
                cell = self.tableWidget.cellWidget(i, col.position)
                result.slice = cell.list.currentIndex()

            col = self.columns.get_column("Event(s)")
            if col is not None and col.visibility:
                cell = self.tableWidget.cellWidget(i, col.position)
                result.dist_to_joc = cell.list.currentIndex()

            col = self.columns.get_column("B.E.C.")
            if col is not None and col.visibility:
                cell = self.tableWidget.cellWidget(i, col.position)
                sub_type = cell.list.currentIndex()
                result.sub_type = sub_type

    def add_group(self):
        """Add a group (the table is rebuilt: keep its current choices)."""
        self.save_hist_prefs()
        widgets_list.widget_results_groups.add_group()

    def short_name_edited(self, row_index, line_edit):
        """A short name was edited: rename the file in the plots."""
        result = shared.exp.results_list[row_index]
        text = line_edit.text().strip()
        if text == experiment.get_short_name(result):
            return
        result.short_name = text
        # Empty: back to the automatic name
        line_edit.setText(experiment.get_short_name(result))
        self.update_MPL()

    def schedule_refresh(self, *args):
        """Recompute the values after a change of the table (the changes
        made quickly one after the other are applied together)."""
        self.refresh_timer.start()

    def button_clicked(self, button):
        """Method called upon a click on a button."""
        if button == "refresh_hist":
            if self.refreshing:
                # The progress bar processes the events: a refresh asked
                # during a refresh is done after it
                self.refresh_timer.start()
                return
            self.refreshing = True
            self.refresh_timer.stop()
            try:
                self.save_hist_prefs()

                stat_tools.get_values()

                self.update_widget()

                stat_tools.fetch_group_data()
                widgets_list.widget_results_groups.update_widget()
                widgets_list.widget_results_experiment.update_widget()
            finally:
                self.refreshing = False

        elif button == "events_per_scan":
            if widgets_list.widget_events_per_scan is None:
                # Create new widget
                EventsPerScanWidget(self)
                widgets_list.widget_events_per_scan.resize(900, 500)
                widgets_list.widget_events_per_scan.show()
            else:
                # Bring to front
                widgets_list.widget_events_per_scan.activateWindow()
                widgets_list.widget_events_per_scan.raise_()

        elif button == "dfs":
            shared.exp.dfs_mode = "single"
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

    def checkbox_clicked(self, name, data_id=None):
        """Called when a checkbox is clicked."""
        if name == "display_hist":
            data_id = int(data_id)

            col = self.columns.get_column("CB").position
            if col is not None:
                cell = self.tableWidget.cellWidget(data_id, col)
                value = cell.checkbox.isChecked()
                shared.exp.results_list[data_id].display = value
                self.check_checkboxes_state()
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
            shared.exp.results_list[int(c.cell_id)].display = headers.select_all
        self.update_widget()

    def update_GUI(self, what):
        """Update some GUI components."""
        if what == "table_labels" or what == "all":
            if shared.single_values is not None:
                for i in range(len(shared.single_values)):
                    # NaN (no values): "-"
                    texts = ["\u2013" if v is None or not numpy.isfinite(v)
                             else str(round(v, 6))
                             for v in shared.single_values[i][:4]]
                    mean, median, std, mode = texts

                    if self.columns.get_column("Mean").visibility:
                        self.means_labels[i].setText(mean)
                    if self.columns.get_column("Median").visibility:
                        self.medians_labels[i].setText(median)
                    if self.columns.get_column("Sd").visibility:
                        self.sds_labels[i].setText(std)
                    if self.columns.get_column("Mode").visibility:
                        self.mode_labels[i].setText(mode)

        if what == "BT_events_per_scan" or what == "all":
            found = False
            for data in shared.exp.list:
                if data.events_calculated:
                    found = True

            if found:
                self.BT_events_per_scan.setEnabled(True)
            else:
                self.BT_events_per_scan.setEnabled(False)

        if what == "BT_dfs" or what == "all":
            found = False
            for data in shared.exp.list:
                if data.loading_rates_calculated:
                    found = True

            if found:
                self.BT_dfs.setEnabled(True)
            else:
                self.BT_dfs.setEnabled(False)

        # if what == "list_values" or what == "all":
        #     if shared.exp.hist_values_or_frequencies == "values":
        #         self.list_values.setCurrentIndex(0)
        #     elif shared.exp.hist_values_or_frequencies == "frequencies":
        #         self.list_values.setCurrentIndex(1)

    def display_selected(self):
        """Display all the selected results."""
        for item in self.tableWidget.selectionModel().selectedRows():
            row_index = item.row()
            col = self.columns.get_column("CB").position
            if col is not None:
                cellwidget = self.tableWidget.cellWidget(row_index, col)
                cellwidget.checkbox.setChecked(True)
                shared.exp.results_list[row_index].display = True

        self.update_MPL()

    def hide_selected(self):
        """Hide all selected results."""
        for item in self.tableWidget.selectionModel().selectedRows():
            row_index = item.row()
            col = self.columns.get_column("CB").position
            if col is not None:
                cell = self.tableWidget.cellWidget(row_index, col)
                cell.checkbox.setChecked(False)
                shared.exp.results_list[row_index].display = False

        self.update_MPL()

    def change_filters(self, name):
        """Change the filters of the selected rows."""
        if name == "slope":
            widget = ChangeFilterWidget(self, "slope_batch")
        else:
            widget = ChangeFilterWidget(self, "dist_batch")
        widget.exec()

    def batch_groups_change(self):
        """Change all the selected groups."""
        # Get the group which was selected
        for i in range(len(shared.exp.groups_list)):
            text = self.submenu_groups.sender().text()
            if shared.exp.groups_list[i].name == text:
                index = i
                break
        # Save the values, change the lists, then recompute the values
        column = self.columns.get_column("Group")
        for item in self.tableWidget.selectionModel().selectedRows():
            row_index = item.row()
            shared.exp.results_list[row_index].group = index
            if column.visibility:
                cell = self.tableWidget.cellWidget(row_index, column.position)
                cell.widget.setCurrentIndex(index)
        self.schedule_refresh()

    def batch_slices_change(self, force_index=None):
        """Allows you to change all the slices menu's in the tableWidget.

        force_index can be used to change to a specific slice programmatically.
        """
        self.logger.debug("Batch slice change called")

        # When using the menu, force_index is the first argument and is set
        # as the QAction by PyQt's signal and slots. Reset the index to None
        # in this case.
        if isinstance(force_index, QtWidgets.QAction):
            force_index = None

        # Get the slice which was selected
        if force_index is None:
            if self.submenu_slices.sender().text() == "All":
                index = "All"  # Will be set to 0
            else:
                for i in range(shared.exp.global_max_indentation_index):
                    if str(i) == self.submenu_slices.sender().text():
                        index = i
                        break
        else:
            if force_index == "All":
                index = force_index
            else:
                index = force_index - 1

        self.logger.debug("Slice selected : %s", index)

        # Save the values, change the lists. You will need to refresh manually
        # to get the correct values / histograms
        for item in self.tableWidget.selectionModel().selectedRows():
            row_index = item.row()
            data_id = shared.exp.results_list[row_index].data_id

            if index == "All":
                newindex = 0
            else:
                if len(shared.exp.list[data_id].indentation_list) == 1:
                    newindex = 0
                else:
                    if shared.exp.list[data_id].max_indentation_index <= index:
                        # Get a new index (Slice 0 = empty slice)
                        dt = shared.exp.list[data_id]
                        newindex = dt.max_indentation_index - 1

                        mess = "Slice nbr " + str(index) + \
                            " doesn't exist for row " + str(row_index + 1) + \
                            ".\r\nIt was replaced by slice nbr " + \
                            str(newindex) + "."

                        if not consts.UNIT_TESTING:
                            QtWidgets.QMessageBox.warning(
                                self,
                                "Information",
                                mess)
                    else:
                        newindex = index + 1

            mess = "Slice change on row %s, slice set to %s"
            self.logger.debug(mess, row_index, newindex)

            col = self.columns.get_column("Slice").position
            if col is not None:
                cell = self.tableWidget.cellWidget(row_index, col)
                cell.list.setCurrentIndex(newindex)

        self.schedule_refresh()

    def batch_rois_change(self):
        """Changes all the selected ROIs to the wanted value."""
        self.logger.debug("Batch change rois")

        # For each file, get the number of ROIs to determine what is the
        # biggest displayable ROI.
        length = 0
        for data in shared.exp.list:
            if len(data.roi_list) > length:
                length = len(data.roi_list)

        # Get the roi which was selected
        if self.submenu_rois.sender().text() == "None":
            index = 0
        else:
            for i in range(1, length + 1):
                if str(i) == self.submenu_rois.sender().text():
                    index = i
                    break

        self.logger.debug("Selected ROI (batch) : %s", index)

        # Save the values, change the lists. You will need to refresh manually
        # to get the correct values / histograms
        for item in self.tableWidget.selectionModel().selectedRows():
            row_index = item.row()
            data_id = shared.exp.results_list[row_index].data_id

            if len(shared.exp.list[data_id].roi_list) < index:
                newindex = len(shared.exp.list[data_id].roi_list)

                if newindex == 0:
                    final = "It was replaced by : None"
                else:
                    final = "It was replaced by ROI nbr " + str(newindex) + "."

                mess = "ROI nbr " + str(index) + " doesn't exist for row " + \
                    str(row_index + 1) + ".\r\n" + final

                QtWidgets.QMessageBox.warning(self, "Information", mess)
            else:
                newindex = index

            shared.exp.results_list[row_index].roi = newindex

            col = self.columns.get_column("ROI").position
            if col is not None:
                cell = self.tableWidget.cellWidget(row_index, col)
                cell.list.setCurrentIndex(newindex)

        self.schedule_refresh()

    def batch_events_change(self):
        """Change selected events."""
        # Get the mode which was selected
        # 0: all the events, 1: the last one
        index = 0 if self.submenu_events.sender().text() == "All" else 1
        # Save the values, change the lists, then recompute the values
        for item in self.tableWidget.selectionModel().selectedRows():
            row_index = item.row()
            shared.exp.results_list[row_index].dist_to_joc = index
            col = self.columns.get_column("Event(s)").position
            if col is not None:
                cell = self.tableWidget.cellWidget(row_index, col)
                cell.list.setCurrentIndex(index)
        self.schedule_refresh()

    def hide_or_display_column(self):
        """Hide or display a column in the tablewidget."""
        # Save prefs before doing this
        self.save_hist_prefs()

        label = self.submenu_columns.sender().text()
        # Find the column
        column = self.columns.get_column(label)
        column.visibility = not column.visibility

        # Reset the positions
        pos = 0
        for column in self.columns.list:
            if column.visibility:
                column.position = pos
                pos += 1

        # Refresh the table
        self.reset_table()

    def popUpMenuTableHeader(self, pos):
        """Pop up menu on displayed when rigth clicking on the table header."""
        # Get the global pos for the menu
        pos = self.tableWidget.mapToGlobal(pos)

        # Create menu
        menu = QtWidgets.QMenu()

        rt = shared.exp.results_type

        # Display selected
        action = QtWidgets.QAction("Check selected", self)
        action.triggered.connect(self.display_selected)
        menu.addAction(action)
        action = QtWidgets.QAction("Uncheck selected", self)
        action.triggered.connect(self.hide_selected)
        menu.addAction(action)
        menu.addSeparator()

        # Hide/display columns
        action = QtWidgets.QAction("Hide/display column", self)
        menu.addAction(action)
        self.submenu_columns = QtWidgets.QMenu()
        for i in range(len(self.columns.list)):
            subaction = QtWidgets.QAction(self.columns.list[i].name, self)
            subaction.setCheckable(True)
            subaction.setChecked(self.columns.list[i].visibility)
            self.submenu_columns.addAction(subaction)
        self.submenu_columns.triggered.connect(self.hide_or_display_column)
        action.setMenu(self.submenu_columns)

        # Groups (as in the grouped plots tab)
        action = QtWidgets.QAction("Add replicate", self)
        action.triggered.connect(self.add_group)
        menu.addAction(action)
        action = QtWidgets.QAction("Rename replicates", self)
        action.triggered.connect(
            widgets_list.widget_results_groups.rename_groups)
        menu.addAction(action)

        # Batch change groups
        action = QtWidgets.QAction("Change selected replicates", self)
        menu.addAction(action)
        self.submenu_groups = QtWidgets.QMenu()
        for group in shared.exp.groups_list:
            name = group.name
            subaction = QtWidgets.QAction(name, self)
            self.submenu_groups.addAction(subaction)
        self.submenu_groups.triggered.connect(self.batch_groups_change)
        action.setMenu(self.submenu_groups)

        # Batch change slices
        if rt == "stiffness":
            action = QtWidgets.QAction("Change selected slices", self)
            menu.addAction(action)
            self.submenu_slices = QtWidgets.QMenu()
            subaction = QtWidgets.QAction("All", self)
            self.submenu_slices.addAction(subaction)
            for slice_nbr in range(shared.exp.global_max_indentation_index):
                subaction = QtWidgets.QAction(str(slice_nbr), self)
                self.submenu_slices.addAction(subaction)
            self.submenu_slices.triggered.connect(self.batch_slices_change)
            action.setMenu(self.submenu_slices)

        # Batch change rois (ROIs are None, 1, 2, ...)
        action = QtWidgets.QAction("Change selected ROI'S", self)
        menu.addAction(action)
        self.submenu_rois = QtWidgets.QMenu()
        length = 0
        for data in shared.exp.list:
            if len(data.roi_list) > length:
                length = len(data.roi_list)

        subaction = QtWidgets.QAction("None", self)
        self.submenu_rois.addAction(subaction)
        for roi in range(length):
            subaction = QtWidgets.QAction(str(roi + 1), self)
            self.submenu_rois.addAction(subaction)
        self.submenu_rois.triggered.connect(self.batch_rois_change)
        action.setMenu(self.submenu_rois)

        # Batch change slopes
        if rt == "events_forces" or rt == "events_per_curve" or \
                rt == "events_rupture_force" or rt == "loading_rates" or \
                rt == "events_distance":
            action = QtWidgets.QAction("Change selected slope filters", self)
            action.triggered.connect(lambda: self.change_filters("slope"))
            menu.addAction(action)
            action = QtWidgets.QAction("Change selected distance filters", self)
            action.triggered.connect(lambda: self.change_filters("dist_filt"))
            menu.addAction(action)

        # Batch change event distance (All or last)
        if rt == "events_distance":
            action = QtWidgets.QAction("Change selected event(s)", self)
            menu.addAction(action)
            self.submenu_events = QtWidgets.QMenu()
            subaction = QtWidgets.QAction("All", self)
            self.submenu_events.addAction(subaction)
            subaction = QtWidgets.QAction("Last", self)
            self.submenu_events.addAction(subaction)
            self.submenu_events.triggered.connect(self.batch_events_change)
            action.setMenu(self.submenu_events)

        # Display menu
        menu.popup(pos, menu.menuAction())
        menu.exec_()
