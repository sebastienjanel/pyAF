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
"""SuperPlots tab: the conditions, with one dot per file (force map, QI
image or single curve), the replicate means and the tests on the replicates
(Lord et al., J Cell Biol 2020)."""

import csv
import logging
import os

import numpy
from qtpy import QtCore, QtWidgets

from .. import shared
from .. import widgets_list
from ..experiment import ConditionGroup
from ..plots.PYAFPlot import PYAFPlot
from ..tools import gui_tools
from ..tools import misc_tools
from ..tools import superplot_tools
from ..tools.gui_tools import PYAFButton, PYAFCheckBox, PYAFComboBox
from ..tools.gui_tools import PYAFTableWidget
from ..tools.PYAFWidget import PYAFWidget
from ..widgets_main.results_single import PLOT_MIN_HEIGHT
from ..widgets_main.subwidgets.result_chooser import ResultsChooserWidget
from ..widgets_small.rename_conditions import RenameConditionsWidget

CONDITION_COLUMNS = ("", "Condition", "Replicates", "Files", "Curves",
                     "Mean", "S.D.", "95% CI")
TEST_COLUMNS = ("Comparison", "Test", "N", "p-value", "",
                "Fold change [95% CI]")


def number(value):
    """Value as text (4 significant digits)."""
    if value is None or not numpy.isfinite(value):
        return "–"
    return "%.4g" % value


def fold_change_text(test):
    """Fold change and its 95 % confidence interval as text."""
    if not numpy.isfinite(test["fc"]):
        return "–"
    if not numpy.isfinite(test["fc_low"]):
        return "%.3g [unbounded]" % test["fc"]
    return "%.3g [%.3g, %.3g]" % (test["fc"], test["fc_low"],
                                  test["fc_high"])


class ResultsExperimentWidget(PYAFWidget):
    """SuperPlots tab.

    Hierarchy: condition > replicate (the groups of the plots tab) > file
    (force map, QI image or single curve) > curves. Each file gives one value
    (median of its curves), each replicate the mean of its files, and the
    tests compare the replicates of
    each condition with the ones of the control (n = number of replicates).
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_results_experiment")

        self.settings = QtCore.QSettings()
        self.logger = logging.getLogger()
        self.parent = parent
        self.canvas_resolution = parent.canvas_resolution
        self.canvas_size = [780, 350]
        sizes = [self.canvas_size[0] / self.canvas_resolution,
                 self.canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]

        # Plot, filling the available space
        self.canva = QtWidgets.QWidget()
        self.MPL_canva = PYAFPlot(self, "results_experiment", self.canva,
                                  sizes)
        self.SA_plot = QtWidgets.QScrollArea()
        self.SA_plot.setWidget(self.MPL_canva)
        self.SA_plot.setWidgetResizable(False)
        self.SA_plot.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarAlwaysOff)
        self.SA_plot.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.SA_plot.setMinimumHeight(PLOT_MIN_HEIGHT)
        self.SA_plot.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        self.plot_resize_timer = QtCore.QTimer(self)
        self.plot_resize_timer.setSingleShot(True)
        self.plot_resize_timer.timeout.connect(self.update_MPL)
        self.SA_plot.viewport().installEventFilter(self)

        # Options
        box = QtWidgets.QGroupBox("SuperPlot")
        form = QtWidgets.QFormLayout()
        self.list_control = PYAFComboBox(self, "control")
        self.list_control.setToolTip(
            "The other conditions are compared to this one.\nNone: all the "
            "pairs of conditions are compared (Tukey's test, or paired "
            "t-tests with the Holm correction)")
        self.list_paired = PYAFComboBox(self, "paired")
        self.list_paired.addItem("Unpaired")
        self.list_paired.addItem("Paired (same exp.)")
        self.list_paired.setToolTip(
            "Unpaired: Welch t-test (one comparison), Dunnett's test "
            "(several) or Tukey's test (all pairs).\nPaired: the "
            "replicates with the same Experiment number (Replicates tab)\n"
            "are compared with paired t-tests (Holm correction for "
            "several).")
        self.list_error = PYAFComboBox(self, "error")
        self.list_error.addItem("SD")
        self.list_error.addItem("95% CI")
        self.list_error.setToolTip(
            "Error bars of the mean of the replicates")
        self.list_paired.fit_to_contents()
        self.list_error.fit_to_contents()
        self.CB_log = PYAFCheckBox(self, "log", "Log scale (plot and tests)")
        self.CB_curves = PYAFCheckBox(self, "curves", "All the curves")
        self.CB_curves.setToolTip(
            "Violin of all the curves of the condition in the background")
        self.CB_white = PYAFCheckBox(self, "white", "White background")
        form.addRow("Control", self.list_control)
        form.addRow("Replicates", self.list_paired)
        form.addRow("Error bars", self.list_error)
        form.addRow(self.CB_log)
        form.addRow(self.CB_curves)
        form.addRow(self.CB_white)
        box.setLayout(form)
        self.widget_options = box

        # Results type and buttons
        self.HL_buttons = QtWidgets.QHBoxLayout()
        name = "widget_results_chooser_experiment"
        self.widget_result_chooser = ResultsChooserWidget(self, name)
        self.BT_export = PYAFButton(self, "export", "Export CSV")
        self.BT_export.setToolTip(
            "One row per file: condition, replicate, experiment, file value, "
            "number of curves, replicate mean.\nAnd the tests: p-values and "
            "fold changes with their 95 % CI (_tests.csv).")
        self.HL_buttons.addWidget(self.widget_result_chooser)
        self.HL_buttons.addStretch(1)
        self.HL_buttons.addWidget(self.BT_export)

        # Conditions and tests
        self.TW_conditions = self.make_table(CONDITION_COLUMNS)
        superplot_tools.set_header_tooltips(self.TW_conditions, "superplots")
        self.TW_conditions.setColumnWidth(0, 35)
        header = self.TW_conditions.horizontalHeader()
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        # Widths of the columns after Condition (it takes the remaining
        # space, at least 150 px: minimum width of the table)
        widths = (80, 50, 70, 80, 80, 80)
        for column, width in enumerate(widths, start=2):
            self.TW_conditions.setColumnWidth(column, width)
        self.TW_conditions.setMinimumWidth(35 + sum(widths) + 150 + 20)
        header.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        header.customContextMenuRequested.connect(self.popUpMenuTableHeader)
        self.TW_tests = self.make_table(TEST_COLUMNS)
        self.TW_tests.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.Stretch)
        # Widths of the columns after Comparison (it takes the remaining
        # space, at least 120 px: minimum width of the table)
        widths = (130, 40, 80, 40, 145)
        for column, width in enumerate(widths, start=1):
            self.TW_tests.setColumnWidth(column, width)
        self.TW_tests.setMinimumWidth(sum(widths) + 120 + 20)
        self.TW_tests.horizontalHeaderItem(5).setToolTip(
            "Ratio of the means of the replicates vs the control (geometric "
            "means with the log scale), with its 95 % confidence interval\n"
            "from the same test as the p-value: the interval excludes 1 when "
            "the difference is significant.\nWith few replicates, the "
            "interval shows how large the effect can be, which the p-value "
            "alone does not.")
        self.HL_tables = QtWidgets.QHBoxLayout()
        self.HL_tables.addWidget(self.TW_conditions, 1)
        self.HL_tables.addWidget(self.TW_tests, 1)

        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.SA_plot, 0, 0)
        layout.addWidget(self.widget_options, 0, 1, QtCore.Qt.AlignTop)
        layout.addLayout(self.HL_buttons, 1, 0, 1, 2)
        layout.addLayout(self.HL_tables, 2, 0, 1, 2)
        layout.setRowStretch(0, 3)
        layout.setRowStretch(2, 2)
        layout.setColumnStretch(0, 1)
        self.setLayout(layout)

        widgets_list.widget_progressbar.update()
        self.update_widget()

    @staticmethod
    def make_table(columns):
        """Table with the given columns."""
        table = PYAFTableWidget(256, 50, 0, len(columns))
        table.setHorizontalHeaderLabels(columns)
        table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        table.horizontalHeader().setHighlightSections(False)
        table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        table.verticalHeader().hide()
        pol = QtWidgets.QSizePolicy.MinimumExpanding
        table.setSizePolicy(QtWidgets.QSizePolicy(pol, pol))
        return table

    def eventFilter(self, obj, event):
        """Redraw the plot when the visible area is resized."""
        if obj is self.SA_plot.viewport() and \
                event.type() == QtCore.QEvent.Resize:
            self.plot_resize_timer.start(150)
        return super().eventFilter(obj, event)

    def update_MPL(self):
        """Update the plot (at the size of the visible area)."""
        viewport = self.SA_plot.viewport().size()
        size = QtCore.QSize(viewport.width(),
                            max(viewport.height(), PLOT_MIN_HEIGHT))
        if self.MPL_canva.size() != size:
            self.MPL_canva.resize(size)
        self.MPL_canva.update_plot()

    def update_widget(self):
        """Update the options, the tables and the plot."""
        widgets_list.widget_results_chooser_experiment.update_widget()
        options = superplot_tools.get_options()

        # Control: one of the conditions, or none (all the pairs)
        names = [c.name for c in shared.exp.conditions_list[1:]]
        ids = [c.condition_id for c in shared.exp.conditions_list[1:]]
        if ids:
            names.append("None: all pairs")
            ids.append(superplot_tools.ALL_PAIRS)
        if options["control"] not in ids and ids:
            options["control"] = ids[0]
        self.list_control.blockSignals(True)
        self.list_control.clear()
        for name in names:
            self.list_control.addItem(name)
        if ids:
            self.list_control.setCurrentIndex(ids.index(options["control"]))
        self.list_control.fit_to_contents()
        self.list_control.blockSignals(False)
        self.list_paired.setCurrentIndex(1 if options["paired"] else 0)
        self.list_error.setCurrentIndex(
            1 if options["error"] == "95% CI" else 0)
        self.CB_log.setChecked(options["log"])
        self.CB_curves.setChecked(options["curves"])
        self.CB_white.setChecked(options["background"] == "White")

        self.update_tables()
        self.update_MPL()

    def update_tables(self):
        """Values of the conditions and results of the tests."""
        options = superplot_tools.get_options()
        data = {c["id"]: c for c in superplot_tools.get_data()}

        table = self.TW_conditions
        table.setRowCount(0)
        for row, condition in enumerate(shared.exp.conditions_list[1:]):
            table.insertRow(row)
            table.setRowHeight(row, 30)
            cb = gui_tools.CenteredCellCheckbox(self, "display", row)
            cb.checkbox.setChecked(condition.display)
            table.setCellWidget(row, 0, cb)
            texts = [condition.name] + ["–"] * 6
            item = data.get(condition.condition_id)
            if item is not None:
                s = superplot_tools.summary(item)
                texts = [condition.name, str(s["n"]), str(s["maps"]),
                         str(s["curves"]), number(s["mean"]),
                         number(s["sd"]), number(s["ci"])]
            for col, text in enumerate(texts, start=1):
                label = QtWidgets.QLabel(text)
                label.setFont(self.parent.smallfont)
                table.setCellWidget(row, col,
                                    gui_tools.CenteredCellwidget(label))

        table = self.TW_tests
        table.setRowCount(0)
        tests = superplot_tools.compare(list(data.values()), options)
        for row, test in enumerate(tests):
            table.insertRow(row)
            table.setRowHeight(row, 30)
            texts = [test["name"], test["test"], str(test["n"]),
                     number(test["p"]), superplot_tools.stars(test["p"]),
                     fold_change_text(test)]
            for col, text in enumerate(texts):
                label = QtWidgets.QLabel(text)
                label.setFont(self.parent.smallfont)
                note = test["fc_note"] if col == 5 else test["note"]
                if note:
                    label.setToolTip(note)
                table.setCellWidget(row, col,
                                    gui_tools.CenteredCellwidget(label))

    def list_updated(self, name):
        """An option changed."""
        options = superplot_tools.get_options()
        if name == "control":
            index = self.list_control.currentIndex()
            conditions = shared.exp.conditions_list[1:]
            if 0 <= index < len(conditions):
                options["control"] = conditions[index].condition_id
            elif index == len(conditions):
                options["control"] = superplot_tools.ALL_PAIRS
        elif name == "paired":
            options["paired"] = self.list_paired.currentIndex() == 1
        elif name == "error":
            options["error"] = self.list_error.currentText()
        self.update_tables()
        self.update_MPL()

    def checkbox_clicked(self, name, row=None):
        """An option or the display of a condition changed."""
        options = superplot_tools.get_options()
        if name == "display":
            row = int(row)
            value = self.TW_conditions.cellWidget(row, 0).checkbox.isChecked()
            shared.exp.conditions_list[row + 1].display = value
        elif name == "log":
            options["log"] = self.CB_log.isChecked()
        elif name == "curves":
            options["curves"] = self.CB_curves.isChecked()
        elif name == "white":
            options["background"] = \
                "White" if self.CB_white.isChecked() else "Default"
        self.update_tables()
        self.update_MPL()

    def button_clicked(self, button):
        """Export of the values."""
        if button == "export":
            self.export_csv()

    def export_csv(self):
        """Saves the values of the files and of the replicates."""
        rows = superplot_tools.export_rows(superplot_tools.get_data())
        if not rows:
            return
        path = misc_tools.get_user_path(self.settings)
        name = "superplot_" + str(shared.exp.results_type) + ".csv"
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export the SuperPlot values", os.path.join(path, name),
            "CSV (*.csv)")
        if not filename:
            return
        if not filename.lower().endswith(".csv"):
            filename += ".csv"
        misc_tools.set_user_path(self.settings, os.path.dirname(filename))
        with open(filename, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        text = "Values exported to " + filename
        # Tests: p-values and fold changes with their confidence intervals
        tests = superplot_tools.test_rows(superplot_tools.get_data(),
                                          superplot_tools.get_options())
        if tests:
            tests_name = filename[:-4] + "_tests.csv"
            with open(tests_name, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(tests[0]))
                writer.writeheader()
                writer.writerows(tests)
            text += "\nTests exported to " + tests_name
        gui_tools.show_success(self, "Export complete", text)

    def add_condition(self):
        """Add a new condition."""
        val = len(shared.exp.conditions_list)
        shared.exp.conditions_list.append(ConditionGroup(val))
        pos = (val - 1) % len(shared.colors_list)
        shared.exp.conditions_list[-1].color = shared.colors_list[pos]
        widgets_list.widget_results_groups.reset_table()
        self.update_widget()

    def rename_conditions(self):
        """Widget with inputs to rename the conditions."""
        if widgets_list.widget_rename_conditions is None:
            RenameConditionsWidget(self)
            widgets_list.widget_rename_conditions.setWindowTitle(
                "Rename Conditions")
            widgets_list.widget_rename_conditions.show()
        widgets_list.widget_rename_conditions.activateWindow()
        widgets_list.widget_rename_conditions.raise_()

    def popUpMenuTableHeader(self, pos):
        """Menu of the header of the conditions table."""
        pos = self.TW_conditions.mapToGlobal(pos)
        menu = QtWidgets.QMenu()
        action = QtWidgets.QAction("Add condition", self)
        action.triggered.connect(self.add_condition)
        menu.addAction(action)
        action = QtWidgets.QAction("Rename conditions", self)
        action.triggered.connect(self.rename_conditions)
        menu.addAction(action)
        menu.exec(pos)
