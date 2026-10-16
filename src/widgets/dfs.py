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

"""Display a list with results of DFS experiments."""
"""Window of the dynamic force spectroscopy (Bell-Evans fits)."""

import numpy
from qtpy import QtCore, QtWidgets
from ..tools import gui_tools
from ..tools import misc_tools
from ..tools import dfs_tools
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFButtonGroup
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFInput
from ..tools.PYAFWidget import PYAFWidget
from .. import widgets_list
from .. import shared

COLUMNS = ["", "Method", "Points", "x [nm]", "koff [s⁻¹]", "R²"]


def number(value, error=None):
    """Value (with its standard error) as text, 3 significant digits."""
    if value is None or not numpy.isfinite(value):
        return "–"
    text = "%.3g" % value
    if error is not None and numpy.isfinite(error):
        text += " ± %.2g" % error
    return text


class DFSWidget(PYAFWidget):
    """Bell-Evans fits of the forces vs the loading rates of the files (plots
    tab) or of the groups (grouped plots tab) displayed in the plot.

    The fit F = a ln(r) + b gives the distance to the barrier x = kBT / a and
    the dissociation rate koff = exp(-b / a) / a.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_dfs")

        self.smallfont = QtWidgets.QApplication.font()
        self.smallfont.setPointSize(12)

        self.VL = QtWidgets.QVBoxLayout()

        self.HL = QtWidgets.QHBoxLayout()
        self.GRP_buttons = PYAFButtonGroup(self, "GRP_buttons")
        self.RBT_single = QtWidgets.QRadioButton("Files (plots tab)")
        self.RBT_groups = QtWidgets.QRadioButton("Groups (grouped plots)")
        self.GRP_buttons.addButton(self.RBT_single, 0)
        self.GRP_buttons.addButton(self.RBT_groups, 1)
        self.HL.addWidget(self.RBT_single)
        self.HL.addWidget(self.RBT_groups)
        self.HL.addStretch(1)

        # Method of the fit
        self.HL_method = QtWidgets.QHBoxLayout()
        self.list_method = PYAFComboBox(self, "method")
        self.list_method.addItem("Medians per file")
        self.list_method.addItem("Binned events")
        self.list_method.setToolTip(
            "Medians per file: one point per displayed file (or group), at "
            "the median loading rate\nand the median force of its events "
            "(several files at different speeds).\n\nBinned events: the "
            "events of the displayed files are put in bins of loading rate\n"
            "with the same number of events; one point per bin, at its "
            "median loading rate and median force\n(uses the spread of the "
            "loading rates of a single speed).")
        self.IN_bins = PYAFInput(self, "bins", "bins", width=40)
        self.IN_bins.input.setValidator(
            misc_tools.validator("RI", bottom=3, top=100))
        self.IN_bins.changeValue(10)
        self.HL_method.addWidget(QtWidgets.QLabel("Fit points:"))
        self.HL_method.addWidget(self.list_method)
        self.HL_method.addWidget(self.IN_bins)
        self.HL_method.addStretch(1)

        self.LB_info = QtWidgets.QLabel()
        self.LB_info.setWordWrap(True)

        self.tableWidget = QtWidgets.QTableWidget(0, len(COLUMNS))
        behavior = QtWidgets.QAbstractItemView.SelectRows
        self.tableWidget.setSelectionBehavior(behavior)
        self.tableWidget.horizontalHeader().setHighlightSections(False)
        self.tableWidget.setHorizontalHeaderLabels(COLUMNS)
        for col, width in enumerate((30, 130, 60, 110, 130, 60)):
            self.tableWidget.setColumnWidth(col, width)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.tableWidget.verticalHeader().hide()

        self.BT_get_dfs = PYAFButton(self, "get_dfs", "Fit")
        self.BT_remove = PYAFButton(self, "remove", "Remove selected")
        self.HL_buttons = QtWidgets.QHBoxLayout()
        self.HL_buttons.addWidget(self.BT_get_dfs)
        self.HL_buttons.addWidget(self.BT_remove)

        self.VL.addLayout(self.HL)
        self.VL.addLayout(self.HL_method)
        self.VL.addWidget(self.LB_info)
        self.VL.addWidget(self.tableWidget)
        self.VL.addLayout(self.HL_buttons)
        self.setLayout(self.VL)

        self.update_widget()

    def input_updated(self, val):
        """Do nothing here (the number of bins is read when fitting)."""

    def list_updated(self, _):
        """Method of the fit changed."""
        self.update_widget()

    def method(self):
        """Selected method."""
        return "binned" if self.list_method.currentIndex() == 1 \
            else "medians"

    def update_plot(self):
        """Redraw the plot of the current tab."""
        if shared.exp.dfs_mode == "single":
            widgets_list.widget_results_single.update_MPL()
        else:
            widgets_list.widget_results_groups.update_MPL()

    def button_clicked(self, button):
        """Called when a button is clicked."""
        mode = shared.exp.dfs_mode

        if button == "get_dfs":
            if shared.exp.results_type != "loading_rates":
                QtWidgets.QMessageBox.warning(
                    self, "DFS", "Display the loading rates in the plots "
                    "tab to fit them.")
                return
            nbr_bins = self.IN_bins.get_int_value() or 10
            result = dfs_tools.new_fit(mode, self.method(), nbr_bins)
            if isinstance(result, str):
                QtWidgets.QMessageBox.warning(self, "DFS", result)
                return
            dfs_tools.get_fits(mode).append(result)
            self.update_widget()
            self.update_plot()

        elif button == "GRP_buttons":
            if self.GRP_buttons.checkedId() == 0:
                shared.exp.dfs_mode = "single"
            elif self.GRP_buttons.checkedId() == 1:
                shared.exp.dfs_mode = "groups"
            self.update_widget()

        elif button == "remove":
            rows = sorted({item.row() for item in
                           self.tableWidget.selectionModel().selectedRows()},
                          reverse=True)
            fits = dfs_tools.get_fits(mode)
            for row in rows:
                del fits[row]
                olds = getattr(shared.exp, "list_dfs_" + mode + "_results")
                if row < len(olds):
                    del olds[row]
            self.update_widget()
            self.update_plot()

    def update_widget(self):
        """Update the GUI."""
        mode = shared.exp.dfs_mode
        if mode == "single":
            self.RBT_single.setChecked(True)
        else:
            self.RBT_groups.setChecked(True)

        what = "group" if mode == "groups" else "file"
        self.list_method.setItemText(0, "Medians per " + what)
        self.IN_bins.setEnabled(self.method() == "binned")
        nbr = len(dfs_tools.get_items(mode))
        self.LB_info.setText(
            "Fit F = (kBT/x) ln(r x / (kBT koff)) on the " + what + "s "
            "displayed in the plot (" + str(nbr) + " " + what +
            ("s" if nbr != 1 else "") + ").")

        fits = dfs_tools.get_fits(mode)
        self.tableWidget.setRowCount(0)
        for i, fit in enumerate(fits):
            self.tableWidget.insertRow(i)
            self.tableWidget.setRowHeight(i, 25)

            cb = gui_tools.CenteredCellCheckbox(self, "display_fit", i)
            self.tableWidget.setCellWidget(i, 0, cb)
            cb.checkbox.setChecked(bool(fit.get("display")))
            cb.setEnabled(fit["method"] != "old")

            if fit["method"] == "old":
                method = "Old fit (not drawn)"
            elif fit["method"] == "binned":
                method = "Binned (" + str(fit["nbr_bins"]) + " bins)"
            else:
                method = "Medians"
            texts = [method, number(fit.get("nbr_points")),
                     number(fit["x"] * 1e9, fit["sx"] * 1e9),
                     number(fit["koff"], fit["skoff"]),
                     number(fit.get("r2"))]
            tooltip = self.tooltip(fit)
            for col, text in enumerate(texts, start=1):
                label = QtWidgets.QLabel(text)
                label.setFont(self.smallfont)
                label.setToolTip(tooltip)
                self.tableWidget.setCellWidget(
                    i, col, gui_tools.CenteredCellwidget(label))

    @staticmethod
    def tooltip(fit):
        """Details of a fit."""
        if fit["method"] == "old":
            return ("Fit of an older version of pyAF (its forces and loading "
                    "rates were sorted separately: the values are not "
                    "reliable).")
        temps = fit["temperatures"]
        temp = number(fit["temperature"] - 273.15) + " °C"
        if temps[1] - temps[0] > 0.01:
            temp += " (mean of " + number(temps[0] - 273.15) + " to " + \
                number(temps[1] - 273.15) + " °C)"
        text = ("F = a ln(r) + b (F in N, r in N/s)\n"
                "a = " + number(fit["a"], fit["sa"]) + " N\n"
                "b = " + number(fit["b"], fit["sb"]) + " N\n"
                "x = kBT / a = " + number(fit["x"] * 1e9, fit["sx"] * 1e9) +
                " nm\nkoff = exp(-b/a) / a = " +
                number(fit["koff"], fit["skoff"]) + " s⁻¹\n"
                "R² = " + number(fit["r2"]) + ", T = " + temp + "\n" +
                str(fit["nbr_points"]) + " points, " +
                str(fit["nbr_events"]) + " events")
        if fit["nbr_negative"]:
            text += "\n" + str(fit["nbr_negative"]) + \
                " events with a loading rate ≤ 0 not used (log)"
        if fit["a"] <= 0:
            text += ("\nThe force does not increase with the loading rate "
                     "(a ≤ 0): no x and koff")
        if fit["nbr_points"] < 3:
            text += "\nNo standard errors with 2 points"
        return text

    def checkbox_clicked(self, checkbox, row):
        """Hides or displays a fit."""
        if checkbox == "display_fit":
            row = int(row)
            value = self.tableWidget.cellWidget(row, 0).checkbox.isChecked()
            dfs_tools.get_fits(shared.exp.dfs_mode)[row]["display"] = value
            self.update_plot()
