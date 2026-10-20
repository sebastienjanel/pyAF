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
"""Window of the dynamic force spectroscopy (Bell-Evans, DHS and
Friddle-Noy-De Yoreo fits, see dfs_tools)."""

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

COLUMNS = ["", "Model", "Points", "x [nm]", "koff [s⁻¹]",
           "\u0394G\u2021 / feq", "R²"]


def number(value, error=None):
    """Value (with its standard error) as text, 3 significant digits."""
    if value is None or not numpy.isfinite(value):
        return "–"
    text = "%.3g" % value
    if error is not None and numpy.isfinite(error):
        text += " ± %.2g" % error
    return text


class DFSWidget(PYAFWidget):
    """Fits of the forces vs the loading rates of the files (plots tab) or
    of the groups (grouped plots tab) displayed in the plot: Bell-Evans (on
    the median forces), DHS and Friddle-Noy-De Yoreo (also on the median
    forces),
    see dfs_tools.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_dfs")

        self.smallfont = QtWidgets.QApplication.font()
        self.smallfont.setPointSize(12)

        self.VL = QtWidgets.QVBoxLayout()

        self.HL = QtWidgets.QHBoxLayout()
        self.GRP_buttons = PYAFButtonGroup(self, "GRP_buttons")
        self.RBT_single = QtWidgets.QRadioButton("Files (plots tab)")
        self.RBT_groups = QtWidgets.QRadioButton("Speeds (replicates tab)")
        self.GRP_buttons.addButton(self.RBT_single, 0)
        self.GRP_buttons.addButton(self.RBT_groups, 1)
        self.HL.addWidget(self.RBT_single)
        self.HL.addWidget(self.RBT_groups)
        self.HL.addStretch(1)

        # Model of the fit
        self.HL_model = QtWidgets.QHBoxLayout()
        self.list_model = PYAFComboBox(self, "model")
        for name, _ in dfs_tools.MODELS.values():
            self.list_model.addItem(name)
        self.list_model.setToolTip(
            "Bell-Evans: F* = kBT/x ln(r x / (kBT koff)), linear in ln(r); "
            "to use unless the spectrum is curved.\nFitted on the median "
            "forces with their exact Bell-Evans expression.\n\nDHS "
            "(Dudko-Hummer-Szabo): curved spectrum (x decreases with the "
            "force), gives also the height\nof the barrier \u0394G\u2021; "
            "\u03bd = 1/2 for a cusp barrier, 2/3 for a linear-cubic one. "
            "Needs a wide range of loading rates.\n\nFriddle-Noy-De "
            "Yoreo: rebinding, plateau at the equilibrium force feq at low "
            "loading rates,\nthen Bell-Evans.\n\nAll the models are "
            "fitted on the median forces (the mean is pulled up by the "
            "multiple bonds\nand the nonspecific adhesions). DHS and "
            "Friddle-Noy-De Yoreo: 3 parameters, at least 4 points.")
        self.list_model.fit_to_contents()
        self.HL_model.addWidget(QtWidgets.QLabel("Model:"))
        self.HL_model.addWidget(self.list_model)
        self.HL_model.addStretch(1)

        # Method of the fit
        self.HL_method = QtWidgets.QHBoxLayout()
        self.list_method = PYAFComboBox(self, "method")
        self.list_method.addItem("Medians per file")
        self.list_method.addItem("Binned events")
        self.list_method.fit_to_contents()
        self.list_method.setToolTip(
            "Medians per file: one point per displayed file (or speed), at "
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
        for col, width in enumerate((30, 175, 50, 110, 130, 120, 60)):
            self.tableWidget.setColumnWidth(col, width)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.tableWidget.verticalHeader().hide()

        self.BT_get_dfs = PYAFButton(self, "get_dfs", "Fit")
        self.BT_remove = PYAFButton(self, "remove", "Remove selected")
        self.HL_buttons = QtWidgets.QHBoxLayout()
        self.HL_buttons.addWidget(self.BT_get_dfs)
        self.HL_buttons.addWidget(self.BT_remove)

        self.VL.addLayout(self.HL)
        self.VL.addLayout(self.HL_model)
        self.VL.addLayout(self.HL_method)
        self.VL.addWidget(self.LB_info)
        self.VL.addWidget(self.tableWidget)
        self.VL.addLayout(self.HL_buttons)
        self.setLayout(self.VL)

        self.update_widget()

    def input_updated(self, val):
        """Do nothing here (the number of bins is read when fitting)."""

    def list_updated(self, _):
        """Model or method of the fit changed."""
        self.update_widget()

    def model(self):
        """Selected model."""
        return list(dfs_tools.MODELS)[self.list_model.currentIndex()]

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
            result = dfs_tools.new_fit(mode, self.method(), nbr_bins,
                                       self.model())
            if isinstance(result, str):
                QtWidgets.QMessageBox.warning(self, "DFS", result)
                return
            dfs_tools.get_fits(mode).append(result)
            self.update_widget()
            self.update_plot()
            warning = self.fit_warning(result)
            if result.get("bell_limit"):
                QtWidgets.QMessageBox.information(
                    self, "DFS", self.bell_limit_text(result))
            elif warning:
                QtWidgets.QMessageBox.warning(self, "DFS", warning)

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

        # For the DFS, the replicates are the speeds
        what = "speed" if mode == "groups" else "file"
        self.list_method.setItemText(0, "Medians per " + what)
        self.list_method.fit_to_contents()
        self.IN_bins.setEnabled(self.method() == "binned")
        nbr = len(dfs_tools.get_items(mode))
        formulas = {
            "bell": "F* = (kBT/x) ln(r x / (kBT koff)), on the median forces",
            "dhs_cusp": "DHS, \u03bd = 1/2, on the median forces",
            "dhs_cubic": "DHS, \u03bd = 2/3, on the median forces",
            "fnd": "Friddle-Noy-De Yoreo, on the median forces"}
        self.LB_info.setText(
            "Fit " + formulas[self.model()] + ", of the " + what + "s "
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
                name = "Old fit (not drawn)"
            else:
                name = dfs_tools.MODELS[fit.get("model", "bell")][0]
                if fit["method"] == "binned":
                    name += ", " + str(fit["nbr_bins"]) + " bins"
                if fit.get("bell_limit"):
                    name += " = Bell-Evans"
                elif self.fit_warning(fit):
                    name = "\u26a0 " + name
            if "dG" in fit:
                third = "\u221e (no curvature)" if fit.get("bell_limit") \
                    else number(fit["dG"], fit["sdG"]) + " kBT"
            elif "feq" in fit:
                third = "0 (no plateau)" if fit.get("bell_limit") \
                    else number(fit["feq"], fit["sfeq"]) + " pN"
            else:
                third = "–"
            texts = [name, number(fit.get("nbr_points")),
                     number(fit["x"] * 1e9, fit["sx"] * 1e9),
                     number(fit["koff"], fit["skoff"]), third,
                     number(fit.get("r2"))]
            tooltip = self.tooltip(fit)
            for col, text in enumerate(texts, start=1):
                label = QtWidgets.QLabel(text)
                label.setFont(self.smallfont)
                label.setToolTip(tooltip)
                self.tableWidget.setCellWidget(
                    i, col, gui_tools.CenteredCellwidget(label))

    @staticmethod
    def bell_limit_text(fit):
        """Explanation of a fit that reduced to Bell-Evans."""
        if fit.get("model") == "fnd":
            return ("No plateau at low loading rates: feq = 0, the "
                    "Friddle-Noy-De Yoreo model is then Bell-Evans (same "
                    "curve).\nNo rebinding is seen in these data: "
                    "Bell-Evans is enough.")
        return ("No curvature: \u0394G\u2021 goes to its maximum, the DHS "
                "model is then Bell-Evans (same curve).\nThe spectrum is "
                "linear in ln(r): Bell-Evans is enough (a wider range of "
                "loading rates may show a curvature).")

    @staticmethod
    def fit_warning(fit):
        """Text if the fit does not describe the data (empty if fine)."""
        if fit.get("method") == "old":
            return ""
        if fit.get("model", "bell") == "bell" and fit.get("a", 1) <= 0:
            return ("The force does not increase with the loading rate: no "
                    "Bell-Evans parameters.")
        if fit.get("at_bound"):
            return ("Parameter at the limit of its range (" +
                    ", ".join(fit["at_bound"]) + "): the model does not "
                    "describe these points (the force should increase with "
                    "the loading rate,\nwith a curvature for DHS or a "
                    "plateau at low loading rates for Friddle-Noy-De "
                    "Yoreo). The curve is not reliable.")
        r2 = fit.get("r2", numpy.nan)
        if numpy.isfinite(r2) and r2 < 0.5:
            return ("The curve describes the points poorly (R\u00b2 = %.2g): "
                    "the model does not fit these data (the force should "
                    "increase with the loading rate).\nThe parameters are "
                    "not reliable." % r2)
        return ""

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
        model = fit.get("model", "bell")
        x = number(fit["x"] * 1e9, fit["sx"] * 1e9)
        koff = number(fit["koff"], fit["skoff"])
        if model == "bell":
            text = ("Bell-Evans: F* = kBT/x ln(r x / (kBT koff))\n"
                    "Linear fit F = a ln(r) + b (F in N, r in N/s) of the "
                    "median forces:\na = " + number(fit["a"], fit["sa"]) +
                    " N, b = " + number(fit["b"], fit["sb"]) + " N\n")
            if fit.get("median_exact"):
                text += ("x and koff fitted with the exact median of the "
                         "Bell-Evans distribution\n(F = kBT/x ln(1 + ln2 "
                         "r x / (kBT koff)); the line gives koff / ln 2)\n")
            elif fit.get("koff_corrected"):
                text += ("Fit of an older version: koff corrected for the "
                         "median forces (x ln 2)\n")
            text += "x = " + x + " nm\nkoff = " + koff + " s⁻¹\n"
        elif model in dfs_tools.NU:
            text = ("DHS, \u03bd = %.3g, fitted on the median forces "
                    "(exact median of the model):\n"
                    "x = " % fit["nu"] + x + " nm\nkoff = " + koff +
                    " s⁻¹\n\u0394G\u2021 = " +
                    number(fit["dG"], fit["sdG"]) + " kBT\n")
        else:
            text = ("Friddle-Noy-De Yoreo, fitted on the median forces:\n"
                    "x = " + x + " nm\nkoff = " + koff + " s⁻¹\n"
                    "feq = " + number(fit["feq"], fit["sfeq"]) + " pN\n"
                    "\u0394Gbu = feq\u00b2 / (2 kc) = " +
                    number(fit["dGbu"]) + " kBT, with the spring constant "
                    "of the cantilever\n(kc = " + number(fit["kc"]) +
                    " N/m; the linker makes the effective kc softer and "
                    "\u0394Gbu larger)\n")
        text += ("R² = " + number(fit["r2"]) + ", T = " + temp + "\n" +
                 str(fit["nbr_points"]) + " points, " +
                 str(fit["nbr_events"]) + " events")
        if fit["nbr_negative"]:
            text += "\n" + str(fit["nbr_negative"]) + \
                " events with a loading rate ≤ 0 not used (log)"
        if model == "bell" and fit["a"] <= 0:
            text += ("\nThe force does not increase with the loading rate "
                     "(a ≤ 0): no x and koff")
        if fit.get("bell_limit"):
            text += "\n" + DFSWidget.bell_limit_text(fit)
        if fit.get("at_bound"):
            text += ("\nParameter at the limit of its range (" +
                     ", ".join(fit["at_bound"]) + "): the data do not "
                     "determine it (curvature or plateau not visible?)")
        if not numpy.isfinite(fit["sx"]):
            text += "\nNo standard errors (too few points)"
        return text

    def checkbox_clicked(self, checkbox, row):
        """Hides or displays a fit."""
        if checkbox == "display_fit":
            row = int(row)
            value = self.tableWidget.cellWidget(row, 0).checkbox.isChecked()
            dfs_tools.get_fits(shared.exp.dfs_mode)[row]["display"] = value
            self.update_plot()
