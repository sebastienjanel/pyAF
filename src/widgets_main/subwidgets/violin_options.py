# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
# Copyright Sébastien Janel (2024- ) sebastien.janel@cnrs.fr
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

"""Options of the violin plots of the plots and grouped plots tabs."""

from qtpy import QtWidgets
from ... import consts
from ... import shared
from ... import widgets_list
from ... import experiment
from ...tools import misc_tools, palettes, stat_tools
from ...tools.gui_tools import PYAFButton, PYAFCheckBox, PYAFComboBox, \
    PYAFInput
from ...tools.PYAFWidget import PYAFWidget

INNER = ("Box", "Notched box", "Quartiles", "None")
POINTS = ("None", "Jitter", "Beeswarm")
DENSITY_NORM = ("Same width", "Same area", "Number of values")

ZEROS_TOOLTIP = (
    "Remove the values lower than or equal to 0 (zeros and negative\n"
    "values), for the plots, the statistics of the table and the exports.\n"
    "A 0 can be a real value (e.g. a curve without events) or a curve\n"
    "without result. For the loading rates, the events with a force <= 0\n"
    "are removed.")


class ViolinOptionsWidget(PYAFWidget):
    """Options of the violin plot of the plots tab (mode "single") or of the
    grouped plots tab (mode "groups").

    The parent has the table of the rows (tableWidget), the plot
    (MPL_canva), update_MPL and update_GUI.
    """

    def __init__(self, parent, name, mode="single"):
        super().__init__(parent, name)

        self.parent = parent
        self.mode = mode

        box = QtWidgets.QGroupBox("Violin plot")
        self.box = box
        form = QtWidgets.QFormLayout()
        form.setLabelAlignment(form.labelAlignment())

        self.list_inner = PYAFComboBox(self, "inner")
        for item in INNER:
            self.list_inner.addItem(item)
        self.list_inner.fit_to_contents()
        self.list_inner.setToolTip(
            "Box: median, quartiles (box) and whiskers (last values within\n"
            "1.5 IQR of the quartiles)\n"
            "Notched box: the notches show the 95% confidence interval of\n"
            "the median (two medians differ if their notches do not overlap)\n"
            "Quartiles: lines at the median and quartiles")
        form.addRow("Inside", self.list_inner)

        self.list_points = PYAFComboBox(self, "points")
        for item in POINTS:
            self.list_points.addItem(item)
        self.list_points.fit_to_contents()
        self.list_points.setToolTip(
            "Display every value on its violin.\n"
            "Jitter: random horizontal shifts (geom_jitter of ggplot2)\n"
            "Beeswarm: shifts following the density of the values, as wide\n"
            "as the violin (quasirandom beeswarm of ggbeeswarm)")
        form.addRow("Values", self.list_points)

        self.list_density_norm = PYAFComboBox(self, "density_norm")
        for item in DENSITY_NORM:
            self.list_density_norm.addItem(item)
        self.list_density_norm.fit_to_contents()
        self.list_density_norm.setToolTip(
            "Same width: all the violins have the same maximum width\n"
            "Same area: all the violins have the same area\n"
            "Number of values: the area is proportional to the number of "
            "values")
        form.addRow("Width", self.list_density_norm)

        self.list_palette = PYAFComboBox(self, "palette")
        for item in palettes.NAMES + [palettes.CUSTOM]:
            self.list_palette.addItem(item)
        self.list_palette.fit_to_contents()
        self.list_palette.setToolTip(
            "Palettes of ggplot2: its default (hue), viridis, Okabe-Ito\n"
            "(colorblind safe), the journals' palettes of ggsci and the\n"
            "ColorBrewer ones. Custom: the colors chosen row by row.")

        self.BT_color = PYAFButton(self, "rows_color", "Rows...", size=60)
        self.BT_color.setToolTip(
            "Choose the color of the rows selected in the table")
        colors = QtWidgets.QHBoxLayout()
        colors.setContentsMargins(0, 0, 0, 0)
        colors.addWidget(self.list_palette)
        colors.addWidget(self.BT_color)
        form.addRow("Colors", colors)

        self.IN_bw_adjust = PYAFInput(self, "bw_adjust", "", width=60)
        self.IN_bw_adjust.input.setValidator(misc_tools.validator("PF"))
        self.IN_bw_adjust.setToolTip(
            "Factor of the bandwidth of the density estimate (smaller: "
            "more details).\nAlso used for the modes of the table (maximum "
            "of the density).")
        form.addRow("Smoothing", self.IN_bw_adjust)

        self.CB_cut = PYAFCheckBox(self, "cut", "Truncate at data range")
        self.CB_cut.setToolTip(
            "Draw the density estimate (KDE) only between the smallest and\n"
            "the largest value. Unticked, it extends 2 bandwidths beyond\n"
            "them, where there is no data.")

        self.CB_log = PYAFCheckBox(self, "log_scale", "Log scale")
        self.CB_log.setToolTip("Log scale (the values <= 0 are not shown)")

        self.CB_counts = PYAFCheckBox(self, "show_counts", "Number of values")
        self.CB_counts.setToolTip(
            "Display the number of values under the names")

        self.CB_zeros = PYAFCheckBox(self, "zeros", "Remove \u2264 0")
        self.CB_zeros.setToolTip(ZEROS_TOOLTIP)

        self.CB_y_manual = PYAFCheckBox(self, "y_manual", "Manual Y limits")
        self.CB_white = PYAFCheckBox(self, "white", "White background")
        self.CB_white.setToolTip(
            "White background and black axes, for the publications (also in\n"
            "the window of Figures > Open results)")

        # Checkboxes on two columns
        checkboxes = QtWidgets.QGridLayout()
        checkboxes.setContentsMargins(0, 0, 0, 0)
        for i, checkbox in enumerate((self.CB_cut, self.CB_log,
                                      self.CB_counts, self.CB_zeros,
                                      self.CB_white, self.CB_y_manual)):
            checkboxes.addWidget(checkbox, i // 2, i % 2)
        form.addRow(checkboxes)

        self.IN_y_min = PYAFInput(self, "y_min", "", width=70)
        self.IN_y_min.input.setValidator(misc_tools.validator("F"))
        self.IN_y_max = PYAFInput(self, "y_max", "", width=70)
        self.IN_y_max.input.setValidator(misc_tools.validator("F"))
        y_range = QtWidgets.QHBoxLayout()
        y_range.setContentsMargins(0, 0, 0, 0)
        y_range.addWidget(self.IN_y_min)
        y_range.addWidget(QtWidgets.QLabel("to"))
        y_range.addWidget(self.IN_y_max)
        y_range.addStretch()
        form.addRow("Y range", y_range)

        box.setLayout(form)

        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(box)
        self.setLayout(layout)

        self.update_widget()

    def update_widget(self):
        """Display the current options."""
        options = self.get_options()
        self.list_inner.setCurrentIndex(INNER.index(options["inner"]))
        self.list_points.setCurrentIndex(POINTS.index(options["points"]))
        self.list_density_norm.setCurrentIndex(
            DENSITY_NORM.index(options["density_norm"]))
        self.IN_bw_adjust.changeValue(misc_tools.format_number(
            options["bw_adjust"]))
        self.list_palette.setCurrentIndex(
            self.list_palette.findText(options["palette"]))
        self.CB_cut.setChecked(options["cut"])
        self.CB_log.setChecked(options["log_scale"])
        self.CB_counts.setChecked(options["show_counts"])
        self.CB_white.setChecked(options["background"] == "White")
        self.CB_zeros.setChecked(shared.exp.hist_remove_zeros)
        manual = options["y_mode"] == "manual"
        self.CB_y_manual.setChecked(manual)
        self.IN_y_min.changeValue(misc_tools.format_number(options["y_min"]))
        self.IN_y_max.changeValue(misc_tools.format_number(options["y_max"]))
        self.IN_y_min.setEnabled(manual)
        self.IN_y_max.setEnabled(manual)

        # The loading rates are plotted as forces vs loading rates (no
        # violins): only the colors, the background and the Y limits apply
        violins = shared.exp.results_type != "loading_rates"
        self.box.setTitle("Violin plot" if violins
                          else "Force vs loading rate plot")
        for widget in (self.list_inner, self.list_points,
                       self.list_density_norm, self.IN_bw_adjust,
                       self.CB_cut, self.CB_log, self.CB_counts):
            widget.setEnabled(violins)

    def list_updated(self, name):
        """Called when a list is changed."""
        options = self.get_options()
        if name == "inner":
            options["inner"] = INNER[self.list_inner.currentIndex()]
        elif name == "points":
            options["points"] = POINTS[self.list_points.currentIndex()]
        elif name == "density_norm":
            options["density_norm"] = \
                DENSITY_NORM[self.list_density_norm.currentIndex()]
        elif name == "palette":
            options["palette"] = self.list_palette.currentText()
            experiment.apply_palette(shared.exp, self.mode)
        self.replot()

    def get_options(self):
        """Options of the violin plot."""
        return experiment.get_violin_options(shared.exp, self.mode)

    def button_clicked(self, name):
        """Called when a button is clicked."""
        if name == "rows_color":
            self.change_rows_color()

    def change_rows_color(self):
        """Ask for a color for the rows selected in the table."""
        table = self.parent.tableWidget
        rows = sorted(index.row()
                      for index in table.selectionModel().selectedRows())
        if not rows:
            QtWidgets.QMessageBox.information(
                self, "Color of the rows",
                "Select the rows in the table first (click on a row, Shift "
                "or Cmd/Ctrl + click for several rows).")
            return
        if consts.UNIT_TESTING:
            return
        items = experiment.get_plotted_items(shared.exp, self.mode)
        color = misc_tools.ask_user_for_color(items[rows[0]].color)
        if color is False:
            return
        for row in rows:
            items[row].color = color
        # The colors do not follow a palette anymore
        self.get_options()["palette"] = palettes.CUSTOM
        self.update_widget()
        self.replot()

    def checkbox_clicked(self, name):
        """Called when a checkbox is clicked."""
        options = self.get_options()
        if name == "cut":
            options["cut"] = self.CB_cut.isChecked()
        elif name == "log_scale":
            options["log_scale"] = self.CB_log.isChecked()
        elif name == "show_counts":
            options["show_counts"] = self.CB_counts.isChecked()
        elif name == "white":
            options["background"] = \
                "White" if self.CB_white.isChecked() else "Theme"
        elif name == "y_manual":
            manual = self.CB_y_manual.isChecked()
            options["y_mode"] = "manual" if manual else "auto"
            if manual:
                # Start from the current limits of the plot
                ymin, ymax = self.parent.MPL_canva.canvas.axes.get_ylim()
                options["y_min"], options["y_max"] = ymin, ymax
            self.update_widget()
        elif name == "zeros":
            # Changes the values: recompute the table and the plot
            shared.exp.hist_remove_zeros = self.CB_zeros.isChecked()
            widgets_list.widget_results_single.button_clicked("refresh_hist")
            return
        self.replot()

    def input_updated(self, name):
        """Called when an input is changed."""
        options = self.get_options()
        if name == "bw_adjust":
            options["bw_adjust"] = self.IN_bw_adjust.get_float_value()
            # The modes of the table use the same density estimate
            if self.mode == "groups":
                stat_tools.update_groups_modes()
            else:
                stat_tools.update_single_modes()
            self.parent.update_GUI("table_labels")
        elif name == "y_min":
            options["y_min"] = self.IN_y_min.get_float_value()
        elif name == "y_max":
            options["y_max"] = self.IN_y_max.get_float_value()
        self.replot()

    def replot(self):
        """Redraw the violins with the new options."""
        self.parent.update_MPL()
