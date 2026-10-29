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

"""Plots profiles."""

from qtpy import QtCore, QtGui, QtWidgets
from ..plots.PYAFPlot import PYAFPlot
from ..tools.gui_tools import PYAFInput, CenteredCellwidget
import numpy
from ..tools import planefit
from ..tools.PYAFWidget import PYAFWidget
from .. import shared


class ProfileWidget(PYAFWidget):
    """The profile widget displays profiles and a tool to measure distances."""

    def __init__(self, parent):
        super().__init__(parent, "widget_profiles")

        self.checkbox_list = []
        self.delta_x_labels = []
        self.delta_z_labels = []
        self.data = shared.exp.current_data
        self.profile_list = shared.exp.current_data.profile_list

        self.canvas_resolution = self.parent.canvas_resolution
        self.canvas_size = [720, 288]

        # Profiles of the values of the map they are drawn on (results tab),
        # computed once here (x and value deltas of the table) and used by
        # the plot
        self.compute_profiles()

        self.VL = QtWidgets.QVBoxLayout()

        # Box with canvas -----------------------------------------------------
        self.box_profiles_plot = QtWidgets.QGroupBox("Profiles")
        self.box_profiles_plot_VL = QtWidgets.QVBoxLayout()
        self.box_profiles_plot_HL = QtWidgets.QHBoxLayout()

        self.W_canvas = QtWidgets.QWidget()
        self.W_canvas.setFixedSize(self.canvas_size[0], self.canvas_size[1])
        sizes = [self.canvas_size[0] / self.canvas_resolution,
                 self.canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]
        self.canvas = PYAFPlot(self, "profiles", self.W_canvas, sizes)
        self.canvas.mpl_connect("button_press_event", self.canvas_press)

        # Bars for calculation
        self.HL_bars = QtWidgets.QHBoxLayout()
        self.input_bar_1 = PYAFInput(
            self, "input_bar_1", ",")
        self.input_bar_2 = PYAFInput(self, "input_bar_2", "Bar positions [nm]")
        self.input_bar_1.changeValue(shared.exp.profiles_bar_1)
        self.input_bar_2.changeValue(shared.exp.profiles_bar_2)
        self.checkbox_show_profiles_bars = QtWidgets.QCheckBox("Show bars")
        self.checkbox_show_profiles_bars.stateChanged.connect(
            lambda: self.checkbox_clicked("checkbox_show_profiles_bars"))

        self.HL_bars.addWidget(self.input_bar_1)
        self.HL_bars.addWidget(self.input_bar_2)
        self.HL_bars.addWidget(self.checkbox_show_profiles_bars)
        self.HL_bars.addStretch(1)

        self.box_profiles_plot_HL.addWidget(self.W_canvas)
        self.box_profiles_plot_VL.addLayout(self.box_profiles_plot_HL)
        self.box_profiles_plot_VL.addLayout(self.HL_bars)
        self.box_profiles_plot.setLayout(self.box_profiles_plot_VL)

        # List of profiles ----------------------------------------------------
        self.box_profiles_list = QtWidgets.QGroupBox("Profiles list")
        self.box_profiles_list_HL = QtWidgets.QHBoxLayout()

        self.tableWidget = QtWidgets.QTableWidget(0, 4)
        self.set_headers()
        # Profiles of the rows of the table (display state of each profile)
        self.table_profiles = []
        self.fill_table()

        # Widths of the columns: the value delta takes the remaining space
        for column, width in enumerate((70, 70, 120)):
            self.tableWidget.setColumnWidth(column, width)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.tableWidget.verticalHeader().hide()

        self.box_profiles_list_HL.addWidget(self.tableWidget)
        self.box_profiles_list.setLayout(self.box_profiles_list_HL)

        self.VL.addWidget(self.box_profiles_plot)
        self.VL.addWidget(self.box_profiles_list)
        self.VL.addStretch(1)

        self.setLayout(self.VL)

        # Check the show colorbar button if needed
        if shared.exp.show_profiles_bars:
            self.checkbox_show_profiles_bars.setChecked(True)

        self.update_MPL("MPL_profiles_canvas")
        self.update_GUI("labels_delta_x")
        self.update_GUI("labels_delta_z")

    def compute_profiles(self):
        """Distances along the profiles [um] and values of the map displayed
        in the results tab at their pixels (in the units of its color bar;
        heights of the topography without map)."""
        dt = shared.exp.current_data
        self.data = dt
        self.profile_list = dt.profile_list
        d_x = (dt.scan_size_x / float(dt.nbr_pixels_x)) / 1000.0
        d_y = (dt.scan_size_y / float(dt.nbr_pixels_y)) / 1000.0
        d_xy = numpy.sqrt((d_x ** 2 + d_y ** 2))

        values, self.value_label = map_values(dt)
        unit = self.value_label.rsplit("[", 1)
        self.value_unit = unit[1].rstrip("]") if len(unit) == 2 else ""

        self.profiles_x = []
        self.profiles_z = []
        pl = self.profile_list
        for profile in range(len(pl)):
            self.profiles_x.append([])
            self.profiles_z.append([])
            for count in range(len(pl[profile][0])):
                if count != 0:
                    if pl[profile][0][count - 1] == pl[profile][0][count]:
                        val = self.profiles_x[profile][-1] + d_x
                    elif pl[profile][1][count - 1] == pl[profile][1][count]:
                        val = self.profiles_x[profile][-1] + d_y
                    else:
                        val = self.profiles_x[profile][-1] + d_xy
                    self.profiles_x[profile].append(val)
                else:
                    self.profiles_x[profile].append(count)
                i = pl[profile][0][count]
                j = pl[profile][1][count]
                self.profiles_z[profile].append(values[i][j])

    def fill_table(self):
        """One row per profile. The profiles keep their display state; at
        the opening only the first one is displayed, the new profiles are
        displayed."""
        opening = not self.table_profiles
        states = {id(profile): checkbox.isChecked() for profile, checkbox
                  in zip(self.table_profiles, self.checkbox_list)}
        self.checkbox_list = []
        self.delta_x_labels = []
        self.delta_z_labels = []
        self.table_profiles = list(self.profile_list)
        self.tableWidget.setRowCount(len(self.table_profiles))
        for i, profile in enumerate(self.table_profiles):
            # Centered values
            item = QtWidgets.QTableWidgetItem(str(i + 1))
            item.setTextAlignment(QtCore.Qt.AlignCenter)
            self.tableWidget.setItem(i, 0, item)

            checkbox = QtWidgets.QCheckBox()
            checkbox.setChecked(states.get(id(profile),
                                           i == 0 if opening else True))
            checkbox.stateChanged.connect(self.list_checkbox_clicked)
            self.checkbox_list.append(checkbox)
            self.tableWidget.setCellWidget(i, 1, CenteredCellwidget(checkbox))

            label_delta_x = QtWidgets.QLabel()
            label_delta_x.setAlignment(QtCore.Qt.AlignCenter)
            self.delta_x_labels.append(label_delta_x)
            self.tableWidget.setCellWidget(i, 2, label_delta_x)
            label_delta_z = QtWidgets.QLabel()
            label_delta_z.setAlignment(QtCore.Qt.AlignCenter)
            self.delta_z_labels.append(label_delta_z)
            self.tableWidget.setCellWidget(i, 3, label_delta_z)

    def set_headers(self):
        """Headers of the table: delta of the values of the map."""
        unit = " [" + self.value_unit + "]" if self.value_unit else ""
        self.tableWidget.setHorizontalHeaderLabels(
            ["Profile", "Display", "Delta X [nm]", "Delta" + unit])
        self.tableWidget.horizontalHeaderItem(3).setToolTip(
            "Difference of the values of the map (" + self.value_label +
            ") between the two bars")

    def update_widget(self):
        """The map of the results tab or the profiles changed: profiles of
        the values of the map."""
        self.compute_profiles()
        self.fill_table()
        self.set_headers()
        self.update_GUI("labels_delta_x")
        self.update_GUI("labels_delta_z")
        self.update_MPL("MPL_profiles_canvas")

    def update_GUI(self, what):
        """Update the GUI."""
        bar1 = shared.exp.profiles_bar_1
        bar2 = shared.exp.profiles_bar_2

        if what == "labels_delta_x":
            for i in range(len(self.profile_list)):
                value = abs(bar1 - bar2)
                self.delta_x_labels[i].setText(str(round(value)))

        elif what == "labels_delta_z":
            # Values of the profiles at the bars (linear interpolation; the
            # distances are in um, the bars in nm)
            for profile in range(len(self.profiles_x)):
                x = numpy.asarray(self.profiles_x[profile], dtype=float)
                z = numpy.asarray(self.profiles_z[profile], dtype=float)
                inside = [x[0] <= bar / 1000.0 <= x[-1]
                          for bar in (bar1, bar2)]
                if len(x) < 2 or not all(inside):
                    # A bar is outside of the profile
                    self.delta_z_labels[profile].setText("\u2013")
                    continue
                z1 = numpy.interp(bar1 / 1000.0, x, z)
                z2 = numpy.interp(bar2 / 1000.0, x, z)
                value = abs(z1 - z2)
                if numpy.isnan(value):
                    # No value at one of the bars (e.g. no POC found)
                    self.delta_z_labels[profile].setText("NaN")
                else:
                    self.delta_z_labels[profile].setText("%.4g" % value)

    def open_single_figure(self):
        """Open the figure in a separate window."""
        PYAFPlot(self, "profiles")

    def list_checkbox_clicked(self):
        """Display or hide a profile."""
        self.update_MPL("MPL_profiles_canvas")

    def update_MPL(self, _):
        """Update the plot."""
        self.canvas.update_plot()

    def input_updated(self, what):
        """Called when an input is updated."""
        if what == "input_bar_1":
            shared.exp.profiles_bar_1 = self.input_bar_1.get_float_value()
        elif what == "input_bar_2":
            shared.exp.profiles_bar_2 = self.input_bar_2.get_float_value()

        self.update_GUI("labels_delta_x")
        self.update_GUI("labels_delta_z")
        self.update_MPL("MPL_profiles_canvas")

    def checkbox_clicked(self, what):
        """Called when a checkbox is clicked."""
        if what == "checkbox_show_profiles_bars":
            state = self.checkbox_show_profiles_bars.isChecked()
            if state:
                shared.exp.show_profiles_bars = True
            else:
                shared.exp.show_profiles_bars = False
            self.update_MPL("MPL_profiles_canvas")

    def canvas_press(self, event):
        """Called when there is a click on the canvas."""
        canv = self.canvas

        if event.button == 3 and event.inaxes == canv.canvas.axes:
            # Position of the mouse on screen (the matplotlib event position
            # is in physical pixels, twice the Qt pixels on Retina screens)
            pos = QtGui.QCursor.pos()

            self.popUpMenu(pos)

    def popUpMenu(self, pos):
        """Menu opened on right click on the figure."""
        # Define actions
        action_open_figure = QtWidgets.QAction("Open figure", self)
        action_open_figure.triggered.connect(self.open_single_figure)

        # Create menu
        menu = QtWidgets.QMenu()
        menu.addAction(action_open_figure)

        # Display action
        menu.popup(pos, menu.menuAction())
        menu.exec_()


def map_values(data):
    """Values of the map displayed in the results tab ([x][y], in the units
    of its color bar, nan without value) and its label. Heights of the
    topography if there is no map."""
    from .. import widgets_list
    canvas = None
    if widgets_list.widget_results is not None:
        canvas = getattr(widgets_list.widget_results.MPL_canvas1, "canvas",
                         None)
    temp = getattr(canvas, "temp_array", None)
    if canvas is not None and temp is not None and \
            getattr(canvas, "data", None) is data and \
            getattr(canvas, "colorbarlabel", None):
        values = numpy.ma.filled(numpy.ma.asarray(temp, dtype=float),
                                 numpy.nan) * canvas.factor
        return values, canvas.colorbarlabel
    return planefit.topo_map(data), "Zero-force topography [nm]"
