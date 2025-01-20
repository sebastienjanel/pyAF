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

"""Plots profiles."""

from PyQt5 import QtCore, QtWidgets
from ..plots.PYAFPlot import PYAFPlot
from ..tools import math_tools
from ..tools.gui_tools import PYAFInput
import numpy
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

        # Prepare the profiles to be displayed. The data preparation is done
        # here because we need these profiles to calculate the x and z deltas
        # in the widget. Furthermore the profiles only need to be calculated
        # once.

        dt = shared.exp.current_data
        d_x = (dt.scan_size_x / float(dt.topography.shape[0])) / 1000.0
        d_y = (dt.scan_size_y / float(dt.topography.shape[1])) / 1000.0
        d_xy = numpy.sqrt((d_x ** 2 + d_y ** 2))

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
                        self.profiles_x[profile].append(val)
                    elif pl[profile][1][count - 1] == pl[profile][1][count]:
                        val = self.profiles_x[profile][-1] + d_y
                        self.profiles_x[profile].append(val)
                    else:
                        val = self.profiles_x[profile][-1] + d_xy
                        self.profiles_x[profile].append(val)
                else:
                    self.profiles_x[profile].append(count)
                i = pl[profile][0][count]
                j = pl[profile][1][count]
                val = dt.topography[i][j] / 1000.0
                self.profiles_z[profile].append(val)

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

        self.tableWidget = QtWidgets.QTableWidget(len(self.profile_list), 4)
        header = ["Profile", "Display", "Delta X [nm]", "Delta Z [nm]"]
        self.tableWidget.setHorizontalHeaderLabels(header)

        # Fill every line of the table
        for i in range(len(self.profile_list)):
            item = QtWidgets.QTableWidgetItem(str(i + 1))
            self.tableWidget.setItem(i, 0, item)

            checkbox = QtWidgets.QCheckBox()
            # We display only the first profile at the beginning
            if i == 0:
                checkbox.setChecked(True)
            checkbox.stateChanged.connect(self.list_checkbox_clicked)
            self.checkbox_list.append(checkbox)
            self.tableWidget.setCellWidget(i, 1, checkbox)

            label_delta_x = QtWidgets.QLabel()
            self.delta_x_labels.append(label_delta_x)
            self.tableWidget.setCellWidget(i, 2, label_delta_x)
            label_delta_z = QtWidgets.QLabel()
            self.delta_z_labels.append(label_delta_z)
            self.tableWidget.setCellWidget(i, 3, label_delta_z)

        self.tableWidget.resizeColumnsToContents()

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

    def update_GUI(self, what):
        """Update the GUI."""
        bar1 = shared.exp.profiles_bar_1
        bar2 = shared.exp.profiles_bar_2

        if what == "labels_delta_x":
            for i in range(len(self.profile_list)):
                value = abs(bar1 - bar2)
                self.delta_x_labels[i].setText(str(round(value)))

        elif what == "labels_delta_z":
            for profile in range(len(self.profiles_x)):
                for i in range(len(self.profiles_x[profile])):
                    x = self.profiles_x[profile][i] * 1000.0
                    if x == shared.exp.profiles_bar_1:
                        z1 = self.profiles_z[profile][i]
                    elif x > shared.exp.profiles_bar_1:
                        x1b = self.profiles_x[profile][i]
                        z1b = self.profiles_z[profile][i]
                        x1a = self.profiles_x[profile][i - 1]
                        z1a = self.profiles_z[profile][i - 1]
                        coeffs, _ = \
                            math_tools.fit_linear([x1a, x1b], [z1a, z1b])
                        z1 = numpy.polyval([coeffs[0], coeffs[1]], bar1)
                        break

                for i in range(len(self.profiles_x[profile])):
                    x = self.profiles_x[profile][i] * 1000.0
                    if x == shared.exp.profiles_bar_2:
                        z2 = self.profiles_z[profile][i]
                    elif x > shared.exp.profiles_bar_2:
                        x2b = self.profiles_x[profile][i]
                        z2b = self.profiles_z[profile][i]
                        x2a = self.profiles_x[profile][i - 1]
                        z2a = self.profiles_z[profile][i - 1]
                        coeffs, _ = \
                            math_tools.fit_linear([x2a, x2b], [z2a, z2b])
                        z2 = numpy.polyval([coeffs[0], coeffs[1]], bar2)
                        break

                value = abs(z1 - z2)
                self.delta_z_labels[profile].setText(str(round(value)))

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
            # Get position in screen coordinates (top left position of canvas)
            globalpos = canv.mapToGlobal(QtCore.QPoint(0, 0))

            # Recalculate positions with canvas position
            x = globalpos.x() + event.x
            y = globalpos.y() + (self.canvas_size[1] - event.y)
            pos = QtCore.QPoint(x, y)

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
