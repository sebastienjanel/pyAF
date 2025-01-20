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

"""Displays the number of events per scan."""

from PyQt5 import QtCore, QtWidgets
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget


class EventsPerScanWidget(PYAFWidget):
    """Widget which displays the events per scan."""

    def __init__(self, parent):
        super().__init__(parent,
                                                  "widget_events_per_scan")

        self.canvas_resolution = self.parent.canvas_resolution
        self.canvas_size = [432, 432]

        VL = QtWidgets.QVBoxLayout()

        self.box = QtWidgets.QGroupBox("Events per scan")
        self.VL_box = QtWidgets.QVBoxLayout()

        self.canvas = QtWidgets.QWidget()
        self.canvas.setFixedSize(self.canvas_size[0], self.canvas_size[1])
        sizes = [
            self.canvas_size[0] / self.canvas_resolution,
            self.canvas_size[1] / self.canvas_resolution,
            self.canvas_resolution]
        self.MPL_canvas = PYAFPlot(
            self, "events_per_scan", self.canvas, sizes)
        self.MPL_canvas.mpl_connect("button_press_event", self.canvas_press)

        self.VL_box.addWidget(self.MPL_canvas)
        self.box.setLayout(self.VL_box)

        VL.addWidget(self.box)
        VL.addStretch(1)
        self.setLayout(VL)

        self.update_MPL()

    def open_single_figure(self):
        """Opens the figure in a separate window."""
        PYAFPlot(self, "events_per_scan")

    def update_MPL(self):
        """Updates the plot."""
        self.MPL_canvas.update_plot()

    def canvas_press(self, event):
        """Called when there is a click on the canvas.

        This opens a drop down menu on right clicking,
        """
        if event.button == 3:
            # Get position in screen coordinates (top left position of canvas)
            globalpos = self.MPL_canvas.mapToGlobal(QtCore.QPoint(0, 0))

            # Recalculate positions with canvas position
            x = globalpos.x() + event.x
            y = globalpos.y() + (self.canvas_size[1] - event.y)
            pos = QtCore.QPoint(x, y)

            # Open the menu
            self.popUpMenu(pos)

    def popUpMenu(self, pos):
        """Opens the menu."""
        # Define actions
        action_open_figure = QtWidgets.QAction("Open figure", self)
        action_open_figure.triggered.connect(self.open_single_figure)

        # Create menu
        menu = QtWidgets.QMenu()
        menu.addAction(action_open_figure)

        # Display action
        menu.popup(pos, menu.menuAction())
        menu.exec_()

    def update_widget(self):
        """Update the GUI."""
        self.update_MPL()
