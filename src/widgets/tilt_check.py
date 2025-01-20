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

"""Widget showing which curve has been corrected for the tilt."""

import logging
from PyQt5 import QtWidgets
from .. import shared
from .. import widgets_list
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFButtonGroup
from ..tools import misc_tools


class CheckTiltWidget(PYAFWidget):
    """Widget with meshgrid."""

    def __init__(self, parent):
        super().__init__(parent, "widget_tilt_check")

        self.logger = logging.getLogger()

        VL = QtWidgets.QVBoxLayout()

        HL_options = QtWidgets.QHBoxLayout()
        self.BTG = PYAFButtonGroup(self, "options")
        self.RBT_nbr_points = QtWidgets.QRadioButton("Nbr of points in fit")
        self.RBT_slopes = QtWidgets.QRadioButton("Fit slopes")
        self.BTG.addButton(self.RBT_nbr_points, 0)
        self.BTG.addButton(self.RBT_slopes, 1)

        HL_options.addWidget(self.RBT_nbr_points)
        HL_options.addWidget(self.RBT_slopes)
        HL_options.addStretch(1)

        self.W_plot = QtWidgets.QWidget()
        self.W_plot.setFixedSize(504, 504)
        self.MPL_canvas = PYAFPlot(
            self, "meshgrid_tilt_check", self.W_plot, [7, 7, 72])
        self.MPL_canvas.mpl_connect(
            "button_press_event", self.mesh_canvas_press)

        LB = QtWidgets.QLabel("Note: Curves with a cross are not corrected at all")

        VL.addLayout(HL_options)
        VL.addWidget(self.W_plot)
        VL.addWidget(LB)

        self.setLayout(VL)

        self.update_widget()

    def update_widget(self):
        """Update the widget."""
        data = shared.exp.current_data
        if data.check_tilt_option == "nbr_points":
            self.RBT_nbr_points.setChecked(True)
        elif data.check_tilt_option == "slopes":
            self.RBT_slopes.setChecked(True)

        self.MPL_canvas.update_plot()

    def button_clicked(self, name):
        """Method called whenever a button is clicked."""
        if name == "options":
            data = shared.exp.current_data

            val = self.BTG.checkedId()
            if val == 1:
                data.check_tilt_option = "slopes"
            else:
                data.check_tilt_option = "nbr_points"

            self.MPL_canvas.update_plot()

    def mesh_canvas_press(self, event):
        """Matplotlib method, called when the user clicks on the meshgrid."""
        if event.button == 1 and event.inaxes == self.MPL_canvas.canvas.axes:
            xpos, ypos = misc_tools.get_position_on_meshgrid(event)

            widgets_list.widget_main.change_curve(xpos, ypos)

            self.MPL_canvas.canvas.update_blit("red_square")
