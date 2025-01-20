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

"""Meshgrid for tilt check"""

import numpy
from .plot_meshgrid import PlotMeshGrid
from ...tools import math_tools
from ...widgets.progressbar import Progressbar
from ... import widgets_list


class PlotTiltCheckMeshGrid(PlotMeshGrid):
    """Meshgrid showing if a curve was corrected for the tilt or not."""

    def __init__(self, parent):
        super().__init__(parent)

        # Create an empty array
        self.temp_array = numpy.empty(
            [self.data.nbr_pixels_x, self.data.nbr_pixels_y])
        self.tilt_mask = numpy.zeros(
            [self.data.nbr_pixels_x, self.data.nbr_pixels_y])

        # Create a progressbar which is displayed during the saving
        Progressbar("Checking data")
        total = self.data.nbr_pixels_x * self.data.nbr_pixels_y
        widgets_list.widget_progressbar.set_range(0, total)
        widgets_list.widget_progressbar.set_label("Tilt correction check")

        # Loop through the curves, apply the tilt correction, and check
        # the returned info dictionnary.
        for i in range(self.data.nbr_pixels_x):
            for j in range(self.data.nbr_pixels_y):
                _, _, res = math_tools.correct_tilt(
                    self.data.curves_approach[i][j],
                    self.data.curves_retraction[i][j],
                    self.data.tilt_limit_1,
                    self.data.tilt_limit_2,
                    self.data.tilt_applied,
                    self.data.approach_positions[i][j],
                    self.data.retraction_positions[i][j])

                if self.data.check_tilt_option == "nbr_points":
                    self.temp_array[i][j] = res["length_of_fit_nbr_points"]
                elif self.data.check_tilt_option == "slopes":
                    self.temp_array[i][j] = res["detected_slope"]
                self.tilt_mask[i][j] = res["is_corrected"]

                widgets_list.widget_progressbar.update()

        widgets_list.widget_progressbar.close()

        self.factor = 1.0
        self.colorbarlabel = None
        self.x_label = None
        self.y_label = None
        self.display_colorbar = False

        # Hide the axes
        self.axes.get_xaxis().set_visible(False)
        self.axes.get_yaxis().set_visible(False)

        # Overwrite colormap to have always a green one
        self.colormap = self.colortables.colortables_mpl_list[6]
        self.colortable_max_value = numpy.amax(self.temp_array)
        self.colortable_min_value = 0

        self.plot_fig()
