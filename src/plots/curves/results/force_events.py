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

"""Plots the force curve with the detected events."""

from ..plot_curve import PlotCurve
from ....tools import math_tools
from ....tools import curve_tools


class ResultsForceEvents(PlotCurve):
    """Plots the force curve with the detected events."""

    def __init__(self, parent):
        super().__init__(parent)

        if self.data.curve_force_units == "nN":
            self.force_factor = 1
        elif self.data.curve_force_units == "pN":
            self.force_factor = 1e3
        if self.data.curve_distance_units == "nm":
            self.x_factor = 1
        elif self.data.curve_distance_units == "um":
            self.x_factor = 1e-3

        # Labels
        if self.data.curve_distance_units == "nm":
            self.xlabel = "Scanner Extension [nm]"
        elif self.data.curve_distance_units == "um":
            self.xlabel = "Scanner Extension [\u03bcm]"
        self.ylabel = "Force [" + self.data.curve_force_units + "]"

        self.curve_retraction = curve_tools.get_force_curves(
            self.retraction[0],
            self.retraction[1],
            self.data.events_jocs2_real[self.xpos][self.ypos],
            self.spring_constant)

        self.joc2_indice = self.data.events_jocs2_indices[self.xpos][self.ypos]

        value = self.data.events_fits_joc[self.xpos][self.ypos]
        params = self.get_params(value, "force")

        # Display jump of contact
        if self.data.events_results_display_joc:
            self.joc2_pos = [0.0, 0.0]

        # Display fit for jump of contact
        if self.data.events_results_display_fit_joc:
            self.setup_lines(params, "ret", "force")

        # Setup the events plotting
        self.setup_events("results")

        # Display events dist filter as a grey square
        if self.data.events_display_results_filter_dist:
            # Get distance between two points on the curve
            xdist = math_tools.get_x_dist(
                self.curve_retraction[0], self.ret_pos[0])

            # xdist can be None for corrupted curves
            if xdist is not None:
                left = \
                    self.data.events_results_filter_dist_left * self.x_factor
                right = \
                    self.data.events_results_filter_dist_right * self.x_factor
                keep_middle = self.data.events_results_filter_dist_keep_middle

                if keep_middle:
                    # First excluded segment
                    # (on the left, at beginning of curve)
                    seg1 = [None, left]

                    # Second excluded segment
                    # (on the right, near the joc2)
                    seg2 = [right, 0]

                else:
                    # Only one excluded segment
                    seg1 = [left, right]
                    seg2 = None

                self.events_filter_dist_square = [seg1, seg2]

        self.define_x_scale()
        self.define_y_scale()

        self.plot_fig()
