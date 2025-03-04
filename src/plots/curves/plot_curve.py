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

"""Plots a curve."""

import numpy
from matplotlib.offsetbox import AnchoredText
from ... import shared
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from ...tools import curve_tools
from ...tools import math_tools
from ..plot_main import MainPlot
from .common import CommonCurveTools


class PlotCurve(MainPlot, CommonCurveTools):
    """Subclass of MainPlot, used to plot curves.

    Note for those not knowing how this works:
    The CommonCurveTools class is the second parent class. You can have two
    parent classes in python. So you can look at the curve_tools module.
    This is a convenient way to share variables and methods between classes.
    If you stumble on a method and you don't know where it is defined, just
    look at the CommonCurveTools class.
    """

    def __init__(self, parent):
        super().__init__(parent)

        # The array of curve which will be finally plotted
        # These arrays are defined in the different subplotting classes
        self.curve_approach = None
        self.curve_retraction = None
        self.curve_pause = None
        self.curve_modulation = None
        self.curve_residuals = None

        # Fit noise values
        self.fit_noise_values = None

        # Approach (default empty values, will not be plotted if not filled)
        self.poc_pos_indice = None
        self.poc_pos = None
        self.fits_lines_poc = None
        self.segments = None
        self.segments_fits = None
        self.chisqr = None
        self.residuals = None
        self.max_indentation = None
        self.min_indentation = None
        self.max_force_range = None
        self.min_force_range = None

        # Retraction (default empty values, will not be plotted if not filled)
        self.joc1_pos = None
        self.joc2_pos = None
        self.joc1_indice = None
        self.joc2_indice = None
        self.fits_lines_joc = None
        self.max_force = None
        self.surface_params = None

        # Events (default empty values, will not be plotted if not filled)
        self.loading_rates_fits = None
        self.events_filter_dist_square = None
        self.display_events_detection_exclude_square = False
        self.current_jump = None
        self.zoom_current_jump = None

        # Fit preview
        self.fit_segment = None
        self.no_fit_found = False
        self.events_preview = None
        self.fits_lines_conv = [None, None]

        # Vertical lines displaying the tilt/stretch limits
        self.curve_mod_limit_1 = None
        self.curve_mod_limit_2 = None
        self.color_curve_mod_limit = None

        # Events
        self.events_ref_segment = None
        self.events_fits = None
        self.force_fits = None
        self.events_force_values = None
        self.events_force_slopes = None
        self.events_fit_vert = None

        # If the data is only made of a single curve, there is no need to
        # have a particular position on the meshgrid, and the only curve to be
        # displayed is stored at the position 0,0 in the arrays.
        if self.data.is_single:
            self.xpos = 0
            self.ypos = 0
        else:
            self.xpos = shared.exp.meshgrid_click_xpos
            self.ypos = shared.exp.meshgrid_click_ypos

        self.app_pos = self.data.approach_positions[self.xpos][self.ypos]
        self.ret_pos = self.data.retraction_positions[self.xpos][self.ypos]

        # Sometimes these positions are stored as floats and cause indexing to fail.
        # Make sure they are int
        self.app_pos = self.app_pos.astype("uint32")
        self.ret_pos = self.ret_pos.astype("uint32")

        # Valid segment of the curve
        self.valid_seg_app = slice(self.app_pos[0], self.app_pos[1])
        self.valid_seg_ret = slice(self.ret_pos[0], self.ret_pos[1])

        # Raw data
        self.display_raw_data = \
            self.settings.value("DisplayRawData", False)

        self.zero_defl = False
        if self.curve_type == "defl_ext" or self.curve_type == "force_ext":
            self.zero_defl = self.settings.value("ZeroDefl", True)

        # Get the curve with corrections
        self.approach, self.retraction, self.pause, self.modulation, self.zero, _ = \
            curve_tools.get_curve(
                self.data, [self.xpos, self.ypos],
                mode=self.plot_type,
                zero_defl=self.zero_defl)

        # Has to be modified manually for each subplot if you want to change
        # this value. 1 is the default.
        self.force_factor = 1
        self.x_factor = 1

        pt = self.plot_type

        if pt == "fit_preview" or pt == "tilt_curve_preview":
            self.spring_constant = self.data.spring_constant * 1e9
            self.deflection_sensitivity = self.data.deflection_sensitivity
        else:
            # In the case of "computed" curves use the stored values
            spring_const = self.data.used_spring_constant
            if spring_const is None:
                # If nothing has been computed, there is no used_spring_const
                # but we still want to display the curve : use the normal one
                spring_const = self.data.spring_constant
            self.spring_constant = spring_const * 1e9
            self.deflection_sensitivity = self.data.used_deflection_sensitivity

    def plot_fig(self):
        """Plots the curve and asked informations."""
        pt = self.plot_type
        ct = shared.exp.compute_type
        dt = self.data

        # Plot the approach curve
        if self.curve_approach is not None:

            if self.curve_type != "indentation" and pt != "curve_dataa" and pt != "curve_datab":
                positions = dt.approach_positions[self.xpos][self.ypos]
                positions = positions.astype('int32')
                curve_x = self.curve_approach[0][positions[0]:positions[1]]
                curve_y = self.curve_approach[1][positions[0]:positions[1]]

            elif pt == "curve_dataa":
                positions = dt.approach_positions[self.xpos][self.ypos]
                positions = positions.astype('int32')
                curve_x = self.curve_approach[2][positions[0]:positions[1]]
                curve_y = self.curve_approach[1][positions[0]:positions[1]]

            elif pt == "curve_datab":
                positions = dt.approach_positions[self.xpos][self.ypos]
                positions = positions.astype('int32')
                curve_x = self.curve_approach[2][positions[0]:positions[1]]
                curve_y = self.curve_approach[0][positions[0]:positions[1]]

            else:
                # In case of the indentation plot, directly plot the curve
                # without removing the corrupted part.
                curve_x = self.curve_approach[0]
                curve_y = self.curve_approach[1]

            self.axes.plot(curve_x, curve_y, "#00FFFF", label="Approach")     # Changed for testing. Implement stylesheet.

            if self.display_raw_data and self.curve_type == "defl_ext" \
                    and pt != "fit_preview" and pt != "tilt_curve_preview":
                # Display raw data only if asked, on the curve on the data tab
                curve_x = self.curve_approach[0][0:positions[0] + 1]
                curve_y = self.curve_approach[1][0:positions[0] + 1]
                self.axes.plot(curve_x, curve_y, "g")
                curve_x = self.curve_approach[0][positions[1] - 1:-1]
                curve_y = self.curve_approach[1][positions[1] - 1:-1]
                self.axes.plot(curve_x, curve_y, "g")

        if self.curve_pause is not None and (pt == "curve_dataa" or pt == "curve_datab"):
            # Display raw data only if asked, on the curve on the data tab
            if pt == "curve_dataa":
                curve_x = self.curve_pause[2]
                curve_y = self.curve_pause[1]

            elif pt == "curve_datab":
                curve_x = self.curve_pause[2]
                curve_y = self.curve_pause[0]

                if self.curve_approach is not None:
                    curve_y += self.curve_approach[0][-1]

            self.axes.plot(curve_x, curve_y, "#FF0000", label="Pause")

        if self.curve_modulation is not None and (pt == "curve_dataa" or pt == "curve_datab"):
            # Display raw data only if asked, on the curve on the data tab
            if pt == "curve_dataa":
                curve_x = self.curve_modulation[2]
                curve_y = self.curve_modulation[1]

                # if self.curve_pause is not None:
                    # curve_y += self.curve_pause[1][-1]

                # else:
                    # curve_y += self.curve_approach[1][-1]

            elif pt == "curve_datab":
                curve_x = self.curve_modulation[2]
                curve_y = self.curve_modulation[0]

                if self.curve_pause is not None:
                    curve_y += self.curve_pause[0][-1]

                else:
                    curve_y += self.curve_approach[0][-1]

            self.axes.plot(curve_x, curve_y, "#00FF80", label="Modulation")

        # Plot the retraction curve
        if self.curve_retraction is not None:

            if pt != "curve_dataa" and pt != "curve_datab":
                positions = dt.retraction_positions[self.xpos][self.ypos]
                positions = positions.astype('int32')
                curve_x = self.curve_retraction[0][positions[0]:positions[1]]
                curve_y = self.curve_retraction[1][positions[0]:positions[1]]

            elif pt == "curve_dataa":
                positions = dt.retraction_positions[self.xpos][self.ypos]
                positions = positions.astype('int32')
                curve_x = self.curve_retraction[2][positions[0]:positions[1]]
                curve_y = numpy.flipud(self.curve_retraction[1][positions[0]:positions[1]])

            elif pt == "curve_datab":
                positions = dt.retraction_positions[self.xpos][self.ypos]
                positions = positions.astype('int32')
                curve_x = self.curve_retraction[2][positions[0]:positions[1]]
                curve_y = numpy.flipud(self.curve_retraction[0][positions[0]:positions[1]])

            self.axes.plot(curve_x, curve_y, "#FFA500", label="Retraction")     # Changed for testing. Implement stylesheet.

            if self.display_raw_data and self.curve_type == "defl_ext" \
                    and pt != "fit_preview" and pt != "tilt_curve_preview":
                # Display raw data only if asked, on the curve on the data tab
                curve_x1 = self.curve_retraction[0][0:positions[0] + 1]
                curve_y1 = self.curve_retraction[1][0:positions[0] + 1]
                self.axes.plot(curve_x1, curve_y1, "g")
                curve_x2 = self.curve_retraction[0][positions[1] - 1:-1]
                curve_y2 = self.curve_retraction[1][positions[1] - 1:-1]
                self.axes.plot(curve_x2, curve_y2, "g")

        # Plot corrupted part on the fit preview window. This part will display
        # a dashed line from the real zero from the curve to the first valid
        # value. We don't want to display the real corrupted data here as this
        # would mess up the plot.
        if pt == "fit_preview" and \
                (ct == "stiffness" or
                 ct == "work_and_rupture_force" or
                 ct == "events" or
                 ct == "loading_rates"):
            curve_x = None
            curve_y = None

            if ct == "stiffness":
                if self.curve_approach is not None:
                    curve_x = self.curve_approach[0]
                    curve_y = self.curve_approach[1]

            elif ct == "work_and_rupture_force" or ct == "events"\
                    or ct == "loading_rates":
                if self.curve_retraction is not None:
                    curve_x = self.curve_retraction[0]
                    curve_y = self.curve_retraction[1]

            # Get distance between two points on the curve
            if curve_x is not None:
                xdist = math_tools.get_x_dist(curve_x, self.ret_pos[0])
            else:
                xdist = None

            if xdist is not None:
                # Estimate the x position of the first value if the curve
                # was not corrupted.
                pos_x = curve_x[positions[0]] - (xdist * positions[0])
                # Check if the curve is long enough to compute a y position for
                # the dashed line
                if len(curve_y) > 40:
                    pos_y = numpy.mean(curve_y[positions[0]:positions[0] + 30])
                else:
                    # Anyway this curve is very short, so we take the last y
                    # position.
                    pos_y = curve_y[-1]
                x = [pos_x, curve_x[positions[0]]]
                y = [pos_y, pos_y]
                dist = abs(pos_x - curve_x[positions[0]])
                # Display the dashed line only in case where there is more than
                # 30 nm of data missing
                if dist > 30:
                    self.axes.plot(x, y, "k--", label="Missing")

        # Plot the labels
        self.plot_labels()

        # Plot the fit noise values
        if self.fit_noise_values is not None and pt == "fit_preview":
            textstr = ""
            for fit_noise_value in self.fit_noise_values:
                textstr += fit_noise_value + "\n"
            # Place a text box in upper left in axes coords
            self.axes.text(
                0.13, 0.97,
                textstr,
                transform=self.axes.transAxes,
                fontsize=10,
                verticalalignment="top")

        # Display or hide the legend (upper left, loc = 2)
        if shared.exp.display_curve_legend and \
                (self.curve_approach is not None or
                 self.curve_retraction is not None):
            if (pt == "curve_dataa" or pt == "curve_datab") and self.modulation is not None:
                self.axes.legend(loc=0, prop=self.font)

            else:
                self.axes.legend(loc=2, prop=self.font)

        # Poc
        if self.poc_pos is not None:
            self.plot_dot(self.poc_pos[0], self.poc_pos[1], "m.")

        # Plot maximum indentation for fit
        if self.max_indentation is not None:
            self.axes.axvline(x=self.max_indentation, color='k', linestyle='--')

        # Plot minimum indentation for fit
        if self.min_indentation is not None:
            self.axes.axvline(x=self.min_indentation, color='k', linestyle='--')

        if self.max_force_range is not None:
            self.axes.axhline(y=self.max_force_range, color='k', linestyle='--')

        # Plot minimum indentation for fit
        if self.min_force_range is not None:
            self.axes.axhline(y=self.min_force_range, color='k', linestyle='--')

        # Fits for the poc
        if self.fits_lines_poc is not None:
            self.plot_fits_lines(
                self.curve_approach[0][self.valid_seg_app],
                self.fits_lines_poc)

        # Segments
        if self.segments is not None:
            for seg in self.segments:
                self.plot_dot(seg[0], seg[1], "r.")

        # Fits for the stiffness
        if self.segments_fits is not None:
            self.plot_stiffness_fits()

        # Plot residuals
        if self.curve_residuals is not None:
            self.plot_fit_residuals()

        # Chisq
        if self.chisqr is not None and self.chisqr != 0:
            anchored_text = AnchoredText("chi-square: {:.5E}".format(self.chisqr),
                                         frameon=False, prop=dict(alpha=None, fontsize='large'),
                                         loc=1)   # Upper right position
            self.axes.add_artist(anchored_text)

        # Plot the point of jump of contact and point of intersection between
        # fit and curve (joc2)
        if self.joc1_pos is not None:
            self.plot_dot(self.joc1_pos[0], self.joc1_pos[1], "k.")
            self.plot_dot(self.joc2_pos[0], self.joc2_pos[1], "k.")

        # Fits for the joc
        if self.fits_lines_joc is not None:
            self.plot_fits_lines(
                self.curve_retraction[0][self.valid_seg_ret],
                self.fits_lines_joc)

        # Plot the maximum force
        if self.max_force is not None:
            self.axes.plot(
                self.max_force[0], self.max_force[1], "r", linewidth=2)

        # Plot the surface
        if self.surface_params is not None:
            self.axes.fill_between(
                self.surface_params[0],
                self.surface_params[1],
                self.surface_params[2],
                color="k")

        # Fits for the joc (events)
        if self.fits_lines_joc is not None:
            self.plot_fits_lines(
                self.curve_retraction[0][self.valid_seg_ret],
                self.fits_lines_joc)

        # Plot the events joc
        if self.joc2_pos is not None:
            self.plot_dot(self.joc2_pos[0], self.joc2_pos[1], "k.")

        # Plot the fits for the loading rates
        if self.loading_rates_fits is not None and \
                self.loading_rates_fits != []:
            for fit in self.loading_rates_fits:
                self.axes.plot(fit[0], fit[1], "k", linewidth=3)

        # Plot the events on results tab (events force curve)
        if self.events_preview is not None and self.events_preview != []:
            # Add annotations if asked
            self.display_events_values()

            for i in range(len(self.events_preview)):
                self.plot_events_force(i)

                self.plot_events_fit_base(i)

        # Plot some data for the fit previews
        if pt == "fit_preview":
            if ct == "stiffness" and self.curve_approach is not None and \
                    self.fit_segment is not None:
                self.axes.plot(self.curve_approach[0][self.fit_segment],
                               self.curve_approach[1][self.fit_segment],
                               "r", linewidth=4.0)

            elif ct == "work_and_rupture_force" and \
                    self.curve_retraction is not None and \
                    self.fit_segment is not None:
                self.axes.plot(self.curve_retraction[0][self.fit_segment],
                               self.curve_retraction[1][self.fit_segment],
                               "r", linewidth=4.0)

            elif ct == "events" and self.curve_retraction is not None:
                ev_prev = self.events_preview

                if ev_prev is not None and ev_prev != []:
                    for i in range(len(self.events_preview)):
                        """
                        positions = self.events_preview[i]
                        for pos in positions:
                            self.plot_dot(
                                self.curve_retraction[0][pos],
                                self.curve_retraction[1][pos],
                                "k.")
                        """

                        if self.data.display_fit_event_seg:
                            self.plot_events_fit(i)

                        if self.data.display_fit_event:
                            self.plot_events_force(i)

            elif ct == "events" and self.curve_retraction is not None:
                if self.events_ref_segment is not None:
                    seg = self.events_ref_segment
                    self.axes.plot(
                        self.curve_retraction[0][seg],
                        self.curve_retraction[1][seg],
                        "k", linewidth=4.0)

            # If no fit found, display a label
            if self.no_fit_found:
                self.axes.annotate(
                    "No fit found",
                    xy=(0.9, 0.92),
                    xycoords = "axes fraction",
                    fontproperties = self.font)

            if pt == "fit_preview" and ct == "events":
                self.plot_fits_lines(
                    self.curve_retraction[0], self.fits_lines_conv)

        # Plot vertical lines for the tilt/stretch limits
        if self.curve_mod_limit_1 is not None:
            col = self.color_curve_mod_limit
            self.axes.axvline(
                self.curve_mod_limit_1,
                linewidth=2,
                color=col,
                linestyle='--')
            self.axes.axvline(
                self.curve_mod_limit_2,
                linewidth=2,
                color=col,
                linestyle='--')

        # Plot a red line for discarded curves
        self.red_line()

        # Set limits for the zooming in
        self.set_lims()

        # On results tab, in detection mode for the events, display a grey
        # square for the excluded area
        # (no events are detected right of the joc)
        if self.display_events_detection_exclude_square:
            size = self.axes.axis()  # (x0, x1, y0, y1)
            width = size[1]
            height = abs(size[3] - size[2])
            rect = Rectangle(
                (0, size[2]),
                width,
                height,
                facecolor="grey",
                alpha=0.4)

            self.axes.add_patch(rect)

        # Plot events dist filter as a grey square
        if self.events_filter_dist_square is not None:
            vals = self.events_filter_dist_square
            size = self.axes.axis()  # (x0, x1, y0, y1)

            if vals[1] is not None:
                seg1 = vals[0]
                seg2 = vals[1]

                height = abs(size[3] - size[2])
                width1 = abs(size[0] - seg1[1])
                width2 = abs(seg2[0])

                # Left rectangle, near beginning of the curve
                rect1 = Rectangle(
                    (size[0], size[2]),
                    width1,
                    height,
                    facecolor="grey",
                    alpha=0.4)

                # Right rectangle, near joc2
                rect2 = Rectangle(
                    (seg2[0], size[2]),
                    width2,
                    height,
                    facecolor="grey",
                    alpha=0.4)

                self.axes.add_patch(rect1)
                self.axes.add_patch(rect2)

            else:
                seg1 = vals[0]

                height = size[3] - size[2]
                width = abs(seg1[1] - seg1[0])

                rect = Rectangle(
                    (seg1[0], size[2]),
                    width,
                    height,
                    facecolor="grey",
                    alpha=0.4)

                self.axes.add_patch(rect)

        # Do the actual plotting, is defined in MainPlot class
        self.start_plot()

    def set_lims(self):
        """Set the x, y limits of the plots.

        For the choice of the fitting parameters there is a possibility
        to zoom in on the detected point of contact (joc or poc)
        Note : control for presence of a poc/joc and fit, it can be that
        there is a corrupted curve. Then it has no meaning to define
        xlim and ylim.
        """
        if self.plot_type == "fit_preview" and \
            (self.poc_pos is not None or self.joc1_pos is not None or
             self.data.auto_zoom_in):
            point = None

            if self.curve_approach is not None:
                point = self.poc_pos_indice
                curve_x = self.curve_approach[0]
                curve_y = self.curve_approach[1]

            elif self.curve_retraction is not None:
                curve_x = self.curve_retraction[0]
                curve_y = self.curve_retraction[1]

                if shared.exp.compute_type == "events":
                    if self.events_preview is not None and \
                            self.events_preview != []:
                        if shared.zoomed_in_element < len(self.events_preview):
                            # We zoom on an event
                            elem = shared.zoomed_in_element
                            point = self.events_preview[elem][1]
                        else:
                            # We want to zoom on the joc
                            point = self.joc2_indice
                    else:
                        # No event, we want to zoom on the joc
                        point = self.joc2_indice

                elif shared.exp.compute_type == "work_and_rupture_force":
                    if shared.zoomed_in_element == 0:
                        point = self.joc1_indice
                    else:
                        point = self.joc2_indice

            if self.no_fit_found:
                # If no fit has been found, do not zoom
                point = None

            # Corrupted curves have sometimes 0 length, in this case don't
            # use the zoom in.
            if len(curve_x) != 0 and point is not None:
                # 10 times zoom on the X axis
                xlength = len(curve_x) / 10.0

                # Get distance between two points on the curve
                if point - 20 < 0:
                    start = 0
                else:
                    start = point - 20
                if point + 10 > len(curve_x):
                    stop = len(curve_x)
                else:
                    stop = point + 10

                seg = slice(start, stop)
                factor = self.data.zoom_fit_preview_factor

                if not factor:
                    return  # Cancel if = 0

                diff = (abs(max(curve_y[seg]) - min(curve_y[seg])) / 2.0)
                ylength = diff * (10.0 / factor)

                if shared.exp.compute_type == "events":
                    # Zoom more for the events
                    ylength / 8.0

                if self.data.auto_zoom_in and point is not None:
                    if point + xlength > len(curve_x):
                        xmax = curve_x[len(curve_x) - 1]
                    else:
                        xmax = curve_x[point + int(xlength)]
                    xmin = curve_x[point - int(xlength)]

                    ymin = numpy.mean(curve_y[seg]) - ylength
                    ymax = numpy.mean(curve_y[seg]) + ylength
                    self.axes.set_xlim([xmin, xmax])
                    self.axes.set_ylim([ymin, ymax])

    def plot_dot(self, pos_x, pos_y, color):
        """Plots a dot on the curve"""
        self.axes.plot(pos_x, pos_y, color, markersize=10.0)

    def plot_fits_lines(self, curve_x, fits):
        """Plot the fitting lines."""
        # The fit
        if fits[0] is not None:
            self.axes.plot(curve_x, fits[0], "r-", linewidth=1)

        # The fit + error
        if fits[1] is not None:
            self.axes.plot(curve_x, fits[1], "g-", linewidth=1)

    def plot_stiffness_fits(self):
        """Plot the fits for each stiffness segment."""
        for seg_fit in self.segments_fits:
            self.axes.plot(seg_fit[0], seg_fit[1], "k-", linewidth=2)

    def plot_fit_residuals(self):
        for seg_res in self.residuals:
            self.axes.plot(seg_res[0], seg_res[1], "g-", linewidth=2)

    def red_line(self):
        """Display red line for discarded curve."""
        if self.data.discarded_curves[self.xpos][self.ypos]:
            values1 = self.axes.transLimits.inverted().transform((0, 0))
            values2 = self.axes.transLimits.inverted().transform((1, 1))
            verts = numpy.array([(values1[0], values1[1]),
                                 (values2[0], values2[1])])
            linecoll = LineCollection(
                [verts],
                linewidths=2,
                colors="r",
                linestyle="solid")
            self.axes.add_collection(linecoll)

    def plot_events_fit_base(self, event):
        return False

        fit = self.events_fits[event]

        # Fitting segments
        lbase = fit["seg_left_fit_base"]
        a = slice(lbase[0], lbase[1] + 1)
        seg_x = self.curve_retraction[0][a]
        seg_y = self.curve_retraction[1][a]
        self.axes.plot(seg_x, seg_y, color="k", linewidth=5)

        rbase = fit["seg_right_fit_base"]
        a = slice(rbase[0], rbase[1] + 1)
        seg_x = self.curve_retraction[0][a]
        seg_y = self.curve_retraction[1][a]
        self.axes.plot(seg_x, seg_y, color="k", linewidth=5)

    def plot_events_fit(self, event):
        fit = self.events_fits[event]

        a = slice(fit["seg_left_fit"][0], fit["seg_left_fit"][1] + 1)
        seg_x = self.curve_retraction[0][a]
        y = numpy.polyval(
            [fit["coeffs_left"][0], fit["coeffs_left"][1]], seg_x)
        self.axes.plot(seg_x, y, color="b", linewidth=1)

        a = slice(fit["seg_right_fit"][0], fit["seg_right_fit"][1] + 1)
        seg_x = self.curve_retraction[0][a]
        y = numpy.polyval(
            [fit["coeffs_right"][0], fit["coeffs_right"][1]], seg_x)
        self.axes.plot(seg_x, y, color="b", linewidth=1)

        fit = self.events_fit_vert[event]
        self.axes.plot(fit[0], fit[1], color="k", linewidth=1)

    def plot_events_force(self, event):
        # The force segment taken into account
        left = self.force_fits[event][0]
        right = self.force_fits[event][1]
        x = [left[0] * self.x_factor, right[0] * self.x_factor]
        y = [left[1] * self.force_factor, right[1] * self.force_factor]
        self.axes.plot(x, y, color="g", linewidth=3)

    def display_events_values(self):
        # Annotate with values (force and slopes)
        if self.events_force_values is not None and self.events_force_slopes is not None:
            y_pos = 0
            min_len = min(len(self.events_force_values), len(self.events_force_slopes))

            for i in range(min_len):  # Use the smaller length
                force = str(
                    round(self.events_force_values[i], 4) * self.force_factor)
                deltaslope = str(self.events_force_slopes[i] * 1000.0)

                unit = self.data.curve_force_units
                text = f"Event {i}, {force} {unit}, dslope {deltaslope}"

                self.axes.annotate(
                    text,
                    xy=(0.015, 0.85 - y_pos),
                    xycoords="axes fraction",
                    fontproperties=self.font)
                y_pos += 0.05
