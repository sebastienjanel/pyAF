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

"""Common module containing utilities for the plotting of AFM curves."""

import math
import numpy
import copy
from ...tools import math_tools
from ...tools import events_tools
from ...tools import curve_tools
from ...compute_tools.stiffness import stiffness_compute


class CommonCurveTools:
    """Common tools used for the plotting of the events."""

    def __init__(self):
        pass

    def setup_events(self, mode):
        """Get the events to display them on the force curve."""
        if mode == "fit_preview":
            fit_size = self.data.events_fit_size

            if self.data.events_algorithm == "kernel":
                # Get threshold
                if self.data.adaptive_threshold_option:
                    threshold = \
                        self.data.events_kernel_adaptive_threshold
                else:
                    threshold = self.data.events_kernel_detection_threshold

                # Parameters for the event detection
                convolution_params = {
                    "kernel_size":
                        self.data.kernel_size,
                    "threshold": threshold,
                    "adaptive_threshold_option":
                        self.data.adaptive_threshold_option,
                    "adaptive_smoothing_window":
                        self.data.adaptive_smoothing_window,
                    "adaptive_smoothing_order":
                        self.data.adaptive_smoothing_order}

                # Get events position
                positions = events_tools.get_events_positions_by_kernel(
                    self.curve_retraction[0],
                    self.curve_retraction[1],
                    self.joc2_indice,
                    self.ret_pos[0],
                    convolution_params)

            elif self.data.events_algorithm == "msf":
                positions = events_tools.get_events_positions_by_msf(
                    self.curve_retraction,
                    self.data.events_msf_detection_threshold,
                    self.data.events_msf_window_size,
                    self.joc2_indice,
                    self.ret_pos[0])

        elif mode == "results":
            fit_size = self.data.used_events_fit_size

            if self.data.events_algorithm == "kernel":
                # Get threshold
                if self.data.used_adaptive_threshold_option:
                    threshold = \
                        self.data.used_events_kernel_adaptive_threshold
                else:
                    threshold = \
                        self.data.used_events_kernel_detection_threshold

                # Parameters for the event detection
                convolution_params = {
                    "kernel_size":
                        self.data.used_kernel_size,
                    "threshold": threshold,
                    "adaptive_threshold_option":
                        self.data.used_adaptive_threshold_option,
                    "adaptive_smoothing_window":
                        self.data.used_adaptive_smoothing_window,
                    "adaptive_smoothing_order":
                        self.data.used_adaptive_smoothing_order}

                # Get events position
                positions = events_tools.get_events_positions_by_kernel(
                    self.curve_retraction[0],
                    self.curve_retraction[1],
                    self.joc2_indice,
                    self.ret_pos[0],
                    convolution_params)

            elif self.data.events_algorithm == "msf":
                positions = events_tools.get_events_positions_by_msf(
                    self.curve_retraction,
                    self.data.used_events_msf_detection_threshold,
                    self.data.used_events_msf_window_size,
                    self.joc2_indice,
                    self.ret_pos[0])

        fits = events_tools.find_fits_around_middle_pos(
            self.curve_retraction, positions, fit_size)

        force_fits, fit_pos_vert, forces = events_tools.find_force_from_fits(
            self.curve_retraction, fits, positions)

        # Plot the events
        self.events_preview = positions
        self.events_fits = fits
        self.force_fits = force_fits
        self.events_fit_vert = fit_pos_vert

        if mode == "results" and self.data.events_results_display_annotations:
            # Used to display the value of the force for the event on the
            # results tab curve plot.
            self.events_force_values = forces
            self.events_force_slopes = \
                self.data.events_slopes[self.xpos][self.ypos]

    def setup_fit_noise(self, curve_type, segment):
        """Display the fit noise values if asked.

        Different coefficients are calculated as there is no clear
        concencuss on what type of noise should be calculated.
        """
        display = self.settings.value("DisplayFitNoise", False)

        if display:
            if curve_type == "approach":
                curve_x = self.curve_approach[0][segment]
                curve_y = self.curve_approach[1][segment]
            elif curve_type == "retraction":
                curve_x = self.curve_retraction[0][segment]
                curve_y = self.curve_retraction[1][segment]

            rms, std, determination, sem = \
                math_tools.compute_noises(curve_x, curve_y)

            self.fit_noise_values = []
            self.fit_noise_values.append(
                "Coefficient of determination: " + str(determination))
            self.fit_noise_values.append(
                "Standard error of the regression: " + str(sem))
            self.fit_noise_values.append("Standard deviation: " + str(std))
            self.fit_noise_values.append("RMS: " + str(rms))

    def get_work_surface(self, params, joc1_indice, joc2_indice):
        """Setup the parameters to plot the surface for the work."""
        [a, b, _] = params

        # Use deepcopy, these are lists, else it will mess up the
        # fits and jocs values ...

        segment = slice(joc1_indice, joc2_indice)
        segx = copy.deepcopy(self.retraction[0][segment])
        segy = copy.deepcopy(self.retraction[1][segment])
        segx[0] = copy.deepcopy(self.joc1_pos[0])
        segy[0] = copy.deepcopy(self.joc1_pos[1])
        segx[len(segx) - 1] = copy.deepcopy(self.joc2_pos[0])
        segy[len(segy) - 1] = copy.deepcopy(self.joc2_pos[1])

        y = numpy.polyval([a, b], segx)

        return [segx, segy, y - self.zero]

    def define_x_scale(self):
        """Set the x scale to microns if needed."""
        if self.plot_type == "fit_preview":
            return False

        if self.curve_approach is not None:
            self.curve_approach[0] = self.curve_approach[0] * self.x_factor
        if self.curve_retraction is not None:
            self.curve_retraction[0] = self.curve_retraction[0] * self.x_factor

    def define_y_scale(self):
        """Set the y scale to pN if needed."""
        if self.plot_type == "fit_preview":
            return False

        if self.curve_approach is not None:
            self.curve_approach[1] = self.curve_approach[1] * self.force_factor
        if self.curve_retraction is not None:
            self.curve_retraction[1] = \
                self.curve_retraction[1] * self.force_factor

    def setup_lines(self, params, app_or_ret, mode):
        """Gets the lines for the fits.

        The fits are defined only on the non-corrupted part of the curve.
        """
        if params is None:
            # If no params are found (no fit is found), then do not display
            # a fit
            return False

        lines = None

        if app_or_ret == "app" and self.curve_approach is not None:
            x = self.curve_approach[0][self.valid_seg_app]
        elif app_or_ret == "ret" and self.curve_retraction is not None:
            x = self.curve_retraction[0][self.valid_seg_ret]
        else:
            x = None

        if x is not None:
            [a, b, berror] = params

            lines = []

            if self.plot_type == "fit_preview" and mode != "force":
                lines.append(numpy.polyval([a, b], x) - self.zero)
                lines.append(numpy.polyval([a, b + berror], x) - self.zero)

            else:
                lines.append(
                    numpy.polyval([a, b], x) * self.force_factor - self.zero)
                if mode == "force":
                    # In this case berror is the value we want
                    lines.append(numpy.polyval(
                        [a, berror], x) * self.force_factor - self.zero)
                else:
                    # Deflexion - extension
                    lines.append(numpy.polyval(
                        [a, b + berror], x) * self.force_factor - self.zero)

        if app_or_ret == "app":
            self.fits_lines_poc = lines
        elif app_or_ret == "ret":
            self.fits_lines_joc = lines

    def get_params(self, params, mode=None):
        """Get the coefficients of the fits.

        If no coefficients are given the method returns None.
        """
        if list(params) != [0, 0, 0] and params is not None:
            if mode == "force":
                # Get the coefficients of the fit
                a = params[0]
                b = params[1]
                berror = params[2]

                # Calculate new coefficients for the fits.
                # berror is the fit passing through the JOC,
                # where indentation = 0 and force = 0, so berror = 0 at the end
                a = (a * self.spring_constant) / (1 - a)
                berror = berror * self.spring_constant
                b = 0 - berror
                berror = 0

                return [a, b, berror]

            else:
                return params

        else:
            return None

    def get_dist_between_jocs(self):
        """Get the distance between joc1 and joc2.

        Used for the plotting of the force curve in the results tab.
        """
        # Find distance between two jocs
        joc1 = self.data.jocs1_real[self.xpos][self.ypos]
        joc2 = self.data.jocs2_real[self.xpos][self.ypos]

        joc1_real_force_x = joc1[0] - joc1[1]
        joc1_real_force_y = joc1[1] * self.spring_constant
        joc2_real_force_x = joc2[0] - joc2[1]
        joc2_real_force_y = joc2[1] * self.spring_constant
        distx = abs(joc1_real_force_x - joc2_real_force_x)
        disty = abs(joc1_real_force_y - joc2_real_force_y)

        # if disty > 0:
        #    disty = -disty

        return [distx, disty]

    def setup_segments(self):
        """Get x,y coordinates of the points segmenting the approach curve.

        Most of the code here is the same as in :
        compute_tools.stiffness.stiffness_compute
        """
        # Display the segments if asked (or if needed for the fits)
        if self.data.display_segments or self.data.display_fits_stiffness:
            curve_x = self.curve_approach[0]
            curve_y = self.curve_approach[1]

            self.segments = []

            parts_x = stiffness_compute.get_x_parts(
                self.data.used_indentation_start,
                self.data.used_indentation_stop,
                self.data.used_indentation_step,
                max(curve_x))

            # The first segment will be between indice 0 and 1
            # else it will be between POC and POC + 1
            if self.curve_type == "indentation":
                newstart = 2
            else:
                # In case of a force_distance curve, the curve_x given has been
                # sliced, so we have to switch back the indices
                # ( == -approach_positions)
                xpos = self.xpos
                ypos = self.ypos
                ind = self.data.pocs_indices[xpos][ypos]
                newstart = int(ind + 1)

            # Look for the segments
            found_last = False
            for x in parts_x:
                found_last = False
                for i in range(newstart, len(curve_x)):
                    if curve_x[i] >= x:
                        coefficients = numpy.polyfit(
                            curve_x[i - 1:i + 1], curve_y[i - 1:i + 1], 1)
                        poly = numpy.poly1d(coefficients)
                        self.segments.append(
                            [x * self.x_factor, poly(x) * self.force_factor])
                        newstart = i
                        found_last = True
                        break
            if not found_last:  # Add last point if didn't find
                self.segments.append([curve_x[-1] * self.x_factor,
                                      curve_y[-1] * self.force_factor])

            if not self.data.used_indentation_start:
                self.segments[0] = [0, 0]  # For PoC being interpolated

    def get_stiffness_fits(self, curve_x):
        """Get the x,y values of the fits for the stiffness, for each segment."""
        dt = self.data

        # Display the fits for the stiffness
        if dt.display_fits_stiffness:
            self.segments_fits = []

            for i in range(len(self.segments) - 1):
                fitcurve_x = []
                positions_x = []
                fitcurve_y = []
                xoffset = 0
                yoffset = 0

                if dt.used_tomography:
                    # Tomography
                    xoffset = self.segments[i][0]
                    yoffset = self.segments[i][1]
                    length = dt.used_indentation_step + 1
                elif dt.used_indentation_step:
                    # Start all segments at zero
                    length = (i + 1) * dt.used_indentation_step + 1
                else:
                    # One segment, to max
                    length = self.segments[-1][0]

                # Get a value each nm
                for j in range(int(length)):
                    fitcurve_x.append(j)
                    posx = xoffset + j
                    positions_x.append(posx)

                # Add the last position
                if dt.used_indentation_step == 0:
                    # Starts anyway at 0
                    fitcurve_x.append(length)
                    posx = xoffset + length
                    positions_x.append(posx)

                E = dt.stiffness_array[self.xpos][self.ypos][i]

                if dt.used_stiffness_model_selected == 0:
                    # Hertz (sphere)
                    a = 1.0 - dt.used_poisson_ratio ** 2
                    b = math.sqrt(dt.used_tip_radius * 10 ** -9)
                    coeff = (3.0 / 4.0) * (a / b)
                elif dt.used_stiffness_model_selected == 1:
                    # Sneddon (cone)
                    a = 2.0 * math.tan(math.radians(dt.used_tip_angle))
                    coeff = (math.pi * (1.0 - dt.used_poisson_ratio ** 2)) / a
                elif dt.used_stiffness_model_selected == 2:
                    # Bilodeau (pyramid)
                    a = math.sqrt(2) * (1.0 - dt.used_poisson_ratio ** 2)
                    coeff = a / math.tan(math.radians(dt.used_tip_angle))
                elif dt.used_stiffness_model_selected == 3:
                    # Slope
                    coeff = 1.0
                if dt.used_stiffness_model_selected == 4:
                    # Flat punch
                    a = 1.0 - dt.used_poisson_ratio ** 2
                    b = dt.used_tip_radius * 10 ** -9
                    coeff = (1 / 2) * (a / b)

                yo = yoffset

                for x in fitcurve_x:
                    if dt.used_stiffness_model_selected == 0:
                        fitcurve_y.append(
                            ((E * (x * 1e-9) ** 1.5) / coeff) * 1e9 + yo)
                    elif dt.used_stiffness_model_selected in (1, 2):
                        fitcurve_y.append(
                            (((E * (x * 1e-9) ** 2) / coeff)) * 1e9 + yo)
                    elif dt.used_stiffness_model_selected in (3,4):
                        fitcurve_y.append(((E * x * 1e-9) / coeff) * 1e9 + yo)

                y = numpy.array(fitcurve_y) * self.force_factor

                self.segments_fits.append(
                    [numpy.array(positions_x) * self.x_factor, y])

    def get_stiffness_lm_fits(self, curve_x=None):
        """Get the x,y values of the fits for the stiffness, for each segment."""

        dt = self.data

        # Display the fits for the stiffness
        if dt.display_fits_stiffness or dt.display_curve_type == "residuals":

            self.segments_fits = []
            self.residuals = []
            self.chisqr = None

            app, _, _, _, _, _ = \
                curve_tools.get_curve(dt, [self.xpos, self.ypos],
                                      mode="curve_results")  # get approach curve data

            poc_pos = dt.pocs_real[self.xpos][self.ypos]
            ind = dt.pocs_indices[self.xpos][self.ypos]

            k = dt.spring_constant * 1e9

            fcapp = curve_tools.get_force_curves(
                app[0], app[1], poc_pos, k)  # Point of contact position is subtracted to set F = 0 and zheight = 0
            app[0] = app[0] * 1e-09  # scaling zheight
            app[1] = app[1] * 1e-09  # scaling deflection
            fcapp[0] = fcapp[0] * 1e-09  # scaling force data
            fcapp[1] = fcapp[1] * 1e-09  # scaling indentation

            if dt.used_stiffness_model_selected == 0:
                # Hertz (sphere)
                fitforce = self.hertz_sphere
                div = numpy.sqrt(dt.used_tip_radius * 1e-9)
                coeff = (4.0 / 3.0) * div / (1.0 - dt.used_poisson_ratio ** 2)
            elif dt.used_stiffness_model_selected == 1:
                # Sneddon (cone)
                fitforce = self.sneddon_cone
                div = numpy.tan(dt.used_tip_angle * numpy.pi / 180)
                coeff = (2.0 / numpy.pi) * div / (1.0 - dt.used_poisson_ratio ** 2)
            elif dt.used_stiffness_model_selected == 2:
                # Bilodeau (pyramid)
                fitforce = self.bilodeau_pyramid
                div = numpy.tan(dt.used_tip_angle * numpy.pi / 180)
                coeff = (1.0 / numpy.sqrt(2)) * div / (1.0 - dt.used_poisson_ratio ** 2)
            if dt.used_stiffness_model_selected == 4:
                # Flat punch
                fitforce = self.flat_punch
                div = dt.used_tip_radius * 1e-9
                coeff = 2 * div / (1.0 - dt.used_poisson_ratio ** 2)

            E = dt.stiffness_array[self.xpos][self.ypos][0]

            delta = fcapp[0][ind::]
            data_to_fit = fcapp[1][ind::]

            # Display elasticity fit within the indentation range selected by user.
            if dt.used_indentation_stop != 0 and dt.fit_range_type == 1:
                start_f = dt.used_indentation_start * 1e-09  # 0 is the default value.
                stop_f = dt.used_indentation_stop * 1e-09

                delta = numpy.clip(delta, a_min=start_f, a_max=stop_f)

            ytemp = fitforce(E, coeff, delta, data_to_fit)

            self.chisqr = self.get_chisqr(ytemp)

            self.residuals.append([delta * 1e09 * self.x_factor, ytemp * 1e09 * self.force_factor])

            y_offset = ytemp[0]
            ytemp = ytemp - y_offset

            yf = data_to_fit - ytemp

            if dt.used_indentation_start == 0:
                x_offset = delta[0]
                delta = delta - x_offset

            self.segments_fits.append(
                [delta * 1e09 * self.x_factor, yf * 1e09 * self.force_factor])


    def hertz_sphere(self, e0, coeff, delta, data):

        y = numpy.zeros(delta.shape)
        for i in range(len(y)):
            y[i] = coeff * e0 * numpy.power((delta[i]), 1.5)

        return data - y

    def sneddon_cone(self, e0, coeff, delta, data):

        y = numpy.zeros(delta.shape)
        for i in range(len(y)):
            y[i] = coeff * e0 * numpy.power((delta[i]), 2)

        return data - y

    def bilodeau_pyramid(self, e0, coeff, delta, data):

        y = numpy.zeros(delta.shape)
        for i in range(len(y)):
            y[i] = coeff * e0 * numpy.power((delta[i]), 2)

        return data - y

    def flat_punch(self, e0, coeff, delta, data):

        y = numpy.zeros(delta.shape)
        for i in range(len(y)):
            y[i] = coeff * e0 * numpy.power((delta[i]), 1)

        return data - y

    def get_chisqr(self, r):
        """
        Chi-square: :math:`\chi^2 = \sum_i^N [{\rm Resid}_i]^2`
        """
        return (r**2).sum()
