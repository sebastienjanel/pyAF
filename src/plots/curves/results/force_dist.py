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

"""Plots the Force - Distance curve."""

import numpy
import copy
from ..plot_curve import PlotCurve
from ....tools import math_tools
from ....tools import breakthrough_tools
from ....tools import curve_tools


class ResultsForceDist(PlotCurve):
    """Plots the Force - Distance curve."""

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

        disp = self.data.display_trace_retrace
        xpos = self.xpos
        ypos = self.ypos

        # Labels
        if self.data.curve_distance_units == "nm":
            self.xlabel = "Distance [nm]"
        elif self.data.curve_distance_units == "um":
            self.xlabel = "Distance [\u03bcm]"
        self.ylabel = "Force [" + self.data.curve_force_units + "]"

        app_calc = self.data.stiffness_calculated
        ret_calc = self.data.work_and_rupture_force1_calculated

        # Display the curve (defined below, can be skipped if there is no
        # data to display).
        skip = False

        if (disp == 0 or disp == 2) and self.approach is not None and app_calc:
            self.curve_approach = curve_tools.get_force_curves(
                self.approach[0],
                self.approach[1],
                self.data.pocs_real[self.xpos][self.ypos],
                self.spring_constant)

        elif (disp == 0 or disp == 2) and self.approach is not None and \
                self.data.breakthrough_calculated:
            # Without the elasticity (no point of contact): distance from
            # the support, force from the baseline (breakthroughs)
            self.curve_approach = list(breakthrough_tools.force_distance(
                self.approach[0], self.approach[1], self.spring_constant))
            self.xlabel = self.xlabel.replace("Distance", "Distance to "
                                              "support")

        if (disp == 1 or disp == 2) and \
                self.retraction is not None and ret_calc:
            self.curve_retraction = curve_tools.get_force_curves(
                self.retraction[0],
                self.retraction[1],
                self.data.jocs1_real[self.xpos][self.ypos],
                self.spring_constant)

        if self.curve_approach is not None and app_calc:
            # The real fit is done on the curve without the corrupted part,
            # so we have to remove it from the curve to display the fit
            tr = [None, None]
            tr[0] = self.curve_approach[0][self.valid_seg_app]
            tr[1] = self.curve_approach[1][self.valid_seg_app]

            # For incomplete files, there is no curve to display
            if len(tr[0]) == 0:
                skip = True

            poc = self.data.pocs_real[xpos][ypos]

            # Display the point of contact (of the fit, on the baseline)
            if self.data.display_poc and skip is False:
                # On the force curve the poc is at 0,0
                self.poc_pos = [0.0, 0.0]
                # Detected point of contact, if the fit moved it
                detected = self.data.poc_detected(xpos, ypos)
                if detected is not None and \
                        abs(detected[0] - float(poc[0])) > 1e-3:
                    x, y = self.force_point(detected, poc)
                    self.poc_detected_pos = [x * self.x_factor,
                                             y * self.force_factor]

            # Display the fit used to determine the point of contact (the
            # baseline and the threshold above its noise)
            if self.data.display_fit_poc and skip is False:
                value = self.data.fits_poc[xpos][ypos]
                params = self.force_line_params(value, poc)
                self.setup_lines(params, "app", "force")

            # Get the fits for each segment
            if skip is False and self.data.used_perform_fit is False:
                self.setup_segments()

            # Get the fits for each segment (pass tr[0] in case there is only
            # one segment, we need the size of it)
            if skip is False:
                if self.data.used_perform_fit:
                    self.get_stiffness_lm_fits(tr[0])
                    self.fits_to_force_distance()
                else:
                    self.get_stiffness_fits(tr[0])

        if self.curve_retraction is not None:
            # Display the point of jump of contact
            if self.data.display_joc:
                self.joc1_pos = [0.0, 0.0]
                self.joc2_pos = self.get_dist_between_jocs()
                self.joc2_pos[0] *= self.x_factor
                self.joc2_pos[1] *= self.force_factor

            # Display the fit used to determine the point of jump of contact
            if self.data.display_fit_joc:
                value = self.data.fits_joc[xpos][ypos]
                params = self.get_params(value, mode="force")
                if params is not None:
                    # Display the noise line below the fit line
                    params[1] = -params[1]

                self.setup_lines(params, "ret", "force")

            # Display the maximum force: largest distance between the curve
            # and the line between the two jumps of contact (the curve is
            # centered on the first one)
            if self.data.display_force:
                value = self.data.fits_joc[xpos][ypos]
                params = self.get_params(value, mode="force")

                if params is not None:
                    joc1_indice = self.data.jocs1_indices[xpos][ypos]
                    joc2_indice = self.data.jocs2_indices[xpos][ypos]
                    segment = slice(joc1_indice, joc2_indice)
                    segx = self.curve_retraction[0][segment]
                    segy = self.curve_retraction[1][segment]
                    [distx, disty] = self.get_dist_between_jocs()
                    a, b = math_tools.find_line_coefs(
                        [0.0, 0.0], [distx, disty])

                    if len(segy) and a is not None:
                        line = numpy.polyval([a, b], segx)
                        k = numpy.argmax(line - segy)
                        self.max_force = [
                            [segx[k] * self.x_factor,
                             segx[k] * self.x_factor],
                            [segy[k] * self.force_factor,
                             line[k] * self.force_factor]]

            # Display the surface used for the work calculations
            if self.data.display_surface:
                value = self.data.fits_joc[xpos][ypos]
                params = self.get_params(value, mode="force")

                if params is not None:
                    # Use deepcopy, these are lists, else it will mess up the
                    # fits and jocs values ...
                    joc1_indice = self.data.jocs1_indices[xpos][ypos]
                    joc2_indice = self.data.jocs2_indices[xpos][ypos]
                    joc1_ind = copy.deepcopy(joc1_indice)
                    joc2_ind = copy.deepcopy(joc2_indice + 1)
                    segment = slice(joc1_ind, joc2_ind)
                    segx = copy.deepcopy(self.curve_retraction[0])
                    segy = copy.deepcopy(self.curve_retraction[1])
                    segx = segx[segment]
                    segy = segy[segment]
                    segx[0] = 0.0
                    segy[0] = 0.0

                    [distx, disty] = self.get_dist_between_jocs()

                    segx[len(segx) - 1] = distx
                    segy[len(segy) - 1] = disty

                    a, b = math_tools.find_line_coefs(
                        [0.0, 0.0], [distx, disty])
                    y = numpy.polyval([a, b], segx)

                    self.surface_params = [
                        segx * self.x_factor,
                        segy * self.force_factor,
                        y * self.force_factor]

        if app_calc and self.curve_approach is not None and \
                self.curve_retraction is not None:
            # Both curves are displayed: use the same origin for both (the
            # point of contact), else the retraction curve (centered on the
            # jump of contact) is offset from the approach curve
            self.center_retraction_on_poc()

        self.define_x_scale()
        self.define_y_scale()

        self.plot_fig()

    def fits_to_force_distance(self):
        """Draw the elasticity fit on the force - distance curve.

        The fit is computed on the indentation curve (baseline slope
        removed, zero of the force on the baseline at the point of contact,
        see stiffness.bec.approach_points), while this curve is the raw
        approach centered on the point of contact: the points of the fit are
        moved back to the deflection - Z curve, then to the force - distance
        curve as the data.
        """
        if not self.segments_fits:
            return
        self.segments_fits = self.to_force_distance(self.segments_fits)
        if self.segments_fits_extra:
            self.segments_fits_extra = self.to_force_distance(
                self.segments_fits_extra)

    def to_force_distance(self, segments):
        """Fit segments of the indentation curve (displayed units) on the
        force - distance curve (see fits_to_force_distance)."""
        dt = self.data
        i, j = self.xpos, self.ypos
        z_poc, d_poc_raw = [float(v) for v in dt.pocs_real[i][j][:2]]
        slope, intercept = [float(v) for v in dt.fits_poc[i][j][:2]]
        if not numpy.isfinite(slope):
            slope = 0.0
        d_poc = d_poc_raw
        if numpy.isfinite(intercept):
            d_poc = slope * z_poc + intercept
        k = dt.used_spring_constant
        if k is None:
            k = dt.spring_constant
        fits = []
        for fit_x, fit_f in segments:
            # Indentation [m] and force [N] of the fit
            delta = numpy.asarray(fit_x) / (1e09 * self.x_factor)
            force = numpy.asarray(fit_f) / (1e09 * self.force_factor)
            # Deflection - Z curve [nm] (inverse of approach_points)
            defl = force / k
            z = z_poc + delta * 1e09 + defl
            d = defl + d_poc + slope * (z - z_poc)
            # Force - distance curve as the data (get_force_curves)
            x, y = curve_tools.get_force_curves(
                z, d, [z_poc, d_poc_raw], self.spring_constant)
            fits.append([x * self.x_factor, y * self.force_factor])
        return fits

    def center_retraction_on_poc(self):
        """Move the retraction curve and its elements to the poc origin.

        The retraction curve and its elements are defined with the jump of
        contact at 0, 0. Shift them so that the point of contact is at 0, 0,
        as for the approach curve (the curves are then at the same relative
        position as on the deflection - Z curve).
        """
        poc = self.data.pocs_real[self.xpos][self.ypos]
        joc = self.data.jocs1_real[self.xpos][self.ypos]

        # Shift in nm and nN (same computation as in get_force_curves)
        dx = (joc[0] - joc[1]) - (poc[0] - poc[1])
        dy = (joc[1] - poc[1]) * self.spring_constant

        self.curve_retraction[0] = self.curve_retraction[0] + dx
        self.curve_retraction[1] = self.curve_retraction[1] + dy

        # The other elements are already in the displayed units
        dx_disp = dx * self.x_factor
        dy_disp = dy * self.force_factor

        if self.joc1_pos is not None:
            self.joc1_pos = [self.joc1_pos[0] + dx_disp,
                             self.joc1_pos[1] + dy_disp]
            self.joc2_pos = [self.joc2_pos[0] + dx_disp,
                             self.joc2_pos[1] + dy_disp]

        if self.fits_lines_joc is not None:
            # The lines are plotted against the x values of the curve
            self.fits_lines_joc = [line + dy_disp
                                   for line in self.fits_lines_joc]

        if self.max_force is not None:
            linex, liney = self.max_force
            self.max_force = [[x + dx_disp for x in linex],
                              [y + dy_disp for y in liney]]

        if self.surface_params is not None:
            segx, segy, y = self.surface_params
            self.surface_params = [segx + dx_disp, segy + dy_disp,
                                   y + dy_disp]
