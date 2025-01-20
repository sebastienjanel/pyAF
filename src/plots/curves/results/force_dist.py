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

"""Plots the Force - Distance curve."""

import numpy
import copy
from ..plot_curve import PlotCurve
from ....tools import math_tools
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

        if (disp == 1 or disp == 2) and \
                self.retraction is not None and ret_calc:
            self.curve_retraction = curve_tools.get_force_curves(
                self.retraction[0],
                self.retraction[1],
                self.data.jocs1_real[self.xpos][self.ypos],
                self.spring_constant)

        if self.curve_approach is not None:
            # The real fit is done on the curve without the corrupted part,
            # so we have to remove it from the curve to display the fit
            tr = [None, None]
            tr[0] = self.curve_approach[0][self.valid_seg_app]
            tr[1] = self.curve_approach[1][self.valid_seg_app]

            # For incomplete files, there is no curve to display
            if len(tr[0]) == 0:
                skip = True

            # Display the point of contact
            if self.data.display_poc and skip is False:
                # On the force curve the poc is at 0,0
                self.poc_pos = [0.0, 0.0]

            # Display the fit used to determine the point of contact
            if self.data.display_fit_poc and skip is False:
                value = self.data.fits_poc[xpos][ypos]
                params = self.get_params(value, mode="force")
                self.setup_lines(params, "app", "force")

            # Get the fits for each segment
            if skip is False and self.data.perform_fit is False:
                self.setup_segments()

            # Get the fits for each segment (pass tr[0] in case there is only
            # one segment, we need the size of it)
            if skip is False:
                if self.data.perform_fit:
                    self.get_stiffness_lm_fits(tr[0])
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

            # Display the maximum force
            if self.data.display_force:
                value = self.data.fits_joc[xpos][ypos]
                params = self.get_params(value, mode="force")

                if params is not None:
                    [a, b, _] = params

                    # Segment curve
                    joc1_indice = self.data.jocs1_indices[xpos][ypos]
                    joc2_indice = self.data.jocs2_indices[xpos][ypos]
                    segment = slice(joc1_indice, joc2_indice)
                    new_curve_y = self.curve_retraction[1][segment]

                    max_force_pos = numpy.argmin(new_curve_y) + joc1_indice

                    xpos0 = self.curve_retraction[0][max_force_pos]
                    ypos0 = self.curve_retraction[1][max_force_pos]

                    # Get intersection with fit
                    ypos1 = xpos0 * a
                    linex = [
                        xpos0 * self.x_factor,
                        xpos0 * self.x_factor]
                    liney = [
                        ypos0 * self.force_factor,
                        ypos1 * self.force_factor]

                    self.max_force = [linex, liney]

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

        self.define_x_scale()
        self.define_y_scale()

        self.plot_fig()
