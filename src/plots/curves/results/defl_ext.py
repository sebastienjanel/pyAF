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

"""Plots the curve with stiffness features."""

import numpy
import copy
from ..plot_curve import PlotCurve


class ResultsDeflExt(PlotCurve):
    """Plots the curve with stiffness features."""

    def __init__(self, parent):
        super().__init__(parent)

        if self.data.curve_distance_units == "nm":
            self.x_factor = 1
        elif self.data.curve_distance_units == "um":
            self.x_factor = 1e-3

        # Labels
        if self.data.curve_distance_units == "nm":
            self.xlabel = "Scanner Extension [nm]"
        elif self.data.curve_distance_units == "um":
            self.xlabel = "Scanner Extension [\u03bcm]"
        self.ylabel = "Deflection [nm]"

        # Show or hide the curves
        if self.data.display_trace_retrace == 0:
            self.curve_approach = self.approach
        elif self.data.display_trace_retrace == 1:
            self.curve_retraction = self.retraction
        elif self.data.display_trace_retrace == 2:
            self.curve_approach = self.approach
            self.curve_retraction = self.retraction

        params_fit_poc = None
        params_fit_joc = None

        if self.data.stiffness_calculated and self.curve_approach is not None:
            params_fit_poc = self.get_params(
                self.data.fits_poc[self.xpos][self.ypos])

        calculated = self.data.work_and_rupture_force1_calculated
        if calculated and self.curve_retraction is not None:
            params_fit_joc = self.get_params(
                self.data.fits_joc[self.xpos][self.ypos])

        if params_fit_poc is not None:
            # Point of contact
            if self.data.display_poc:
                self.poc_pos_indice = \
                    self.data.pocs_indices[self.xpos][self.ypos]
                self.poc_pos = self.data.pocs_real[self.xpos][self.ypos]
                self.poc_pos[1] = self.poc_pos[1] - self.zero
                self.poc_pos[0] *= self.x_factor

            # Display the fit used to determine the point of contact
            if self.data.display_fit_poc:
                self.setup_lines(params_fit_poc, "app", "defl_ext")

        if params_fit_joc is not None:
            joc1_pos = self.data.jocs1_real[self.xpos][self.ypos]
            joc2_pos = self.data.jocs2_real[self.xpos][self.ypos]
            # Correct position of the joc
            joc1_pos[1] = joc1_pos[1] - self.zero
            joc2_pos[1] = joc2_pos[1] - self.zero

            # Display the point of jump of contact (work and rupture force)
            # A joc is also needed when displaying the surface
            if self.data.display_joc:
                self.joc1_pos = joc1_pos
                self.joc2_pos = joc2_pos
                self.joc1_pos[0] *= self.x_factor
                self.joc2_pos[0] *= self.x_factor

            # Display the fit used to determine the point of jump of contact
            # (work and rupture force)

            if self.data.display_fit_joc:
                self.setup_lines(params_fit_joc, "ret", "defl_ext")

            # Display the surface used for the work calculations
            if self.data.display_surface:
                joc1_ind = self.data.jocs1_indices[self.xpos][self.ypos]
                joc2_ind = self.data.jocs2_indices[self.xpos][self.ypos]

                [a, b, _] = params_fit_joc

                a /= self.x_factor

                # Use deepcopy, these are lists, else it will mess up the
                # fits and jocs values ...
                segment = slice(joc1_ind, joc2_ind)
                segx = copy.deepcopy(
                    self.retraction[0][segment]) * self.x_factor
                segy = copy.deepcopy(self.retraction[1][segment])
                segx[0] = copy.deepcopy(joc1_pos[0])
                segy[0] = copy.deepcopy(joc1_pos[1])
                segx[len(segx) - 1] = copy.deepcopy(joc2_pos[0])
                segy[len(segy) - 1] = copy.deepcopy(joc2_pos[1])

                y = numpy.polyval([a, b], segx)

                self.surface_params = [segx, segy, y - self.zero]

        # Set the x axis to microns
        self.define_x_scale()

        self.plot_fig()
