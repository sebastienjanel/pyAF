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

"""Plots the indentation curve."""

from ..plot_curve import PlotCurve
from ....tools import curve_tools


class ResultsIndentation(PlotCurve):
    """Plots the Force - Indentation curve."""

    def __init__(self, parent):
        super().__init__(parent)

        # Labels
        if self.data.curve_distance_units == "nm":
            self.xlabel = "Indentation [nm]"
        elif self.data.curve_distance_units == "um":
            self.xlabel = "Indentation [\u03bcm]"
        self.ylabel = "Force [" + self.data.curve_force_units + "]"

        if self.data.curve_force_units == "nN":
            self.force_factor = 1
        elif self.data.curve_force_units == "pN":
            self.force_factor = 1e3
        if self.data.curve_distance_units == "nm":
            self.x_factor = 1
        elif self.data.curve_distance_units == "um":
            self.x_factor = 1e-3

        if self.approach is not None:
            force_curve = curve_tools.get_force_curves(
                self.approach[0],
                self.approach[1],
                self.data.pocs_real[self.xpos][self.ypos],
                self.spring_constant)

            pos = self.data.pocs_indices[self.xpos][self.ypos]
            segment = slice(pos, len(force_curve[0]))

            # poc_indice is only the indice of the value before the real poc,
            # so we have to replace the first value of the curve with the real
            # poc value (0, 0)
            curve_x = force_curve[0][segment]
            curve_y = force_curve[1][segment]
            curve_x[0] = 0.0
            curve_y[0] = 0.0

            self.curve_approach = [curve_x, curve_y]

            # Get the fits for each segment
            if self.data.perform_fit is False:
                self.setup_segments()

            # Get the fits for each segment (pass curve_x in case there is only
            # one segment, we need the size of it)

            if self.data.perform_fit:
                self.get_stiffness_lm_fits(curve_x)

                if self.data.fit_range_type == 0:
                    setpoint = self.data.trig_threshold * self.data.deflection_sensitivity * self.data.spring_constant * 1E9
                    self.max_force_range = self.data.force_stop/100 * setpoint * self.force_factor
                    self.min_force_range = self.data.force_start/100 * setpoint * self.force_factor

                elif self.data.fit_range_type == 1:

                    if self.data.indentation_stop > max(curve_x):
                        self.max_indentation = max(curve_x)

                    else:
                        self.max_indentation = self.data.indentation_stop * self.x_factor

                    if self.data.indentation_start != 0:
                        self.min_indentation = self.data.indentation_start * self.x_factor

            else:
                self.get_stiffness_fits(curve_x)

        self.define_x_scale()
        self.define_y_scale()

        self.plot_fig()
