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

"""Plots the force - extension curve"""

from ..plot_curve import PlotCurve


class ResultsForceExt(PlotCurve):
    """Plots the force - extension curve.

    This plot is not much used, so it displays only the curve for the moment.
    It could be extended to display more stuff (poc, joc, ...).
    """

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

        fac = self.spring_constant * self.force_factor
        disp = self.data.display_trace_retrace

        # Show or hide the curves
        if disp == 0 or disp == 2:
            self.curve_approach = [self.approach[0], self.approach[1] * fac]

        if disp == 1 or disp == 2:
            ret = self.retraction
            self.curve_retraction = [ret[0], ret[1] * fac]

        # Set the x axis to microns
        self.define_x_scale()

        self.plot_fig()
