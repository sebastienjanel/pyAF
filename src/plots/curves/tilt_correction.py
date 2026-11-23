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

"""Plots the curve for the tilt correction."""

import numpy

from .plot_curve import PlotCurve
from ...tools import misc_tools
from ...tools import tilt_tools
from ... import shared


class PlotTiltPreview(PlotCurve):
    """Preview for the tilt correction."""

    def __init__(self, parent):
        super().__init__(parent)

        # Labels
        self.xlabel = misc_tools.z_label(self.data)
        self.ylabel = "Deflection [nm]"

        self.curve_approach = self.approach
        self.curve_retraction = self.retraction

        # Limits for the tilt correction; in automatic mode, the baseline
        # found on the curve displayed: on the curve the correction is
        # fitted on, before: on the retraction for the adhesion and the
        # events, otherwise on the approach
        if not self.data.tilt_auto:
            self.curve_mod_limit_1 = self.data.tilt_limit_1
            self.curve_mod_limit_2 = self.data.tilt_limit_2
            self.color_curve_mod_limit = "r"
        else:
            retraction = shared.exp.compute_type in (
                "work_and_rupture_force", "events")
            if self.data.tilt_applied:
                retraction = self.data.tilt_applied == "retrace"
            if retraction:
                curve, start = self.curve_retraction, self.ret_pos[0]
            else:
                curve, start = self.curve_approach, self.app_pos[0]
            found = None
            if curve is not None:
                # Corrected curve: its baseline is flat
                slope = 0.0 if self.data.tilt_applied else None
                found = tilt_tools.baseline(curve[0][start:],
                                            curve[1][start:], slope)
            if found is not None:
                z = curve[0][start:]
                z = z[numpy.isfinite(z)]
                self.curve_mod_limit_1 = z[found[2]]
                self.curve_mod_limit_2 = z[found[3] - 1]
                self.color_curve_mod_limit = "g"

        self.plot_fig()
