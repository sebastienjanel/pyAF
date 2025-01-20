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

"""Plots the force curve."""

from ..plot_curve import PlotCurve
from ....tools import events_tools
from ....tools import search_joc2_events
from ....tools import curve_tools


class ComputeEventsDetectionMSF(PlotCurve):
    """Show the detection with the MSF algorithm."""

    def __init__(self, parent):
        super().__init__(parent)

        # Labels
        self.xlabel = "Distance [nm]"
        self.ylabel = "MSF"

        curve_x = self.retraction[0]
        curve_y = self.retraction[1]

        self.joc2_indice, self.joc2_real, _, _ = \
            search_joc2_events.get_JOC_events(
                self.retraction,
                self.ret_pos[0],
                self.data.fit_params_joc_events)

        force_curve = curve_tools.get_force_curves(
            curve_x, curve_y, self.joc2_real, self.spring_constant)

        # Remove the corrupted part from the curve
        curve_x = force_curve[0][self.valid_seg_ret]
        curve_y = force_curve[1][self.valid_seg_ret]

        window_size = self.data.events_msf_window_size

        msf_curve_y = events_tools.apply_msf(curve_x, curve_y, window_size)
        w2 = window_size / 2

        self.curve_retraction = [curve_x[w2:-w2], msf_curve_y]

        threshold = self.data.events_msf_detection_threshold

        # Plot threshold as green line
        self.fits_lines_conv[1] = [threshold] * len(curve_x[w2:-w2])

        self.axes.set_ylim([-threshold * 5, threshold * 5])

        self.display_events_detection_exclude_square = True

        self.plot_fig()
