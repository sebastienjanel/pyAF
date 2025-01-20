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

from ..plot_curve import PlotCurve
from ....compute_tools import loading_rates
from ....tools import curve_tools


class ComputeLoadingRates(PlotCurve):
    """Plots the Force - Distance curve with the loading rates for the compute
    tab.
    """

    def __init__(self, parent):
        super().__init__(parent)

        # Labels
        self.xlabel = "Distance [\u03bcm]"
        self.ylabel = "Force [nN]"

        if self.data.events_calculated:
            # Use the same spring constant as the one for the events
            # computation.
            self.curve_retraction = curve_tools.get_force_curves(
                self.retraction[0],
                self.retraction[1],
                self.data.events_jocs2_real[self.xpos][self.ypos],
                self.data.used_spring_constant)

            # Force in nN
            self.curve_retraction[1] = self.curve_retraction[1] * 1e9

            # Get the fit for the loading rates. The speed is set to 1, as
            # we are not interested in the real value here. We want only the
            # fit values.
            _, fit_boundaries = \
                loading_rates.get_loading_rate_for_curve(
                    self.curve_retraction,
                    self.data.events_positions_start[self.xpos][self.ypos],
                    self.data.events_positions_stop[self.xpos][self.ypos],
                    1,
                    self.data.lr_coef)

            self.loading_rates_fits = fit_boundaries

        self.plot_fig()
