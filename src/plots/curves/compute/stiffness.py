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

"""Plots the stiffness fit preview."""

from ..plot_curve import PlotCurve
from ....tools import search_poc


class ComputeStiffness(PlotCurve):
    """Plots the stiffness fit preview."""

    def __init__(self, parent):
        super().__init__(parent)

        # Labels
        self.xlabel = "Scanner Extension [nm]"
        self.ylabel = "Deflection [nm]"

        self.curve_approach = self.approach

        # Get a fit and a poc
        poc_indice, poc, fit, fit_segment = search_poc.get_POC(
            self.curve_approach, self.app_pos, self.data.fit_params_poc)
        params = self.get_params(fit)

        if params is not None:
            # Set the point of contact
            self.poc_pos_indice = poc_indice
            self.poc_pos = poc

            self.setup_lines(params, "app", "fit_preview")

            # Segment used for the fit
            self.fit_segment = slice(fit_segment[0], fit_segment[1] + 1)

        else:
            self.no_fit_found = True

        # Setup the values for the noise of the fits
        self.setup_fit_noise("approach", self.fit_segment)

        self.plot_fig()
