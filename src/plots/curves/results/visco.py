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

"""Plots the viscoelastic fit (Ting model) of a curve."""

import numpy

from ..plot_curve import PlotCurve
from .... import shared
from ....compute_tools.viscoelasticity import main as visco
from ....compute_tools.viscoelasticity import ting
from ....gui_styles.theme_handler import curve_color


class ResultsVisco(PlotCurve):
    """Force vs indentation (distance from the point of contact of the
    elasticity fit) of the approach and the retraction (fitted part,
    downsampled as for the fit), with the fitted Ting model (computed in
    time, displayed vs the indentation).

    The curve is fitted again with the parameters used for the map: the
    values are the ones of the map. In the compute tab (fit preview), the
    curve is fitted with the current parameters.
    """

    def __init__(self, parent):
        super().__init__(parent)

        if self.data.curve_distance_units == "um":
            self.xlabel = "Indentation [\u03bcm]"
            self.x_factor = 1e6
        else:
            self.xlabel = "Indentation [nm]"
            self.x_factor = 1e9
        self.ylabel = "Force [" + self.data.curve_force_units + "]"
        self.force_factor = 1e9 if self.data.curve_force_units == "nN" \
            else 1e12

        self.used = self.plot_type != "fit_preview"
        self.fit = None
        self.message = "No viscoelastic fit for this curve"
        if not self.used:
            reason = visco.check_file(self.data)
            if reason is None:
                self.fit = visco.fit_curve(self.data, self.xpos, self.ypos,
                                           used=False)
            else:
                self.message = "Preview not available:\n" + reason
        elif self.data.visco_calculated:
            self.fit = visco.fit_curve(self.data, self.xpos, self.ypos)

        self.plot_fig()

    def plot_fig(self):
        """Plots the force, the fit and the fitted values."""
        if self.fit is None:
            self.axes.text(0.5, 0.5, self.message,
                           transform=self.axes.transAxes,
                           ha="center", va="center")
        else:
            inputs, result = self.fit
            time, delta, force, idx_tm, tc_init = inputs[2:7]
            has_pause, height = inputs[8], inputs[10]
            x = delta * self.x_factor
            f = force * self.force_factor
            self.axes.plot(x[:idx_tm + 1], f[:idx_tm + 1],
                           color=curve_color("approach"), label="Approach",
                           **self.curve_style())
            self.axes.plot(x[idx_tm:], f[idx_tm:],
                           color=curve_color("retraction"), label="Retract",
                           **self.curve_style())
            if result["success"]:
                settings = visco.get_settings(self.data, used=self.used)
                common = visco.get_model(self.data, settings)
                solver = "numerical" if common["solver"] == 1 or has_pause \
                    else "analytical"
                model = ting.TingModel(
                    common["model"], common["coeff"], common["exponent"],
                    solver, tip_radius=common["tip_radius"],
                    tip_angle=common["tip_angle"],
                    height=height if numpy.isfinite(height) else None)
                curve = ting.TingCurve(time, delta, force, idx_tm)
                velocities = ting.get_velocities(curve, tc_init) \
                    if solver == "analytical" else None
                fitted = model.evaluate(curve, result["e0"], result["tc"],
                                        result["beta"], result["f0"],
                                        velocities)
                self.axes.plot(x, fitted * self.force_factor, "r",
                               label="Ting fit")
                # Point of contact of the Ting fit
                x_tc = numpy.interp(result["tc"], time, delta) * \
                    self.x_factor
                self.axes.axvline(x_tc, color="grey", lw=0.8, ls="--")
                text = "E0 = " + format(result["e0"] * 1e-3, ".3g") + \
                    " kPa\nFluidity = " + format(result["beta"], ".3f") + \
                    "\ntan δ = " + format(result["tan_delta"], ".3f") + \
                    "\nR² = " + format(result["r2"], ".3f")
                if numpy.isfinite(height):
                    text += "\nB.E.C. h = " + format(height * 1e9, ".0f") + \
                        " nm"
                if shared.exp.display_curve_legend:
                    # Values in the legend (no overlap with it)
                    for line in text.split("\n"):
                        self.axes.plot([], [], " ", label=line)
                else:
                    self.axes.text(0.02, 0.97, text,
                                   transform=self.axes.transAxes,
                                   ha="left", va="top", fontsize=9)
            else:
                self.axes.text(0.98, 0.97, "The fit failed",
                               transform=self.axes.transAxes,
                               ha="right", va="top", fontsize=9)
            if shared.exp.display_curve_legend:
                self.axes.legend(loc=2, prop=self.font)

        self.plot_labels()
        self.red_line()
        self.start_plot()
