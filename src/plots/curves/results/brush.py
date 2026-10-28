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


"""Plots the brush model of a curve."""

import numpy

from ..plot_curve import PlotCurve
from .... import shared
from ....compute_tools.stiffness import brush
from ....gui_styles.theme_handler import curve_color


class ResultsBrush(PlotCurve):
    """Force - distance curve of the approach (zero of the force on the
    baseline of the point of contact) with the two fits of the brush model:
    the cell body (Hertz type model) fitted above the threshold, and the
    brush, its exponential force F = A exp(-2 pi h / L) at the distance h
    from the body, drawn vs the distance.

    The curve is analysed again with the parameters used for the maps: the
    values are the ones of the maps.
    """

    def __init__(self, parent):
        super().__init__(parent)

        if self.data.curve_distance_units == "um":
            self.xlabel = "Distance [\u03bcm]"
            self.x_factor = 1e6
        else:
            self.xlabel = "Distance [nm]"
            self.x_factor = 1e9
        self.ylabel = "Force [" + self.data.curve_force_units + "]"
        self.force_factor = 1e9 if self.data.curve_force_units == "nN" \
            else 1e12

        self.result = None
        if self.data.brush_calculated:
            self.result = brush.curve(self.data, self.xpos, self.ypos,
                                      full=True)

        self.plot_fig()

    def plot_fig(self):
        """Plots the force - distance curve and the fits."""
        result = self.result
        if result is None or "delta" not in result:
            self.axes.text(0.5, 0.5, "No brush model for this curve",
                           transform=self.axes.transAxes, ha="center",
                           va="center")
        else:
            xf, ff = self.x_factor, self.force_factor
            delta, force = result["delta"], result["force"]
            idx = result["brush"]
            model, radius, poisson_ratio, _ = brush.settings(self.data)
            threshold = result["threshold"]
            modulus = result["modulus"]
            surface = result["surface"] * 1e-9

            self.axes.plot(delta * xf, force * ff,
                           color=curve_color("approach"), label="Approach",
                           **self.curve_style())
            # Brush points outside the window of validity of the brush force
            # (0.1 < h / L < 0.8, not fitted), and the fitted ones
            others = numpy.setdiff1d(result["candidates"], idx)
            if others.size:
                self.axes.plot(delta[others] * xf, force[others] * ff, ".",
                               ms=3, color="#9a9a9a",
                               label="Brush, outside the window")
            if idx.size:
                self.axes.plot(delta[idx] * xf, force[idx] * ff, ".", ms=3,
                               color=curve_color("retraction"),
                               label="Brush, fitted")

            # Cell body: from its surface to the end of the approach
            top = numpy.nanmax(delta)
            body_x = numpy.linspace(surface, top, 200)
            body_f = modulus * brush.body_unit_force(
                model, radius, poisson_ratio, body_x - surface)
            self.axes.plot(body_x * xf, body_f * ff, "--",
                           color="grey", label="Cell body fit")

            # Brush: distance of the probe for the forces of the brush
            # points (indentation of the body minus the distance h)
            if "slope" in result:
                # Range of the fitted points (0.1 < h / L < 0.8)
                forces = numpy.linspace(numpy.min(force[idx]),
                                        numpy.max(force[idx]), 200)
                h = (numpy.log(forces) - result["intercept"]) / \
                    result["slope"]
                brush_x = surface + brush.body_indentation(
                    model, radius, poisson_ratio, modulus, forces) - h
                self.axes.plot(brush_x * xf, forces * ff, "r-",
                               label="Brush fit")

            self.axes.axhline(threshold * ff, ls=":", lw=1,
                              color=curve_color("breakthrough"),
                              label="Threshold")

            text = ["E = " + format(modulus * 1e-3, ".3g") +
                    " kPa (cell body)"]
            if "slope" in result:
                if numpy.isfinite(result["length"]):
                    text += ["L = " + format(result["length"], ".0f") +
                             " nm",
                             "N = " + format(result["density"], ".0f") +
                             " \u00b5m\u207b\u00b2"]
                else:
                    text.append("Brush fit rejected (R\u00b2 < " +
                                format(brush.MIN_R2, "g") + ")")
                text.append("R\u00b2 = " + format(result["r2"], ".3f"))
                text.append("Fit window 0.1 < h/L < 0.8: " +
                            format(numpy.min(force[idx]) * ff, ".2g") +
                            "\u2013" +
                            format(numpy.max(force[idx]) * ff, ".2g") +
                            " " + self.data.curve_force_units)
            else:
                text.append("No brush fitted")
            if shared.exp.display_curve_legend:
                for line in text:
                    self.axes.plot([], [], " ", label=line)
                self.axes.legend(loc=2, prop=self.font)
            else:
                self.axes.text(0.02, 0.97, "\n".join(text),
                               transform=self.axes.transAxes, ha="left",
                               va="top", fontsize=9)

        self.plot_labels()
        self.red_line()
        self.start_plot()
