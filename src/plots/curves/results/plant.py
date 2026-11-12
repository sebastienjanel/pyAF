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


"""Plots the apparent stiffness of a plant cell."""

import numpy

from ..plot_curve import PlotCurve
from .... import shared
from ....compute_tools import plant
from ....gui_styles.theme_handler import curve_color


class ResultsPlant(PlotCurve):
    """Apparent stiffness of a curve (see compute_tools.plant), approach and
    retraction from the point of contact of the elasticity: derivative of
    the force vs the indentation with the fits (Yamasaki), or force vs the
    indentation with the linear fits (Beauzamy), and the values of the
    cell.

    The curve is analysed again with the parameters used for the maps: the
    values are the ones of the maps. In the compute tab (fit preview), with
    the current parameters.
    """

    def __init__(self, parent):
        super().__init__(parent)

        if self.data.curve_distance_units == "um":
            self.xlabel = "Indentation [μm]"
            self.x_factor = 1e6
        else:
            self.xlabel = "Indentation [nm]"
            self.x_factor = 1e9

        self.used = self.plot_type != "fit_preview"
        self.result = None
        self.message = "No apparent stiffness for this curve"
        if not self.used:
            reason = plant.check_file(self.data)
            if reason is None:
                self.result = plant.curve(self.data, self.xpos, self.ypos,
                                          used=False)
            else:
                self.message = "Preview not available:\n" + reason
        elif self.data.plant_calculated:
            self.result = plant.curve(self.data, self.xpos, self.ypos)
        method = plant.settings(self.data, self.used)["method"]
        self.derivative = method == 0
        if self.derivative:
            self.ylabel = "dF/dδ [N/m]"
            self.force_factor = 1.0
        else:
            self.ylabel = "Force [" + self.data.curve_force_units + "]"
            self.force_factor = 1e9 if self.data.curve_force_units == "nN" \
                else 1e12

        self.plot_fig()

    def plot_fig(self):
        """Plots the curves, the fits and the values."""
        result = self.result
        if result is None:
            self.axes.text(0.5, 0.5, self.message,
                           transform=self.axes.transAxes, ha="center",
                           va="center")
        else:
            xf, yf = self.x_factor, self.force_factor
            curves = result["curves"]
            for index, (name, label) in enumerate(
                    (("approach", "Approach"), ("retraction", "Retract"))):
                color = curve_color(name)
                part = result[name]
                if self.derivative:
                    if part is None:
                        continue
                    self.axes.plot(part["delta"] * xf, part["k"] * yf,
                                   color=color, label=label,
                                   **self.curve_style())
                    self.axes.plot(part["delta"] * xf, part["fit"] * yf,
                                   "-", color="red", lw=1.2)
                    self.axes.axhline(part["k_as"] * yf, ls=":", lw=1,
                                      color=color)
                else:
                    delta, force = curves[index]
                    self.axes.plot(delta * xf, force * yf, color=color,
                                   label=label, **self.curve_style())
                    if part is not None:
                        self.axes.plot(part["delta"] * xf, part["fit"] * yf,
                                       "-", color="red", lw=1.2)
            self.show_values()

        self.plot_labels()
        self.red_line()
        self.start_plot()

    def show_values(self):
        """Values of the curve and of its cell, in the legend or in a box."""
        values = self.result["values"]
        cell = self.result["cell"]

        def number(value, factor=1.0, digits=3):
            if not numpy.isfinite(value):
                return "-"
            return format(value * factor, "." + str(digits) + "g")

        text = ["$k_\\mathrm{as}$ = " + number(values["k_as"]) +
                " N/m (approach " +
                number(values["k_approach"]) + ", retract " +
                number(values["k_retraction"]) + ")",
                "R² = " + number(values["r2_approach"], digits=3) +
                ", " + number(values["r2_retraction"], digits=3),
                "E Hertz = " + number(values["modulus_hertz"], 1e-6) +
                " MPa"]
        if cell is None:
            text.append("Not in a cell (ROI of a map): $k_\\mathrm{as}$ "
                        "only")
        else:
            text += ["Cell: ROI " + str(int(values["roi"])) + " of " +
                     str(cell[1].filename) + ", f = " +
                     number(values["f"]),
                     "Yamasaki: E = " + number(values["modulus"], 1e-6) +
                     " MPa, P = " + number(values["pressure"], 1e-6) +
                     " MPa, τ = " + number(values["tau"], digits=2),
                     "E Hertz: P = " +
                     number(values["pressure_hertz"], 1e-6) +
                     " MPa, τ = " +
                     number(values["tau_hertz"], digits=2)]
            if not cell[0]["dome"]:
                text.append("Saddle: the surface is not a dome here")
        if shared.exp.display_curve_legend:
            for line in text:
                self.axes.plot([], [], " ", label=line)
            self.axes.legend(loc=4 if self.derivative else 2,
                             prop=self.font)
        else:
            self.axes.text(0.98 if self.derivative else 0.02, 0.03 if
                           self.derivative else 0.97, "\n".join(text),
                           transform=self.axes.transAxes,
                           ha="right" if self.derivative else "left",
                           va="bottom" if self.derivative else "top",
                           fontsize=9)
