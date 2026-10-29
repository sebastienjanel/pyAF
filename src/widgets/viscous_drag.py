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

"""Dialog with the viscous drag of the cantilever measured at several
heights above the surface, and its extrapolation to the surface b(0)."""

import os

import numpy
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from qtpy import QtWidgets

from ..compute_tools import viscous_drag
from ..gui_styles.theme_handler import curve_color


class ViscousDragDialog(QtWidgets.QDialog):
    """b(h) of a file (median of the curves, standard deviation as error
    bars) with the fit b(h) = A / (h + h_eff); "Use b(0)" accepts the
    dialog (the value is in self.result["b0"])."""

    def __init__(self, parent, data, result, sample=None):
        super().__init__(parent)
        self.result = result
        self.setWindowTitle("Viscous drag")

        figure = Figure(figsize=(5.5, 3.6), tight_layout=True)
        canvas = FigureCanvasQTAgg(figure)
        axes = figure.add_subplot(111)
        heights = result["heights"] * 1e9
        axes.errorbar(heights, result["b"] * 1e6, yerr=result["b_std"] * 1e6,
                      fmt="o", color=curve_color("approach"),
                      label="Measured (median of the curves)")
        b0, h_eff = result["b0"], result["h_eff"]
        if numpy.isfinite(b0):
            x = numpy.linspace(0, heights.max(), 200)
            axes.plot(x, viscous_drag.model(x * 1e-9, b0, h_eff) * 1e6, "r-",
                      label="b(h) = A / (h + h_eff)")
            axes.plot([0], [b0 * 1e6], "rs", label="b(0)")
        axes.set_xlabel("Height above the surface [nm]")
        axes.set_ylabel("b [µN·s/m]")
        axes.set_xlim(left=0)
        axes.set_ylim(bottom=0)
        axes.legend(fontsize=8)

        name = os.path.basename(str(data.filename))
        if numpy.isfinite(b0):
            text = ("b(0) = <b>" + format(b0 * 1e6, ".3g") +
                    " µN·s/m</b>, h_eff = " +
                    format(h_eff * 1e9, ".0f") + " nm, R² = " +
                    format(result["r2"], ".3f") + "<br>" + name + ": " +
                    str(result["nbr"]) + " curves, " + str(len(heights)) +
                    " heights from " + format(heights.min(), ".0f") +
                    " nm<br><small>Extrapolated below the lowest height: "
                    "check that the fit follows the points close to the "
                    "surface.</small>")
        else:
            text = name + ": the fit failed (b(h) does not decrease with " \
                "the height)."
        # Dates of the acquisitions (same cantilever, same day expected)
        if sample is not None and sample is not data and data.date and \
                sample.date:
            text += "<br><small>Acquired " + str(data.date) + \
                ", sample " + os.path.basename(str(sample.filename)) + \
                " acquired " + str(sample.date) + "</small>"
        # Piezo correction applied to F / delta (set before the drag)
        if sample is not None and sample.microrheo_piezo:
            text += "<br><small>With the piezo correction " + \
                sample.microrheo_piezo["file"] + " [" + ", ".join(
                    str(k) for k in sample.microrheo_piezo["curve"]) + \
                "] (" + ("phase and amplitude" if
                         sample.microrheo_piezo_amplitude else "phase") + \
                ")</small>"
        else:
            text += "<br><small>Without piezo correction (set it before "\
                "the drag to apply it to the drag measurement)</small>"
        warning = cantilever_warning(data, sample)
        if warning:
            text += "<br><b>" + warning + "</b>"
        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)

        buttons = QtWidgets.QDialogButtonBox()
        self.BT_use = buttons.addButton(
            "Use b(0)", QtWidgets.QDialogButtonBox.AcceptRole)
        self.BT_use.setEnabled(bool(numpy.isfinite(b0)))
        buttons.addButton(QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(canvas)
        layout.addWidget(label)
        layout.addWidget(buttons)
        self.setLayout(layout)


def cantilever_warning(data, sample, tolerance=0.01):
    """Warning if the drag measurement and the sample file were not
    acquired with the same cantilever (spring constant or deflection
    sensitivity different by more than the tolerance), None otherwise."""
    if sample is None or sample is data:
        return None
    differences = []
    for name, attribute, unit, factor in (
            ("spring constant", "spring_constant", "N/m", 1e9),
            ("deflection sensitivity", "deflection_sensitivity", "nm/V",
             1)):
        a, b = getattr(data, attribute), getattr(sample, attribute)
        if a and b and abs(a - b) > tolerance * abs(b):
            differences.append(name + " " + format(a * factor, ".4g") +
                               " vs " + format(b * factor, ".4g") + " " +
                               unit)
    if not differences:
        return None
    return "Not the cantilever of " + os.path.basename(
        str(sample.filename)) + ": " + ", ".join(differences) + "."
