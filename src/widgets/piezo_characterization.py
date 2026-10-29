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


"""Dialog with the piezo characterization (lag and gain of the piezo vs the
frequency) measured on a curve of a rigid sample, for the microrheology."""

import numpy
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from qtpy import QtWidgets

from ..gui_styles.theme_handler import curve_color

# Below this amplitude quotient at the lowest frequency, the deflection
# does not follow the piezo: not a rigid sample
MIN_QUOTIENT = 0.8


class PiezoCharacterizationDialog(QtWidgets.QDialog):
    """Phase and amplitude quotient of the deflection over the height vs
    the frequency; "Use" accepts the dialog."""

    def __init__(self, parent, result):
        super().__init__(parent)
        self.setWindowTitle("Piezo characterization")

        freqs = numpy.array(result["frequencies"])
        figure = Figure(figsize=(5.5, 3.6), tight_layout=True)
        canvas = FigureCanvasQTAgg(figure)
        axes = figure.add_subplot(111)
        axes.plot(freqs, result["phase"], "o-", color=curve_color("approach"),
                  label="Phase (piezo lag)")
        axes.axhline(0, color="grey", lw=0.8, ls="--")
        axes.set_xscale("log")
        axes.set_xlabel("Frequency [Hz]")
        axes.set_ylabel("Phase [°]", color=curve_color("approach"))
        axes2 = axes.twinx()
        axes2.plot(freqs, result["quotient"], "s-",
                   color=curve_color("retraction"), label="Amplitude quotient")
        axes2.axhline(1, color=curve_color("retraction"), lw=0.8, ls=":")
        axes2.set_ylabel("Amplitude quotient", color=curve_color("retraction"))
        lines = axes.get_lines()[:1] + axes2.get_lines()[:1]
        axes.legend(lines, [line.get_label() for line in lines], fontsize=8)

        rigid = result["quotient"][0] >= MIN_QUOTIENT
        text = (result["file"] + ", curve [" +
                ", ".join(str(k) for k in result["curve"]) + "]: " +
                str(len(freqs)) + " frequencies.<br><small>The height is "
                "corrected by the lag and the gain of the piezo measured on "
                "a rigid sample (deflection following the piezo).</small>")
        if not rigid:
            text += ("<br><b>The deflection does not fully follow the "
                     "piezo</b> (amplitude quotient " +
                     format(result["quotient"][0], ".2f") + " at " +
                     format(freqs[0], "g") + " Hz, 1 on a rigid sample): "
                     "check that the curve is on the substrate.")
        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)

        buttons = QtWidgets.QDialogButtonBox()
        buttons.addButton("Use", QtWidgets.QDialogButtonBox.AcceptRole)
        buttons.addButton(QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(canvas)
        layout.addWidget(label)
        layout.addWidget(buttons)
        self.setLayout(layout)
