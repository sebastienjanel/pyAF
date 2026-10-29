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


"""Plots the microrheology of a curve."""

import numpy

from ..plot_curve import PlotCurve
from .... import shared
from ....compute_tools import microrheology
from ....gui_styles.theme_handler import curve_color
from ....tools import curve_tools


class ResultsMicrorheology(PlotCurve):
    """Storage and loss moduli G'(f) and G''(f) of the curve (log - log),
    one point per modulation segment.

    The curve is analysed again with the parameters used for the maps: the
    values are the ones of the maps. In the compute tab (preview), with the
    current parameters, or the modulation segments vs time if the
    microrheology can not be computed (no elasticity fit).
    """

    def __init__(self, parent):
        super().__init__(parent)

        self.used = self.plot_type != "fit_preview"
        self.results = None
        self.segments = None
        # Results tab, map average: median and quartiles of the curves
        self.statistics = None
        self.message = "No microrheology for this curve"
        if self.used and \
                self.data.display_curve_type == "microrheology_map":
            self.statistics = microrheology.map_statistics(self.data)
            if self.statistics is None:
                self.message = "No microrheology for this file"
        elif not self.used and \
                getattr(self.data, "microrheo_preview_raw", False) and \
                microrheology.available(self.data):
            # Compute tab: the oscillations (raw data)
            self.segments = curve_tools.get_modulation_segments(
                self.data, self.xpos, self.ypos)
            self.message = ""
        elif not self.used:
            reason = microrheology.check_file(self.data)
            if reason is None:
                self.results = microrheology.curve(
                    self.data, self.xpos, self.ypos, used=False)
                skip = microrheology.skip_reason(self.data, self.xpos,
                                                 self.ypos)
                if skip is not None:
                    self.message = "Curve not analysed:\n" + skip
            elif microrheology.available(self.data):
                self.segments = curve_tools.get_modulation_segments(
                    self.data, self.xpos, self.ypos)
                self.message = "G′, G″ not available: " + reason
            else:
                self.message = "Not available:\n" + reason
        elif self.data.microrheo_calculated:
            self.results = microrheology.curve(self.data, self.xpos,
                                               self.ypos)
            skip = microrheology.skip_reason(self.data, self.xpos, self.ypos)
            if skip is not None:
                self.message = "Curve not analysed:\n" + skip

        if self.results is not None or self.statistics is not None:
            self.xlabel = "Frequency [Hz]"
            self.ylabel = "G′, G″ [Pa]"
        elif self.segments is not None:
            self.xlabel = "Time [s]"
            self.ylabel = "Deflection [nm]"
        self.plot_fig()

    def plot_fig(self):
        """Plots G' and G'' vs the frequency (or the modulation segments)."""
        if self.statistics is not None:
            self.plot_statistics()
        elif self.results is not None:
            freq = numpy.array([r["frequency"] for r in self.results])
            storage = numpy.array([r["storage"] for r in self.results])
            loss = numpy.array([r["loss"] for r in self.results])
            # Log scales: the negative values (not physical: piezo lag or
            # viscous drag not corrected) are drawn as open symbols at their
            # absolute value
            for values, name, color, marker in (
                    (storage, "G′ (storage)", curve_color("approach"),
                     "o"),
                    (loss, "G″ (loss)", curve_color("retraction"), "s")):
                pos = values > 0
                self.axes.plot(freq[pos], values[pos], marker=marker,
                               color=color, ls="-", ms=5, label=name)
                neg = values < 0
                if numpy.any(neg):
                    self.axes.plot(freq[neg], -values[neg], marker=marker,
                                   color=color, ls="none", ms=6,
                                   mfc="none", label=name + " < 0")
            self.axes.set_xscale("log")
            if numpy.any(numpy.abs(numpy.r_[storage, loss]) > 0):
                self.axes.set_yscale("log")
            # Crossover frequency G' = G''
            cross = microrheology.crossover(freq, storage, loss)
            if numpy.isfinite(cross):
                self.axes.axvline(cross, color="grey", lw=0.8, ls="--",
                                  label="Crossover " + format(cross, ".3g") +
                                  " Hz")
            delta0 = self.results[0]["indentation"]
            text = "Working indentation " + format(delta0 * 1e9, ".0f") + \
                " nm"
            if shared.exp.display_curve_legend:
                self.axes.plot([], [], " ", label=text)
                self.axes.legend(loc=2, prop=self.font)
            else:
                self.axes.text(0.02, 0.97, text,
                               transform=self.axes.transAxes, ha="left",
                               va="top", fontsize=9)
        elif self.segments is not None:
            for frequency, _, defl, time in self.segments:
                label = format(frequency, "g") + " Hz"
                self.axes.plot(time, defl, lw=0.6, label=label)
            if self.data.auto_zoom_in:
                self.zoom_oscillations()
            if shared.exp.display_curve_legend:
                self.axes.legend(loc=2, prop=self.font, ncol=3)
            self.axes.text(0.98, 0.03, self.message,
                           transform=self.axes.transAxes, ha="right",
                           va="bottom", fontsize=9)
        else:
            self.axes.text(0.5, 0.5, self.message,
                           transform=self.axes.transAxes, ha="center",
                           va="center")

        self.plot_labels()
        if self.statistics is None:
            self.red_line()
        self.start_plot()

    def zoom_oscillations(self):
        """Zoom in (compute tab toggle): y axis on the oscillations of a
        typical segment (median over the segments of their 5th and 95th
        percentiles, without the transients at their start and end, 50 %
        margin); the large oscillations at high frequency go out of the
        plot."""
        lows, highs = [], []
        for _, _, defl, _ in self.segments:
            defl = numpy.asarray(defl, float)
            cut = max(1, len(defl) // 50)
            defl = defl[cut:-cut] if len(defl) > 2 * cut else defl
            defl = defl[numpy.isfinite(defl)]
            if defl.size >= 2:
                low, high = numpy.percentile(defl, [5, 95])
                lows.append(low)
                highs.append(high)
        if not lows:
            return
        low, high = float(numpy.median(lows)), float(numpy.median(highs))
        margin = 0.5 * (high - low) if high > low else 1.0
        self.axes.set_ylim(low - margin, high + margin)

    def plot_statistics(self):
        """Median of the curves of the map with the interquartile range
        (band) for G' and G'', and the crossover of the medians."""
        stats = self.statistics
        freq = stats["frequencies"]
        for name, label, color, marker in (
                ("storage", "G\u2032 (storage)", curve_color("approach"),
                 "o"),
                ("loss", "G\u2033 (loss)", curve_color("retraction"), "s")):
            low, median, high = stats[name].T
            ok = median > 0
            self.axes.plot(freq[ok], median[ok], marker=marker, color=color,
                           ls="-", ms=5, label=label + ", median")
            band = ok & (low > 0) & (high > 0)
            self.axes.fill_between(freq[band], low[band], high[band],
                                   color=color, alpha=0.2, lw=0)
        self.axes.set_xscale("log")
        if numpy.any(numpy.r_[stats["storage"][:, 1],
                              stats["loss"][:, 1]] > 0):
            self.axes.set_yscale("log")
        cross = microrheology.crossover(freq, stats["storage"][:, 1],
                                        stats["loss"][:, 1])
        if numpy.isfinite(cross):
            self.axes.axvline(cross, color="grey", lw=0.8, ls="--",
                              label="Crossover of the medians " +
                              format(cross, ".3g") + " Hz")
        text = "Average: " + str(max(stats["nbr"]) if stats["nbr"] else 0) + \
            " curves (band: 25\u201375 %)"
        if shared.exp.display_curve_legend:
            self.axes.plot([], [], " ", label=text)
            self.axes.legend(loc=2, prop=self.font)
        else:
            self.axes.text(0.02, 0.97, text, transform=self.axes.transAxes,
                           ha="left", va="top", fontsize=9)
