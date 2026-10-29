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


"""SuperPlot of the conditions (Lord et al., J Cell Biol 2020)."""

import numpy

from .plot_results import PlotResults
from .single import set_white_background, quasirandom_offsets
from ... import shared
from ...gui_styles import theme_handler
from ...tools import palettes
from ...tools import superplot_tools

# Width of the beeswarm of the maps (the space of a condition is 1)
SWARM_WIDTH = 0.3


def experiment_colors(conditions):
    """Color of each experiment number (Okabe-Ito, colorblind safe)."""
    numbers = sorted({r["experiment"] for c in conditions
                      for r in c["replicates"]})
    colors = palettes.get_colors("Okabe-Ito", len(numbers))
    return dict(zip(numbers, colors))


class PlotSuperPlot(PlotResults):
    """SuperPlot: the files (small dots, one value per file) colored by
    experiment, the replicate means (large dots), the mean of the replicates
    of each condition with its SD or 95% CI, and the tests against the
    control (on the replicates)."""

    def __init__(self, parent):
        super().__init__(parent)

        # The label of the values goes on the y axis
        self.ylabel = self.xlabel
        self.xlabel = None

        self.plot_fig()

    def plot_fig(self):
        """Draw the SuperPlot."""
        ax = self.axes
        opts = superplot_tools.get_options()
        white = opts["background"] == "White"
        self.fig.set_frameon(self.parent.single or white)
        conditions = superplot_tools.get_data()

        if not conditions:
            if shared.exp.results_type == "loading_rates":
                text = ("No SuperPlot for the loading rates "
                        "(see the DFS of the Replicates tab)")
            else:
                text = ("Put the files in replicates (Plots tab), and the "
                        "replicates in conditions (Replicates tab)")
            ax.text(0.5, 0.5, text, ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            if white:
                set_white_background(self.fig, ax)
            self.start_plot()
            return

        log = opts["log"]
        colors = experiment_colors(conditions)
        dark = theme_handler.current_theme == "Dark" and not white
        ink = "0.9" if dark else "black"

        # All the curves in the background
        if opts["curves"]:
            curves = [c["curves"][c["curves"] > 0] if log else c["curves"]
                      for c in conditions]
            positions = [i for i, v in enumerate(curves) if len(v)]
            parts = ax.violinplot([curves[i] for i in positions],
                                  positions=positions, widths=0.8,
                                  showextrema=False)
            for body in parts["bodies"]:
                body.set_facecolor("0.6")
                body.set_edgecolor("none")
                body.set_alpha(0.25)

        # Maps (small dots), beeswarm in each condition
        map_values = [numpy.array([m["value"] for r in c["replicates"]
                                   for m in r["maps"]]) for c in conditions]
        map_colors = [[colors[r["experiment"]] for r in c["replicates"]
                       for _ in r["maps"]] for c in conditions]
        if log:
            # No log of the values <= 0
            keep = [v > 0 for v in map_values]
            map_values = [v[k] for v, k in zip(map_values, keep)]
            map_colors = [[col for col, ok in zip(cols, k) if ok]
                          for cols, k in zip(map_colors, keep)]
        offsets = quasirandom_offsets(map_values, 1.0, "width", log)
        for i, (vals, offs, cols) in enumerate(
                zip(map_values, offsets, map_colors)):
            if len(vals):
                ax.scatter(i + offs * SWARM_WIDTH / 0.4, vals, s=22,
                           color=cols, alpha=0.6, linewidths=0, zorder=3)

        # Replicate means (large dots), spread in the condition
        positions = {}
        for i, c in enumerate(conditions):
            reps = c["replicates"]
            spread = numpy.linspace(-0.12, 0.12, len(reps)) \
                if len(reps) > 1 else [0.0]
            for x, rep in zip(spread, reps):
                if log and rep["value"] <= 0:
                    continue
                ax.scatter(i + x, rep["value"], s=110,
                           color=colors[rep["experiment"]],
                           edgecolors=ink, linewidths=1, zorder=5)
                positions[(c["id"], rep["experiment"])] = (i + x,
                                                           rep["value"])

        # Lines between the replicates of the same experiment
        if opts["paired"]:
            for number in colors:
                points = [positions[(c["id"], number)] for c in conditions
                          if (c["id"], number) in positions]
                if len(points) > 1:
                    xs, ys = zip(*points)
                    ax.plot(xs, ys, "-", color=colors[number], linewidth=1,
                            alpha=0.7, zorder=4)

        # Mean of the replicates with its SD or 95% CI
        for i, c in enumerate(conditions):
            s = superplot_tools.summary(c)
            error = s["ci"] if opts["error"] == "95% CI" else s["sd"]
            ax.plot([i - 0.25, i + 0.25], [s["mean"]] * 2, "-",
                    color=ink, linewidth=2, zorder=6)
            if numpy.isfinite(error):
                ax.errorbar(i, s["mean"], yerr=error, fmt="none",
                            ecolor=ink, elinewidth=1.5, capsize=8, zorder=6)

        # Tests against the control
        results = superplot_tools.compare(conditions, opts)
        tests = {t["name"]: t for t in results}
        control = next((c["name"] for c in conditions
                        if c["id"] == opts["control"]), None)
        if log:
            ax.set_yscale("log")
        bottom, top = ax.get_ylim()
        if opts["control"] == superplot_tools.ALL_PAIRS:
            # All pairs: brackets above the significant pairs, the shortest
            # ones first
            significant = sorted(
                (t for t in results if numpy.isfinite(t["p"]) and
                 t["p"] < 0.05),
                key=lambda t: (t["pair"][1] - t["pair"][0], t["pair"][0]))
            level = top
            for test in significant:
                i, j = test["pair"]
                if log:
                    y, tick = level * 1.12, level * 1.04
                else:
                    step = 0.07 * (top - bottom)
                    y, tick = level + step, level + 0.35 * step
                ax.plot([i, i, j, j], [tick, y, y, tick], "-", color=ink,
                        linewidth=1, zorder=6)
                ax.text((i + j) / 2, y, superplot_tools.stars(test["p"]),
                        ha="center", va="bottom", color=ink,
                        fontproperties=self.font)
                level = y * 1.06 if log else y + 0.04 * (top - bottom)
            if significant:
                ax.set_ylim(bottom, level * 1.1 if log else
                            level + 0.06 * (top - bottom))
        for i, c in enumerate(conditions):
            test = tests.get(c["name"] + " vs " + str(control))
            if test is not None:
                ax.text(i, top, superplot_tools.stars(test["p"]),
                        ha="center", va="bottom", color=ink,
                        fontproperties=self.font)

        # Names of the conditions, with the numbers of replicates and maps
        ax.set_xticks(range(len(conditions)))
        labels = []
        for c in conditions:
            s = superplot_tools.summary(c)
            label = c["name"]
            if c["name"] == control:
                label += " (control)"
            labels.append(label + "\nN=%d, n=%d" % (s["n"], s["maps"]))
        ax.set_xticklabels(labels)
        ax.set_xlim(-0.6, len(conditions) - 0.4)

        # Legend of the experiments
        for number, color in colors.items():
            ax.scatter([], [], s=60, color=color, edgecolors=ink,
                       label="Experiment " + str(number))
        ax.legend(fontsize="small", loc="upper left",
                  bbox_to_anchor=(1.0, 1.0))

        ax.set_ylabel(self.ylabel or "", fontproperties=self.font)
        self.ylabel = None
        if white:
            set_white_background(self.fig, ax)
        self.start_plot()
