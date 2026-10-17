# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
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

"""Plots the results of the single files (one violin per file) or of the
groups (one violin per group)."""

import matplotlib
import numpy
import pandas
import seaborn
from scipy import stats

from .plot_results import PlotResults
from ... import shared
from ...gui_styles import theme_handler
from ...experiment import get_short_name, get_violin_options
from ...tools import dfs_tools

# Maximum horizontal shift of the values with the jitter (the space of a
# file is 1)
JITTER = 0.15


# Maximum space given to each file on the x axis, in pixels: with few files
# the violins would otherwise fill the whole width of the plot
MAX_SLOT_WIDTH = 200

# Width of the box inside the violins (the space of a file is 1)
BOX_WIDTH = 0.12


class PlotSingleResults(PlotResults):
    """Plots the results of the displayed files, one violin per file."""

    # "single" (plots tab) or "groups" (grouped plots tab)
    MODE = "single"

    def __init__(self, parent):
        super().__init__(parent)

        # The label of the values (defined for the histograms in
        # PlotResults) goes on the y axis
        self.ylabel = self.xlabel
        self.xlabel = None

        self.names, self.values, self.colors = \
            get_displayed_values(self.MODE)

        self.plot_fig()

    def plot_fig(self):
        """Draw the violins."""
        opts = get_violin_options(shared.exp, self.MODE)
        ax = self.axes

        if shared.exp.results_type == "loading_rates":
            self.plot_loading_rates(opts)
            return

        if self.parent.single and self.names:
            # Separate window: wide enough for the files (1.3 inch each)
            self.fig.set_size_inches(
                max(6.4, 1.6 + 1.3 * len(self.names)), 4.8)

        # The figures in the tabs have no background (the one of the
        # window is visible), except for the white background
        self.fig.set_frameon(self.parent.single or
                             opts["background"] == "White")

        if not self.names:
            if self.MODE == "groups":
                text = "Put files in replicates (Plots tab) and tick the replicates"
            else:
                text = "Tick some files in the table to display them"
            ax.text(0.5, 0.5, text,
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            if opts["background"] == "White":
                set_white_background(self.fig, ax)
            self.start_plot()
            return

        names = self.names

        # Unique labels (two files can have the same name)
        labels = []
        for i, name in enumerate(names):
            label = name
            if opts["show_counts"]:
                label += "\n(n=" + str(len(self.values[i])) + ")"
            while label in labels:
                label += " "
            labels.append(label)

        log_scale = opts["log_scale"]
        values = self.values
        if log_scale:
            # No log of the values <= 0
            values = [v[v > 0] for v in values]

        frame = pandas.DataFrame({
            "file": numpy.repeat(labels, [len(v) for v in values]),
            "value": numpy.concatenate(values)})

        # The boxes are drawn by matplotlib (larger than the ones of
        # seaborn, and with notches)
        inner = {"Box": None, "Notched box": None, "Quartiles": "quart",
                 "None": None}[opts["inner"]]
        density_norm = {"Same area": "area", "Same width": "width",
                        "Number of values": "count"}[opts["density_norm"]]

        seaborn.violinplot(
            data=frame, x="file", y="value", hue="file", order=labels,
            hue_order=labels, palette=self.colors, legend=False,
            inner=inner, bw_adjust=opts["bw_adjust"],
            cut=0 if opts["cut"] else 2, density_norm=density_norm,
            # Each file is a hue: normalize all the violins together, else
            # "area" and "count" would be the same as "width"
            common_norm=True,
            log_scale=log_scale, linewidth=1, ax=ax)

        if opts["inner"] in ("Box", "Notched box"):
            self.plot_boxes(ax, values, opts["inner"] == "Notched box")

        # Individual values on top of the violins
        points = opts["points"]
        if points != "None":
            dark = theme_handler.current_theme == "Dark" and \
                opts["background"] != "White"
            point_color = "0.85" if dark else "0.15"
            if points == "Beeswarm":
                offsets = quasirandom_offsets(
                    values, opts["bw_adjust"], density_norm, log_scale)
            else:
                rng = numpy.random.default_rng(0)
                offsets = [rng.uniform(-JITTER, JITTER, len(v))
                           for v in values]
            for i, (v, offset) in enumerate(zip(values, offsets)):
                ax.scatter(i + offset, v, s=4, color=point_color,
                           alpha=0.5, linewidths=0, zorder=4)

        # Limit the width of the violins when there are few files (centered)
        # (in screen pixels: the dpi of the figure is doubled on a retina
        # screen)
        ratio = getattr(self.fig.canvas, "device_pixel_ratio", 1) or 1
        width_px = self.fig.get_figwidth() * self.fig.dpi / ratio
        slots = max(len(labels), width_px / MAX_SLOT_WIDTH)
        center = (len(labels) - 1) / 2
        ax.set_xlim(center - slots / 2, center + slots / 2)

        if opts["y_mode"] == "manual":
            ax.set_ylim(opts["y_min"], opts["y_max"])

        # Labels (seaborn sets the names of the columns)
        ax.set_ylabel(self.ylabel or "", fontproperties=self.font)
        self.ylabel = None
        ax.set_xlabel("")
        if opts["background"] == "White":
            set_white_background(self.fig, ax)

        # Rotate the names when they are wider than the space of a file
        # (approximate width of the characters: 0.6 em)
        size = matplotlib.font_manager.FontProperties(
            size=matplotlib.rcParams["xtick.labelsize"]).get_size_in_points()
        longest = max(len(line) for name in names for line in name.split())
        name_width = longest * 0.6 * size / 72
        slot_width = self.fig.get_figwidth() * ax.get_position().width / slots
        if name_width > slot_width:
            for tick in ax.get_xticklabels():
                tick.set_rotation(30)
                tick.set_horizontalalignment("right")
                tick.set_rotation_mode("anchor")

        self.start_plot()

    def plot_loading_rates(self, opts):
        """Force vs loading rate (log scale), for the dynamic force
        spectroscopy: the events of each file (or group), their median point
        with the quartiles, and the Bell-Evans fits of the DFS window."""
        ax = self.axes
        self.fig.set_frameon(self.parent.single or
                             opts["background"] == "White")
        items = dfs_tools.get_items(self.MODE)

        if not items:
            if self.MODE == "groups":
                text = "Put files in replicates (Plots tab) and tick the replicates"
            else:
                text = "Tick some files in the table to display them"
            ax.text(0.5, 0.5, text,
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            for name, forces, lrs, color, _ in items:
                positive = lrs > 0
                label = name + " (n=" + str(forces.size)
                if not positive.all():
                    # No place for them on the log scale
                    label += ", " + str(int((~positive).sum())) + \
                        " with r ≤ 0 not shown"
                label += ")"
                ax.scatter(lrs[positive], forces[positive], s=6,
                           color=color, alpha=0.35, linewidths=0, zorder=2)
                r, f, r_q, f_q = dfs_tools.median_point(forces, lrs)
                if r > 0:
                    ax.errorbar(
                        r, f, xerr=[[r - r_q[0]], [r_q[1] - r]],
                        yerr=[[f - f_q[0]], [f_q[1] - f]], fmt="o",
                        markersize=8, color=color, markeredgecolor="black",
                        ecolor="black", elinewidth=1, capsize=3, zorder=4,
                        label=label)
                else:
                    ax.scatter([], [], color=color, label=label)

            self.plot_dfs_fits(ax)
            ax.set_xscale("log")
            ax.legend(fontsize="small", loc="upper left")

        if opts["y_mode"] == "manual":
            ax.set_ylim(opts["y_min"], opts["y_max"])
        ax.set_xlabel("Loading rate [pN/s]", fontproperties=self.font)
        ax.set_ylabel("Force [pN]", fontproperties=self.font)
        self.xlabel = None
        self.ylabel = None
        if opts["background"] == "White":
            set_white_background(self.fig, ax)
        self.start_plot()

    def plot_dfs_fits(self, ax):
        """Lines of the displayed Bell-Evans fits, with the points of the
        binned fits."""
        number = 0
        for fit in dfs_tools.get_fits(self.MODE):
            number += 1
            if not fit.get("display") or fit["method"] == "old":
                continue
            r = numpy.array(fit["points_r"])
            line_r = numpy.logspace(numpy.log10(r.min()),
                                    numpy.log10(r.max()), 50)
            # Fit in SI units, plot in pN and pN/s
            line_f = (fit["a"] * numpy.log(line_r * 1e-12) + fit["b"]) * 1e12
            if numpy.isfinite(fit["x"]):
                label = ("Fit " + str(number) + ": x = " +
                         "%.3g nm, koff = %.3g s$^{-1}$" % (fit["x"] * 1e9,
                                                           fit["koff"]))
            else:
                label = "Fit " + str(number) + ": slope ≤ 0"
            ax.plot(line_r, line_f, "-", color="red", linewidth=2, zorder=5,
                    label=label)
            if fit["method"] == "binned":
                ax.plot(r, fit["points_f"], "s", markersize=6,
                        markerfacecolor="none", markeredgecolor="red",
                        zorder=5)

    @staticmethod
    def plot_boxes(ax, values, notch):
        """Box plots inside the violins: median, quartiles and whiskers at
        the last values within 1.5 IQR of the quartiles. The notches show
        the 95% confidence interval of the median (1.57 IQR / sqrt(n))."""
        positions = [i for i, v in enumerate(values) if len(v) > 0]
        values = [numpy.asarray(v) for v in values if len(v) > 0]
        if not values:
            return
        line = {"color": "black", "linewidth": 1}
        ax.boxplot(values, positions=positions, widths=BOX_WIDTH,
                   notch=notch, patch_artist=True, showfliers=False,
                   manage_ticks=False, showcaps=False, zorder=3,
                   boxprops={"facecolor": "white", "edgecolor": "black",
                             "linewidth": 1},
                   medianprops={"color": "black", "linewidth": 2},
                   whiskerprops=line)


def set_white_background(fig, ax):
    """White background and black axes and texts (for the publications),
    whatever the theme."""
    fig.patch.set_facecolor("white")
    # Saved images (Figures > Open results) with the color of the figure,
    # not the one of the theme
    matplotlib.rcParams["savefig.facecolor"] = "auto"
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_edgecolor("black")
    ax.tick_params(which="both", colors="black", labelcolor="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    for text in ax.texts:
        text.set_color("black")


def quasirandom_offsets(values, bw_adjust, density_norm, log_scale):
    """Horizontal shifts of the values inside the violins (quasirandom
    beeswarm of ggbeeswarm): proportional to the density at each value,
    distributed with a van der Corput sequence along the values. Fast, and
    all the values are drawn (unlike a swarm)."""
    densities = []
    for v in values:
        x = numpy.log10(v) if log_scale else numpy.asarray(v, dtype=float)
        if len(x) < 2 or numpy.ptp(x) == 0:
            densities.append(numpy.ones(len(x)))
            continue
        kde = stats.gaussian_kde(x)
        kde.set_bandwidth(kde.factor * bw_adjust)
        # Density on a grid, interpolated at the values (fast)
        grid = numpy.linspace(x.min(), x.max(), 256)
        densities.append(numpy.interp(x, grid, kde(grid)))

    # Same scaling as the violins (see seaborn.violinplot)
    peaks = [d.max() if len(d) else 1.0 for d in densities]
    max_peak = max(peaks)
    max_count = max(len(v) for v in values)
    half_width = 0.4  # Violins of width 0.8

    offsets = []
    for v, d, peak in zip(values, densities, peaks):
        if density_norm == "area":
            scale = d / max_peak
        elif density_norm == "count":
            scale = d / peak * len(v) / max_count
        else:
            scale = d / peak
        # van der Corput sequence (base 2), following the order of the
        # values
        order = numpy.argsort(numpy.argsort(v, kind="stable"))
        sequence = numpy.array([van_der_corput(k + 1) for k in order])
        # Slightly inside the outline of the violin
        offsets.append((sequence - 0.5) * 2 * half_width * 0.9 * scale)
    return offsets


def van_der_corput(n, base=2):
    """n-th number of the van der Corput sequence (between 0 and 1)."""
    value, denominator = 0.0, 1.0
    while n:
        n, remainder = divmod(n, base)
        denominator *= base
        value += remainder / denominator
    return value


def get_displayed_values(mode="single"):
    """Values of the files (mode "single") or of the groups (mode "groups")
    ticked in their table.

    Returns the names, the values and the colors. For the loading rates,
    the values are the loading rates (the forces are not used).
    """
    names = []
    values = []
    colors = []

    if mode == "groups":
        all_data = shared.groups_data
        # The first group is the None group (the files without group)
        items = shared.exp.groups_list
        first = 1
    else:
        all_data = shared.single_data
        items = shared.exp.results_list
        first = 0

    if all_data is None:
        return names, values, colors

    for i in range(first, min(len(all_data), len(items))):
        data = all_data[i]
        item = items[i]
        if not item.display:
            continue

        if shared.exp.results_type == "loading_rates":
            if len(data) != 2:
                continue
            data = data[1]

        data = numpy.asarray(data, dtype=float).ravel()
        data = data[numpy.isfinite(data)]
        if data.size == 0:
            continue

        if mode == "groups":
            names.append(item.name)
        else:
            names.append(get_short_name(item))
        values.append(data)
        colors.append(item.color)

    return names, values, colors
