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

"""Main results plotting module."""
import random

import numpy
import textwrap

import pandas as pd
from matplotlib.pyplot import xlabel, ylabel, margins, xticks
from seaborn import axes_style
import numpy as np
import matplotlib.pyplot as pyplot
import seaborn as sns

from ... import shared
from ..plot_main import MainPlot


class PlotResults(MainPlot):
    """Parent class of the plotting of the single and grouped results."""

    def __init__(self, parent):
        super().__init__(parent)

        self.data_x = None
        self.data_y = None
        self.lr_fits = None
        self.lr_error_x = None
        self.lr_error_y = None
        self.pdfs = None
        self.colors = []

        self.hist_x_min = None
        self.hist_x_max = None
        self.hist_y_max = None

        self.display_legend = False

        rt = shared.exp.results_type
        pt = self.plot_type

        self.plot_selected = shared.exp.hist_single_plot_selected

        # Common labels for the histograms
        if rt == "stiffness" or rt == "stiffness_corr":
            # Labels for stiffness or slope (can be normal or Log10)
            # Check if we have some slopes, in this case the label is different
            found_slope = False
            found_stiffness = False
            if pt == "results_single":
                for res in shared.exp.results_list:
                    data = shared.exp.list[res.data_id]
                    model = data.used_stiffness_model_selected
                    displayed = res.display
                    if model == 3 and displayed:
                        # model 3 = slope
                        found_slope = True
                    if model != 3 and displayed:
                        found_stiffness = True

            elif pt == "results_groups":
                for res in shared.exp.results_list:
                    data = shared.exp.list[res.data_id]
                    model = data.used_stiffness_model_selected
                    group = res.group  # Group 0 = no group
                    if model == 3 and group != 0:
                        # model 3 = slope
                        found_slope = True
                    if model != 3 and group != 0:
                        found_stiffness = True

            elif pt == "results_experiment":
                for res in shared.exp.results_list:
                    data = shared.exp.list[res.data_id]
                    model = data.used_stiffness_model_selected
                    condition = res.condition  # Condition 0 = no condition
                    if model == 3 and condition != 0:
                        # model 3 = slope
                        found_slope = True
                    if model != 3 and condition != 0:
                        found_stiffness = True

            if self.plot_selected == "box":
                self.xlabel = ""
                pt = None
                if found_slope and found_stiffness:
                    if shared.exp.hist_log is False:
                        self.ylabel = "E [kPa] and Apparent stiffness [N/m]"
                    else:
                        self.ylabel = "Log(E [Pa]) and Log([N/m])"
                elif found_slope and found_stiffness is False:
                    if shared.exp.hist_log is False:
                        self.ylabel = "Apparent stiffness [N/m]"
                    else:
                        self.ylabel = "Log([N/m])"
                elif found_slope is False and found_stiffness:
                        if shared.exp.hist_log is False:
                            self.ylabel = "E [kPa]"
                        else:
                            self.ylabel = "Log(E [Pa])"
                else:
                    self.ylabel = ""

            else:
                if found_slope and found_stiffness:
                    if shared.exp.hist_log is False:
                        self.xlabel = "E [kPa] and Stiffness [N/m]"
                    else:
                        self.xlabel = "Log(E [Pa]) and Log([N/m])"
                elif found_slope and found_stiffness is False:
                    if shared.exp.hist_log is False:
                        self.xlabel = "Stiffness [N/m]"
                    else:
                        self.xlabel = "Log([N/m])"
                elif found_slope is False and found_stiffness:
                        if shared.exp.hist_log is False:
                            self.xlabel = "E [kPa]"
                        else:
                            self.xlabel = "Log(E [Pa])"
                else:
                    self.xlabel = ""

        elif rt == "work":
            if self.plot_selected == "box":
                self.xlabel = ""
            else:
                self.xlabel = "Detachment work [fJ]"
        elif rt == "rupture_force":
            if self.plot_selected == "box":
                self.xlabel = ""
            else:
                self.xlabel = "Detachment force [pN]"
        elif rt == "events_forces":
            if self.plot_selected == "box":
                self.xlabel = ""
            else:
                self.xlabel = "Force [pN]"
        elif rt == "events_per_curves":
            if self.plot_selected == "box":
                self.xlabel = ""
            else:
                self.xlabel = "Events per curve"
        elif rt == "events_rupture_force":
            if self.plot_selected == "box":
                self.xlabel = ""
            else:
                self.xlabel = "Detachment force [pN]"
        elif rt == "events_distance":
            if self.plot_selected == "box":
                self.xlabel = ""
            else:
                self.xlabel = "Distance [nm]"
        elif rt == "loading_rates":
            if self.plot_selected == "box":
                self.xlabel = ""
            else:
                self.xlabel = "Loading rate [pN/s]"

        if pt == "results_single":
            if self.plot_selected == "box":
                if rt == "work":
                    self.ylabel = "Detachment work [fJ]"
                elif rt == "rupture_force":
                    self.ylabel = "Detachment force [pN]"
                elif rt == "events_forces":
                    self.ylabel = "Force [pN]"
                elif rt == "events_per_curves":
                    self.ylabel = "Events per curve"
                elif rt == "events_rupture_force":
                    self.ylabel = "Detachment force [pN]"
                elif rt == "events_distance":
                    self.ylabel = "Distance [nm]"
                elif rt == "loading_rates":
                    self.ylabel = "Loading rate [pN/s]"
            else:
                if shared.exp.norm_hist_single:
                    self.ylabel = r"Probability density"
                else:
                    self.ylabel = r"Frequency"

        elif pt == "results_groups":
            if shared.exp.norm_hist_groups:
                self.ylabel = r"Probability density"
            else:
                self.ylabel = r"Frequency"

        elif pt == "results_experiment":
            if shared.exp.norm_hist_experiment:
                self.ylabel = r"Probability density"
            else:
                self.ylabel = r"Frequency"

    def plot_fig(self):
        """Plots the results and asked informations."""
        pt = self.plot_type
        sm = shared.exp.hist_lr_single_display_mode
        gm = shared.exp.hist_lr_groups_display_mode
        sf = shared.exp.hist_single_display_pdf_fit
        gf = shared.exp.hist_groups_display_pdf_fit

        # Leave some space on the right for the legend
        box = self.axes.get_position()
        self.axes.set_position([box.x0, box.y0, box.width * 0.75, box.height])

        box_axes = None

        if self.plot_selected == "box":
            box_axes = self.fig.add_subplot(111)

        if shared.exp.results_type != "loading_rates"\
                and self.data_x is not None:
            for i in range(len(self.data_x)):
                hist = self.data_x[i][0]
                bins = self.data_x[i][1]
                width = (bins[1] - bins[0])
                center = (bins[:-1] + bins[1:]) / 2

                if self.hist_labels[i] is not None:
                    text = "\n".join(textwrap.wrap(self.hist_labels[i], 24))
                else:
                    text = ""

                # Conditions for display another plot according to the combobox "plot_selector"
                # For a histogram
                if self.plot_selected == "hist":
                    self.axes.bar(
                        center,
                        hist,
                        align="center",
                        width=width,
                        alpha=0.5,
                        color=self.colors[i],
                        label=text,
                        edgecolor="white"
                    )

                # For a box plot
                if self.plot_selected == "box":
                    # Adjust the data format for the Seaborn boxplot in a Dataframe object
                    dataframe = pd.DataFrame(data=hist, columns=['values']).astype(float)
                    dataframe['file_index'] = "#" + str(i+1)

                    box = sns.boxplot(
                        data=dataframe,
                        x='file_index',
                        y="values",
                        color=self.colors[i],
                        ax=box_axes,
                        width=0.5,
                        linecolor="white",
                        orient="v",
                        showfliers=False,
                    )

                    # Remove tick labels on the main axes
                    self.axes.set(xticklabels=[])
                    self.axes.set(yticklabels=[])

                    # Remove ticks on the main axes
                    self.axes.tick_params(bottom=False)
                    self.axes.tick_params(left=False)

                    # Extend the space between the Y-axe and his label
                    self.axes.tick_params(axis='y', pad=30)

                    # Remove the label on the box plot's axes
                    box.set(xlabel=None)
                    box.set(ylabel=None)

                    #print(dataframe)

        elif shared.exp.results_type == "loading_rates"\
                and self.data_x is not None:
            for i in range(len(self.data_x)):
                values_x = self.data_x[i]
                values_y = self.data_y[i]

                if self.hist_labels[i] is not None:
                    text = "\n".join(textwrap.wrap(self.hist_labels[i], 26))
                else:
                    text = ""

                self.axes.semilogx(
                    values_x,
                    values_y,
                    "o",
                    markersize=5,
                    c=self.colors[i],
                    label=text)

                # Error bars
                if (pt == "results_single" and sm == "points") or \
                        (pt == "results_groups" and gm == "points"):
                    self.axes.errorbar(values_x, values_y,
                                       xerr=self.lr_error_x[i],
                                       yerr=self.lr_error_y[i],
                                       color="k", label="_nolegend_'")

            self.plot_loading_rate_fits()

        # X and Y scales
        if self.hist_x_max is not None and self.hist_x_min is not None:
            self.axes.set_xlim(self.hist_x_min, self.hist_x_max)
        if self.hist_y_max is not None and self.hist_y_min is not None:
            self.axes.set_ylim(self.hist_y_min, self.hist_y_max)

        if self.display_legend:
            self.axes.legend(
                loc="center left", bbox_to_anchor=(1.0, 0.5), prop = self.font)

        # Plot the density function
        if (pt == "results_single" and sf) or (pt == "results_groups" and gf):  # Implement for experiment
            if self.data_x is not None and self.pdfs is not None and \
                    self.pdfs != []:
                for pdf_x, pdf_y, factor in self.pdfs:
                    if isinstance(pdf_x, numpy.ndarray) and pdf_x.size > 0:
                        self.axes.plot(
                            pdf_x, pdf_y * factor, "r-", linewidth=3)

        # Plot the labels
        self.plot_labels()

        # Do the actual plotting, is defined in MainPlot class
        self.start_plot()

    def plot_loading_rate_fits(self):
        """Plot the fit for the DFS experiment."""
        if self.lr_fits is not None and self.lr_fits != []:
            for fit in self.lr_fits:
                if fit[4]:
                    x = [fit[2], fit[3]]
                    y = [fit[0] * fit[2] + fit[1], fit[0] * fit[3] + fit[1]]
                    self.axes.plot(x, y, "r", linewidth=2)

    def prepare_min_max_x(self, ranges):
        """Method used to prepare the min max values for the x axis."""
        if ranges != []:
            ranges_min = []
            ranges_max = []
            for r in ranges:
                ranges_min.append(r[0])
                ranges_max.append(r[1])
        else:
            ranges_min = 0
            ranges_max = 1

        self.hist_x_min = numpy.amin(ranges_min)
        self.hist_x_max = numpy.amax(ranges_max)
