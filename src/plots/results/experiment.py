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

"""Plots the experiment results."""

import numpy
from ... import shared
from .plot_results import PlotResults
from ...tools import stat_tools


class PlotExperimentResults(PlotResults):
    """Plots the experiment results."""

    def __init__(self, parent):
        super().__init__(parent)

        found_one_label_experiment = False
        self.lr_fits = shared.exp.list_dfs_groups_fits  # To be implemented for experiment plots

        self.hist_labels = []
        datas_x = []
        datas_y = []

        ranges = []

        themax = 1

        # Sort the data
        if shared.conditions_data is not None:
            if shared.exp.results_type != "loading_rates":
                self.pdfs = []
                self.data_x = []

                for i in range(len(shared.conditions_data)):
                    result = shared.exp.conditions_list[i]

                    if shared.conditions_factors is None:
                        factor = 1.0
                    elif shared.exp.norm_hist_groups:
                        factor = 1.0
                    elif shared.conditions_factors[i]:
                        factor = shared.conditions_factors[i]

                    if len(shared.conditions_data[i]) > 0 and result.display:
                        found_one_label_experiment = True
                        datas_x.append(shared.conditions_data[i])
                        ma = numpy.amax(shared.conditions_data[i])
                        if ma > themax:
                            themax = ma

                        r = stat_tools.get_range(
                            "experiment", themax, shared.conditions_data[i])
                        ranges.append(r)

                        self.hist_labels.append(result.name)
                        self.colors.append(result.color)

                        # Density functions
                        if shared.conditions_pdfs_x is not None:
                            self.pdfs.append([
                                shared.conditions_pdfs_x[i],
                                shared.conditions_pdfs_y[i],
                                factor])

            else:
                self.data_x = []
                self.data_y = []

                for i in range(len(shared.conditions_data)):
                    result = shared.exp.conditions_list[i]

                    if shared.conditions_data[i] and result.display:
                        found_one_label_experiment = True

                        ma = numpy.amax(shared.conditions_data[i][1])
                        if ma > themax:
                            themax = ma

                        r = stat_tools.get_range(
                            "experiment", themax, shared.conditions_data[i][1])
                        ranges.append(r)

                        # Loading rates
                        datas_x.append(shared.conditions_data[i][1])
                        # Forces
                        datas_y.append(shared.conditions_data[i][0])

                        self.hist_labels.append(result.name)
                        self.colors.append(result.color)

        self.prepare_min_max_x(ranges)

        # Y limits
        if shared.exp.hist_experiment_y_mode == "manual":   # Implement for experiment
            self.hist_y_min = shared.exp.hist_experiment_min_y
            self.hist_y_max = shared.exp.hist_experiment_max_y

        if shared.exp.results_type != "loading_rates":
            if datas_x:
                hist_experiment_bins = stat_tools.get_bins("experiment")
                norm = shared.exp.norm_hist_groups  # Implement for experiment

                for data in datas_x:
                    hist, bins = numpy.histogram(
                        data,
                        bins=hist_experiment_bins,
                        range=r,
                        density=norm)

                    self.data_x.append([hist, bins])

        elif shared.exp.results_type == "loading_rates":
            self.ylabel = "Force [pN]"

            self.lr_error_x = []
            self.lr_error_y = []

            if datas_x:
                for i in range(len(datas_x)):
                    data_x = datas_x[i]
                    data_y = datas_y[i]

                    if shared.exp.hist_lr_experiment_display_mode == "points":
                        self.lr_error_x.append(numpy.std(data_x) / 2.0)
                        self.lr_error_y.append(numpy.std(data_y) / 2.0)
                        self.data_x.append(numpy.median(data_x))
                        self.data_y.append(numpy.median(data_y))
                    else:
                        self.data_x.append(data_x)
                        self.data_y.append(data_y)

        # Disable legend if none
        self.display_legend = shared.exp.display_legend_hist_experiment
        if not found_one_label_experiment:
            self.display_legend = False

        self.plot_fig()
