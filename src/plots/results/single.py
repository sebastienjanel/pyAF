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

"""Plots the single results."""

import numpy

from .plot_results import PlotResults
from ... import shared
from ...tools import stat_tools


class PlotSingleResults(PlotResults):
    """Plots the single results."""

    def __init__(self, parent):
        super().__init__(parent)

        found_one_label_single = False
        self.lr_fits = shared.exp.list_dfs_single_fits

        self.hist_labels = []
        datas_x = []
        datas_y = []

        ranges = []

        themax = 1

        # Sort the data
        if shared.single_data is not None:
            if shared.exp.results_type != "loading_rates":
                self.pdfs = []
                self.data_x = []

                for i in range(len(shared.single_data)):
                    result = shared.exp.results_list[i]

                    if shared.single_factors is None:
                        factor = 1.0
                    elif shared.exp.norm_hist_single:
                        factor = 1.0
                    elif shared.single_factors[i]:
                        factor = shared.single_factors[i]

                    if len(shared.single_data[i]) > 0 and result.display:
                        found_one_label_single = True
                        datas_x.append(shared.single_data[i])
                        ma = numpy.amax(shared.single_data[i])
                        if ma > themax:
                            themax = ma

                        r = stat_tools.get_range(
                            "single", themax, data=shared.single_data[i])
                        ranges.append(r)

                        self.hist_labels.append(result.name)
                        self.colors.append(result.color)

                        # Density functions
                        if shared.single_pdfs_x is not None:
                            if i >= len(shared.single_pdfs_x) or i >= len(shared.single_pdfs_y):
                                # print(f"Skipping index {i} because PDFs are missing.")
                                continue  # Avoid accessing out-of-range indexes

                            self.pdfs.append([
                                shared.single_pdfs_x[i],
                                shared.single_pdfs_y[i],
                                factor])

            else:
                self.data_x = []
                self.data_y = []

                for i in range(len(shared.single_data)):
                    result = shared.exp.results_list[i]

                    if shared.single_data[i] != [[], []] and result.display:
                        found_one_label_single = True

                        ma = numpy.amax(shared.single_data[i][1])
                        if ma > themax:
                            themax = ma

                        r = stat_tools.get_range(
                            "single", themax, data=shared.single_data[i][1])
                        ranges.append(r)

                        # Loading rates
                        datas_x.append(shared.single_data[i][1])
                        # Forces
                        datas_y.append(shared.single_data[i][0])

                        self.hist_labels.append(result.name)
                        self.colors.append(result.color)

        self.prepare_min_max_x(ranges)

        # Y limits
        if shared.exp.hist_single_y_mode == "manual":
            self.hist_y_min = shared.exp.hist_single_min_y
            self.hist_y_max = shared.exp.hist_single_max_y

        if shared.exp.results_type != "loading_rates":
            if datas_x:
                hist_single_bins = stat_tools.get_bins("single")
                norm = shared.exp.norm_hist_single

                for data in datas_x:
                    hist, bins = numpy.histogram(
                        data,
                        bins=hist_single_bins,
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

                    if shared.exp.hist_lr_single_display_mode == "points":
                        self.lr_error_x.append(numpy.std(data_x) / 2.0)
                        self.lr_error_y.append(numpy.std(data_y) / 2.0)
                        self.data_x.append(numpy.median(data_x))
                        self.data_y.append(numpy.median(data_y))
                    else:
                        self.data_x.append(data_x)
                        self.data_y.append(data_y)

        # Disable legend if none
        self.display_legend = shared.exp.display_legend_hist_single
        if not found_one_label_single:
            self.display_legend = False

        self.plot_fig()
