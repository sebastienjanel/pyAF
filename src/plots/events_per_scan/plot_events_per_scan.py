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

"""Plots events per scan."""

from ...tools import results_sorting
from ... import shared
from ..plot_main import MainPlot


class PlotEventsPerScan(MainPlot):
    """Plots events per scan.

    For each event a red dot is plotted, so that the user can compare the
    number of events per scan. The plotted events are coming from the
    results single table, and are therefore filtered by the user, so that
    only the wanted events are taken into account.
    """

    def __init__(self, parent):
        super().__init__(parent)

        ev_per_scan = []

        res_util = results_sorting.GetResults(force_type="events_per_curve")

        for i in range(len(shared.exp.results_list)):
            data_id = shared.exp.results_list[i].data_id
            if shared.exp.results_list[i].display:
                res_util.load_data(data_id=data_id)
                nbr = res_util.count_events_per_scan()
                ev_per_scan.append([i, nbr, shared.exp.list[data_id].filename])

        xlabels_indexes = []
        xlabels_names = []

        max_val = 0
        for i in range(len(ev_per_scan)):
            self.axes.plot(i + 1, ev_per_scan[i][1], "r.", markersize=15.0)

            xlabels_indexes.append(i + 1)
            xlabels_names.append(ev_per_scan[i][2])

            if ev_per_scan[i][1] > max_val:
                max_val = ev_per_scan[i][1]

        self.axes.set_xlim([0, len(ev_per_scan) + 1])
        self.axes.set_ylim([0, max_val + 30])
        self.axes.set_xticks(xlabels_indexes)
        self.axes.set_xticklabels(xlabels_names, rotation="vertical",
                                  fontproperties=self.font)
        self.axes.set_ylabel("Events per scan", fontproperties=self.font)

        # Do the actual plotting, is defined in MainPlot class
        self.start_plot()
