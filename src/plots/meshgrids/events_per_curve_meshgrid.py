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

"""Events per curve meshgrid."""

from .plot_meshgrid import PlotMeshGrid


class PlotEventsPerCurveMeshGrid(PlotMeshGrid):
    """Plot a meshgrid with each pixel showing the number of events per curve.

    The events per curve can be filtered by the user. The colorbar's label
    displays the currently selected slope and distance limit.
    """

    def __init__(self, parent):
        super().__init__(parent)

        self.temp_array[:] = self.data.events_per_curve

        self.factor = 1

        left = self.data.events_results_filter_dist_left
        right = self.data.events_results_filter_dist_right
        middle = self.data.events_results_filter_dist_keep_middle

        if middle:
            text = "(" + str(left) + " < event < " + str(right) + ")"
        else:
            text = "(" + str(left) + " > event > " + str(right) + ")"

        # print(left, right, middle, text)

        if left != 0 and right != 0:
            self.colorbarlabel = "Events per curve" + text
        else:
            self.colorbarlabel = "Events per curve"

        self.plot_fig()
