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

"""Plots the breakthroughs of supported lipid bilayers."""

import numpy
from .plot_meshgrid import PlotMeshGrid
from ...tools import map_values


class PlotBreakthroughMeshGrid(PlotMeshGrid):
    """Plots a map of the breakthroughs (first breakthrough of each curve):
    force, jump of the separation, or number per curve.

    The breakthroughs are computed in the Breakthrough tab of the compute
    tab.
    """

    def __init__(self, parent):
        super().__init__(parent)

        mt = self.mesh_type
        self.factor = map_values.MAPS[mt][2]
        self.colorbarlabel = map_values.MAPS[mt][1]
        values = getattr(self.data, mt)

        # The curves without breakthrough have a different color
        self.temp_array = numpy.ma.masked_invalid(
            numpy.asarray(values, dtype=float))

        self.plot_fig()
