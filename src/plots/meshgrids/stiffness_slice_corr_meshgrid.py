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

"""Plots the meshgrid with the corrected stiffness (z slice)."""

import numpy
from .plot_meshgrid import PlotMeshGrid


class PlotStiffnessSliceCorrMeshGrid(PlotMeshGrid):
    """Plots the meshgrid with the corrected stiffness (z slice).

    Will search for all the values found at a defined height.
    """

    def __init__(self, parent):
        super().__init__(parent)

        # Stiffness (at a defined height) + Corrected
        height = self.data.stiffness_slice_depth  # nm

        for i in range(self.nbr_pixels_x):
            for j in range(self.nbr_pixels_y):
                dist = self.data.topography[i][j] - height
                if dist > 0 and int(dist / self.data.indentation_step) < \
                        len(self.data.stiffness_corrected_array[i][j]):
                    depth = int(dist / self.data.indentation_step)
                    self.temp_array[i][j] = \
                        self.data.stiffness_corrected_array[i][j][depth]
                else:
                    self.temp_array[i][j] = numpy.nan

        self.temp_array = numpy.array(self.temp_array)
        self.temp_array = numpy.ma.masked_where(
            numpy.isnan(self.temp_array), self.temp_array)

        self.factor = 1e-3
        self.colorbarlabel = "E [kPa]"

        self.plot_fig()
