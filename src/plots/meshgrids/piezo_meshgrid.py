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

"""Plots a piezo extension meshgrid."""

from ...tools import math_tools
from .plot_meshgrid import PlotMeshGrid


class PlotPiezoMeshGrid(PlotMeshGrid):
    """Plots a piezo extension meshgrid.

    This meshgrid corresponds to the last position of the fully extended piezo
    This is called the Piezo height by Bruker.
    """

    def __init__(self, parent):
        super().__init__(parent)

        if self.data.flatten_applied and self.plot_type == "meshgrid":
            # Only apply the flatten for the meshgrid on the results tab
            self.temp_array[:] = math_tools.flatten(
                self.data.piezo_image,
                self.data.applied_flatten_order)
        else:
            self.temp_array[:] = self.data.piezo_image

        if self.plot_type == "meshgrid_data":
            # Overwrite colormap to have always a copper colormap
            self.colormap = self.colortables.colortables_mpl_list[0]
            self.colortable_max_value = self.data.max_piezo
            self.colortable_min_value = 0

        self.factor = 1.0
        self.colorbarlabel = "Piezo height [nm]"

        self.plot_fig()
