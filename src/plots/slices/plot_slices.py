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

"""Plots a slice in the stiffness tomography."""

import numpy
import copy
import math
from itertools import zip_longest
from matplotlib.collections import PolyCollection, LineCollection
from ...tools.colortables import ColorTables
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ... import shared
from ..plot_main import MainPlot



class PlotSlices(MainPlot):
    """Plots a slice in the stiffness tomography.

    A matplotlib PolyCollection is used to have fast plotting. All the sizes of
    the different voxels are taken into account (for example for diagonal slices
    you have three voxel sizes : x, y, and xydiagonal).
    """

    def __init__(self, parent):
        super().__init__(parent)

        # Force the aspect of the scales depending on the choice of the
        # user.
        if self.data.slice_aspect_ratio:
            self.axes.set_aspect("equal")
        else:
            self.axes.set_aspect("auto")

        self.display_colorbar = True

        # Colormap
        self.colortables = ColorTables()
        self.colortableid = self.data.color_opts_stiffness[0]
        color_table_list = self.colortables.colortables_mpl_list
        self.colormap = color_table_list[self.colortableid]
        self.color_saturation = self.data.color_saturation
        self.color_negative = self.data.color_negative
        self.color_nan = self.data.color_nan
        self.colortable_min_value = self.data.color_opts_stiffness[1]
        self.colortable_middle_value = self.data.color_opts_stiffness[3]
        self.colortable_max_value = self.data.color_opts_stiffness[2]

        self.color_bar = None
        self.colortable = None

        if self.data.used_stiffness_model_selected == 3:
            # Slope
            self.factor = 1.0
            self.colorbarlabel = "Slope [N/m]"
        else:
            self.factor = 1e-3
            self.colorbarlabel = "E [kPa]"
        self.xlabel = "X [\u03bcm]"
        self.ylabel = "Z [\u03bcm]"

        self.plot_colorbar = True

        self.slices_z_mode = "auto"
        if self.data.slices_z_mode == "manual":
            self.slices_z_mode = "manual"
            self.z_min = self.data.slices_z_min / 1000.0
            self.z_max = self.data.slices_z_max / 1000.0

        # The data
        E_array = []
        X_array = []
        Z_array = []
        nbr_pixels_x = self.data.nbr_pixels_x
        nbr_pixels_y = self.data.nbr_pixels_y
        self.scan_size_x = self.data.scan_size_x
        self.scan_size_y = self.data.scan_size_y
        x_size = self.scan_size_x / nbr_pixels_x
        y_size = self.scan_size_y / nbr_pixels_y

        old_pos_x = 0
        old_pos_y = 0

        z_size = self.data.indentation_step
        self.custom_max_x = 0  # For custom profile

        # Prepare the arrays with the data, depending on the direction

        if self.data.slice_view_from_direction == 0:
            pos_top = self.data.slice_position_top

            title = copy.deepcopy(self.data.title_slice_top)
            self.title = title.replace("%i", str(pos_top + 1))

            for i in range(nbr_pixels_x):
                if self.data.display_slice_with_topo:
                    ztop = self.data.topography[i][pos_top]
                else:
                    ztop = 0

                e = []
                z = []

                for k in range(len(self.data.stiffness_array[i][pos_top])):
                    e.append(self.data.stiffness_array[i][pos_top][k])
                    z.append(ztop - k * z_size)

                X_array.append(i * x_size)
                E_array.append(numpy.array(e))
                Z_array.append(numpy.array(z))

            X_array.append(X_array[-1] + x_size)

        elif self.data.slice_view_from_direction == 1:
            pos_bot = self.data.slice_position_bottom

            title = copy.deepcopy(self.data.title_slice_bottom)
            self.title = title.replace("%i", str(pos_bot + 1))

            for i in range(nbr_pixels_x):
                if self.data.display_slice_with_topo:
                    ztop = self.data.topography[i][pos_bot]
                else:
                    ztop = 0

                e = []
                z = []

                for k in range(len(self.data.stiffness_array[i][pos_bot])):
                    e.append(self.data.stiffness_array[i][pos_bot][k])
                    z.append(ztop - k * z_size)

                X_array.append(i * x_size)
                E_array.append(numpy.array(e))
                Z_array.append(numpy.array(z))

            X_array.append(X_array[-1] + x_size)

        elif self.data.slice_view_from_direction == 2:
            pos_left = self.data.slice_position_left

            title = copy.deepcopy(self.data.title_slice_left)
            self.title = title.replace("%i", str(pos_left + 1))

            for i in range(nbr_pixels_y):
                if self.data.display_slice_with_topo:
                    ztop = self.data.topography[pos_left][i]
                else:
                    ztop = 0

                e = []
                z = []

                for k in range(len(self.data.stiffness_array[pos_left][i])):
                    e.append(self.data.stiffness_array[pos_left][i][k])
                    z.append(ztop - k * z_size)

                X_array.append(i * y_size)
                E_array.append(numpy.array(e))
                Z_array.append(numpy.array(z))

            X_array.append(X_array[-1] + y_size)

        elif self.data.slice_view_from_direction == 3:
            pos_right = self.data.slice_position_right

            title = copy.deepcopy(self.data.title_slice_right)
            self.title = title.replace("%i", str(pos_right + 1))

            for i in range(nbr_pixels_y):
                if self.data.display_slice_with_topo:
                    ztop = self.data.topography[pos_right][i]
                else:
                    ztop = 0

                e = []
                z = []

                for k in range(len(self.data.stiffness_array[pos_right][i])):
                    e.append(self.data.stiffness_array[pos_right][i][k])
                    z.append(ztop - k * z_size)

                X_array.append(i * y_size)
                E_array.append(numpy.array(e))
                Z_array.append(numpy.array(z))

            X_array.append(X_array[-1] + y_size)

        elif self.data.slice_view_from_direction == 4:
            self.title = self.data.title_slice_profile
            profile_x = self.data.profile_list[0][0]
            profile_y = self.data.profile_list[0][1]
            nbr_pixels_x = len(profile_x)

            for i in range(len(profile_x)):
                pos_x = int(profile_x[i])
                pos_y = int(profile_y[i])
                if i != 0:
                    # Diagonal
                    d = math.sqrt((pos_x * x_size - old_pos_x * x_size) ** 2 +
                                  (pos_y * y_size - old_pos_y * y_size) ** 2)

                    self.custom_max_x += d

                old_pos_x = pos_x
                old_pos_y = pos_y

                if self.data.display_slice_with_topo:
                    ztop = self.data.topography[pos_x][pos_y]
                else:
                    ztop = 0
                values = self.data.stiffness_array[pos_x][pos_y]

                e = []
                z = []

                for k in range(len(values) - 2):
                    e.append(values[k])
                    z.append(ztop - k * z_size)

                E_array.append(numpy.array(e))
                Z_array.append(numpy.array(z))
                X_array.append(self.custom_max_x)

        elif self.data.slice_view_from_direction == 5:
            self.title = self.data.title_slice_profile
            profile_x = self.data.profile_list[0][0]
            profile_y = self.data.profile_list[0][1]
            nbr_pixels_x = len(profile_x)

            for i in range(len(profile_x)):
                pos_x = int(profile_x[i])
                pos_y = int(profile_y[i])

                if i != 0:
                    d = math.sqrt((old_pos_x * x_size - pos_x * x_size) ** 2 +
                                  (old_pos_y * y_size - pos_y * y_size) ** 2)

                    self.custom_max_x += d

                old_pos_x = pos_x
                old_pos_y = pos_y

                if self.data.display_slice_with_topo:
                    ztop = self.data.topography[pos_x][pos_y]
                else:
                    ztop = 0
                values = self.data.stiffness_array[pos_x][pos_y]

                e = []
                z = []

                for k in range(len(values) - 2):
                    e.append(values[k])
                    z.append(ztop - k * z_size)

                E_array.append(numpy.array(e))
                Z_array.append(numpy.array(z))
                X_array.append(self.custom_max_x)

            # Invert positions in this case
            X_array = X_array[::-1]

        # Pad shorter lists with NaN
        E_array_padded = list(zip_longest(*E_array, fillvalue=numpy.nan))
        E_array = numpy.array(E_array_padded).T  # Transpose back
        Z_array_padded = list(zip_longest(*Z_array, fillvalue=numpy.nan))
        Z_array = numpy.array(Z_array_padded).T  # Transpose back

        self.x_size = x_size
        self.z_size = z_size
        self.nbr_pixels_x = nbr_pixels_x
        self.nbr_pixels_y = nbr_pixels_y
        self.y_height = numpy.amax(self.data.topography)
        self.data_stiffness_slices = [numpy.array(E_array),
                                      numpy.array(X_array),
                                      numpy.array(Z_array)]

        self.plot_fig()

        self.plot_labels()

        self.axes.set_title(self.title, fontproperties=self.font)

        # Do the actual plotting, is defined in MainPlot class
        self.start_plot()

    def plot_fig(self):
        """Plot the figure."""
        self.colortable = ColorTables(
            self.data.colortableid,
            self.data.color_saturation,
            self.data.color_negative,
            self.data.colortable_max_value *
            self.factor,
            self.data.colortable_middle_value *
            self.factor,
            mode="scalarMap",
            color_nan=self.data.color_nan)

        elasticity = self.data_stiffness_slices[0] * self.factor
        X = numpy.array(self.data_stiffness_slices[1]) / 1000.0
        Z = numpy.array(self.data_stiffness_slices[2]) / 1000.0
        z_size = self.z_size / 1000.0
        z_height = math.ceil(self.y_height / 1000.0)

        direct = self.data.slice_view_from_direction
        if direct == 0 or direct == 1:
            xmax = self.data.scan_size_x / 1000.0
        elif direct == 2 or direct == 3:
            xmax = self.data.scan_size_y / 1000.0
        else:
            xmax = self.custom_max_x / 1000.0
        xmin = 0

        if self.slices_z_mode == "manual":
            self.axes.axis([xmin, xmax, self.z_min, self.z_max])
        else:
            if self.data.display_slice_with_topo is False:
                self.axes.axis([xmin, xmax, -z_height, 0])
            else:
                self.axes.axis([xmin, xmax, 0, z_height])

        # Scales are inverted when view is from top or left
        if direct == 0:
            self.axes.set_xlim(self.data.scan_size_x / 1000.0, 0)
        elif direct == 2:
            self.axes.set_xlim(self.data.scan_size_y / 1000.0, 0)

        # Custom profile, define a special max
        if direct == 4 or direct == 5:
            self.axes.set_xlim(xmin, xmax)
            nbr_polys = len(Z) - 1
        else:
            nbr_polys = len(Z)

        verts = []
        colors = []
        j = 0

        for i in range(nbr_polys):
            for j in range(len(Z[i])):
                z = Z[i][j]
                e = elasticity[i][j]
                x = X[i]

                dist = X[i + 1] - X[i]

                tri = ((x, z), (x, z - z_size), (x + dist, z))
                verts.append(tri)
                tri = ((x, z - z_size), (x + dist, z), (x + dist, z - z_size))
                verts.append(tri)
                colors.append(self.colortable.get_color_as_list(e))
                colors.append(self.colortable.get_color_as_list(e))

        verts_np = numpy.array(verts)
        pcoll_np = PolyCollection(verts_np, closed=True,
                                  facecolor=colors, edgecolor=colors)
        self.axes.add_collection(pcoll_np, autolim=True)

        # Plot the color bar
        if self.plot_colorbar:
            m = self.colortable.scalarMap
            m.set_array(elasticity)

            # Option for colorbar alignement to the meshgrid. See :
            # http://matplotlib.sourceforge.net/users/tight_layout_guide.html
            # (part on the colorbar caveat)
            divider = make_axes_locatable(self.fig.gca())
            cax = divider.append_axes("right", "5%", pad="3%")
            self.color_bar = self.fig.colorbar(m, cax=cax)
            # Color bar label
            if self.colorbarlabel is not None:
                self.color_bar.set_label(self.colorbarlabel,
                                         fontproperties=self.font)

            # Pretty plot the colorbar (mainly for export)
            # http://stackoverflow.com/questions/15003353/
            self.color_bar.solids.set_edgecolor("face")

        # Plot a red line at a defined height (slice)
        if shared.exp.display_line_stiffness_slice_at_height:
            verts = []
            x1 = 0
            x2 = self.x_size * 1000.0
            z1 = self.data.stiffness_slice_depth / 1000.0
            z2 = self.data.stiffness_slice_depth / 1000.0

            verts.append(((x1, z1), (x2, z2)))
            verts = numpy.array(verts)
            coll = LineCollection(verts, linewidths=2,
                                  colors="r", linestyle="solid")
            self.axes.add_collection(coll)
