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

"""Main meshgrid plotting module."""

import numpy
import matplotlib
from matplotlib.colors import Normalize
from ... import widgets_list
from ... import shared
from PyQt5 import QtCore
import matplotlib.font_manager as fm
from matplotlib.collections import PolyCollection, LineCollection
from ...tools.colortables import ColorTables
from ...tools import misc_tools
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..plot_main import MainPlot
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar


class PlotMeshGrid(MainPlot):
    """Plots a meshgrid. Subclasses MainPlot."""

    def __init__(self, parent):
        super().__init__(parent)

        self.display_colorbar = True
        self.color_bar = None
        self.data_meshgrid = None
        self.colortable = None

        # For the tilt check meshgrid
        self.tilt_mask = None

        # Load the settings
        self.settings = QtCore.QSettings()

        # Colormap
        self.colortables = ColorTables()
        if self.plot_type == "multimeshgrid":
            colors = misc_tools.get_color_opts_by_meshgrid_type(
                self.data, shared.exp.multi_meshgrid_type)
            self.colormap = self.colortables.colortables_mpl_list[colors[0]]
            self.colortable_min_value = colors[1]
            self.colortable_middle_value = colors[3]
            self.colortable_max_value = colors[2]
            self.color_saturation = self.data.color_saturation
            self.color_negative = self.data.color_negative
            self.color_nan = self.data.color_nan
        else:
            self.colormap = \
                self.colortables.colortables_mpl_list[self.data.colortableid]
            self.color_saturation = self.data.color_saturation
            self.color_negative = self.data.color_negative
            self.color_nan = self.data.color_nan
            self.colortable_min_value = self.data.colortable_min_value
            self.colortable_middle_value = self.data.colortable_middle_value
            self.colortable_max_value = self.data.colortable_max_value

        self.nbr_pixels_x = self.data.nbr_pixels_x
        self.nbr_pixels_y = self.data.nbr_pixels_y

        # Get the sizes of the different elements and axes
        if self.nbr_pixels_x != 0 and self.data.scan_size_x != 0:
            if self.data.meshgrid_units == "um":
                meshfactor = 1000.0
            elif self.data.meshgrid_units == "nm":
                meshfactor = 1.0

            self.scan_size_x = self.data.scan_size_x / meshfactor
            self.scan_size_y = self.data.scan_size_y / meshfactor
            self.x_size = self.scan_size_x / float(self.nbr_pixels_x)
            self.y_size = self.scan_size_y / float(self.nbr_pixels_y)

            self.x_range = numpy.linspace(0, self.scan_size_x,
                                          num=self.data.nbr_pixels_x + 1)
            self.y_range = numpy.linspace(0, self.scan_size_y,
                                          num=self.data.nbr_pixels_y + 1)
        else:
            self.scan_size_x = self.nbr_pixels_x
            self.scan_size_y = self.nbr_pixels_y
            self.x_size = 1
            self.y_size = 1
            self.x_range = numpy.linspace(0, self.data.nbr_pixels_x,
                                          num=self.data.nbr_pixels_x + 1)
            self.y_range = numpy.linspace(0, self.data.nbr_pixels_y,
                                          num=self.data.nbr_pixels_y + 1)

        self.temp_array = numpy.zeros([self.nbr_pixels_x, self.nbr_pixels_y])

        # Set the X,Y labels
        if self.data.scan_size_x != 0:
            if self.data.meshgrid_units == "um":
                self.xlabel = "X [\u03bcm]"
                self.ylabel = "Y [\u03bcm]"
            elif self.data.meshgrid_units == "nm":
                self.xlabel = "X [nm]"
                self.ylabel = "Y [nm]"
        else:
            self.xlabel = "Positions"
            self.ylabel = "Positions"

        # The meshgrid is always a square, so we set the aspect.
        # This should also be kept as is for rectangular maps. If we don't
        # do this the pixels will get rectangular, and we don't want this to
        # happen.
        self.axes.set_aspect("equal", "box")

        # Multimeshgrids
        if self.plot_type == "multimeshgrid":
            self.x_label = None
            self.y_label = None
            self.display_colorbar = False

            # Hide the axes
            self.axes.get_xaxis().set_visible(False)
            self.axes.get_yaxis().set_visible(False)

            # The title will be plotted in a qwidget in the
            # widgets/multi_meshgrid.py file

    def plot_fig(self):
        """Do the actual plotting."""

        # Transpose self.data for matplotlib to have it the right way
        self.data_meshgrid = None
        if self.temp_array is not None:
            self.data_meshgrid = self.temp_array.transpose() * self.factor

        norm = None

        # Hacks for the displaying of the meshgrid
        if not self.parent.single:
            edgecolor = "None"
            linewidth = None
        else:
            # Draw edges with a certain width. 1/72 = 1 pdf unit
            # It seems that this value (30/72) is ideal for adobe illustrator
            # on OS X with matplotlib eps exporting.
            edgecolor = "face"
            linewidth = 30.0 / 72.0

        # Plot the meshgrid
        if self.colormap != "Rainbow 1" and self.colormap != "Cool-Warm" and \
                self.data_meshgrid is not None:
            if self.mesh_type == "events_per_curve":
                # Discrete self.colormap for the events
                self.colormap = matplotlib.cm.get_cmap(self.colormap)

                # Adjust the number of bins to fit the colormap's size
                ncolors = self.colormap.N
                bounds = list(range(int(self.colortable_min_value * self.factor),
                                    int(self.colortable_max_value * self.factor)))

                # Ensure bounds don't exceed colormap's size
                if len(bounds) > ncolors:
                    # If too many bins, reduce the range or step size
                    step_size = max(1, len(bounds) // ncolors)
                    bounds = bounds[::step_size]

                norm = matplotlib.colors.BoundaryNorm(bounds, self.colormap.N)
            else:
                # Normal self.colormap for the other meshgrid types
                self.colormap = matplotlib.cm.get_cmap(self.colormap)
                if self.color_saturation is not None:
                    self.colormap.set_over(self.color_saturation)
                    self.colormap.set_under(self.color_negative)

            # If some values are masked, define their color
            if self.color_nan is not None:
                self.colormap.set_bad(self.color_nan)

            # Define norm with vmin and vmax
            norm = Normalize(
                vmin=self.colortable_min_value * self.factor,
                vmax=self.colortable_max_value * self.factor
            )

            mesh = self.axes.pcolormesh(
                self.x_range,
                self.y_range,
                self.data_meshgrid,
                cmap=self.colormap,
                norm=norm,
                edgecolor=edgecolor,
                linewidth=linewidth)

        elif (self.colormap in ("Rainbow 1", "Cool-Warm")
              and self.data_meshgrid is not None):
            # Rainbow 1 and Cool-Warm
            self.colortable = ColorTables(
                self.data.colortableid,
                self.data.color_saturation,
                self.data.color_negative,
                self.data.colortable_max_value *
                self.factor,
                self.data.colortable_middle_value *
                self.factor,
                color_nan=self.color_nan)

            if self.mesh_type == "events_per_curve":
                # Discrete self.colormap for the events
                diff = (self.colortable_max_value -
                        self.colortable_min_value) * self.factor
                self.colormap, norm = self.colortable.get_color_as_cdict(
                    freq=int(diff) - 1)
                bounds = list(range(int(self.colortable_min_value * self.factor),
                               int(self.colortable_max_value * self.factor)))
                norm = matplotlib.colors.BoundaryNorm(bounds, self.colormap.N)
            else:
                # Normal self.colormap for the other meshgrid types
                self.colormap, norm = self.colortable.get_color_as_cdict()

            mesh = self.axes.pcolormesh(self.x_range, self.y_range,
                                        self.data_meshgrid,
                                        cmap=self.colormap,
                                        norm=norm,
                                        edgecolor=edgecolor,
                                        linewidth=linewidth)

        # Set the ranges
        self.axes.axis([0, self.x_range[-1], 0, self.y_range[-1]])

        # Plot the color bar
        if self.display_colorbar and self.noborder is False:
            # Option for colorbar alignment to the meshgrid
            # http://matplotlib.sourceforge.net/users/tight_layout_guide.html
            # (part on the colorbar caveat)
            divider = make_axes_locatable(self.fig.gca())
            cax = divider.append_axes("right", "5%", pad="3%")
            if norm is None:
                self.color_bar = self.fig.colorbar(mesh,
                                                   cax=cax, cmap=self.colormap)
            else:
                self.color_bar = self.fig.colorbar(
                    mesh,
                    cax=cax,
                    cmap=self.colormap,
                    norm=norm)
            # Color bar label
            self.color_bar.set_label(self.colorbarlabel,
                                     fontproperties=self.font)

            # Pretty plot the colorbar (mainly for export)
            # http://stackoverflow.com/questions/15003353/
            self.color_bar.solids.set_edgecolor("face")

            # Set the font
            for tick in self.color_bar.ax.get_yticklabels():
                tick.set_fontproperties(self.font)

        if self.settings.value("ClassicMeshgridStyle", True):
            self.plot_labels()

        pt = self.plot_type
        tp = (pt == "meshgrid" or pt == "meshgrid_data")
        st = self.settings.value("ClassicMeshgridStyle", True)

        if tp and st is False:
            vert = self.scan_size_y / 100.0

            # Remove x, y ticks and their labels
            self.axes.get_xaxis().set_ticklabels([])
            self.axes.get_xaxis().set_ticks([])
            self.axes.get_yaxis().set_ticklabels([])
            self.axes.get_yaxis().set_ticks([])

            # Add scale bar
            # Font (use the font I ship)
            path = misc_tools.get_app_path() + "/fonts/HelveticaNowText-Regular.ttf"
            ft = fm.FontProperties(fname=path, size=10)

            # label = u"1 \u03bcm" Label deactivated for the moment
            label = ""
            scalebar = AnchoredSizeBar(
                self.axes.transData,
                1,
                label,
                pad=1,
                loc=4,
                sep=0,
                borderpad=0,
                frameon=True,
                color="white",
                size_vertical=vert,
                fontproperties=ft)
            scalebar.patch.set(alpha=0)
            self.axes.add_artist(scalebar)

        # Do the actual plotting, is defined in MainPlot class
        self.start_plot()

    def draw_check_tilt_squares(self):
        """Show tilted correction.

        Draw squares on the meshgrid for pixels were there was no tilt
        correction.
    """

        if self.tilt_mask is None:
            return False

        all_verts = []
        for i in range(self.data.nbr_pixels_x):
            for j in range(self.data.nbr_pixels_y):
                if self.tilt_mask[i][j] == 0:
                    x1 = i * self.x_size + self.x_size
                    y1 = j * self.y_size + self.y_size
                    x2 = i * self.x_size + self.x_size
                    y2 = j * self.y_size
                    x3 = i * self.x_size
                    y3 = j * self.y_size
                    x4 = i * self.x_size
                    y4 = j * self.y_size + self.y_size
                    all_verts.append(list(zip([x1, x2, x3, x4], [y1, y2, y3, y4])))

                    verts = numpy.array([(x1, y1), (x3, y3)])
                    coll = LineCollection(
                        [verts],
                        linewidths=1,
                        colors="k",
                        linestyle="solid")
                    self.axes.add_collection(coll)

                    if self.parent.canvas is not None:
                        self.axes.draw_artist(coll)

                    verts = numpy.array([(x2, y2), (x4, y4)])
                    coll = LineCollection(
                        [verts],
                        linewidths=1,
                        colors="k",
                        linestyle="solid")
                    self.axes.add_collection(coll)

                    if self.parent.canvas is not None:
                        self.axes.draw_artist(coll)

    def draw_red_square(self):
        """Draws a red square for the selected curve."""
        pos = [shared.exp.meshgrid_click_xpos, shared.exp.meshgrid_click_ypos]
        if (shared.exp.meshgrid_display_red_square or
                self.plot_type == "small_meshgrid_compute"):
            x1 = pos[0] * self.x_size + self.x_size
            y1 = pos[1] * self.y_size + self.y_size
            x2 = pos[0] * self.x_size + self.x_size
            y2 = pos[1] * self.y_size
            x3 = pos[0] * self.x_size
            y3 = pos[1] * self.y_size
            x4 = pos[0] * self.x_size
            y4 = pos[1] * self.y_size + self.y_size

            verts = numpy.array([(x1, y1), (x2, y2), (x3, y3), (x4, y4)])
            coll = PolyCollection([verts], facecolors="none",
                                  edgecolors="r", linewidths=2)
            self.axes.add_collection(coll)

            if self.parent.canvas is not None:
                self.axes.draw_artist(coll)

    def draw_profiles(self, mode, profile=None,
                      profile_id=None, closing=False):
        """Draw all the profiles on the meshgrid."""
        if mode == "profile_preview":
            if profile_id is None:
                count = self.data.profile_count
                diff = int(count / len(self.data.profiles_color_names))
            else:
                count = profile_id
                diff = int(count / len(self.data.profiles_color_names))

            cur_id = count - diff * len(self.data.profiles_color_names)

            color = self.data.profiles_color_names[cur_id]
            self.draw_single_profile(profile, color)

        elif mode == "profile_preview_erase":
            self.draw_single_profile(profile, "k")

        elif mode == "all":
            if self.data.profile_list != [] and \
                    self.data.profile_list is not None:
                i = 0
                for profile in self.data.profile_list:
                    color = self.data.profiles_color_names[i]
                    self.draw_single_profile(profile, color)
                    i = i + 1
                    if i == len(self.data.profiles_color_names):
                        i = 0

            if widgets_list.widget_slices is not None and closing is False:
                color = "red"
                profile = [[], []]
                if self.data.slice_view_from_direction == 0:
                    for i in range(self.data.nbr_pixels_x):
                        profile[0].append(i)
                        profile[1].append(self.data.slice_position_top)
                    self.draw_single_profile(profile, color)
                elif self.data.slice_view_from_direction == 1:
                    for i in range(self.data.nbr_pixels_x):
                        profile[0].append(i)
                        profile[1].append(self.data.slice_position_bottom)
                    self.draw_single_profile(profile, color)
                elif self.data.slice_view_from_direction == 2:
                    for i in range(self.data.nbr_pixels_y):
                        profile[0].append(self.data.slice_position_left)
                        profile[1].append(i)
                    self.draw_single_profile(profile, color)
                elif self.data.slice_view_from_direction == 3:
                    for i in range(self.data.nbr_pixels_y):
                        profile[0].append(self.data.slice_position_right)
                        profile[1].append(i)
                    self.draw_single_profile(profile, color)

    def draw_single_profile(self, profile, color):
        """Draw a profile on the meshgrid."""
        verts = []
        for i in range(len(profile[0]) - 1):
            x1 = profile[0][i] * self.x_size + self.x_size / 2.0
            y1 = profile[1][i] * self.y_size + self.y_size / 2.0
            x2 = profile[0][i + 1] * self.x_size + self.x_size / 2.0
            y2 = profile[1][i + 1] * self.y_size + self.y_size / 2.0
            verts.append(((x1, y1), (x2, y2)))
        verts = numpy.array(verts)

        coll = LineCollection(verts, linewidths=2,
                              colors=color, linestyle="solid")

        self.axes.add_collection(coll)
        if self.parent.canvas is not None:
            self.axes.draw_artist(coll)

    def draw_roi(self, mode, roi=None, roi_mask=None):
        """Draw the ROIS on a meshgrid."""
        if mode == "all":
            if self.data.roi_list is not None:
                for roi in self.data.roi_list:
                    if roi.display and roi.values != []:
                        color = self.data.roi_color_names[roi.color]
                        self.draw_single_roi(roi.values, color)

        elif mode == "local_add" or mode == "preview_add":
            roi_id = self.data.roi_selected_row
            color = self.data.roi_color_names[self.data.roi_list[roi_id].color]
            self.draw_single_roi(roi, color, roi_mask)

        elif mode == "local_erase" or mode == "preview_erase":
            # --- Replot background ---
            self.colortable = ColorTables(
                self.data.colortableid,
                self.data.color_saturation,
                self.data.color_negative,
                self.data.colortable_max_value *
                self.factor,
                self.data.colortable_middle_value *
                self.factor,
                mode="scalarMap")
            colors = []

            # Get colors
            for i in range(len(roi)):
                if roi_mask[i]:
                    value = self.data_meshgrid[roi[i][1]][roi[i][0]]
                    color = self.colortable.get_color_as_list(value)
                    colors.append(color)
            self.draw_single_roi(roi, colors, roi_mask, alpha=1.0)

    def draw_single_roi(self, roi, color, roi_mask=None, alpha=0.5):
        """Draws a ROI on the meshgrid."""
        all_verts = []
        for i in range(len(roi)):
            square = roi[i]
            if square is not None and roi_mask is None or roi_mask[i]:
                x1 = square[0] * self.x_size + self.x_size
                y1 = square[1] * self.y_size + self.y_size
                x2 = square[0] * self.x_size + self.x_size
                y2 = square[1] * self.y_size
                x3 = square[0] * self.x_size
                y3 = square[1] * self.y_size
                x4 = square[0] * self.x_size
                y4 = square[1] * self.y_size + self.y_size
                all_verts.append(list(zip([x1, x2, x3, x4], [y1, y2, y3, y4])))
        all_verts = numpy.array(all_verts)
        pcoll_np = PolyCollection(all_verts, facecolors=color,
                                  edgecolors=color, linewidth=1)
        pcoll_np.set_alpha(alpha)
        self.axes.add_collection(pcoll_np)
        if self.parent.canvas is not None:
            self.axes.draw_artist(pcoll_np)

    def draw_corrupted(self):
        """Draws a sign on the meshgrid when a curve is corrupted.

        (Grey square and a red cross)
        """
        if (self.data.discarded_curves is not None
                and shared.exp.meshgrid_display_discarded):
            for i in range(self.nbr_pixels_x):
                for j in range(self.nbr_pixels_y):
                    if self.data.discarded_curves[i][j]:
                        x1 = i * self.x_size + self.x_size
                        y1 = j * self.y_size + self.y_size
                        x2 = i * self.x_size + self.x_size
                        y2 = j * self.y_size
                        x3 = i * self.x_size
                        y3 = j * self.y_size
                        x4 = i * self.x_size
                        y4 = j * self.y_size + self.y_size

                        verts = numpy.array([(x1, y1), (x2, y2),
                                             (x3, y3), (x4, y4)])
                        coll = PolyCollection([verts],
                                              facecolors="w", edgecolors="w")
                        self.axes.add_collection(coll)

                        if self.parent.canvas is not None:
                            self.axes.draw_artist(coll)

                        verts = numpy.array([(x1, y1), (x3, y3)])
                        coll = LineCollection(
                            [verts],
                            linewidths=1,
                            colors="r",
                            linestyle="solid")
                        self.axes.add_collection(coll)

                        if self.parent.canvas is not None:
                            self.axes.draw_artist(coll)

                        verts = numpy.array([(x2, y2), (x4, y4)])
                        coll = LineCollection(
                            [verts],
                            linewidths=1,
                            colors="r",
                            linestyle="solid")
                        self.axes.add_collection(coll)

                        if self.parent.canvas is not None:
                            self.axes.draw_artist(coll)

    def draw_missing_z(self):
        """Draws a sign on the meshgrid when a z value in a piezo image is missing.

        (Blue square)
        """
        if (shared.exp.missing_z_positions is not None \
                and shared.exp.meshgrid_display_missing_z):
            for x, y in shared.exp.missing_z_positions:
                x1 = x * self.x_size + self.x_size
                y1 = y * self.y_size + self.y_size
                x2 = x * self.x_size + self.x_size
                y2 = y * self.y_size
                x3 = x * self.x_size
                y3 = y * self.y_size
                x4 = x * self.x_size
                y4 = y * self.y_size + self.y_size

                verts = numpy.array([(x1, y1), (x2, y2),
                                     (x3, y3), (x4, y4)])
                coll = PolyCollection([verts],
                                      facecolors="b", edgecolors="b")
                self.axes.add_collection(coll)

                if self.parent.canvas is not None:
                    self.axes.draw_artist(coll)

    def update_blit(self, mode, roi=None, roi_mask=None, profile=None,
                    profile_id=None, closing=False):
        """Update the blitting for the meshgrids.

        http://matplotlib.org/examples/animation/animation_blit_qt4.html
        """
        if mode == "all":
            self.parent.restore_region(self.empty_bbox)

        elif mode == "preview_add" or mode == "preview_erase" or \
                mode == "profile_preview" or \
                mode == "profile_preview_erase" or \
                mode == "red_square" or mode == "discarded":
            if self.current_bbox is not None:
                # On Fedora it sometimes happens that self.current_bbox is
                # None (I don't know why). In this case don't restore the
                # region. This does not break anything so it's okay to skip
                # this if self.current_bbox is None.
                self.parent.restore_region(self.current_bbox)

        elif mode == "local_add" or mode == "local_erase":
            if self.first_blit:
                self.parent.restore_region(self.current_bbox)
                self.first_blit = False
            else:
                bbox = self.fig.canvas.copy_from_bbox(self.axes.bbox)
                self.parent.restore_region(bbox)

        if (self.plot_type in ("meshgrid", "meshgrid_data",
                               "meshgrid_tilt_check")):
            if self.plot_type != "meshgrid_tilt_check":
                # Draw ROI first
                if roi is not None or mode == "all":
                    self.draw_roi(mode, roi, roi_mask)

                # Draw profiles
                if profile is not None or mode == "all":
                    self.draw_profiles(mode, profile, profile_id, closing)

            if mode == "all" and self.parent is not None:
                self.current_bbox = \
                    self.parent.copy_from_bbox(self.axes.bbox)

            if self.plot_type != "meshgrid_tilt_check":
                # Draw corrupted curves if needed
                self.draw_corrupted()

        if self.mesh_type in ("piezo", "topo") or self.plot_type == "small_meshgrid_compute":
            self.draw_missing_z()

        # Always draw the red square on top
        self.draw_red_square()

        self.parent.blit(self.axes.bbox)
