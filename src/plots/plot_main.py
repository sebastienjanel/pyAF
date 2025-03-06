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

"""Main plotting module."""

from .. import shared
from ..tools import misc_tools
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.font_manager as fm
import platform
from PyQt5 import QtCore


class MainPlot:
    """Main plotting class.

    To be inherited.
    """

    def __init__(self, parent):
        self.parent = parent
        self.plot_type = self.parent.plot_type
        self.curve_type = self.parent.curve_type
        self.mesh_type = self.parent.mesh_type
        self.save_path = self.parent.save_path
        self.noborder = self.parent.noborder
        self.data = self.parent.data
        self.first_blit = True  # For the meshgrids

        self.empty_bbox = None
        self.current_bbox = None

        self.fig = self.parent.fig
        self.window_color_level = self.parent.window_color_level

        # Load the settings
        self.settings = QtCore.QSettings()

        self.xlabel = None
        self.ylabel = None

        if self.save_path is not None:
            # Default savefig directory
            matplotlib.rcParams["savefig.directory"] = self.save_path

        # Font (use the Lucida Grande font I ship)
        path = misc_tools.get_app_path() + "/fonts/HelveticaNowText-Regular.ttf"
        font_size = shared.exp.mpl_labels_font_size
        try:
            with open(path):
                self.font = fm.FontProperties(fname=path, size=font_size)
        except IOError:
            self.font = fm.FontProperties(size=font_size)

        matplotlib.rcParams["font.size"] = shared.exp.mpl_labels_font_size
        # Adjust distance between xticks and plot, prevents overlapping
        # of xlabels and ylabels at origin
        matplotlib.rcParams["xtick.major.pad"] = "12"
        # Smaller font for tick labels
        tick_labels_font_size = shared.exp.mpl_tick_labels_font_size
        matplotlib.rcParams["xtick.labelsize"] = tick_labels_font_size
        matplotlib.rcParams["ytick.labelsize"] = tick_labels_font_size

        if self.parent.single:
            # Create a new empty figure
            self.fig = pyplot.figure()
            self.fig.clf()
        else:
            # Manually set the background color
            rect = self.fig.patch

            if platform.uname()[0] == "Linux":
                if self.window_color_level == 0:
                    rect.set_facecolor("#f2f1f0")
                elif self.window_color_level == 1 or \
                        self.window_color_level == 2:
                    rect.set_facecolor("#edeceb")
                else:
                    rect.set_facecolor("#ffffff")  # white
            elif platform.uname()[0] == "Darwin":
                # OS X background color
                if self.window_color_level == 0:
                    rect.set_facecolor("#e5e5e5")
                elif self.window_color_level == 1:
                    rect.set_facecolor("#ededed")
                elif self.window_color_level == 2:
                    rect.set_facecolor("#e4e4e4")
                else:
                    rect.set_facecolor("#ffffff")  # white
            elif platform.uname()[0] == "Windows":
                rect.set_facecolor("white")  # Windows 7 (Not tested)

        self.axes = self.fig.add_subplot(111)

    def plot_labels(self):
        """Plots the x and y lables."""
        if self.xlabel is not None:
            self.axes.set_xlabel(self.xlabel, fontproperties=self.font)
        if self.ylabel is not None:
            self.axes.set_ylabel(self.ylabel, fontproperties=self.font)

    def start_plot(self):
        """Do the actual plotting.

        Can plot the data in an embedded canvas for the GUI, or in a
        separate matplotlib window.
        """
        p_type = self.plot_type

        if self.parent.single:
            pyplot.draw()

            if p_type == "meshgrid_tilt_check":
                self.draw_check_tilt_squares()

            if p_type in ("meshgrid", "meshgrid_data"):
                # Save the current and empty bboxes for blitting
                cv = self.parent
                self.empty_bbox = cv.copy_from_bbox(self.axes.bbox)
                self.current_bbox = cv.copy_from_bbox(self.axes.bbox)

                self.draw_roi("all")
                self.draw_profiles("all")
                self.draw_corrupted()
                self.draw_red_square()

                # If it's for a save (meshgrid)
                if self.save_path is not None:
                    if self.noborder:
                        # Remove borders
                        self.fig.set_size_inches(8, 8)
                        self.fig.set_dpi(100)
                        pyplot.axis("off")
                        pyplot.subplots_adjust(0, 0, 1, 1, None, None)
                        pyplot.savefig(self.save_path, bbox_inches="tight",
                                       pad_inches=0)
                    else:
                        pyplot.savefig(self.save_path)

                    # Important to use close, because event if this class is
                    # garbage collected at close, the ref count in matplotlib
                    # for the figure stays > 0, and the figure is not freed
                    # from memory, which gives you a nice memory leak.
                    pyplot.close()

            # Open the curve in a new window
            if self.save_path is None:
                pyplot.show()

        else:
            if p_type in ("results_single", "results_groups", "results_experiment", "results_stats"):
                # In the case of the results plot, due to the legend which
                # is placed outside of the box with a small hack, we can not
                # use the tight layout. In this case we define the borders
                # manually, by changing the bottom value to 0.2 (to leave
                # enough space for the x label), and the right value to 0.7
                # (to leave enough space for the legend).
                self.fig.subplots_adjust(left=0.125, right=0.7,
                                         bottom=0.2, top=0.9)

            elif p_type == "small_meshgrid" or p_type == "multimeshgrid":
                self.axes.set_axis_off()  # Remove black border
                self.fig.subplots_adjust(0, 0, 1, 1, None, None)  # Adjust

            else:
                # Use tight layout for all the other figures
                if p_type != "meshgrid" and \
                        p_type != "small_meshgrid_compute" and \
                        p_type != "meshgrid_data":
                    self.fig.tight_layout()

            # Draw in the canvas
            self.fig.canvas.draw()

            # Draw ROI for multimeshgrid
            display_roi = shared.exp.display_roi_in_multimeshgrid
            if p_type == "multimeshgrid" and display_roi:
                self.draw_roi("all")

            if p_type == "meshgrid_tilt_check":
                self.draw_check_tilt_squares()
