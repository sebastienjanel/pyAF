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

"""Menu displayed when right-clicking on the meshgrids."""

import os
from PyQt5 import QtWidgets

from .. import shared
from .. import widgets_list
from ..plots.PYAFPlot import PYAFPlot
from ..utils import discarding
from ..widgets.plot_options import PlotOptionsWidget
from ..widgets.progressbar import Progressbar
from ..widgets_small.indentation import WidgetIndentation
from ..widgets_small.slice_depth import ChangeSliceDepthWidget


def create_meshgrid_menu(parent, widget, pos):
    """Create the menu."""
    data = shared.exp.current_data

    # Create menu
    menu = QtWidgets.QMenu()

    # Units for the distances
    if shared.exp.current_data.meshgrid_units == "nm":
        txt = "um"
    elif shared.exp.current_data.meshgrid_units == "um":
        txt = "nm"

    action = QtWidgets.QAction("Set distance units to " + txt, parent)
    action.triggered.connect(change_mesh_units)
    menu.addAction(action)

    menu.addSeparator()

    # Show options widget
    if widget != "data":
        action = QtWidgets.QAction("Color Options", parent)
        action.triggered.connect(lambda: open_color_options(parent))
        menu.addAction(action)

        menu.addSeparator()

    discarding.add_discard_actions_to_menu(parent, menu)

    # Display/hide red square
    if shared.exp.meshgrid_display_red_square:
        action = QtWidgets.QAction("Hide red square", parent)
    else:
        action = QtWidgets.QAction("Display red square", parent)
    action.triggered.connect(display_square)
    menu.addAction(action)
    # Display/hide missing z
    if shared.exp.meshgrid_display_missing_z:
        action = QtWidgets.QAction("Hide missing z values", parent)
    else:
        action = QtWidgets.QAction("Display missing z values", parent)
    action.triggered.connect(display_missing_z)
    menu.addAction(action)
    # Display/hide discarded
    if shared.exp.meshgrid_display_discarded:
        action = QtWidgets.QAction("Hide discarded", parent)
    else:
        action = QtWidgets.QAction("Display discarded", parent)
    action.triggered.connect(discarding.display_discarded)
    menu.addAction(action)

    # Change indentation depth or heigth
    if widget == "results":
        if data.meshgrid_type == "stiffness" or \
                data.meshgrid_type == "stiffness_corr":
            action = QtWidgets.QAction("Indentation depth", parent)
            action.triggered.connect(lambda: open_indentation_widget(parent))
            menu.addAction(action)
        elif data.meshgrid_type == "stiffness_slice" or \
                data.meshgrid_type == "stiffness_corr_slice":
            action = QtWidgets.QAction("Slice height", parent)
            action.triggered.connect(lambda: change_slice_heigth(parent))
            menu.addAction(action)

        menu.addSeparator()

        if data.roi_selected_row is not None:
            # Add ROI
            text = "Draw to ROI " + str(data.roi_selected_row + 1)
            action = QtWidgets.QAction(text, parent)
            action.triggered.connect(
                lambda: parent.change_cursor_mode("roi_add"))
            menu.addAction(action)

            # Erase from ROI
            text = "Erase from ROI " + str(data.roi_selected_row + 1)
            action = QtWidgets.QAction(text, parent)
            action.triggered.connect(
                lambda: parent.change_cursor_mode("roi_erase"))
            menu.addAction(action)

        # Add a profile
        action_plot_profile = QtWidgets.QAction("Plot profile", parent)
        action_plot_profile.triggered.connect(
            lambda: parent.change_cursor_mode("profile_add"))
        if data.stiffness_array is None:
            action_plot_profile.setEnabled(False)
        menu.addAction(action_plot_profile)

        # Remove a profile
        action_remove_profile = QtWidgets.QAction("Remove profile", parent)
        action_remove_profile.triggered.connect(
            lambda: parent.change_cursor_mode("profile_erase"))
        if len(data.profile_list) == 0:
            action_remove_profile.setEnabled(False)
        menu.addAction(action_remove_profile)
        menu.addSeparator()

        # Add action to save all figures
        if data.meshgrid_type == "stiffness" or \
                data.meshgrid_type == "stiffness_corr":
            menu.addSeparator()
            action = QtWidgets.QAction("Save slices as", parent)
            action.triggered.connect(lambda: save_slices_as(parent))
            menu.addAction(action)

        # Add action to save figure without border
        if data.stiffness_calculated:
            action = QtWidgets.QAction("Save without border", parent)
            action.triggered.connect(lambda: save_without_border(parent))
            menu.addAction(action)

        # Change slope for events per curve meshgrid
        if data.meshgrid_type == "events_per_curve" or \
                data.meshgrid_type == "events_rupture_force":
            menu.addSeparator()
            action = QtWidgets.QAction("Change slope filter", parent)
            action.triggered.connect(lambda: parent.change_filter("slope"))
            menu.addAction(action)
            action = QtWidgets.QAction("Change dist filter", parent)
            action.triggered.connect(lambda: parent.change_filter("dist"))
            menu.addAction(action)

    # Display menu
    menu.popup(pos, menu.menuAction())
    menu.exec_()


def open_color_options(parent):
    """Open the color options widget.

    Open the widget that is used to change the color scales of the
    meshgrids.
    """
    if widgets_list.widget_meshgrid_options is None:
        # Create new widget
        PlotOptionsWidget(parent)
        widgets_list.widget_meshgrid_options.resize(250, 250)
        widgets_list.widget_meshgrid_options.show()
    else:
        # Bring to front
        widgets_list.widget_meshgrid_options.update_widget()
        widgets_list.widget_meshgrid_options.activateWindow()
        widgets_list.widget_meshgrid_options.raise_()


def change_mesh_units():
    """Change the units of the distances in the meshgrid plot.

    The units can be either nm or um.
    """
    if shared.exp.current_data.meshgrid_units == "nm":
        shared.exp.current_data.meshgrid_units = "um"
    elif shared.exp.current_data.meshgrid_units == "um":
        shared.exp.current_data.meshgrid_units = "nm"

    refresh_meshgrids()


def open_indentation_widget(parent):
    """Open the widget used to change the indentation depth."""
    if widgets_list.widget_indentation is None:
        WidgetIndentation(parent)
        widgets_list.widget_indentation.setWindowTitle("Indendation depth")
        widgets_list.widget_indentation.resize(300, 100)
        widgets_list.widget_indentation.activateWindow()
        widgets_list.widget_indentation.show()
    else:
        widgets_list.widget_indentation.activateWindow()
        widgets_list.widget_indentation.raise_()


def change_slice_heigth(parent):
    """Opens a small widget to change the height of the slice."""
    widget = ChangeSliceDepthWidget(parent)
    widget.setWindowTitle("Change height (nm)")
    widget.resize(200, 100)
    widget.activateWindow()
    widget.show()

    refresh_meshgrids()


def display_square():
    """Display or hide the red selection square."""
    if shared.exp.meshgrid_display_red_square:
        shared.exp.meshgrid_display_red_square = False
    else:
        shared.exp.meshgrid_display_red_square = True
    refresh_meshgrids()


def display_missing_z():
    """Display or hide the blue squares corresponding to missing z values in the piezo image."""
    if shared.exp.meshgrid_display_missing_z:
        shared.exp.meshgrid_display_missing_z = False
    else:
        shared.exp.meshgrid_display_missing_z = True
    refresh_meshgrids()


def refresh_meshgrids():
    """Refresh the meshgrids."""
    widgets_list.widget_data.meshgrid.update_plot()
    widgets_list.widget_results.MPL_canvas1.update_plot()


def save_without_border(parent):
    """Save the meshgrid without borders."""
    ffilename, _ = QtWidgets.QFileDialog.getSaveFileName(parent,
        "Choose a file name", ".", "TIFF (*.tiff)")

    if not filename:
        return
    else:
        orig_save_path = str(filename).split(".")[0]
        ext = str(filename).split(".")[1]

        save_path = orig_save_path + "." + ext
        PYAFPlot(parent, "meshgrid", save_path=save_path, noborder=True)


def save_slices_as(parent):
    """Save all the slices."""
    data = shared.exp.current_data
    filename, _ = QtWidgets.QFileDialog.getSaveFileName(parent,
        "Choose a file name", ".", "TIFF (*.tiff)")
    if not filename:
        return
    else:
        orig_save_path = str(filename).split(".")[0]
        ext = str(filename).split(".")[1]

    # Save the oldposition to reset it afterwards
    old_pos = data.stiffness_depth_view

    # Create a progressbar
    themax = shared.exp.global_max_indentation_index
    Progressbar()
    widgets_list.widget_progressbar.set_label("Saving ...")
    widgets_list.widget_progressbar.set_range(0, themax)

    # Loop through the slices
    for i in range(shared.exp.global_max_indentation_index):
        data.stiffness_depth_view = i
        save_path = orig_save_path + "_" + str(i) + "." + ext
        text = "Saving " + os.path.basename(save_path) + " ..."
        widgets_list.widget_progressbar.set_label(text)
        PYAFPlot(parent, "meshgrid", save_path=save_path)

        widgets_list.widget_progressbar.update()

    # Reset the position
    data.stiffness_depth_view = old_pos

    # Close the progressbar
    widgets_list.widget_progressbar.close()
