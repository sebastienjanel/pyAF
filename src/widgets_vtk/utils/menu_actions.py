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

"""Actions called by the menu. Used by the VTK widget."""

from ... import widgets_list
from ... import shared


def remove_layer():
    """Can be called from the menu to delete a layer from VTK."""
    # Get the current layer
    layer_id = widgets_list.widget_vtk.get_current_layer()

    # Remove from 3D
    shared.layer_list[layer_id].delete()

    # Delete
    del shared.layer_list[layer_id]

    # Highlight the picked actor with a red box (layer 0 will be selected)
    layer = shared.layer_list[0]
    widgets_list.widget_vtk.canvas.style.HighlightProp(layer.actor)

    # Reset the list in the menu
    widgets_list.widget_vtk.menu_layer_list.reset_table()


def ask_user_for_afm_layers():
    """Opens a widget asking the user which AFM tomograpy he wants to add.

    The user can chose tomographies coming from the result's list. It is
    even possible to load multiple times the same tomography, if you want
    to show the difference between two different color scales for example.
    """
    widgets_list.widget_results.open_vtk_entry_widget("add_more")


def copy_pos():
    """Called by the menu or by Ctrl-C. Copies the layer's position.

    The position is copied to an internal clipboard in shared.
    """
    # Get the current layer
    layer_id = widgets_list.widget_vtk.get_current_layer()
    layer = shared.layer_list[layer_id]

    # Store the values in a dict
    vals = {"id": layer_id,
            "type": layer.type,
            "filename": layer.filename,
            "x": layer.pos_x,
            "y": layer.pos_y,
            "z": layer.pos_z,
            "angle": layer.angle,
            "size_x": layer.size_x,
            "size_y": layer.size_y,
            "size_z": layer.size_z,
            "res_x": layer.res_x,
            "res_y": layer.res_y,
            "res_z": layer.res_z}

    # Save the dict in the clipboard
    shared.clipboard_layer_positions = vals


def paste_pos():
    """Called by the menu or by Ctrl-V. Pastes the layer's position.

    The position is taken from the internal clipboard in shared.
    """
    # Get the values from the clipboard
    vals = shared.clipboard_layer_positions

    # Get the current layer
    layer_id = widgets_list.widget_vtk.get_current_layer()
    layer = shared.layer_list[layer_id]

    old_angle = layer.angle

    # Paste the values
    layer.pos_x = vals["x"]
    layer.pos_y = vals["y"]
    layer.pos_z = vals["z"]
    layer.angle = vals["angle"]

    # Update the sizes only for tiff layers
    if vals["type"] == "optics":
        layer.size_x = vals["size_x"]
        layer.size_y = vals["size_y"]
        layer.size_z = vals["size_z"]
        layer.res_x = vals["res_x"]
        layer.res_y = vals["res_y"]
        layer.res_z = vals["res_z"]

    if vals["type"] == "afm":
        # For AFM layers, just change the position
        layer.translate(0, 0, 0)
        layer.rotate(None, None, increment=(old_angle - layer.angle))

    if vals["type"] == "optics":
        # For tiff layers, reinitialize the size
        layer.update_actor()
        widgets_list.widget_vtk.menu.tab_Optics.update_widget()
