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

"""Actions called by the menu. Used by the VTK widget."""

from ... import widgets_list
from ... import shared


def remove_layer():
    """Can be called from the menu to delete a layer from VTK."""
    if widgets_list.widget_vtk is None:
        # No 3D window
        return

    # Get the current layer (the last layer can not be removed)
    layer_id = widgets_list.widget_vtk.get_current_layer()
    if layer_id is None or len(shared.layer_list) < 2:
        return

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
