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

"""Import saved layers in VTK."""

import pickle
from qtpy import QtCore, QtWidgets

from ... import shared
from ... import widgets_list
from ...tools import misc_tools
from ...widgets_vtk.utils.create_layer import create_layer


def load_layers(parent, path=None):
    """Import saved layers in VTK."""
    # Get the saved save path from the settings
    settings = QtCore.QSettings()
    settings_path = misc_tools.get_user_path(settings)

    # Open dialog to save the file
    if path is None:
        path = QtWidgets.QFileDialog.getExistingDirectory(
            parent, "Load layers", settings_path)

    # Save the file
    if path:
        # Save current path
        misc_tools.set_user_path(settings, str(path))

        # Get number of layers to load
        global_file = open(str(path) + "/global.txt", "r")
        nbr_layers = int(global_file.read().split("=")[1])
        global_file.close()

        # Gloabl parameters
        global_params_file = open(path + "/global_params.txt", "rb")
        global_params = pickle.loads(global_params_file.read())
        global_params_file.close()
        camera = widgets_list.widget_vtk.canvas.renderer.GetActiveCamera()
        camera.SetPosition(global_params["camera_pos"])
        widgets_list.widget_vtk.canvas.renderer.ResetCameraClippingRange()
        widgets_list.widget_vtk.scales_color = global_params["scales_color"]
        # Label of the scales button
        widgets_list.widget_vtk.menu_layer_list.update_color_buttons()

        # The saved AFM layers are applied to the AFM layers of the window,
        # in the same order; the other layers are created
        afm_layers = [layer for layer in shared.layer_list
                      if layer.type == "afm"]
        nbr_afm = 0
        # Index of each saved layer in the window
        new_index = {}

        # Read each layer file
        for i in range(nbr_layers):
            with open(str(path) + "/layer_" + str(i) + ".txt", "rb") as fi:
                layer_prefs = pickle.loads(fi.read())

            layer_type = layer_prefs["type"]
            nbr_before = len(shared.layer_list)

            if layer_type == "afm":
                if nbr_afm >= len(afm_layers):
                    # No AFM layer left in the window for this one
                    continue
                layer = afm_layers[nbr_afm]
                nbr_afm += 1

                for key in layer_prefs:
                    # Do not update the pixel sizes for AFM
                    if key not in ("nbr_pixels_x", "nbr_pixels_y",
                                   "nbr_pixels_z"):
                        misc_tools.setattr_special(
                            layer, key, layer_prefs[key])

                layer.draw_at_position()
                new_index[i] = shared.layer_list.index(layer)
                continue

            if layer_type in ("single_image", "stack"):
                new_path = str(path) + "/layer_" + str(i) + ".tiff"
                create_layer(layer_type, new_path)

            elif layer_type == "isosurface":
                stack_id = layer_prefs["isosurf_orig_layer_id"]
                if stack_id not in new_index:
                    # The stack was not loaded
                    continue
                new_path = str(path) + "/layer_" + str(stack_id) + ".tiff"
                create_layer(
                    layer_type,
                    isosurf_stack_id=new_index[stack_id],
                    iso_file=new_path)

            elif layer_type == "glass" and "afm_id" in layer_prefs:
                widgets_list.widget_vtk.canvas.addGlassLayer(
                    len(shared.layer_list), layer_prefs["afm_id"],
                    layer_prefs["roi_id"])
                newid = len(shared.layer_list) - 1
                shared.layer_list[-1].layer_id = newid
                menu = widgets_list.widget_vtk.menu_layer_list
                menu.tableWidget.insertRow(newid)
                menu.add_line_in_tablewidget(newid, "glass")

            if len(shared.layer_list) == nbr_before:
                # The layer could not be created
                continue

            # Update the values in the new layer
            layer = shared.layer_list[-1]
            new_index[i] = len(shared.layer_list) - 1
            for key in layer_prefs:
                if key not in ("afm_id", "roi_id", "isosurf_orig_layer_id"):
                    misc_tools.setattr_special(layer, key, layer_prefs[key])

            if layer_type != "glass":
                # Use the new path of the file as image path
                layer.filepath = new_path

            layer.update_actor()

        # Update the GUI
        widgets_list.widget_vtk.menu_layer_list.reset_table()
        widgets_list.widget_vtk.menu_layer_list.update_widget()
        widgets_list.widget_vtk.menu_options.update_widget()
