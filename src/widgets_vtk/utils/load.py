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
from PyQt5 import QtCore, QtWidgets

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
        global_file = open(str(path) + "/global.txt", "rb")
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

        # Read each layer file
        for i in range(nbr_layers):
            layer_file = open(
                str(path) +
                "/layer_" +
                str(i) +
                ".txt",
                "rb")

            layer_prefs = pickle.loads(layer_file.read())

            layer_type = layer_prefs["type"]

            if layer_type == "single_image":
                new_path = str(path) + "/layer_" + str(i) + ".tiff"

                create_layer(layer_type, new_path)

                # Update the values in the layer
                layer = shared.layer_list[i]
                for key in layer_prefs:
                    misc_tools.setattr_special(layer, key, layer_prefs[key])

                # Use the new path of the file as image path
                layer.filepath = new_path

                layer.update_actor()

            elif layer_type == "afm":
                # For AFM layers, just update
                layer = shared.layer_list[i]

                for key in layer_prefs:
                    # Do not update the pixel sizes for AFM
                    if key != "nbr_pixels_x" and key != "nbr_pixels_y" and \
                            key != "nbr_pixels_z":
                        misc_tools.setattr_special(
                            layer, key, layer_prefs[key])

                layer.draw_at_position()

            elif layer_type == "stack":
                new_path = str(path) + "/layer_" + str(i) + ".tiff"

                create_layer(layer_type, new_path)

                # Update the values in the layer
                layer = shared.layer_list[i]
                for key in layer_prefs:
                    misc_tools.setattr_special(layer, key, layer_prefs[key])

                # Use the new path of the file as image path
                layer.filepath = new_path

                layer.update_actor()

            elif layer_type == "isosurface":
                stack_id = layer_prefs["isosurf_orig_layer_id"]
                new_path = str(path) + "/layer_" + str(i-1) + ".tiff"
                create_layer(
                    layer_type,
                    isosurf_stack_id=stack_id,
                    iso_file=new_path)

                # Update the values in the layer
                layer = shared.layer_list[i]
                for key in layer_prefs:
                    misc_tools.setattr_special(layer, key, layer_prefs[key])

                # Use the new path of the file as image path
                layer.filepath = new_path

                layer.update_actor()

            elif layer_type == "glass":
                print("No substrate loading")

            layer = shared.layer_list[i]

        # Update the GUI
        widgets_list.widget_vtk.menu_layer_list.reset_table()
        widgets_list.widget_vtk.menu_layer_list.update_widget()
        widgets_list.widget_vtk.menu_options.update_widget()
