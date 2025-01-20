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

"""Module used to save the 3D layers."""

import pickle
import shutil
from PyQt5 import QtCore, QtWidgets

from ... import consts
from ... import shared
from ... import widgets_list
from ...tools import misc_tools
from ...widgets.progressbar import Progressbar

if consts.ALLOW_ITK:
    import itk


def save_layers(parent):
    """Module used to save the 3D layers.

    All the images and the positions are saved in a folder.
    """
    # Get the saved save path from the settings
    settings = QtCore.QSettings()
    settings_path = misc_tools.get_user_path(settings)

    # Open dialog to save the file
    folder = QtWidgets.QFileDialog.getExistingDirectory(
        parent, "Save layers", settings_path)

    # Save the file
    if folder:
        # Create a progressbar which is displayed during the saving
        Progressbar("Preparing afm data")
        widgets_list.widget_progressbar.set_range(0, len(shared.layer_list))

        # Save current path
        misc_tools.set_user_path(settings, str(folder))

        # Save a global file
        nbr_layers = len(shared.layer_list)
        global_file = open(folder + "/global.txt", "wb")
        global_file.write("nbr_layers=" + str(nbr_layers))
        global_file.close()

        # Save global parameters
        to_save = {}

        camera = widgets_list.widget_vtk.canvas.renderer.GetActiveCamera()
        campos = camera.GetPosition()
        to_save["camera_pos"] = campos
        to_save["scales_color"] = widgets_list.widget_vtk.scales_color

        global_params = open(folder + "/global_params.txt", "wb")
        pickle.dump(to_save, global_params, pickle.HIGHEST_PROTOCOL)
        global_params.close()

        # Save a file for each layer
        for i in range(len(shared.layer_list)):
            layer = shared.layer_list[i]

            text = "Saving " + str(layer.filename) + "..."
            widgets_list.widget_progressbar.set_label(text)

            if layer.type == "single_image":
                shutil.copy(
                    layer.filepath,
                    folder + "/layer_" + str(i) + ".tiff")

            if layer.type == "stack":
                image_type = itk.Image[itk.UC, 3]
                writer = itk.ImageFileWriter[image_type].New()
                fname = str(folder + "/layer_" + str(i) + ".tiff")
                writer.SetFileName(fname)
                writer.SetInput(layer.input_image)
                writer.Update()

            to_save = {}

            # to_save["layer_id"] = layer.layer_id
            to_save["type"] = layer.type
            to_save["filename"] = layer.filename
            to_save["pos_x"] = layer.pos_x
            to_save["pos_y"] = layer.pos_y
            to_save["pos_z"] = layer.pos_z
            to_save["angle"] = layer.angle
            to_save["visible"] = layer.visible
            to_save["axes_visible"] = layer.axes_visible
            to_save["opacity"] = layer.opacity
            # to_save["fiducials"] = layer.fiducials # Contains vtk objet ...
            to_save["use_lut"] = layer.use_lut
            to_save["lut_threshold"] = layer.lut_threshold
            to_save["lut_invert"] = layer.lut_invert
            to_save["flip_x"] = layer.flip_x
            to_save["flip_y"] = layer.flip_y
            to_save["display_z_scale"] = layer.display_z_scale

            to_save["lower_thresh"] = layer.lower_thresh
            to_save["upper_thresh"] = layer.upper_thresh
            to_save["lower_thresh_value"] = layer.lower_thresh_value
            to_save["upper_thresh_value"] = layer.upper_thresh_value
            to_save["use_binary_thresholding"] = layer.use_binary_thresholding

            to_save["lut_min_color"] = layer.lut_min_color
            to_save["lut_max_color"] = layer.lut_max_color
            to_save["bind"] = layer.bind

            # Smoothing of images
            to_save["smooth_gaussian"] = layer.smooth_gaussian

            # Color of the isosurface
            to_save["isosurface_color"] = layer.isosurface_color

            to_save["nbr_pixels_x"] = layer.nbr_pixels_x
            to_save["nbr_pixels_y"] = layer.nbr_pixels_y
            to_save["nbr_pixels_z"] = layer.nbr_pixels_z
            to_save["size_x"] = layer.size_x
            to_save["size_y"] = layer.size_y
            to_save["size_z"] = layer.size_z
            to_save["res_x"] = layer.res_x
            to_save["res_y"] = layer.res_y
            to_save["res_z"] = layer.res_z

            if layer.type == "single_image":
                to_save["filepath"] = layer.filepath

            if layer.type == "stack":
                to_save["slide_pos"] = layer.slide_pos
                to_save["single_opacity"] = layer.single_opacity
                to_save["others_opacity"] = layer.others_opacity

            if layer.type == "isosurface":
                to_save["isosurf_orig_layer_id"] = layer.isosurf_orig_layer_id

            # Save layer's parameters
            path = folder + "/layer_" + str(i) + ".txt"
            pickle.dump(to_save, open(path, "wb"), pickle.HIGHEST_PROTOCOL)

            widgets_list.widget_progressbar.update()

        # Close the progressbar
        widgets_list.widget_progressbar.close()
