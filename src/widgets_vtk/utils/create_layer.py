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

"""Common tools."""

import os
import sys
from PyQt5 import QtCore, QtWidgets

from ... import consts
from ... import shared
from ... import widgets_list
from ...tools import misc_tools


def create_layer(layer_type, path=None, isosurf_stack_id=None, iso_file=None):
    """Creates a new empty layer.

    Defines also the sizes of the layers.
    """
    if layer_type in ("single_image", "stack"):
        filenames = ask_user_for_files(path)
    elif layer_type == "isosurface":
        filenames = [iso_file]

    print(filenames)

    wg = widgets_list.widget_vtk

    for path in filenames:
        # Add the layer
        if layer_type == "single_image":
            ok = wg.canvas.addSingleImageLayer(path)
        elif layer_type == "HR":
            ok = wg.canvas.addOpticsHRLayer(path)
        elif layer_type == "stack":
            ok = wg.canvas.addStackLayer(path)
        elif layer_type == "isosurface":
            ok = wg.canvas.addIsoSurfaceLayer(isosurf_stack_id, path)

        if ok:
            # Add line to tablewidget
            newid = len(shared.layer_list) - 1
            wg.menu_layer_list.tableWidget.insertRow(newid)
            wg.menu_layer_list.add_line_in_tablewidget(newid, layer_type)

            # Select the added line and select image tab
            wg.menu_layer_list.tableWidget.setCurrentCell(newid, 0)
            wg.menu_layer_list.tableWidget_clicked(newid, 0)

            # Add id
            layer = shared.layer_list[newid]
            layer.layer_id = newid

            # Change values
            if layer_type == "single_image" or layer_type == "stack":
                tab = wg.menu_options.tab_Optics
                tab.IN_x_size.changeValue(layer.size_x / 1000.0)
                tab.IN_y_size.changeValue(layer.size_y / 1000.0)
                tab.IN_x_res.changeValue(layer.res_x / 1000.0)
                tab.IN_y_res.changeValue(layer.res_y / 1000.0)
        else:
            # The loading failed, tell the user
            name = os.path.basename(path)

            if layer_type == "single_image":
                txt = ("The file {name} did not load... Only 2D 8 "
                       "bit tiff files are accepted !".format(name=name))
            elif layer_type == "stack":
                txt = ("The file {name} is not an image stack..."
                       "".format(name=name))
            QtWidgets.QMessageBox.warning(wg, "Warning", txt)


def ask_user_for_files(path=None):
    """Ask the user to load files.

    Return the list of files selected, which may be empty.
    """
    # Load the settings
    settings = QtCore.QSettings()

    # Get saved path
    settings_path = misc_tools.get_user_path(settings)

    if path is not None:
        fnames = [path]
    elif consts.AUTO_LOAD_TIFF_IN_3D:
        # Will autoload the tiff file, is used for debugging purposes
        tiff_file = "/../examples/data/tem.tiff"
        fnames = [sys.path[0] + tiff_file]
    elif consts.AUTO_LOAD_STACK_IN_3D:
        # Will autoload the tiff stack, is used for debugging purposes
        tiff_stack = "/../examples/data/conf_stack.tiff"
        fnames = [sys.path[0] + tiff_stack]
    else:
        fnames, _ = QtWidgets.QFileDialog.getOpenFileNames(
            directory=settings_path)

    if fnames:
        dirpath = os.path.dirname(fnames[-1])
        # Save path for next use
        misc_tools.set_user_path(settings, dirpath)

    return fnames
