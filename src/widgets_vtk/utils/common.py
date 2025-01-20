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

import vtk


def center_data(data, size_x, size_y, off=None):
    """Centers the polydata to the 0, 0 position.

    An supplementary offset parameter can be passed as argument.
    """
    off_x = 0
    off_y = 0
    if off is not None:
        off_x = off[0]
        off_y = off[1]

    trans = vtk.vtkTransform()
    trans.Translate(-size_x / 2.0 + off_x, -size_y / 2.0 + off_y, 0.0)
    trans_filter = vtk.vtkTransformPolyDataFilter()
    trans_filter.SetTransform(trans)

    if vtk.VTK_MAJOR_VERSION <= 5:
        trans_filter.SetInput(data)
    else:
        trans_filter.SetInputData(data)
        trans_filter.Update()

    return trans_filter.GetOutput()


def fill_image_sizes(parent, size, cal):
    """Saves the sizes of the images in the different variables of the layer."""
    parent.nbr_pixels_x = size[0]
    parent.nbr_pixels_y = size[1]
    if len(size) > 2:
        parent.nbr_pixels_z = size[2]
    else:
        parent.nbr_pixels_z = 1.0

    if cal[0] is None:
        # No calibration was found, set a default value.
        parent.size_x = 50000.0
        parent.size_y = 50000.0
        parent.size_z = 1000.0
        parent.res_x = parent.size_x / parent.nbr_pixels_x
        parent.res_y = parent.size_y / parent.nbr_pixels_y
        parent.res_z = parent.size_z / parent.nbr_pixels_z

    else:
        # A calibration was found
        parent.res_x = cal[0] * 1e9
        parent.res_y = cal[1] * 1e9
        parent.res_z = cal[2] * 1e9
        parent.size_x = parent.res_x * parent.nbr_pixels_x
        parent.size_y = parent.res_y * parent.nbr_pixels_y
        parent.size_z = parent.res_z * parent.nbr_pixels_z

    # Store original values
    parent.original_nbr_pixels_x = parent.nbr_pixels_x
    parent.original_nbr_pixels_y = parent.nbr_pixels_y
    parent.original_nbr_pixels_z = parent.nbr_pixels_z
    parent.original_res_x = parent.res_x
    parent.original_res_y = parent.res_y
    parent.original_res_z = parent.res_z
    parent.original_size_x = parent.size_x
    parent.original_size_y = parent.size_y
    parent.original_size_z = parent.size_z
