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

"""Resampling of an image layer on the area of an AFM map."""

import numpy
from PIL import Image


def afm_footprint(afm):
    """Corner (x, y) and sizes of the area of the AFM map, in the
    coordinates of its actor [nm] (the surface is shifted by half a pixel,
    see Layer.set_axes_bounds)."""
    x0 = -afm.size_x / 2.0 - afm.res_x / 2.0
    y0 = -afm.size_y / 2.0 - afm.res_y / 2.0
    return x0, y0, afm.size_x, afm.size_y


def image_array(layer):
    """Pixels of an image layer as displayed (gray or RGB): row j of the
    array is the row j of the texture (bottom at y = -size_y / 2), with the
    flips."""
    array = layer.pixels
    if layer.flip_x:
        array = array[:, ::-1]
    if layer.flip_y:
        array = array[::-1, :]
    return array


def matrix(layer):
    """Matrix (4 x 4) of the actor of a layer: actor -> world."""
    mat = layer.actor.GetMatrix()
    return numpy.array([[mat.GetElement(i, j) for j in range(4)]
                        for i in range(4)])


def resample_images(image, afm, path):
    """Cuts the image layer on the area of the AFM layer (with the
    positions and angles of both layers), with the pixel size of the image,
    and saves it as a tiff in path.

    Returns the position (x, y) of the center of the new image and its
    angle, to put it on the AFM map when it is loaded."""
    array = image_array(image)
    ny_img, nx_img = array.shape[:2]

    # Pixels of the new image, about the pixel size of the image
    x0, y0, size_x, size_y = afm_footprint(afm)
    nx = max(1, int(round(size_x / image.res_x)))
    ny = max(1, int(round(size_y / image.res_y)))
    res_x, res_y = size_x / nx, size_y / ny

    # Centers of the pixels: AFM actor -> world -> image actor
    xs = x0 + (numpy.arange(nx) + 0.5) * res_x
    ys = y0 + (numpy.arange(ny) + 0.5) * res_y
    gx, gy = numpy.meshgrid(xs, ys)
    points = numpy.stack([gx.ravel(), gy.ravel(), numpy.zeros(gx.size),
                          numpy.ones(gx.size)])
    to_image = numpy.linalg.inv(matrix(image)) @ matrix(afm)
    local = to_image @ points

    # Pixels of the image (nearest), 0 outside of the image
    col = numpy.floor((local[0] / image.size_x + 0.5) * nx_img).astype(int)
    row = numpy.floor((local[1] / image.size_y + 0.5) * ny_img).astype(int)
    inside = (col >= 0) & (col < nx_img) & (row >= 0) & (row < ny_img)
    values = numpy.zeros((gx.size,) + array.shape[2:], dtype=numpy.uint8)
    values[inside] = array[row[inside], col[inside]]

    write_to_tiff(values.reshape((ny, nx) + array.shape[2:]),
                  [res_x / 1000.0, res_y / 1000.0], path)

    # Center of the area, in the world
    center = matrix(afm) @ [x0 + size_x / 2.0, y0 + size_y / 2.0, 0.0, 1.0]
    return center[0], center[1], afm.angle


def write_to_tiff(array, spacing, path):
    """Write the tiff to disk, with the pixel size [um] in the description
    (as in the tiffs exported by ImageJ, read back by pyAF). The row j of
    the array is the row j of the file (read back in the same order); gray
    or RGB."""
    description = ("Voxel_size_X: " + repr(float(spacing[0])) + " \xb5m\n"
                   "Voxel_size_Y: " + repr(float(spacing[1])) + " \xb5m\n"
                   "Voxel_size_Z: 0.0 \xb5m\n")
    Image.fromarray(numpy.ascontiguousarray(array, dtype=numpy.uint8)).save(
        path, tiffinfo={270: description.encode("latin-1")})
