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

"""Plane with image"""

import matplotlib
import numpy
import vtk
from PIL import Image
from vtkmodules.util import numpy_support

# LUTs made from matplotlib colormaps (the others go from the min to the
# max color)
COLORMAPS = {"Fire": "afmhot", "Viridis": "viridis"}


def frame_values(im):
    """Values of a frame of an image: 8 bits gray (rows, cols) or RGB
    (rows, cols, 3), or float for more bits (gray). None if not
    supported."""
    if im.mode in ("1", "L", "LA"):
        return numpy.array(im.convert("L"))
    if im.mode in ("P", "PA", "RGB", "RGBA", "CMYK", "YCbCr", "LAB", "HSV"):
        return numpy.array(im.convert("RGB"))
    if im.mode.startswith(("I", "F")):
        return numpy.array(im, dtype=float)
    return None


def to_8bits(values):
    """Values on 8 bits: more bits from the min to the max of all the
    frames."""
    if values.dtype == numpy.uint8:
        return values
    low, high = numpy.nanmin(values), numpy.nanmax(values)
    scale = 255.0 / (high - low) if high > low else 0.0
    return numpy.nan_to_num((values - low) * scale).astype(numpy.uint8)


def read_frames(path):
    """Pixels of all the frames of an image (tiff, 8 or 16 bits, gray or
    color), the row 0 being the first row of the file (as read by the
    vtkTIFFReader).

    Returns the pixels to display (8 bits, gray or RGB, (frames, rows,
    cols[, 3])) and their gray levels for the LUTs (8 bits; brightest
    channel of the color images, e.g. a green fluorescence keeps its
    intensity). None if the file is not an image."""
    try:
        im = Image.open(path)
        frames = []
        for index in range(getattr(im, "n_frames", 1)):
            im.seek(index)
            values = frame_values(im)
            if values is None:
                return None
            frames.append(values)
        values = numpy.stack(frames)
    except (IOError, ValueError):
        return None

    values = to_8bits(values)
    color = values.ndim == 4
    if color and (values[..., 0] == values[..., 1]).all() and \
            (values[..., 0] == values[..., 2]).all():
        # Gray levels (e.g. palette of grays)
        values, color = values[..., 0], False
    gray = values.max(axis=-1) if color else values
    return numpy.ascontiguousarray(values), numpy.ascontiguousarray(gray)


def read_image(path):
    """Pixels and gray levels (see read_frames) of a 2D image, (rows,
    cols[, 3]). None if the file is not a 2D image (or a stack)."""
    try:
        if getattr(Image.open(path), "n_frames", 1) > 1:
            return None
    except (IOError, ValueError):
        return None
    images = read_frames(path)
    if images is None or images[0].shape[0] > 1:
        return None
    return images[0][0], images[1][0]


def image_source(pixels):
    """VTK source (output port) of the pixels (8 bits, gray or RGB)."""
    rows, cols = pixels.shape[:2]
    data = vtk.vtkImageData()
    data.SetDimensions(cols, rows, 1)
    scalars = numpy_support.numpy_to_vtk(
        pixels.reshape(rows * cols, -1), deep=True,
        array_type=vtk.VTK_UNSIGNED_CHAR)
    data.GetPointData().SetScalars(scalars)
    source = vtk.vtkTrivialProducer()
    source.SetOutput(data)
    return source


def get_actor_plane(image, sizes, pos_z, lut_options, flip_options):
    """Plane with image."""
    # Create a plane
    plane = vtk.vtkPlaneSource()
    plane.SetNormal(0.0, 0.0, 1.0)

    size_x = sizes[0]
    size_y = sizes[1]

    plane.SetOrigin(size_x/2, size_y/2, pos_z)
    plane.SetPoint1(size_x/2, -size_y/2, pos_z)
    plane.SetPoint2(-size_x/2, size_y/2, pos_z)

    # Flip image if needed
    if flip_options["flip_x"]:
        flipX = vtk.vtkImageFlip()
        flipX.SetInputConnection(image.GetOutputPort())
        flipX.SetFilteredAxis(0)
        image = flipX
    if flip_options["flip_y"]:
        flipY = vtk.vtkImageFlip()
        flipY.SetInputConnection(image.GetOutputPort())
        flipY.SetFilteredAxis(1)
        image = flipY

    if lut_options["use"]:
        # Create a lookup table to map cell data to colors
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(256)
        lut.SetTableRange(0, 255)

        # Linear ramp from the min to the max color, or colormap (reversed
        # if inverted). The background is transparent: the gray
        # levels below the threshold (above 255 - threshold if inverted,
        # for a bright background)
        low = lut_options["min_color"]
        high = lut_options["max_color"]
        threshold = int(lut_options["threshold"])
        colormap = COLORMAPS.get(lut_options.get("name"))
        if colormap is not None:
            colormap = matplotlib.colormaps[colormap]
        for i in range(256):
            t = i / 255.0
            if lut_options["invert"]:
                t = 1.0 - t
                transparent = i > 255 - threshold
            else:
                transparent = i < threshold
            if colormap is not None:
                rgb = colormap(t)[:3]
            else:
                rgb = [low[k] + (high[k] - low[k]) * t for k in range(3)]
            lut.SetTableValue(i, *rgb, 0.0 if transparent else 1.0)

        # Lut build() has to be called
        lut.Build()

        # Map the image through the lookup table
        color = vtk.vtkImageMapToColors()
        color.SetLookupTable(lut)
        color.SetInputConnection(image.GetOutputPort())

        # Apply the texture
        texture = vtk.vtkTexture()
        texture.SetInputConnection(color.GetOutputPort())

    else:
        # Apply the texture
        texture = vtk.vtkTexture()
        texture.SetInputConnection(image.GetOutputPort())

    texture_plane = vtk.vtkTextureMapToPlane()
    texture_plane.SetInputConnection(plane.GetOutputPort())

    # Get the mapper for the actor
    plane_mapper = vtk.vtkPolyDataMapper()
    plane_mapper.SetInputConnection(texture_plane.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(plane_mapper)
    actor.SetTexture(texture)

    return actor
