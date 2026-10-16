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

"""Uses tiffs to make a tiff stack."""

import math
import itk
import numpy
from PIL import Image
from vtkmodules.util import numpy_support
from ... import shared
from ... import widgets_list


def resample_images(layer_id, path):
    """Cuts the image of a layer along the drawn square (with the position,
    angle and pixel size of the layer) and saves it as a tiff in path.

    Used to make a tiff stack out of single tiffs of different sizes."""
    layer_list = shared.layer_list
    layer = layer_list[layer_id]
    wl = widgets_list.widget_vtk

    img = get_itk_image_from_vtk_tiff_2D(layer.tiff)
    img.SetSpacing([layer.res_x / 1000.0, layer.res_y / 1000.0])

    img_type = itk.Image[itk.UC, 2]

    # Flip image (for photoshop tiffs, should be made conditionnal ?)
    flip_filter = itk.FlipImageFilter[img_type].New()
    flip_filter.SetFlipAxes([0, 1])
    flip_filter.SetInput(img)

    itk_img = flip_filter.GetOutput()
    itk_img.Update()
    itk_img.SetOrigin([0, 0])  # Reset origin AFTER the update

    rec = wl.rectangle_pos

    x = abs(rec[2] - rec[0])
    y = abs(rec[3] - rec[1])

    wanted_size = [x, y]

    # Keep the angle between -90 and 90 degrees
    angle = layer.angle
    while angle > 90:
        angle -= 180

    while angle < -90:
        angle += 180

    orig_size = [layer.size_x, layer.size_y]

    # Base postion = lower left ...
    #  + 90 degrees = upper left position of the initial image
    ang1 = (orig_size[0] / 2.0) * math.cos(math.radians(angle + 90))
    ang2 = (orig_size[0] / 2.0) * math.sin(math.radians(angle + 90))
    ang3 = (orig_size[1] / 2.0) * math.sin(math.radians(angle + 90))
    ang4 = (orig_size[1] / 2.0) * math.cos(math.radians(angle + 90))

    act_x0 = ang1 - ang3
    act_y0 = ang2 + ang4

    # New position
    pos = [rec[0] - act_x0, -(rec[1] - act_y0)]

    resample_filter, spacing = resample(
        itk_img, angle, wanted_size, pos, orig_size)
    resample_filter.Update()

    write_to_tiff(resample_filter.GetOutput(), spacing, path)


def write_to_tiff(img, spacing, path):
    """Write the tiff to disk, with the pixel size [um] in the description
    (as in the tiffs exported by ImageJ, read back by pyAF)."""
    array = itk.array_from_image(img).astype(numpy.uint8)
    description = ("Voxel_size_X: " + repr(float(spacing[0])) + " \xb5m\n"
                   "Voxel_size_Y: " + repr(float(spacing[1])) + " \xb5m\n"
                   "Voxel_size_Z: 0.0 \xb5m\n")
    Image.fromarray(array).save(
        path, tiffinfo={270: description.encode("latin-1")})


def get_itk_image_from_vtk_tiff_2D(vtk_tiff):
    """Transfom VTk tiff in ITK image."""
    vtk_tiff.Update()
    vtk_tiff_imageData = vtk_tiff.GetOutput()

    # The 2D tiff is read as a volume with a single slice
    dims = vtk_tiff_imageData.GetDimensions()
    scalars = vtk_tiff_imageData.GetPointData().GetScalars()
    array = numpy_support.vtk_to_numpy(scalars).reshape(dims[1], dims[0], -1)
    return itk.image_from_array(numpy.ascontiguousarray(
        array[:, :, 0], dtype=numpy.uint8))


def resample(itk_image, angle, wanted_size, position, orig_size):
    """Do the actual resampling here."""
    img_type = itk.Image[itk.UC, 2]

    # Get the angle in radians
    angle = math.radians(angle)

    # Get the new and old sizes
    old_pos1 = [-orig_size[0], -orig_size[1]]
    new_pos1 = [0, 0]
    new_pos1[0] = math.cos(angle) * old_pos1[0] - math.sin(angle) * old_pos1[1]
    new_pos1[1] = math.sin(angle) * old_pos1[0] + math.cos(angle) * old_pos1[1]
    old_pos2 = orig_size
    new_pos2 = [0, 0]
    new_pos2[0] = math.cos(angle) * old_pos2[0] - math.sin(angle) * old_pos2[1]
    new_pos2[1] = math.sin(angle) * old_pos2[0] + math.cos(angle) * old_pos2[1]
    val = max(abs(new_pos2[0]) + abs(new_pos1[0]),
              abs(new_pos2[1]) + abs(new_pos1[1]))

    old_size = orig_size
    new_size = [val, val]

    # Get the new spacing
    old_nbr_pix = itk_image.GetLargestPossibleRegion().GetSize()

    new_nbr_pix = [0, 0]
    new_nbr_pix[0] = int(old_nbr_pix[0] * new_size[0] / old_size[0])
    new_nbr_pix[1] = int(old_nbr_pix[1] * new_size[1] / old_size[1])
    new_spacing = [0, 0]
    new_spacing[0] = (new_size[0] / new_nbr_pix[0]) / 1000.0
    new_spacing[1] = (new_size[1] / new_nbr_pix[1]) / 1000.0

    # The final image is always bigger, add more pixels
    new_nbr_pix[0] = int(wanted_size[0] / (new_spacing[0] * 1000.0))
    new_nbr_pix[1] = int(wanted_size[1] / (new_spacing[1] * 1000.0))

    # New origin of the image
    origin = [position[0] / 1000.0, position[1] / 1000.0]

    # The transformation
    transform = itk.Rigid2DTransform[itk.D].New()
    transform.SetAngle(angle)

    # Use an interpolator
    interp = itk.NearestNeighborInterpolateImageFunction[img_type, itk.D].New()

    # Resample the image
    resample_filter = itk.ResampleImageFilter[img_type, img_type].New()
    resample_filter.SetInterpolator(interp)
    resample_filter.SetDefaultPixelValue(0)  # The new pixels color
    resample_filter.SetOutputSpacing(new_spacing)
    resample_filter.SetOutputOrigin(origin)
    resample_filter.SetSize(new_nbr_pix)
    resample_filter.SetTransform(transform)
    resample_filter.SetInput(itk_image)

    return resample_filter, new_spacing
