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

"""Stack of images."""

import vtk
import itk


def get_stack_images(input_image):
    """Create a list of images."""
    image_type = itk.Image[itk.UC, 3]

    # Convert from ITK to VTK image
    itk_vtk_converter = itk.ImageToVTKImageFilter[image_type].New()
    itk_vtk_converter.SetInput(input_image)
    itk_vtk_converter.Update()
    image_data = vtk.vtkImageData()
    image_data.DeepCopy(itk_vtk_converter.GetOutput())

    # Number of slices in the stack
    nbr_slices = image_data.GetDimensions()[2]

    list_of_images = []

    for slice_pos in range(nbr_slices):
        reslice_filter = reslice(image_data, slice_pos)
        list_of_images.append(reslice_filter)

    return list_of_images


def reslice(image_data, slice_pos):
    """Extract a slice from the stack."""
    # Get some sizes
    (xMin, xMax, yMin, yMax, _, _) = image_data.GetExtent()
    (xSpacing, ySpacing, _) = image_data.GetSpacing()
    (x0, y0, _) = image_data.GetOrigin()

    # Get the center
    center = [x0 + xSpacing * 0.5 * (xMin + xMax),
              y0 + ySpacing * 0.5 * (yMin + yMax),
              slice_pos]

    axial = vtk.vtkMatrix4x4()
    axial.DeepCopy((1, 0, 0, center[0],
                    0, 1, 0, center[1],
                    0, 0, 1, center[2],
                    0, 0, 0, 1))

    reslice_filter = vtk.vtkImageReslice()
    reslice_filter.SetOutputDimensionality(2)
    reslice_filter.SetResliceAxes(axial)
    reslice_filter.SetInputData(image_data)

    return reslice_filter
