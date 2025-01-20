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

"""Makes an isosurface from a stack."""

import vtk
import itk
from .common import center_data


def get_isosurface(input_image, sizes, res, color, flip):
    """Uses Marching cubes."""
    image_type = itk.Image[itk.UC, 3]

    # Flip the image if needed
    if flip != [0, 0, 0]:
        flip_filter = itk.FlipImageFilter[image_type].New()
        flip_filter.SetInput(input_image)
        flip_filter.SetFlipAxes(flip)
        image = flip_filter.GetOutput()
    else:
        image = input_image

    itk_vtk_converter = itk.ImageToVTKImageFilter[image_type].New()
    itk_vtk_converter.SetInput(image)
    itk_vtk_converter.Update()

    imageData = vtk.vtkImageData()
    imageData.DeepCopy(itk_vtk_converter.GetOutput())
    # Spacing = 1, 1, 1 by default
    imageData.SetSpacing(res[0], res[1], res[2])

    # Create a 3D model using marching cubes
    mc = vtk.vtkMarchingCubes()
    mc.SetInputData(imageData)
    mc.ComputeNormalsOn()
    mc.ComputeGradientsOn()
    mc.SetValue(0, 1)

    # Get all the regions
    confilter = vtk.vtkPolyDataConnectivityFilter()
    confilter.SetInputConnection(mc.GetOutputPort())
    confilter.SetExtractionModeToAllRegions()

    # Can be useful sometimes but mostly buggy for normal confocal stacks
    # fillHolesFilter = vtk.vtkFillHolesFilter()
    # fillHolesFilter.SetInputConnection(confilter.GetOutputPort())
    # fillHolesFilter.SetHoleSize(1000000.0)

    # Make the triangle order consistent (else, the order is not the same
    # for the filled triangles, resuling on bad colors)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(confilter.GetOutputPort())
    normals.ConsistencyOn()
    normals.SplittingOff()
    normals.Update()

    # Center the data
    poly_data = center_data(normals.GetOutput(), sizes[0], sizes[1])

    mapper = vtk.vtkPolyDataMapper()
    actor = vtk.vtkActor()

    actor.SetMapper(mapper)
    mapper.SetInputData(poly_data)
    mapper.ScalarVisibilityOff()

    prop = vtk.vtkProperty()
    prop.SetAmbient(0.1)
    prop.SetDiffuse(0.1)
    prop.SetSpecular(0.5)
    prop.SetColor(color[0], color[1], color[2])
    prop.SetLineWidth(1.0)
    prop.SetRepresentationToSurface()

    actor.SetProperty(prop)

    return actor
