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

import vtk


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

        dx = 255 - 0

        dy = lut_options["max_color"][0] - lut_options["min_color"][0]
        slope1 = dy/dx
        dy = lut_options["max_color"][1] - lut_options["min_color"][1]
        slope2 = dy/dx
        dy = lut_options["max_color"][2] - lut_options["min_color"][2]
        slope3 = dy/dx

        # Define the linear LUT
        for i in range(256):
            lut.SetTableValue(i, i*slope1, i*slope2, i*slope3, 1.0)

        # Set values to transparent (depending on threshold value and
        # if LUT is inverted or not)
        if lut_options["invert"]:
            start = 255
            end = 255 - lut_options["threshold"]
            step = -1
        else:
            start = 0
            end = lut_options["threshold"]
            step = 1

        for i in range(start, end, step):
            lut.SetTableValue(i, 0, 0, 0, 0)

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
