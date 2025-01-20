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

"""Creates the vtk actors for the tomographies."""

import vtk
import numpy
from ... import widgets_list
from ...tools.colortables import ColorTables
from .common import center_data


def get_surface(parent, data, mode):
    """Get the actor for the surface and the bottom of the tomography.

    The tomography has 6 faces : surface, bottom, and 4 sides

    For the smoothing, one has to ensure that the tomography is a closed
    volume with no redundant triangles. For this, the surface and bottom
    arrays are shifted by a z_offset value of 0.001. In fact they shoud be
    at the same position when the indenation is maximal, but having the offset
    ensures that those are two separate surfaces (surface and bottom).
    """
    colortable = ColorTables(data.colortableid, data.color_saturation,
                             data.color_negative, data.colortable_max_value,
                             data.colortable_middle_value, mode="scalarMap")

    # Get the array with the stiffness
    if data.stiffness_corrected:
        array = data.stiffness_corrected_array
    else:
        array = data.stiffness_array

    widgets_list.widget_progressbar.set_label("Prepare data (" + mode + ")")

    # Fill the vals array containing the stiffness values
    vals = numpy.zeros([data.nbr_pixels_x, data.nbr_pixels_y])
    nans = numpy.zeros([data.nbr_pixels_x, data.nbr_pixels_y])
    if mode == "top":
        tp = data.opengl_surf_type

        for i in range(data.nbr_pixels_x):
            for j in range(data.nbr_pixels_y):
                if tp == "Topo_surf":
                    if len(array[i][j]) > data.stiffness_depth_view:
                        vals[i][j] = data.topography[i][j]
                    else:
                        nans[i][j] = True
                else:
                    if len(array[i][j]) > data.stiffness_depth_view:
                        # Take the stiffness value of the slice
                        vals[i][j] = array[i][j][data.stiffness_depth_view]
                    else:
                        # Display as nan
                        nans[i][j] = True

            widgets_list.widget_progressbar.update()

    elif mode == "bottom":
        for i in range(data.nbr_pixels_x):
            for j in range(data.nbr_pixels_y):
                if len(array[i][j]) > data.stiffness_depth_view:
                    # Display the last aviable position
                    vals[i][j] = array[i][j][-1]
                else:
                    # Display as nan
                    nans[i][j] = True

            widgets_list.widget_progressbar.update()

    widgets_list.widget_progressbar.set_label("Assign colors (" + mode + ")")

    # Convert the stiffness values to color values with matplotlib
    carray = numpy.array(colortable.scalarMap.to_rgba(vals, bytes=True))

    # scalarMap.to_rgba can not replace nans. In this case we just
    # post-process the color array (this is very fast) and replace the
    # nans by a transparent color (or the color chosen by the user)
    if data.nan_color_transparent:
        color = [0, 0, 0, 0]
    else:
        color = numpy.array(data.color_nan) * 255.0

    # Replace the nans by the nan color
    for i in range(data.nbr_pixels_x):
        for j in range(data.nbr_pixels_y):
            if nans[i][j]:
                carray[i][j] = color

    # Make a z_array containing the height for each pixel
    dpt = data.stiffness_depth_view * data.used_indentation_step
    z_array = numpy.zeros([data.nbr_pixels_x, data.nbr_pixels_y])
    bottom_array = numpy.zeros([data.nbr_pixels_x, data.nbr_pixels_y])

    # Z Offset for the bottom of the tomography. Ensures that those triangles
    # are not at the exact same position as the topograhpy : this makes the
    # smoothing algorithm segfault.
    z_offset = 0.001

    used_st = data.used_indentation_step
    for i in range(data.opengl_slice_bottom, data.opengl_slice_top + 1, 1):
        for j in range(data.opengl_slice_left, data.opengl_slice_right + 1, 1):
            if mode == "top":
                if data.opengl_afm_flat_view:
                    z_array[i][j] = parent.pos_z
                else:
                    z_array[i][j] = data.topography[i][j] - dpt
            if mode == "top" or mode == "bottom":
                val = (len(array[i][j]) - 1) * used_st - z_offset
                bottom_array[i][j] = data.topography[i][j] - val

    # The z positions which will be used are the ones from the z_array
    # So in case of the bottom array redefine the z_array
    if mode == "bottom":
        z_array = bottom_array

    widgets_list.widget_progressbar.set_label("Create vertices (" + mode + ")")

    # Create the triangles
    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(4)

    points = vtk.vtkPoints()

    triangles = vtk.vtkCellArray()

    count = 0

    # Store all the z values in a new array, we want to know the height of the
    # tomography for the z scale bar
    tmp_z = []

    for i in range(data.opengl_slice_bottom, data.opengl_slice_top, 1):
        for j in range(data.opengl_slice_left, data.opengl_slice_right, 1):
            z1 = z_array[i][j]
            z2 = z_array[i][j + 1]
            z3 = z_array[i + 1][j]

            # Storing z1 is enough
            tmp_z.append(z1)

            # In this case, the point is below the maximum depth, so use
            # just the maximum depth for the point
            # Reshift at + z_offset to be at the right position
            if data.opengl_afm_flat_view is False and mode == "top":
                if z1 < bottom_array[i][j]:
                    z1 = bottom_array[i][j] + z_offset
                if z2 < bottom_array[i][j + 1]:
                    z2 = bottom_array[i][j + 1] + z_offset
                if z3 < bottom_array[i + 1][j]:
                    z3 = bottom_array[i + 1][j] + z_offset

            points.InsertNextPoint(i * data.x_size, j * data.y_size, z1)
            points.InsertNextPoint(i * data.x_size,
                                   j * data.y_size + data.y_size, z2)
            points.InsertNextPoint(i * data.x_size + data.x_size,
                                   j * data.y_size, z3)

            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, count)
            triangle.GetPointIds().SetId(1, count + 1)
            triangle.GetPointIds().SetId(2, count + 2)

            triangles.InsertNextCell(triangle)

            z1 = z_array[i][j + 1]
            z2 = z_array[i + 1][j + 1]
            z3 = z_array[i + 1][j]

            # Storing z1 is enough
            tmp_z.append(z1)

            # In this case, the point is below the maximum depth, so use
            # just the maximum depth for the point
            # Reshift at + z_offset to be at the right position
            if data.opengl_afm_flat_view is False and mode == "top":
                if z1 < bottom_array[i][j + 1]:
                    z1 = bottom_array[i][j + 1] + z_offset
                if z2 < bottom_array[i + 1][j + 1]:
                    z2 = bottom_array[i + 1][j + 1] + z_offset
                if z3 < bottom_array[i + 1][j]:
                    z3 = bottom_array[i + 1][j] + z_offset

            points.InsertNextPoint(i * data.x_size,
                                   j * data.y_size + data.y_size, z1)
            points.InsertNextPoint(i * data.x_size + data.x_size,
                                   j * data.y_size + data.y_size, z2)
            points.InsertNextPoint(i * data.x_size + data.x_size,
                                   j * data.y_size, z3)

            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, count + 3)
            triangle.GetPointIds().SetId(1, count + 4)
            triangle.GetPointIds().SetId(2, count + 5)

            count += 6

            triangles.InsertNextCell(triangle)
            colors.InsertNextTuple(carray[i][j])
            colors.InsertNextTuple(carray[i][j + 1])
            colors.InsertNextTuple(carray[i + 1][j])
            colors.InsertNextTuple(carray[i][j + 1])
            colors.InsertNextTuple(carray[i + 1][j + 1])
            colors.InsertNextTuple(carray[i + 1][j])

            widgets_list.widget_progressbar.update()

    # Create a polydata object
    polydata = vtk.vtkPolyData()

    # Add the geometry and topology to the polydata
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(colors)
    polydata.SetPolys(triangles)

    # Translate to 0,0,0
    polydata = center_data(polydata, data.scan_size_x, data.scan_size_y)

    # Get the z scale
    tmp_z = numpy.array(tmp_z)
    if mode == "bottom":
        z_scale = numpy.amin(tmp_z)
    elif mode == "top":
        z_scale = numpy.amax(tmp_z)

    return polydata, z_scale


def get_side(data, side):
    """Method to build the 4 actors for the sides of the tomography."""
    widgets_list.widget_progressbar.set_label("Prepare data (" + side + ")")

    colortable = ColorTables(data.colortableid, data.color_saturation,
                             data.color_negative, data.colortable_max_value,
                             data.colortable_middle_value, mode="scalarMap")

    if data.stiffness_corrected:
        array = data.stiffness_corrected_array
    else:
        array = data.stiffness_array

    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(4)

    points = vtk.vtkPoints()

    triangles = vtk.vtkCellArray()

    step = data.indentation_step
    depth = data.stiffness_depth_view

    # Get the maximum indentation depth for a side, and the number of triangles
    # to draw for this side (=length).
    maxdepth = 0
    if side == "bottom":
        j = data.opengl_slice_left  # One position stays fixed

        for i in range(data.opengl_slice_bottom, data.opengl_slice_top + 1, 1):
            if maxdepth < len(array[i][j]):
                maxdepth = len(array[i][j])

        length = data.opengl_slice_top - data.opengl_slice_bottom + 1

    elif side == "top":
        j = data.opengl_slice_right  # One position stays fixed

        for i in range(data.opengl_slice_bottom, data.opengl_slice_top + 1, 1):
            if maxdepth < len(array[i][j]):
                maxdepth = len(array[i][j])

        length = data.opengl_slice_top - data.opengl_slice_bottom + 1

    elif side == "left":
        i = data.opengl_slice_bottom  # One position stays fixed

        for j in range(data.opengl_slice_left, data.opengl_slice_right + 1, 1):
            if maxdepth < len(array[i][j]):
                maxdepth = len(array[i][j])

        length = data.opengl_slice_right - data.opengl_slice_left + 1

    elif side == "right":
        i = data.opengl_slice_top  # One position stays fixed

        for j in range(data.opengl_slice_left, data.opengl_slice_right + 1, 1):
            if maxdepth < len(array[i][j]):
                maxdepth = len(array[i][j])

        length = data.opengl_slice_right - data.opengl_slice_left + 1

    # Fill an empty array. An array needs to be rectangular at least, so some
    # cells will stay empty, but those will never be called, so it's not
    # important.
    # Note : at the beggining I used a numpy.empty array, but the unchanged
    # values would then be numpy.nan or val*1e307 ... This was not important
    # because these values were never used. But scalarMap.to_rgba will
    # complain with a warning because it actually parses these values.
    # So it is better to use a numpy.zeros array here.
    vals = numpy.zeros((length, maxdepth))

    # Fill the vals array with stiffness values comming from array.
    # The loop is going from 0 to lenght, but the values to be fetched are
    # offseted by the slice position (a = i + offset) or (a = j + offset)

    if side == "bottom":
        for i in range(length):
            a = i + data.opengl_slice_bottom
            lencurrent = len(array[a][j])

            for z in range(depth, lencurrent, 1):
                vals[i][z] = array[a][j][z]

    elif side == "top":
        for i in range(length):
            a = i + data.opengl_slice_bottom
            lencurrent = len(array[a][j])

            for z in range(depth, lencurrent, 1):
                vals[i][z] = array[a][j][z]

    elif side == "left":
        for j in range(length):
            a = j + data.opengl_slice_left
            lencurrent = len(array[i][a])

            for z in range(depth, lencurrent, 1):
                vals[j][z] = array[i][a][z]

    elif side == "right":
        for j in range(length):
            a = j + data.opengl_slice_left
            lencurrent = len(array[i][a])

            for z in range(depth, lencurrent, 1):
                vals[j][z] = array[i][a][z]

    widgets_list.widget_progressbar.set_label("Assign colors (" + side + ")")

    # Fill the color_array with values. Matplotlib is used to fill in the
    # values, which is the fastest solution I found for the moment.
    carray = numpy.array(colortable.scalarMap.to_rgba(vals, bytes=True))

    widgets_list.widget_progressbar.set_label("Create vertices (" + side + ")")

    # Create the triangles
    if side == "bottom" or side == "top":
        start = 0
        stop = data.opengl_slice_top - data.opengl_slice_bottom
        b = j  # Fixed position
        count = 0

        for a in range(start, stop, 1):
            # a goes from start to stop == from 0 to the size of the slice
            # in fact a corresponds to the indices used in the color array
            # which has a length of the size of the slice
            # c is the real indice (in real wold coordinates), so it needs to
            # be shifted
            c = a + data.opengl_slice_bottom

            lencurrent = len(array[c][b])
            lennext = len(array[c + 1][b])
            finallen = max(lencurrent, lennext) - 1

            for z in range(depth, finallen, 1):
                # Break out of the loop and do not draw anything. Else the
                # triangles will get out of bounds
                if (finallen > lennext - 2) and z >= (lencurrent - 1):
                    break

                z1 = data.topography[c][b] - z * step
                if z >= lennext:
                    z2 = data.topography[c + 1][b] - (lennext - 1) * step
                else:
                    z2 = data.topography[c + 1][b] - z * step
                z3 = data.topography[c][b] - z * step - step

                points.InsertNextPoint(c * data.x_size, b * data.y_size, z1)
                points.InsertNextPoint(c * data.x_size + data.x_size,
                                       b * data.y_size, z2)
                points.InsertNextPoint(c * data.x_size, b * data.y_size, z3)

                triangle = vtk.vtkTriangle()
                triangle.GetPointIds().SetId(0, count)
                triangle.GetPointIds().SetId(1, count + 1)
                triangle.GetPointIds().SetId(2, count + 2)

                triangles.InsertNextCell(triangle)

                count += 3

                colors.InsertNextTuple(carray[a][z])
                colors.InsertNextTuple(carray[a + 1][z])
                colors.InsertNextTuple(carray[a][z + 1])

            c += 1
            a += 1

            lencurrent = len(array[c][b])
            lenbefore = len(array[c - 1][b])
            finallen = max(lencurrent, lenbefore) - 1

            for z in range(depth, finallen, 1):
                # Break out of the loop and do not draw anything. Else the
                # triangles will get out of bounds
                if (finallen > lenbefore - 2) and z >= (lencurrent - 1):
                    break

                z1 = data.topography[c][b] - z * step
                if z >= lenbefore - 1:
                    z2 = data.topography[
                        c - 1][b] - lenbefore * step + step
                else:
                    z2 = data.topography[c - 1][b] - z * step - step

                z3 = data.topography[c][b] - z * step - step

                points.InsertNextPoint(c * data.x_size, b * data.y_size, z1)
                points.InsertNextPoint(c * data.x_size - data.x_size,
                                       b * data.y_size, z2)
                points.InsertNextPoint(c * data.x_size, b * data.y_size, z3)

                triangle = vtk.vtkTriangle()
                triangle.GetPointIds().SetId(0, count)
                triangle.GetPointIds().SetId(1, count + 1)
                triangle.GetPointIds().SetId(2, count + 2)

                triangles.InsertNextCell(triangle)

                count += 3

                colors.InsertNextTuple(carray[a][z])
                colors.InsertNextTuple(carray[a - 1][z + 1])
                colors.InsertNextTuple(carray[a][z + 1])

            c -= 1
            a -= 1

            widgets_list.widget_progressbar.update()

    elif side == "left" or side == "right":
        start = 0
        stop = data.opengl_slice_right - data.opengl_slice_left
        a = i  # Fixed position
        count = 0

        for b in range(start, stop, 1):
            # b goes from start to stop == from 0 to the size of the slice
            # in fact a corresponds to the indices used in the color array
            # which has a length of the size of the slice
            # c is the real indice (in real wold coordinates), so it needs to
            # be shifted
            c = b + data.opengl_slice_left

            lencurrent = len(array[a][c])
            lennext = len(array[a][c + 1])
            finallen = max(lencurrent, lennext) - 1

            for z in range(depth, finallen, 1):
                # Break out of the loop and do not draw anything. Else the
                # triangles will get out of bounds
                if (finallen > lennext - 2) and z >= (lencurrent - 1):
                    break

                z1 = data.topography[a][c] - z * step
                if z >= lennext:
                    z2 = data.topography[a][c + 1] - (lennext - 1) * step
                else:
                    z2 = data.topography[a][c + 1] - z * step
                z3 = data.topography[a][c] - z * step - step

                points.InsertNextPoint(a * data.x_size, c * data.y_size, z1)
                points.InsertNextPoint(a * data.x_size,
                                       c * data.y_size + data.y_size, z2)
                points.InsertNextPoint(a * data.x_size, c * data.y_size, z3)

                triangle = vtk.vtkTriangle()
                triangle.GetPointIds().SetId(0, count)
                triangle.GetPointIds().SetId(1, count + 1)
                triangle.GetPointIds().SetId(2, count + 2)

                triangles.InsertNextCell(triangle)

                count += 3

                colors.InsertNextTuple(carray[b][z])
                colors.InsertNextTuple(carray[b + 1][z])
                colors.InsertNextTuple(carray[b][z + 1])

            c += 1
            b += 1

            lencurrent = len(array[a][c])
            lenbefore = len(array[a][c - 1])
            finallen = max(lencurrent, lenbefore) - 1

            for z in range(depth, finallen, 1):
                # Break out of the loop and do not draw anything. Else the
                # triangles will get out of bounds
                if (finallen > lenbefore - 2) and z >= (lencurrent - 1):
                    break

                z1 = data.topography[a][c] - z * step
                if z >= lenbefore - 1:
                    z2 = data.topography[a][
                        c - 1] - lenbefore * step + step
                else:
                    z2 = data.topography[a][c - 1] - z * step - step

                z3 = data.topography[a][c] - z * step - step

                points.InsertNextPoint(a * data.x_size, c * data.y_size, z1)
                points.InsertNextPoint(a * data.x_size,
                                       c * data.y_size - data.y_size, z2)
                points.InsertNextPoint(a * data.x_size, c * data.y_size, z3)

                triangle = vtk.vtkTriangle()
                triangle.GetPointIds().SetId(0, count)
                triangle.GetPointIds().SetId(1, count + 1)
                triangle.GetPointIds().SetId(2, count + 2)

                triangles.InsertNextCell(triangle)

                count += 3

                colors.InsertNextTuple(carray[b][z])
                colors.InsertNextTuple(carray[b - 1][z + 1])
                colors.InsertNextTuple(carray[b][z + 1])

            c -= 1
            b -= 1
            widgets_list.widget_progressbar.update()

    # Create a polydata object
    polydata = vtk.vtkPolyData()

    # Add the geometry and topology to the polydata
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(colors)
    polydata.SetPolys(triangles)

    # Translate to 0,0,0
    polydata = center_data(polydata, data.scan_size_x, data.scan_size_y)

    return polydata


def get_surface_as_squares(data, off):
    """Uses triangles to draw rectangles. No color interpolation is done.

    The offset parameter is used to center the data like the triangle surface.
    """
    colortable = ColorTables(data.colortableid, data.color_saturation,
                             data.color_negative, data.colortable_max_value,
                             data.colortable_middle_value, mode="scalarMap")

    # Get the array with the stiffness
    if data.stiffness_corrected:
        array = data.stiffness_corrected_array
    else:
        array = data.stiffness_array

    # Fill the vals array containing the stiffness values
    vals = numpy.zeros([data.nbr_pixels_x, data.nbr_pixels_y])
    nans = numpy.zeros([data.nbr_pixels_x, data.nbr_pixels_y])
    depth = data.stiffness_depth_view
    for i in range(data.nbr_pixels_x):
        for j in range(data.nbr_pixels_y):
            if len(array[i][j]) > data.stiffness_depth_view:
                vals[i][j] = data.stiffness_array[i][j][depth]
            else:
                nans[i][j] = True

    # Convert the stiffness values to color values with matplotlib
    carray = numpy.array(colortable.scalarMap.to_rgba(vals, bytes=True))

    # scalarMap.to_rgba can not replace nans. In this case we just
    # post-process the color array (this is very fast) and replace the
    # nans by a transparent color (or the color chosen by the user)
    if data.nan_color_transparent:
        color = [0, 0, 0, 0]
    else:
        color = numpy.array(data.color_nan) * 255.0

    # Replace the nans by the nan color
    for i in range(data.nbr_pixels_x):
        for j in range(data.nbr_pixels_y):
            if nans[i][j]:
                carray[i][j] = color

    # Create the triangles
    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(4)

    points = vtk.vtkPoints()

    triangles = vtk.vtkCellArray()

    count = 0

    for i in range(data.opengl_slice_bottom, data.opengl_slice_top + 1, 1):
        for j in range(data.opengl_slice_left, data.opengl_slice_right + 1, 1):
            points.InsertNextPoint(i * data.x_size, j * data.y_size, 0)
            points.InsertNextPoint(i * data.x_size,
                                   j * data.y_size + data.y_size, 0)
            points.InsertNextPoint(i * data.x_size + data.x_size,
                                   j * data.y_size, 0)

            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, count)
            triangle.GetPointIds().SetId(1, count + 1)
            triangle.GetPointIds().SetId(2, count + 2)

            triangles.InsertNextCell(triangle)

            points.InsertNextPoint(i * data.x_size,
                                   j * data.y_size + data.y_size, 0)
            points.InsertNextPoint(i * data.x_size + data.x_size,
                                   j * data.y_size + data.y_size, 0)
            points.InsertNextPoint(i * data.x_size + data.x_size,
                                   j * data.y_size, 0)

            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, count + 3)
            triangle.GetPointIds().SetId(1, count + 4)
            triangle.GetPointIds().SetId(2, count + 5)

            count += 6

            triangles.InsertNextCell(triangle)
            colors.InsertNextTuple(carray[i][j])
            colors.InsertNextTuple(carray[i][j])
            colors.InsertNextTuple(carray[i][j])
            colors.InsertNextTuple(carray[i][j])
            colors.InsertNextTuple(carray[i][j])
            colors.InsertNextTuple(carray[i][j])

    # Create a polydata object
    polydata = vtk.vtkPolyData()

    # Add the geometry and topology to the polydata
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(colors)
    polydata.SetPolys(triangles)

    # Translate to 0,0,0
    polydata = center_data(polydata, data.scan_size_x, data.scan_size_y, off)

    return polydata, 0


def get_actor_surface_as_cubes(self, data, mode):
    """Not working for the moment"""
    colortable = ColorTables(data.colortableid, data.color_saturation,
                             data.color_negative, data.colortable_max_value,
                             data.colortable_middle_value, mode="scalarMap",
                             color_nan=data.color_nan)

    if data.stiffness_corrected:
        array = data.stiffness_corrected_array
    else:
        array = data.stiffness_array

    points = vtk.vtkPoints()
    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(4)
    val_E = numpy.nan
    for i in range(data.opengl_slice_bottom, data.opengl_slice_top + 1, 1):
        for j in range(data.opengl_slice_left, data.opengl_slice_right + 1, 1):
            if mode == "top":
                if data.opengl_afm_flat_view:
                    z = self.pos_z
                else:
                    z = data.topography[i][j] - \
                        data.stiffness_depth_view * data.indentation_step
            elif mode == "bottom":
                z = data.topography[i][j] - \
                    (len(array[i][j]) - 1) * data.indentation_step

            points.InsertNextPoint(i * data.x_size, j * data.y_size, z)

            if len(array[i][j]) > data.stiffness_depth_view:
                val_E = array[i][j][data.stiffness_depth_view]
            else:
                val_E = numpy.nan

            color_st = colortable.get_color_as_list(val_E)
            color = [int(color_st[0] * 255),
                     int(color_st[1] * 255),
                     int(color_st[2] * 255),
                     int(color_st[3] * 255)]
            colors.InsertNextTuple(color)

    # Combine into a polydata
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(colors)

    # Translate to 0,0,0
    polydata = center_data(polydata, data.scan_size_x, data.scan_size_y)

    # Create anything you want here, we will use a cube for the demo.
    cubeSource = vtk.vtkCubeSource()
    cubeSource.SetXLength(data.x_size)
    cubeSource.SetYLength(data.y_size)
    cubeSource.SetZLength(data.indentation_step)

    glyph3D = vtk.vtkGlyph3D()
    glyph3D.SetColorModeToColorByScalar()
    glyph3D.SetSourceConnection(cubeSource.GetOutputPort())
    glyph3D.SetInput(polydata)
    glyph3D.ScalingOff()
    glyph3D.Update()

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph3D.GetOutputPort())

    actorsurface = vtk.vtkActor()
    actorsurface.SetMapper(mapper)

    return actorsurface


def get_surface_as_points(self, data, mode):
    """Topography as single dots."""
    colortable = ColorTables(data.colortableid, data.color_saturation,
                             data.color_negative, data.colortable_max_value,
                             data.colortable_middle_value, mode="scalarMap",
                             color_nan=data.color_nan)

    points = vtk.vtkPoints()
    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(4)
    val_E = numpy.nan
    for i in range(data.opengl_slice_bottom, data.opengl_slice_top + 1, 1):
        for j in range(data.opengl_slice_left, data.opengl_slice_right + 1, 1):
            if mode == "top":
                if data.opengl_afm_flat_view:
                    z = self.pos_z
                else:
                    z = data.topography[i][j]

            points.InsertNextPoint(i * data.x_size, j * data.y_size, z)

            if len(data.stiffness_array[i][j]) > data.stiffness_depth_view:
                val_E = data.topography[i][j]
            else:
                val_E = numpy.nan

            color_st = colortable.get_color_as_list(val_E)
            color = [int(color_st[0] * 255),
                     int(color_st[1] * 255),
                     int(color_st[2] * 255),
                     int(color_st[3] * 255)]
            colors.InsertNextTuple(color)

    # Combine into a polydata
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(colors)

    # Translate to 0,0,0
    polydata = center_data(polydata, data.scan_size_x, data.scan_size_y)

    # Create anything you want here, we will use a cube for the demo.
    cubeSource = vtk.vtkDiskSource()
    cubeSource.SetInnerRadius(0)
    cubeSource.SetOuterRadius(50)

    glyph3D = vtk.vtkGlyph3D()
    glyph3D.SetColorModeToColorByScalar()
    glyph3D.SetSourceConnection(cubeSource.GetOutputPort())
    glyph3D.ScalingOff()
    glyph3D.SetInputData(polydata)
    glyph3D.Update()

    return glyph3D


def get_actor_profile(sizes, trans, z_scale, profile, polydata_tomo):
    """Build a profile"""
    x_size = sizes["x_size"]
    y_size = sizes["y_size"]
    trans_x = trans["trans_x"]
    trans_y = trans["trans_y"]
    trans_z = trans["trans_z"]
    z_min = z_scale[0]
    z_max = z_scale[1]

    x = profile[0]
    y = profile[1]
    # Start position
    x1 = (x[0])*x_size - trans_x
    y1 = (y[0])*y_size - trans_y
    # End position
    x2 = (x[-1])*x_size - trans_x
    y2 = (y[-1])*y_size - trans_y

    # Define a cellLocator to be able to compute intersections between lines
    # and the surface
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(polydata_tomo.GetOutput())
    locator.BuildLocator()

    maxloop = 1000
    tolerance = 0.001

    # Make a list of points. Each point is the intersection of a vertical line
    # defined by p1 and p2 and the surface.
    distx = x2 - x1
    disty = y2 - y1
    resx = distx/maxloop
    resy = disty/maxloop

    points = vtk.vtkPoints()
    for i in range(maxloop):
        # For z, use the the z min and z max values, and add 100 as security
        p1 = [x1, y1, trans_z + z_min - 100]
        p2 = [x1+resx, y1+resy, trans_z + z_max + 100]
        x1 += resx
        y1 += resy

        # Outputs (we need only pos which is the x, y, z position
        # of the intersection)
        t = vtk.mutable(0)
        pos = [0.0, 0.0, 0.0]
        pcoords = [0.0, 0.0, 0.0]
        subId = vtk.mutable(0)
        locator.IntersectWithLine(p1, p2, tolerance, t, pos, pcoords, subId)

        # Add a slight offset in z
        pos[2] += 0.1
        # Add the x, y, z position of the intersection
        points.InsertNextPoint(pos)

    # Define a line
    lines = vtk.vtkCellArray()
    lines.InsertNextCell(maxloop)
    for i in range(maxloop):
        lines.InsertCellPoint(i)

    # Create a polydata for the tube
    tube_data = vtk.vtkPolyData()
    tube_data.SetPoints(points)
    tube_data.SetLines(lines)

    # Define a tube radius (constant radius for each position)
    tubeRadius = vtk.vtkDoubleArray()
    tubeRadius.SetName("TubeRadius")
    tubeRadius.SetNumberOfTuples(maxloop)
    for i in range(maxloop):
        tubeRadius.SetTuple1(i, 20)
    tube_data.GetPointData().AddArray(tubeRadius)
    tube_data.GetPointData().SetActiveScalars("TubeRadius")

    # Color
    colors = vtk.vtkUnsignedCharArray()
    colors.SetName("Colors")
    colors.SetNumberOfComponents(4)
    colors.SetNumberOfTuples(maxloop)
    for i in range(maxloop):
        colors.InsertTuple4(i, 255, 0, 0, 255)
    tube_data.GetPointData().AddArray(colors)

    # Create the tube
    tube = vtk.vtkTubeFilter()
    tube.SetInputData(tube_data)
    tube.SetNumberOfSides(20)
    tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    # Close the tube at the ends
    tube.CappingOn()

    # Map the spline
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    # Activate the color
    mapper.ScalarVisibilityOn()
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray("Colors")

    # Define the actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor
