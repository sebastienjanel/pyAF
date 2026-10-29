# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
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

"""Layer for the AFM tomography."""

from .layer import Layer
import numpy
import vtk
from vtkmodules.util import numpy_support
from ... import widgets_list
from ... import shared
from ..utils import vtk_afm
from ...tools import map_values
from ...widgets.progressbar import Progressbar


class LayerAFM(Layer):
    """Layer for the AFM tomography."""

    def __init__(self, parent, data_id):
        super().__init__()

        self.parent = parent
        self.afm_id = data_id
        self.data = shared.exp.list[self.afm_id]
        self.type = "afm"
        self.filename = self.data.filename
        self.append_filter = None
        self.surface_polydata = None

        # Map drawn: one of the computed maps (elasticity by default)
        if self.data.opengl_map_type not in \
                map_values.available_maps(self.data):
            self.data.opengl_map_type = "stiffness"

        # The tomography actor
        self.actor = None
        # The profiles
        self.actor_profiles = []

        self.size_x = self.data.scan_size_x
        self.size_y = self.data.scan_size_y
        self.nbr_pixels_x = self.data.nbr_pixels_x
        self.nbr_pixels_y = self.data.nbr_pixels_y
        self.res_x = self.data.x_size
        self.res_y = self.data.y_size

        self.init_axes()

    def update_tomography(self):
        """Removes the tomography and recreates a new one."""
        # Create a progressbar which is displayed during the computations
        Progressbar("Preparing afm data")

        dt = self.data

        poly_data_bottom = None
        z_min = 0

        value = dt.nbr_pixels_x * dt.nbr_pixels_y + dt.nbr_pixels_x

        if dt.opengl_afm_flat_view is False:
            value = 2 * value
            value += abs(dt.opengl_slice_bottom - dt.opengl_slice_top)
            value += abs(dt.opengl_slice_left - dt.opengl_slice_right)

        widgets_list.widget_progressbar.set_range(0, value)

        # Remove old actor (not the first time)
        if shared.VTK_first is False:
            self.parent.renderer.RemoveActor(self.actor)

        # Create new actors: surface of the selected map, and the bottom and
        # sides of the elasticity volume (Tomography)
        poly_data_top, z_max = vtk_afm.get_surface(self, dt, "top")

        if vtk_afm.is_volume(dt):
            poly_data_bottom, z_min = vtk_afm.get_surface(self, dt, "bottom")
            sides = [vtk_afm.get_side(dt, side)
                     for side in ("bottom", "top", "left", "right")]
            if dt.vtk_smoothing_type is None:
                # Colors per triangle, as the surface without smoothing (the
                # appended polydata keep only the arrays common to all)
                sides = [point_to_cell_colors(side) for side in sides]
            poly_data_sidebottom, poly_data_sidetop, poly_data_sideleft, \
                poly_data_sideright = sides

        self.append_filter = vtk.vtkAppendPolyData()
        self.append_filter.AddInputData(poly_data_top)
        if poly_data_bottom is not None:
            # Not the empty sides (e.g. along a substrate without
            # indentation: no colors, which disturbs the append)
            for part in (poly_data_bottom, poly_data_sidebottom,
                         poly_data_sidetop, poly_data_sideleft,
                         poly_data_sideright):
                if part.GetNumberOfCells() > 0:
                    self.append_filter.AddInputData(part)

        self.append_filter.Update()
        self.poly_data_tomo = self.append_filter

        self.actor = self.get_smoothed_actor()

        self.z_scale = [z_min, z_max]

        self.finish_update_tomography()

    def finish_update_tomography(self):
        """Display axes and render."""
        # Add actor to renderer
        self.parent.renderer.AddActor(self.actor)

        self.draw_at_position()

        # Reinit the axes (if slicing was applied, will update the axes)
        self.init_axes()

        # Hide or display layer, mainly to hide the axes if we are repainting
        # and the axes were hidden
        self.hide_or_display_layer()

        # Close the progressbar
        if widgets_list.widget_progressbar is not None:
            widgets_list.widget_progressbar.close()

    def get_smoothed_actor(self):
        """Smooth the tomography actor."""
        # Use a filter to smooth the data (will add triangles and smooth),
        # off by default
        smooth = None
        if self.data.vtk_smoothing_type is not None:
            smooth = vtk.vtkLoopSubdivisionFilter()

        if smooth is not None:
            # Clean the polydata so that the edges are shared !
            clean_polydata = vtk.vtkCleanPolyData()
            tomo = self.poly_data_tomo.GetOutputPort()
            clean_polydata.SetInputConnection(tomo)

            smooth.SetNumberOfSubdivisions(self.data.vtk_smoothing_iterations)
            smooth.SetInputConnection(clean_polydata.GetOutputPort())

            decimate = vtk.vtkDecimatePro()
            decimate.SetInputConnection(smooth.GetOutputPort())
            decimate.SetTargetReduction(
                self.data.vtk_decimate_target_reduction)
            decimate.SetFeatureAngle(80)  # Angle above 80 > == edge
            decimate.PreserveTopologyOff()
            decimate.SplittingOn()

        # Create a mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        if smooth is not None:
            self.poly_data_tomo = decimate
            mapper.SetInputConnection(self.poly_data_tomo.GetOutputPort())
        else:
            clean_polydata = vtk.vtkCleanPolyData()
            tomo = self.poly_data_tomo.GetOutputPort()
            clean_polydata.SetInputConnection(tomo)
            self.poly_data_tomo = clean_polydata
            mapper.SetInputConnection(self.poly_data_tomo.GetOutputPort())
        mapper.Update()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        return actor

    def draw_at_position(self):
        """Some drawing updates and render call."""
        self.translate(0, 0, 0)
        self.rotate(None, None)

        # Update status
        self.hide_or_display_layer()

        self.set_opacity()

        if shared.VTK_first is False:
            widgets_list.widget_vtk.update_scalar_bar()

            # Update the rendering
            self.parent.interactor.GetRenderWindow().Render()

    def update_colors(self):
        """Update the colorscale of the tomography.

        Should only change colors, reupdates everything for the moment (slow).
        """
        self.update_tomography()

    def set_opacity(self, value=None):
        """Sets the opacity of the actors."""
        if value is None:
            value = self.opacity

        self.actor.GetProperty().SetOpacity(value / 100.0)

        # The profile's opacity also needs to be changed
        if self.actor_profiles:
            for actor in self.actor_profiles:
                actor.GetProperty().SetOpacity(value / 100.0)

        # Call parent's method set the axes opacity
        super().set_opacity(value)

    def update_profiles(self):
        """Hides or displays profiles on the tomography.

        Goes through the list of profiles and displays them on the surface.
        """
        if self.actor_profiles:
            # Remove the actors
            for actor in self.actor_profiles:
                self.parent.renderer.RemoveActor(actor)
            # Reset the profile actor list
            self.actor_profiles = []

        display = widgets_list.widget_vtk.display_profiles
        data = shared.exp.list[self.afm_id]

        if display and data.profile_list:
            sizes = {
                "x_size": data.x_size,
                "y_size": data.y_size}
            trans = {
                "trans_x": self.size_x / 2.0,
                "trans_y": self.size_y / 2.0,
                "trans_z": 0}

            for profile in data.profile_list:
                actor = vtk_afm.get_actor_profile(sizes,
                                                  trans,
                                                  self.z_scale,
                                                  profile,
                                                  self.poly_data_tomo)
                self.actor_profiles.append(actor)
                self.parent.renderer.AddActor(actor)

        # Update the rendering
        self.draw_at_position()


def point_to_cell_colors(polydata):
    """Polydata with the colors of its points moved to its triangles (mean of
    the three points of each triangle)."""
    colors = polydata.GetPointData().GetScalars()
    if colors is None or polydata.GetNumberOfCells() == 0:
        return polydata
    point_colors = numpy_support.vtk_to_numpy(colors).astype(float)
    # Connectivity of the triangles: 3 point ids per cell
    ids = numpy_support.vtk_to_numpy(
        polydata.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    cell_colors = numpy.round(point_colors[ids].mean(axis=1)).astype(
        numpy.uint8)
    array = numpy_support.numpy_to_vtk(cell_colors, deep=True)
    output = vtk.vtkPolyData()
    output.ShallowCopy(polydata)
    output.GetPointData().Initialize()
    output.GetCellData().SetScalars(array)
    return output
