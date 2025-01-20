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

"""Layer (substrate)"""

import vtk
from ... import shared
from .layer import Layer


class LayerGlass(Layer):
    """Layer for the glass."""

    def __init__(self, parent, layer_id, afm_id, roi_id):
        super().__init__()

        self.parent = parent
        self.layer_id = layer_id
        self.afm_id = afm_id
        self.roi_id = roi_id
        self.type = "glass"
        self.axes_visible = True
        self.bind = True

        self.size_x = 20000.0
        self.size_y = 20000.0
        self.original_size_x = 20000.0
        self.original_size_y = 20000.0
        self.xdist = 0.0
        self.ydist = 0.0

        self.actor = self.get_actor_glass()

        self.init_axes()
        self.addToRenderer()

    def addToRenderer(self):
        """Add the actor to the renderer."""
        self.parent.renderer.AddActor(self.actor)

    def get_actor_glass(self):
        """Create the glass actor."""
        data = shared.exp.list[self.afm_id]

        d, a, b, _ = data.roi_list[self.roi_id].glass_coeffs
        h1 = d
        h2 = a * self.size_x + d
        h3 = a * self.size_x + b * self.size_x + d
        h4 = b * self.size_x + d

        # Setup four points
        points = vtk.vtkPoints()
        points.InsertNextPoint(-self.size_x / 2, -self.size_y / 2, h1)
        points.InsertNextPoint(self.size_x / 2, -self.size_y / 2, h4)
        points.InsertNextPoint(self.size_x / 2, self.size_y / 2, h3)
        points.InsertNextPoint(-self.size_x / 2, self.size_y / 2, h2)

        # Create the polygon
        polygon = vtk.vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(4)
        polygon.GetPointIds().SetId(0, 0)
        polygon.GetPointIds().SetId(1, 1)
        polygon.GetPointIds().SetId(2, 2)
        polygon.GetPointIds().SetId(3, 3)

        # Add the polygon to a list of polygons
        polygons = vtk.vtkCellArray()
        polygons.InsertNextCell(polygon)

        # Create a PolyData
        polygonPolyData = vtk.vtkPolyData()
        polygonPolyData.SetPoints(points)
        polygonPolyData.SetPolys(polygons)

        # Create a mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.SetInput(polygonPolyData)
        else:
            mapper.SetInputData(polygonPolyData)
            mapper.Update()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        return actor

    def remove_actors(self):
        """Remove the glass actor from the scenery."""
        self.parent.renderer.RemoveActor(self.actor)

    def update_actor(self):
        """Update the actor."""
        self.remove_actors()
        self.actor = self.get_actor_glass()

        self.init_axes()
        self.addToRenderer()

        self.rotate(None, None, new_angle=self.angle, new_actor=True)
        self.translate(0, 0, 0)

        self.hide_or_display_layer()

        self.set_opacity()

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def set_opacity(self, value=None):
        """Sets the opacity of the actors."""
        if value is None:
            value = self.opacity

        self.actor.GetProperty().SetOpacity(value / 100.0)

        # Call parent's method set the axes opacity
        super().set_opacity(value)
