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

"""High resolution optics layer."""

from .layer import Layer
from ... import shared
import vtk


class LayerOpticsHR(Layer):
    """High resolution optics layer."""

    def __init__(self, parent, path):
        super().__init__()

        self.parent = parent
        self.layer_id = len(shared.layer_list)
        self.type = "HR"
        self.filename = path.split("/")[-1]

        self.original_size_x = 52100
        self.original_size_y = 52100
        self.size_x = self.original_size_x
        self.size_y = self.original_size_y
        self.xdist = 0
        self.ydist = 0
        self.color = [0, 1, 0, 1]

        self.palm_resolution = 70.0

        self.palm_data = []
        afile = open(path, "rb")
        count = 0
        for line in afile:
            if line.split("\x00")[0] == "":
                # Found last line
                break
            if count != 0 and line.split("\x00")[0] != "":
                aline = line.split("\t")
                val = [float(aline[5]), float(aline[6]), float(aline[7])]
                self.palm_data.append(val)
            count += 1

        afile.close()

        self.actor = None

        self.create_new_actor()
        self.init_axes()

        self.addToRenderer()

    def addToRenderer(self):
        """Add actor to renderer."""
        self.parent.renderer.AddActor(self.actor)

    def create_new_actor(self):
        """Create new actor."""
        append_polydata = vtk.vtkAppendPolyData()

        for j in range(len(self.palm_data)):
            line = self.palm_data[j]
            res = float(line[2])
            if res < self.palm_resolution and res > 20.0:
                polygonSource = vtk.vtkRegularPolygonSource()
                polygonSource.SetNumberOfSides(10)
                polygonSource.SetRadius(res / 2.0)
                ct = [line[0] - self.original_size_x / 2.0,
                      line[1] - self.original_size_y / 2.0, 0.0]
                polygonSource.SetCenter(ct)
                sc = polygonSource.GetOutputPort()
                append_polydata.AddInputConnection(sc)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(append_polydata.GetOutputPort())

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)
        col = [self.color[0], self.color[1], self.color[2]]
        self.actor.GetProperty().SetColor(col)

    def remove_actors(self):
        """Remove the actor from the layer."""
        self.parent.renderer.RemoveActor(self.actor)

    def update_actor(self):
        """Update the actor."""
        self.remove_actors()
        self.create_new_actor()

        self.init_axes()
        self.addToRenderer()

        self.rotate(None, None, None)
        self.translate(0, 0, 0)

        self.hide_or_display_layer()

        self.set_opacity(self.opacity)

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
