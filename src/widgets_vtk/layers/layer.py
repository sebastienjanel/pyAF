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

"""Main vtk layer."""

import vtk
from ... import shared
from ... import widgets_list


class Layer:
    """A simple wrapper for the layers.

    The layer object is a generic layer, which is the parent class for
    different types of layers (AFM, optics, ...). It contains the common
    methods and elements for the different types of layers, for example the
    scales, the methods for translation and rotation, and so on.
    """

    def __init__(self):
        super().__init__()

        self.afm_id = None
        self.layer_id = None
        self.filename = None
        self.actor_axes = None
        self.pos_x = 0
        self.pos_y = 0
        self.pos_z = 0
        self.angle = 0
        self.visible = True
        self.axes_visible = True
        self.opacity = 100
        self.fiducials = Fiducials(self)

        self.load_error = False

        self.poly_data_tomo = None

        self.nbr_pixels_x = None
        self.nbr_pixels_y = None
        self.nbr_pixels_z = None
        self.size_x = None
        self.size_y = None
        self.size_z = None
        self.res_x = None
        self.res_y = None
        self.res_z = None
        self.original_nbr_pixels_x = None
        self.original_nbr_pixels_y = None
        self.original_nbr_pixels_z = None
        self.original_res_x = None
        self.original_res_y = None
        self.original_res_z = None
        self.original_size_x = None
        self.original_size_y = None
        self.original_size_z = None
        self.z_scale = None
        self.bind = None

        self.glass_color = [207 / 255.0, 207 / 255.0, 207 / 255.0, 1.0]

        # For the layer optics mainly
        self.filepath = None
        self.use_lut = False
        self.lut_threshold = 0
        self.lut_invert = False
        self.lut_min_color = [0.0, 0.0, 0.0, 1.0]  # Black RGBA
        self.lut_max_color = [1.0, 1.0, 1.0, 1.0]  # White RGBA

        # Thresholding
        self.use_binary_thresholding = False
        self.lower_thresh = 50
        self.upper_thresh = 255
        self.lower_thresh_value = 0
        self.upper_thresh_value = 255

        # For the stack layers
        self.slide_pos = 0
        self.single_opacity = 100
        self.others_opacity = 0

        # Flip the images
        self.flip_x = False
        self.flip_y = False

        # Display or hide the z scale
        self.display_z_scale = True

        # Smoothing (0 = No smoothing)
        self.smooth_gaussian = 0

        # Isosurface
        self.isosurface_color = [1.0, 0.0, 0.0, 1.0]
        self.isosurf_orig_layer_id = None

    @property
    def sizes(self):
        """Getter for the sizes of the actor."""
        return [self.size_x, self.size_y, self.size_z]

    @property
    def res(self):
        """Getter for the resolutions of the actor."""
        return [self.res_x, self.res_y, self.res_z]

    def init_axes(self):
        """Creates a vtkCubeAxesActor."""
        # The first time there is no self.actor_axes
        if self.actor_axes is not None:
            self.parent.renderer.RemoveActor(self.actor_axes)

        self.actor_axes = vtk.vtkCubeAxesActor()

        # Position of the axes, x and y will be in the front
        # z on the left or the right
        self.actor_axes.SetFlyModeToOuterEdges()

        # The the size of the axes
        self.set_axes_bounds()

        val_x1 = 0
        val_x2 = self.size_x
        val_y1 = 0
        val_y2 = self.size_y

        self.actor_axes.SetXAxisRange(val_x1 / 1000.0, val_x2 / 1000.0)
        self.actor_axes.SetYAxisRange(val_y1 / 1000.0, val_y2 / 1000.0)

        # Use z axis only for AFM layers (for the moment)
        if self.type != "afm":
            self.actor_axes.SetZAxisVisibility(False)
        else:
            # Disable z axis if asked
            if self.display_z_scale:
                self.actor_axes.SetZAxisRange(0, self.data.max_topo / 1000.0)
            else:
                self.actor_axes.SetZAxisRange(0, 0)
            self.actor_axes.SetZAxisVisibility(self.display_z_scale)
            self.actor_axes.SetZAxisLabelVisibility(self.display_z_scale)

        # No minor ticks (could be made optional)
        self.actor_axes.XAxisMinorTickVisibilityOff()
        self.actor_axes.YAxisMinorTickVisibilityOff()
        if self.type == "afm":
            self.actor_axes.ZAxisMinorTickVisibilityOff()

        self.set_axes_colors()

        # Add the axes to the renderer
        self.parent.renderer.AddActor(self.actor_axes)

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def set_axes_bounds(self):
        """Sets the size of the axes in nm.

        For tomographies this is special. If we have 16 values to plot in
        one line of 10 um, on a meshgrid you would plot 16 squares of
        10000/16 = 625 nm border. Plotting it with triangles in 3D is
        different, you plot only 15 squares, because the values are set on the
        border, and there is a value missing to have a line of 10 um. In fact,
        here you plot only a line of 9375 nm. So we shift the actor_axes by
        size/2 to have it still display 10 um, but at the right place.
        This is done only for surfaces made of triangles.
        """
        if self.type == "afm":
            # Full bounds
            x1 = - self.size_x / 2.0 - self.res_x / 2.0 + self.pos_x
            x2 = self.size_x / 2.0 - self.res_x / 2.0 + self.pos_x
            y1 = - self.size_y / 2.0 - self.res_y / 2.0 + self.pos_y
            y2 = self.size_y / 2.0 - self.res_y / 2.0 + self.pos_y

            # Slice the bounds
            x1 = x1 + self.data.opengl_slice_bottom * self.res_x
            diffx = self.nbr_pixels_x - self.data.opengl_slice_top - 1
            x2 = x2 - (diffx) * self.res_x
            y1 = y1 + self.data.opengl_slice_left * self.res_y
            diffy = self.nbr_pixels_y - self.data.opengl_slice_right - 1
            y2 = y2 - (diffy) * self.res_y

            if self.z_scale is not None:
                z1 = self.z_scale[0] + self.pos_z
                z2 = self.z_scale[1] + self.pos_z
            else:
                # During initialisation we have no values
                z1 = 0
                z2 = 0

        else:
            x1 = - self.size_x / 2.0 + self.pos_x
            x2 = self.size_x / 2.0 + self.pos_x
            y1 = - self.size_y / 2.0 + self.pos_y
            y2 = self.size_y / 2.0 + self.pos_y
            z1 = self.pos_z
            z2 = self.pos_z

        # self.actor_axes.SetUseAxisOrigin(1)
        # self.actor_axes.SetAxisOrigin(x2, y2, z1)
        self.actor_axes.SetBounds(x1, x2, y1, y2, z1, z2)

        self.actor_axes.SetCamera(self.parent.renderer.GetActiveCamera())

    def set_axes_colors(self):
        """Sets the color of the axes, axe's labels and titles."""
        # Get the color
        r = widgets_list.widget_vtk.scales_color[0]
        g = widgets_list.widget_vtk.scales_color[1]
        b = widgets_list.widget_vtk.scales_color[2]

        # Set the color of the lines
        self.actor_axes.GetXAxesLinesProperty().SetColor(r, g, b)
        self.actor_axes.GetYAxesLinesProperty().SetColor(r, g, b)
        self.actor_axes.GetZAxesLinesProperty().SetColor(r, g, b)
        # Set the color of the labels
        self.actor_axes.GetLabelTextProperty(0).SetColor(r, g, b)
        self.actor_axes.GetLabelTextProperty(1).SetColor(r, g, b)
        self.actor_axes.GetLabelTextProperty(2).SetColor(r, g, b)
        # Set the color of the titles
        self.actor_axes.GetTitleTextProperty(0).SetColor(r, g, b)
        self.actor_axes.GetTitleTextProperty(1).SetColor(r, g, b)
        self.actor_axes.GetTitleTextProperty(2).SetColor(r, g, b)

    def updateAxesColors(self):
        """Changes the color of the axes."""
        # Set the colors
        self.set_axes_colors()

        # Needed to refresh title color, force a refresh of the titles
        # (Waiting for patch in vtk library, normally not needed)
        self.actor_axes.SetXTitle("")
        self.actor_axes.SetXTitle("X-Axis")

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def translate(self, x, y, z):
        """Method used to translate the layer in the z plane."""
        # Store the position of the layer
        self.pos_x = self.pos_x + x
        self.pos_y = self.pos_y + y
        self.pos_z = self.pos_z + z

        # Translate the actor
        self.actor.SetPosition(self.pos_x, self.pos_y, self.pos_z)

        # For AFM layers, we can have profiles to move
        if self.type == "afm":
            if self.actor_profiles:
                for actor in self.actor_profiles:
                    actor.SetPosition(self.pos_x, self.pos_y, self.pos_z)

        # Update the values in the inputs
        widgets_list.widget_vtk.menu_options.tab_options.update_widget()

        # Translate the axes
        self.set_axes_bounds()

        # Translate the fiducials
        self.fiducials.update_fiducials()

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def rotate(self, cur_pos, last_pos, increment=None, new_angle=None,
               new_actor=False):
        """
        Method used to rotate the layer around z axis.

        Could be simplified / cleaned up ?

        """

        if cur_pos is not None:
            # Angle increment for rotation with the mouse
            new_angle = cur_pos.x() - last_pos.x()
            self.angle = self.angle + new_angle
        else:
            if increment is not None:
                # Add an angle increment (when updating angle from the inputs
                # in the menu)
                new_angle = increment
                self.angle = self.angle + increment
            elif new_angle is not None:
                # Rotate back
                if not new_actor:
                    self.actor.RotateZ(-self.angle)

                self.angle = new_angle

            else:
                new_angle = 0

        # Rotate
        self.actor.RotateZ(new_angle)

        # For AFM layers, we can have profiles to rotate
        if self.type == "afm":
            if self.actor_profiles:
                for actor in self.actor_profiles:
                    actor.RotateZ(new_angle)

        # Update the values in the inputs
        widgets_list.widget_vtk.menu_options.tab_options.update_widget()

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def hide_or_display_layer(self):
        """Sets the visibility to True or False."""
        self.actor.SetVisibility(self.visible)

        value = self.axes_visible
        if self.visible is False:
            value = False
        self.actor_axes.SetVisibility(value)
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def set_opacity(self, value):
        """Sets the axe's opacity. (Not working for the moment)

        To change the opacity of the actors belonging to this layer, you have
        to add a set_opacity method to your specific layer class.
        """
        # Save the opacity
        self.opacity = value

        # Needed ?
        self.hide_or_display_layer()

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def delete(self):
        """Can be called to delete the layer.

        The actors will be removed from the renderer.
        """
        self.parent.renderer.RemoveActor(self.actor_axes)
        self.parent.renderer.RemoveActor(self.actor)

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def get_current_lut_options(self):
        """Returns a new lut dictionnary with the current options."""
        lut = {
            "use": self.use_lut,
            "invert": self.lut_invert,
            "threshold": self.lut_threshold,
            "min_color": self.lut_min_color,
            "max_color": self.lut_max_color}

        return lut

    def get_current_flip_options(self):
        """Returns a dictionnary with the current flip options."""
        flip_options = {"flip_x": self.flip_x, "flip_y": self.flip_y}

        return flip_options


class Fiducials:
    """Object containing the fiducials and some options."""

    def __init__(self, parent):
        super().__init__()

        self.parent = parent
        self.positions = []
        self.color = [1.0, 0.0, 0.0]
        self.display_list = []
        self.actors_list = []
        self.radius = 200  # nm
        self.current_fiducial = 0  # The currently selected fiducial

    def add_to_fiducials(self, positions, display):
        """Add new fiducials."""
        self.positions.extend(positions)
        self.display_list.extend(display)
        self.update_fiducials()

    def add_empty_fiducial(self):
        """Add a signle empty fiducial."""
        self.positions.append([0.0, 0.0, 0.0])
        self.display_list.append(True)
        self.update_fiducials()

    def remove_fiducial(self, i):
        """Remove a fiducial."""
        del self.positions[i]
        del self.display_list[i]
        actor = self.actors_list[i]
        self.parent.parent.renderer.RemoveActor(actor)
        del self.actors_list[i]
        self.update_fiducials()

    def update_fiducials(self):
        """Draws all the fiducials."""
        # Remove old actors
        for actor in self.actors_list:
            self.parent.parent.renderer.RemoveActor(actor)
        self.actors_list = []

        for i in range(len(self.positions)):
            x, y, z = self.positions[i]

            # Create source
            source = vtk.vtkSphereSource()
            # The positions are relative to the layer positions
            x = x + self.parent.pos_x
            y = y + self.parent.pos_y
            z = z + self.parent.pos_z
            source.SetCenter(x, y, z)
            source.SetRadius(self.radius)
            # Default resolutions is 8, 20 is to make it look like a real
            # sphere
            source.SetThetaResolution(20)
            source.SetPhiResolution(20)

            # Mapper
            mapper = vtk.vtkPolyDataMapper()
            if vtk.VTK_MAJOR_VERSION <= 5:
                mapper.SetInput(source.GetOutput())
            else:
                mapper.SetInputConnection(source.GetOutputPort())
                mapper.Update()

            # Actor
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            if i == self.current_fiducial:
                # Use a white fiducial for the current one
                actor.GetProperty().SetColor([1.0, 1.0, 1.0])
            else:
                # Set the defined color for the other fiducials
                actor.GetProperty().SetColor(self.color)

            self.actors_list.append(actor)

            self.parent.parent.renderer.AddActor(actor)

        # Hide or display (will also render)
        self.hide_or_display()

    def hide_or_display(self):
        """Hides or displays a fiducial."""
        for i in range(len(self.actors_list)):
            actor = self.actors_list[i]
            actor.SetVisibility(self.display_list[i])

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.parent.interactor.GetRenderWindow().Render()

    def set_current_fiducial(self, i):
        """Changes the color of the current fiducial"""
        self.current_fiducial = i
        self.update_fiducials()

    def center_fiducials(self):
        """Center the fiducials around the parents actor."""
        for position in self.positions:
            position[0] = position[0] - self.parent.size_x / 2
            position[1] = position[1] - self.parent.size_y / 2

        self.update_fiducials()
