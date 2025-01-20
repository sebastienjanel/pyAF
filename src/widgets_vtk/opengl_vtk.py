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

"""The 3D canvas."""

import vtk
import math
from .. import shared
from .. import consts
from PyQt5 import QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from .layers.layer_afm import LayerAFM
from .layers.layer_single_image import LayerSingleImage
from .layers.layer_optics_HR import LayerOpticsHR
from .layers.layer_glass import LayerGlass
from .. import widgets_list

if consts.ALLOW_ITK:
    from .layers.layer_stack import LayerStack
    from .layers.layer_isosurface import LayerIsoSurface


class VTKWidget(QtWidgets.QWidget):
    """QFrame containing the 3D scenery."""

    def __init__(self, load_list):
        super().__init__()
        policy = QtWidgets.QSizePolicy.Expanding
        self.setSizePolicy(policy, policy)

        VL = QtWidgets.QVBoxLayout()
        self.frame = QtWidgets.QFrame()

        VL_frame = QtWidgets.QVBoxLayout()
        VL_frame.setContentsMargins(0, 0, 0, 0)

        # See PYAFInteractor class defintion at end of file.
        # It's just a subclass of a QVTKRenderWindowInteractor.
        self.interactor = PYAFInteractor(self.frame)
        VL_frame.addWidget(self.interactor)

        self.renderer = vtk.vtkRenderer()
        self.interactor.GetRenderWindow().AddRenderer(self.renderer)

        self.style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor._Iren.SetInteractorStyle(self.style)
        self.style.SetCurrentRenderer(self.renderer)

        self.frame.setLayout(VL_frame)
        VL.addWidget(self.frame)
        self.setLayout(VL)

        # Add light
        self.light_kit = vtk.vtkLightKit()
        self.light_kit.SetKeyLightIntensity(1.0)
        self.light_kit.SetKeyLightWarmth(0.5)
        self.light_kit.SetFillLightWarmth(0.5)
        self.light_kit.SetHeadLightWarmth(0.5)
        self.light_kit.SetBackLightWarmth(0.5)
        self.light_kit.AddLightsToRenderer(self.renderer)

        # Sign up to receive TimerEvent
        self.anim = Animation(24, self)

        # Create the AFM layers
        for i in load_list:
            shared.layer_list.append(LayerAFM(self, i))
            shared.layer_list[-1].layer_id = len(shared.layer_list) - 1

        camera = self.renderer.GetActiveCamera()
        camera.SetViewUp(0, 0, 1)  # Set up vector to z axis
        camera.SetFocalPoint(0, 0, 0)

        size = 0
        # Get the biggest size to set the camera's position
        for layer in shared.layer_list:
            if layer.size_x > size:
                size = layer.size_x
            if layer.size_y > size:
                size = layer.size_y

        camera.SetPosition(-size * 2, -size * 2, 30000)

        # http://www.vtk.org/pipermail/vtkusers/2011-July/117968.html
        # Add renderer.ResetCameraClippingRange() after setting up the
        # camera.  The interactor calls this automatically for certain kinds of
        # interactions, but when you call Render() manually, you have to call
        # ResetCameraClippingRange() yourself when the object moves out
        # of the clipping range. It is called automatically for the first
        # render, but not for subsequent renders.

        self.renderer.ResetCameraClippingRange()

        self.picker = vtk.vtkPropPicker()
        self.interactor.SetPicker(self.picker)
        # Register the pick callback
        self.picker.AddObserver("PickEvent", self.pick)

    def sizeHint(self):
        """Redefine the size of the canvas, to force it's expansion.

        Use 4k resolution, to be sure that the canvas expands to it's
        maximum on any screen.
        """
        # pylint: disable=R0201
        # R0201 Method could be a function : it's not the case, we are
        # just redefining sizeHint here !

        return QtCore.QSize(3840, 2160)

    def pick(self, obj, event):
        """When clicking on a actor, draw a red border around the actor."""
        # Do nothing with these arguments :
        _ = obj
        _ = event

        menu = widgets_list.widget_vtk.menu_layer_list

        for layer in shared.layer_list:
            if layer.type == "afm":
                actor = self.picker.GetAssembly()
            else:
                actor = self.picker.GetActor()
            if layer.actor == actor:
                # We are in this layer, select it
                menu.tableWidget_clicked(layer.layer_id, None)
                break

    def add_more_afm_layers(self, load_list):
        """Once the user choosed the tomographies to load, add them to the 3D."""
        for result_id in load_list:
            shared.layer_list.append(LayerAFM(self, result_id))
            shared.layer_list[-1].update_tomography()
            shared.layer_list[-1].layer_id = len(shared.layer_list) - 1

        # Reset the list in the menu
        widgets_list.widget_vtk.menu_layer_list.reset_table()

    def addSingleImageLayer(self, path):
        """Adds an optics layer."""
        layer = LayerSingleImage(self, path)
        if layer.load_error is False:
            shared.layer_list.append(layer)
            if shared.VTK_first is False:
                self.renderer.Render()
            return True
        else:
            return False

    def addOpticsHRLayer(self, path):
        """Adds an high resolution optics layer"""
        shared.layer_list.append(LayerOpticsHR(self, path))
        if shared.VTK_first is False:
            self.renderer.Render()
        return True

    def addGlassLayer(self, layer_id, afm_id, roi_id):
        """Adds glass layers."""
        # Create a new glass layer
        shared.layer_list.append(LayerGlass(self, layer_id, afm_id, roi_id))

        if shared.VTK_first is False:
            self.renderer.Render()
        return True

    def addStackLayer(self, path):
        """Adds a stack layer (3D tiffs)"""
        layer = LayerStack(self, path)
        if not layer.load_error:
            shared.layer_list.append(layer)
            if not shared.VTK_first:
                self.renderer.Render()
            return True
        else:
            return False

    def addIsoSurfaceLayer(self, isosurf_stack_id=None, iso_file=None):
        """Add an isosurface."""
        index = widgets_list.widget_vtk.get_current_layer()

        if isosurf_stack_id is not None:
            index = isosurf_stack_id

        shared.layer_list.append(LayerIsoSurface(self, index, iso_file))
        if shared.VTK_first is False:
            self.renderer.Render()
        return True

    def TakeScreenshot(self, filename, return_img=False, magn=None):
        """Take a screenshot.

        A magnification can be used to increase the size of the resultion
        picture.

        The return_img argument can be used so that the method returns the
        image without writing it to the disk. This allows some post-treatment
        on the file.
        """
        # Either use the new magnification or the stored one
        if magn is None:
            magn = shared.exp.opengl_screenshot_magn
        else:
            shared.exp.opengl_screenshot_magn = magn

        # Save old camera position
        camera = self.renderer.GetActiveCamera()
        pos1 = camera.GetPosition()
        pos2 = camera.GetFocalPoint()
        pos3 = camera.GetViewAngle()
        pos4 = camera.GetViewUp()
        pos5 = camera.GetClippingRange()
        pos6 = camera.GetParallelScale()

        windowToImageFilter = vtk.vtkWindowToImageFilter()
        windowToImageFilter.SetInput(self.interactor.GetRenderWindow())
        windowToImageFilter.SetScale(magn, magn)
        windowToImageFilter.SetInputBufferTypeToRGB()
        windowToImageFilter.Update()

        # Write file
        if not return_img:
            writer = vtk.vtkPNGWriter()
            # Encode to utf-8 if there are some special characters in the path.
            # The writer needs an 8 bit encoding.
            writer.SetFileName(filename.encode("utf-8"))
            writer.SetInputConnection(windowToImageFilter.GetOutputPort())
            writer.Write()

        # Reset old position
        camera.SetPosition(pos1)
        camera.SetFocalPoint(pos2)
        camera.SetViewAngle(pos3)
        camera.SetViewUp(pos4)
        camera.SetClippingRange(pos5)
        camera.SetParallelScale(pos6)

        # Render
        if shared.VTK_first is False:
            self.interactor.GetRenderWindow().Render()

        # Return image for further usage
        if return_img:
            return windowToImageFilter

    def lock_top(self):
        """Sets the view above the actors.

        Used for doing screenshots or during the translation of the actors.
        """
        camera = self.renderer.GetActiveCamera()
        camera.SetPosition(0, 0, 50000)
        camera.SetFocalPoint(0, 0, 0)
        camera.SetViewUp(0, 1, 0)  # Set up vector to y axis
        self.renderer.ResetCamera()
        if shared.VTK_first is False:
            self.interactor.GetRenderWindow().Render()

    def update_rectangle_actor(self):
        """Build a rectangle and draw it."""
        if widgets_list.widget_vtk.rectangle_actor is not None:
            self.renderer.RemoveActor(widgets_list.widget_vtk.rectangle_actor)

        pos = widgets_list.widget_vtk.rectangle_pos

        positions = []
        positions.append([pos[0], pos[1], 0.0])
        positions.append([pos[2], pos[1], 0.0])
        positions.append([pos[2], pos[3], 0.0])
        positions.append([pos[0], pos[3], 0.0])

        points = vtk.vtkPoints()
        for pos in positions:
            points.InsertNextPoint(pos)

        lines = vtk.vtkCellArray()
        lines.InsertNextCell(5)
        lines.InsertCellPoint(0)
        lines.InsertCellPoint(1)
        lines.InsertCellPoint(2)
        lines.InsertCellPoint(3)
        lines.InsertCellPoint(0)

        polygon = vtk.vtkPolyData()
        polygon.SetPoints(points)
        polygon.SetLines(lines)

        polygonMapper = vtk.vtkPolyDataMapper()
        polygonMapper.SetInputData(polygon)

        widgets_list.widget_vtk.rectangle_actor = vtk.vtkActor()
        widgets_list.widget_vtk.rectangle_actor.SetMapper(polygonMapper)

        self.renderer.AddActor(widgets_list.widget_vtk.rectangle_actor)

        self.interactor.GetRenderWindow().Render()


class Animation:
    """Used to let the 3D scenery turn."""

    def __init__(self, fps, parent):
        self.fps = fps
        self.msec = int(1000 / fps)
        self.current_frame = 0
        self.start_angle = 0
        self.radius = 0
        self.parent = parent
        self.running = False

        # Define a timer
        self.timerObj = QtCore.QTimer()
        self.timerObj.timeout.connect(self.animate)

    def play(self):
        """Start the animation."""
        # Reset current frame
        self.current_frame = 0

        # Get the camera
        camera = self.parent.renderer.GetActiveCamera()

        # Get the distance between the camera and the center
        x = camera.GetPosition()[0]
        y = camera.GetPosition()[1]
        self.radius = math.sqrt(x ** 2 + y ** 2)

        # Get the angle between the x, y position and the 0, 0 vector
        # See : http://fr.wikipedia.org/wiki/Atan2
        # Atan2 is needed to get the right angle
        # (self.start_angle is in radians)
        self.start_angle = math.atan2(float(y), float(x))

        # Start the animation
        self.running = True
        self.timerObj.start(self.msec)

    def stop(self):
        """Stop the animation."""
        self.running = False
        self.timerObj.stop()

    def animate(self):
        """Method defining how the scenery turns.

        Will update the positon by one degree.
        """
        # Get the current camera
        camera = self.parent.renderer.GetActiveCamera()

        # Turn around at the current height
        z = camera.GetPosition()[2]

        # Get new x and y positions for the camera, take into account the
        # current angle (start_angle).
        delta = math.radians(self.current_frame)
        x = self.radius * math.cos(self.start_angle + delta)
        y = self.radius * math.sin(self.start_angle + delta)

        # Make sure the camera looks always up
        camera.SetViewUp(0, 0, 1)

        # Set the new camera position
        camera.SetPosition(x, y, z)

        # Render
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

        # Update the angle for the next animation step
        self.current_frame += 1


class PYAFInteractor(QVTKRenderWindowInteractor):
    """Subclass of QVTKRenderWindowInteractor, fetches mouse events.

    This class allows to fetch mouse movements and mouse clicks. Depending on
    the selected action mode, different thimgs are done, like translating or
    rotating the layers. If we are in "normal" camera mode, the event is
    re-sent to the parent QVTKRenderWindowInteractor class.
    """

    def __init__(self, frame):
        super().__init__(frame)

        self.started_drawing_square = False

    def mousePressEvent(self, event):
        """Fetch mouse press events."""
        if widgets_list.widget_vtk.action_mode == "translate":
            widgets_list.widget_vtk.last_pos = event.pos()

        elif widgets_list.widget_vtk.action_mode == "rotate":
            widgets_list.widget_vtk.last_pos = event.pos()

        elif widgets_list.widget_vtk.action_mode == "draw_square":
            # Reset pos
            widgets_list.widget_vtk.rectangle_pos = [0, 0, 0, 0]

            vals = self.get_world_pos(event.x(), event.y())

            widgets_list.widget_vtk.rectangle_pos[0] = vals[0]
            widgets_list.widget_vtk.rectangle_pos[1] = vals[1]
            widgets_list.widget_vtk.canvas.update_rectangle_actor()

            widgets_list.widget_vtk.last_pos = event.pos()
            self.started_drawing_square = True

        else:
            # Normal camera mode, send the event to the parent class.
            super().mousePressEvent(event)

    def mouseReleaseEvent(self, event):
        """Fetch mouse release events."""
        if widgets_list.widget_vtk.action_mode == "translate":
            widgets_list.widget_vtk.last_pos = event.pos()

        elif widgets_list.widget_vtk.action_mode == "draw_square":
            if self.started_drawing_square:
                widgets_list.widget_vtk.last_pos = event.pos()
                self.started_drawing_square = False
                widg = widgets_list.widget_vtk
                widg.action_mode = "camera"
                widg.menu_layer_list.BT_resample.setEnabled(True)

        else:
            # Normal camera mode, send the event to the parent class.
            super().mouseReleaseEvent(event)

    def mouseMoveEvent(self, event):
        """Fetch mouse move events."""
        # Chose the actor to be translated
        index = widgets_list.widget_vtk.get_current_layer()

        if widgets_list.widget_vtk.action_mode == "translate":
            if event.buttons() == QtCore.Qt.LeftButton:
                layer = shared.layer_list[index]

                canvas = widgets_list.widget_vtk.canvas
                camera = canvas.renderer.GetActiveCamera()
                campos = camera.GetPosition()

                # Factor campos (z)/1000.0 to make the movement's amplitude
                # depend on the cameras z position.
                x = (event.x() - widgets_list.widget_vtk.last_pos.x()) * \
                    (campos[2] / 1000.0)
                y = -(event.y() - widgets_list.widget_vtk.last_pos.y()) * \
                    (campos[2] / 1000.0)

                # Translate the whole layer (no z movement)
                layer.translate(x, y, 0)

        elif widgets_list.widget_vtk.action_mode == "rotate":
            if event.buttons() == QtCore.Qt.LeftButton:
                # Rotate the whole layer
                shared.layer_list[index].rotate(
                    event,
                    widgets_list.widget_vtk.last_pos)

        elif widgets_list.widget_vtk.action_mode == "draw_square":
            if self.started_drawing_square:
                vals = self.get_world_pos(event.x(), event.y())

                widgets_list.widget_vtk.rectangle_pos[2] = vals[0]
                widgets_list.widget_vtk.rectangle_pos[3] = vals[1]

                widgets_list.widget_vtk.canvas.update_rectangle_actor()

        else:
            # Normal camera mode, send the event to the parent class.
            super().mouseMoveEvent(event)

        # Update last position
        widgets_list.widget_vtk.last_pos = event.pos()

    def get_world_pos(self, x, y):
        """Converts mouse position in world coordinates.

        z = 0 near plane
        z = 1 far plane
        Perhaps z should be calculated but Z=0.5 seems to work.
        There is still a small offset depending on the zooming factor ...
        """
        z = 0.5

        size = self.GetSize()
        canvas = widgets_list.widget_vtk.canvas
        canvas.renderer.SetDisplayPoint([x, size[1] - y, z])
        canvas.renderer.DisplayToWorld()
        vals = canvas.renderer.GetWorldPoint()

        return vals
