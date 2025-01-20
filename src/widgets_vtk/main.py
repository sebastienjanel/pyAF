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

"""Main VTK widget."""

from PyQt5 import QtCore, QtWidgets
from .menu_options import MenuOptionsWidget
from .menu_layer_list import MenuLayerListWidget
from .opengl_vtk import VTKWidget
from .utils.load import load_layers
from ..tools.PYAFWidget import PYAFWidget
from .. import widgets_list
from .. import shared
import vtk
from ..tools.colortables import ColorTables
from .. import consts


class OpenglWidget(PYAFWidget):
    """The OpenglWidget displays the 3D scenery and a menu.

    No rendering is allow before the self.show(). But Render() calls are needed
    in the methods which can be called by the GUI. During the creation of the
    VTK widget, these Render() calls are not triggered with the help of the
    shared.VTK_first which is set to True by default. Once show() has been
    called shared.VTK_first is set to False to allow the Render() calls.

    Note : the Render() calls are mandatory, else if you change something in
    the GUI the rendering will only update when you move your mouse over the
    scenery ... and we want it to update immediately.
    """

    def __init__(self, parent, load_list):
        super().__init__(parent, "widget_vtk")

        self.setGeometry(0, 0, 1200, 700)

        shared.layer_list = []
        self.action_mode = "camera"
        self.last_pos = QtCore.QPoint(0, 0)
        self.rectangle_pos = [0, 0, 0, 0]  # x1, y1, x2, y2
        self.rectangle_actor = None
        self.scales_color = [1.0, 1.0, 1.0, 1.0]
        self.display_profiles = False

        HL = QtWidgets.QHBoxLayout()

        self.canvas = VTKWidget(load_list)  # Layers are created here

        self.menu_layer_list = MenuLayerListWidget(self, load_list)
        self.menu_options = MenuOptionsWidget()
        self.menu_layer_list.update_widget(index=0)

        HL.addWidget(self.menu_layer_list)
        HL.addWidget(self.canvas)
        HL.addWidget(self.menu_options)

        self.setLayout(HL)

        for i in range(len(shared.layer_list)):
            if shared.layer_list[i].type == "afm":
                shared.layer_list[i].update_tomography()
            shared.layer_list[i].hide_or_display_layer()

        # Pick first layer
        self.canvas.style.HighlightProp(shared.layer_list[0].actor)

        self.menu_options.tab_afm.update_widget()

        # Add colorbar to the renderer
        self.scalar_bar = get_color_bar()
        self.canvas.renderer.AddActor2D(self.scalar_bar)

        # Automatically load saved layers, is used for debugging
        if consts.AUTO_LOAD_LAYER_IN_3D:
            path = consts.TEST_PATH + "layer"
            load_layers(self, path)

        self.show()
        self.canvas.interactor.Initialize()
        # Allow calls to Render() now that show() was called
        shared.VTK_first = False

    def get_current_layer(self):
        """Get the index of the currently selected layer."""
        index = None
        rows = self.menu_layer_list.tableWidget.selectionModel().selectedRows()
        for item in rows:
            index = item.row()

        return index

    def closeEvent(self, event):
        """Cleanup when closing the 3D widget.

        Close the color options widget (to be sure it's closed)
        Close also the indentation widget, then close OpenglWidget
        Reset also VTK_first to True for the next opening
        """
        shared.VTK_first = True
        widgets_list.widget_meshgrid_options.close()
        widgets_list.widget_indentation.close()
        super().closeEvent(event)

    def update_AFM_colors(self):
        """Updates the colorscale of the AFM tomographies."""
        if shared.exp.apply_to_all_data is False:
            tw = self.menu_layer_list.tableWidget
            for item in tw.selectionModel().selectedRows():
                index = item.row()
            shared.layer_list[index].update_colors()
        else:
            # Update all
            for actor in shared.layer_list:
                # Only update the colors if it's an AFM layer
                if actor.type == "afm":
                    actor.update_colors()
        if shared.VTK_first is False:
            self.canvas.interactor.GetRenderWindow().Render()

    def update_indentation(self):
        """Changes the indentation depth of the AFM tomographies."""
        if shared.exp.apply_to_all_data is False:
            tw = self.menu_layer_list.tableWidget
            for item in tw.selectionModel().selectedRows():
                index = item.row()
            shared.layer_list[index].update_tomography()
        else:
            # Update all
            for actor in shared.layer_list:
                # Only update the indentation if it's an AFM layer
                if actor.type == "afm":
                    actor.update_tomography()
        if shared.VTK_first is False:
            self.canvas.interactor.GetRenderWindow().Render()

    def update_scalar_bar(self):
        """Update the scalar bar."""
        # Remove old scalar_bar
        self.canvas.renderer.RemoveActor(self.scalar_bar)

        # Get a new scalar bar and add it to the renderer
        self.scalar_bar = get_color_bar()
        self.canvas.renderer.AddActor2D(self.scalar_bar)

        # Render
        if shared.VTK_first is False:
            self.canvas.interactor.GetRenderWindow().Render()


def get_color_bar():
    """Create a lookup table for the colorbar

    (Defined with the values of the first layer)
    """
    data = shared.exp.list[0]

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfColors(4)
    lut.SetNumberOfTableValues(512)
    tp = data.opengl_surf_type
    if tp == "Stiffness" or tp == "Stiffness_squares":
        themax = data.colortable_max_value
        lut.SetTableRange(0, themax / 1000.0)
    elif tp == "Topo_dots" or tp == "Topo_surf":
        themax = data.colortable_max_value
        lut.SetTableRange(0, themax)
    colortable = ColorTables(data.colortableid,
                             data.color_saturation,
                             data.color_negative,
                             data.colortable_max_value,
                             data.colortable_middle_value,
                             mode="scalarMap",
                             color_nan=data.color_nan)
    step = themax / 512.0

    for i in range(512):
        color_st = colortable.get_color_as_list(i * step)
        lut.SetTableValue(i, [color_st[0], color_st[1], color_st[2], 1.0])
    lut.Build()

    # Create a scalar bar
    scalar_bar = vtk.vtkScalarBarActor()
    # Set the LUT
    scalar_bar.SetLookupTable(lut)

    # Some cosmetic options
    scalar_bar.SetTextPosition(0)  # Labels on the left
    scalar_bar.GetLabelTextProperty().SetFontSize(10)
    scalar_bar.SetLabelFormat("%.0f")
    scalar_bar.SetPosition(0.9, 0.1)
    scalar_bar.SetOrientationToVertical()
    scalar_bar.SetWidth(0.07)
    scalar_bar.SetHeight(0.9)
    scalar_bar.GetTitleTextProperty().ItalicOff()
    scalar_bar.GetLabelTextProperty().ItalicOff()

    # Add the title
    if tp == "Stiffness":
        scalar_bar.SetTitle("E (kPa)")
    elif tp == "Topo_dots" or tp == "Topo_surf":
        scalar_bar.SetTitle("Height \n (nm)  ")

    # Default value is 64, we want a smooth scalarbar
    scalar_bar.SetMaximumNumberOfColors(512)

    # In case it should be hidden, just change the visibility
    if data.opengl_display_colorbar is False:
        scalar_bar.SetVisibility(False)

    return scalar_bar
