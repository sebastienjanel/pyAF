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

from qtpy import QtCore, QtWidgets
from .menu_options import MenuOptionsWidget
from .menu_layer_list import MenuLayerListWidget
from .opengl_vtk import VTKWidget
from .utils.load import load_layers
from ..tools.PYAFWidget import PYAFWidget
from .. import widgets_list
from .. import shared
import vtk
from ..tools import map_values
from .utils import vtk_afm
from .utils.common import use_app_font
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
        self.scales_color = [1.0, 1.0, 1.0, 1.0]
        self.display_profiles = False

        HL = QtWidgets.QHBoxLayout()

        self.canvas = VTKWidget(load_list)  # Layers are created here

        self.menu_layer_list = MenuLayerListWidget(self, load_list)
        self.menu_options = MenuOptionsWidget()
        self.menu_layer_list.add_position_box(self.menu_options.tab_options)
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

        # Background of the experiment (black or white)
        color = shared.exp.opengl_background_color
        self.canvas.renderer.SetBackground(color[0], color[1], color[2])

        # Add colorbar to the renderer
        self.scalar_bar = get_color_bar()
        self.canvas.renderer.AddViewProp(self.scalar_bar)

        # Automatically load saved layers, is used for debugging
        if consts.AUTO_LOAD_LAYER_IN_3D:
            path = consts.TEST_PATH + "layer"
            load_layers(self, path)

        self.show()
        self.canvas.interactor.Initialize()
        # Allow calls to Render() now that show() was called
        shared.VTK_first = False

    def get_current_layer(self):
        """Get the index of the currently selected layer (the current row,
        or the first one, if none is selected)."""
        table = self.menu_layer_list.tableWidget
        index = None
        rows = table.selectionModel().selectedRows()
        for item in rows:
            index = item.row()

        if index is None and table.rowCount() > 0:
            index = max(table.currentRow(), 0)
        return index

    def closeEvent(self, event):
        """Cleanup when closing the 3D widget.

        Close the color options widget (to be sure it's closed)
        Close also the indentation widget, then close OpenglWidget
        Reset also VTK_first to True for the next opening
        """
        shared.VTK_first = True
        for widget in (widgets_list.widget_meshgrid_options,
                       widgets_list.widget_indentation):
            if widget is not None:
                widget.close()
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
        self.canvas.renderer.AddViewProp(self.scalar_bar)

        # Render
        if shared.VTK_first is False:
            self.canvas.interactor.GetRenderWindow().Render()


def colorbar_data():
    """File whose color scale is displayed: the selected AFM layer, else the
    first AFM layer of the window."""
    afm_layers = [layer for layer in shared.layer_list if layer.type == "afm"]
    widget = widgets_list.widget_vtk
    if widget is not None:
        index = widget.get_current_layer()
        if index is not None and shared.layer_list[index].type == "afm":
            return shared.exp.list[shared.layer_list[index].afm_id]
    if afm_layers:
        return shared.exp.list[afm_layers[0].afm_id]
    return shared.exp.list[0]


def get_color_bar():
    """Create a lookup table for the colorbar

    (Defined with the values of the selected AFM layer, or of the first AFM
    layer of the window)
    """
    data = colorbar_data()
    mt = vtk_afm.map_type(data)
    opts = data.get_color_opts(mt)
    fac = map_values.factor(data, mt)

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfColors(4)
    lut.SetNumberOfTableValues(512)
    # Same range as the colors of the surface (color options of the map):
    # from the min value, except for Rainbow 1 and Cool-Warm which start at 0
    # (as in the 2D meshgrids)
    themin = 0 if opts[0] in (1, 6) else opts[1]
    themax = opts[2]
    if themax <= themin:
        themax = themin + 1
    lut.SetTableRange(themin * fac, themax * fac)
    colortable = vtk_afm.get_colortable(data, mt)
    step = (themax - themin) / 512.0

    for i in range(512):
        color_st = colortable.get_color_as_list(themin + i * step)
        lut.SetTableValue(i, [color_st[0], color_st[1], color_st[2], 1.0])
    lut.Build()

    # Create a scalar bar
    scalar_bar = vtk.vtkScalarBarActor()
    # Set the LUT
    scalar_bar.SetLookupTable(lut)

    # Some cosmetic options
    scalar_bar.SetTextPosition(0)  # Labels on the left
    scalar_bar.GetLabelTextProperty().SetFontSize(10)
    scalar_bar.SetLabelFormat("%.3g")
    scalar_bar.SetPosition(0.9, 0.1)
    scalar_bar.SetOrientationToVertical()
    scalar_bar.SetWidth(0.07)
    scalar_bar.SetHeight(0.9)
    scalar_bar.GetTitleTextProperty().ItalicOff()
    scalar_bar.GetLabelTextProperty().ItalicOff()
    use_app_font(scalar_bar.GetTitleTextProperty())
    use_app_font(scalar_bar.GetLabelTextProperty())
    # Text in the color of the scales (black or white)
    widget = widgets_list.widget_vtk
    color = widget.scales_color[:3] if widget is not None else [1, 1, 1]
    scalar_bar.GetTitleTextProperty().SetColor(*color)
    scalar_bar.GetLabelTextProperty().SetColor(*color)

    # Title: label of the map
    scalar_bar.SetTitle(map_values.label(data, mt) + "  ")

    # Default value is 64, we want a smooth scalarbar
    scalar_bar.SetMaximumNumberOfColors(512)

    # In case it should be hidden, just change the visibility
    if data.opengl_display_colorbar is False:
        scalar_bar.SetVisibility(False)

    return scalar_bar
