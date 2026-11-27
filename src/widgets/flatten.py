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

"""Widget used to flatten the piezo image."""

import numpy
from qtpy import QtWidgets
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFToggle
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from . import roi_manager
from .. import shared


class FlattenWidget(PYAFWidget):
    """Allows to flatten the piezo image with order 0, 1, 2."""

    def __init__(self, parent):
        super().__init__(parent, "widget_flatten")

        # Order None: no flatten (only the z motor steps levelled)
        if shared.exp.current_data.flatten_applied:
            flatten_order = shared.exp.current_data.applied_flatten_order
        else:
            flatten_order = 1
        shared.exp.selected_flatten_order = flatten_order
        if shared.exp.current_data.flatten_applied:
            shared.exp.selected_flatten_motor = \
                shared.exp.current_data.applied_flatten_motor

        self.box_flatten = QtWidgets.QGroupBox("Flatten")
        self.VL_box_flatten = QtWidgets.QVBoxLayout()

        # The maps fill the free space of the window
        self.HL_flatten_canvas = QtWidgets.QHBoxLayout()
        self.W_original = self.canvas_widget()
        mesh_type = "small_meshgrid_for_correction_original"
        self.MPL_canvas_original = \
            PYAFPlot(self, mesh_type, self.W_original, [3, 3, 72])
        self.W_original.layout().addWidget(self.MPL_canvas_original)
        self.canvas_flatten_widget = self.canvas_widget()
        mesh_type = "small_meshgrid_for_correction"
        self.MPL_canvas_flatten = \
            PYAFPlot(self, mesh_type, self.canvas_flatten_widget, [3, 3, 72])
        self.canvas_flatten_widget.layout().addWidget(self.MPL_canvas_flatten)
        if shared.exp.current_data.is_single is False:
            self.update_MPL("MPL_canvas_flatten")

        self.HL_flatten_canvas.addWidget(self.W_original)
        self.HL_flatten_canvas.addWidget(self.canvas_flatten_widget)

        # Rectangles left out of the fits (high features), drawn on the
        # original map
        if not shared.exp.current_data.flatten_applied:
            shared.exp.current_data.flatten_excluded = []
        else:
            shared.exp.current_data.flatten_excluded = list(
                shared.exp.current_data.applied_flatten_excluded)
        self.drag_start = None
        for event, method in (("button_press_event", self.mouse_press),
                              ("motion_notify_event", self.mouse_move),
                              ("button_release_event", self.mouse_release)):
            self.MPL_canvas_original.mpl_connect(event, method)
        self.MPL_canvas_original.setToolTip(
            "Drag rectangles on the high features (cells, debris): their "
            "pixels are left out of the fits of the lines\n(still "
            "flattened with their line)")
        self.BT_clear_excluded = PYAFButton(self, "clear_excluded",
                                            "Clear rectangles", size=0)
        self.BT_clear_excluded.setToolTip(
            "All the pixels used for the fits of the lines")

        self.flatten_buttons_HL = QtWidgets.QHBoxLayout()
        self.BT_apply_flatten = PYAFButton(self, "apply_flatten", "Apply")
        self.BT_remove_flatten = \
            PYAFButton(self, "remove_flatten", "Reset")
        self.BT_remove_flatten.setToolTip("Back to the original heights")
        self.update_GUI("buttons_flatten")

        self.TG_apply_flatten_on_all = \
            PYAFToggle(self, "apply_flatten_on_all", "Apply on all")

        if len(shared.exp.list) < 2:
            self.TG_apply_flatten_on_all.setEnabled(False)
        if shared.exp.flatten_correction_all:
            self.TG_apply_flatten_on_all.setChecked(True)

        # Offsets at the steps of the z motor (JPK Smart Mapping)
        self.TG_flatten_motor = \
            PYAFToggle(self, "flatten_motor", "Level the z motor steps")
        self.TG_flatten_motor.setToolTip(
            "Removes the height offsets at each move of the z motor "
            "(maps higher than the range of the piezo), before the flatten")
        self.update_GUI("checkbox_flatten_motor")

        self.label_radiobuttons_order = QtWidgets.QLabel()
        self.label_radiobuttons_order.setText("Order : ")
        self.RV_flatten_none = QtWidgets.QRadioButton("None")
        self.RV_flatten_order0 = QtWidgets.QRadioButton("0")
        self.RV_flatten_order1 = QtWidgets.QRadioButton("1")
        self.RV_flatten_order2 = QtWidgets.QRadioButton("2")
        self.RV_flatten_none.setToolTip(
            "No flatten of the lines (e.g. only the z motor steps levelled)")
        self.RV_flatten_order0.setToolTip("Offset of each line")
        self.RV_flatten_order1.setToolTip("Offset and slope of each line")
        self.RV_flatten_order2.setToolTip("Offset, slope and curvature of each line")
        self.update_GUI("radiobuttons_flatten_order")
        for name, button in (("none", self.RV_flatten_none),
                             ("order0", self.RV_flatten_order0),
                             ("order1", self.RV_flatten_order1),
                             ("order2", self.RV_flatten_order2)):
            button.clicked.connect(
                lambda _=False, name=name:
                self.button_clicked("radiobutton_flatten_" + name))

        # Order and rectangles; z motor steps; apply
        HL_order = QtWidgets.QHBoxLayout()
        HL_order.addWidget(self.label_radiobuttons_order)
        for button in (self.RV_flatten_none, self.RV_flatten_order0,
                       self.RV_flatten_order1, self.RV_flatten_order2):
            HL_order.addWidget(button)
        HL_order.addStretch(1)
        HL_order.addWidget(self.BT_clear_excluded)

        self.flatten_buttons_HL.addWidget(self.TG_apply_flatten_on_all)
        self.flatten_buttons_HL.addStretch(1)
        self.flatten_buttons_HL.addWidget(self.BT_apply_flatten)
        self.flatten_buttons_HL.addWidget(self.BT_remove_flatten)

        self.VL_box_flatten.addLayout(self.HL_flatten_canvas, 1)
        self.VL_box_flatten.addLayout(HL_order)
        self.VL_box_flatten.addWidget(self.TG_flatten_motor)
        self.VL_box_flatten.addLayout(self.flatten_buttons_HL)
        self.box_flatten.setLayout(self.VL_box_flatten)

        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.addWidget(self.box_flatten, 0, 0)
        self.setLayout(self.gridLayout)

    @staticmethod
    def canvas_widget():
        """Widget holding a map, which grows with the window."""
        widget = QtWidgets.QWidget()
        widget.setMinimumSize(216, 216)
        widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                             QtWidgets.QSizePolicy.Expanding)
        layout = QtWidgets.QVBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        return widget

    def button_clicked(self, name):
        """Called when a button is clicked."""
        data = shared.exp.current_data

        if name == "apply_flatten":
            if shared.exp.flatten_correction_all:
                for i in range(len(shared.exp.list)):
                    data = shared.exp.list[i]
                    data.applied_flatten_order = \
                        shared.exp.selected_flatten_order
                    data.applied_flatten_motor = \
                        shared.exp.selected_flatten_motor and \
                        data.motor_parts is not None
                    # Rectangles of each file (drawn on the current one)
                    data.applied_flatten_excluded = tuple(
                        data.flatten_excluded)
                    data.flatten_applied = True

            else:
                data.applied_flatten_order = shared.exp.selected_flatten_order
                data.applied_flatten_motor = \
                    shared.exp.selected_flatten_motor and \
                    data.motor_parts is not None
                data.applied_flatten_excluded = tuple(data.flatten_excluded)
                data.flatten_applied = True
            self.heights_changed()
            self.BT_remove_flatten.setEnabled(True)
            self.BT_apply_flatten.setEnabled(False)

        elif name == "remove_flatten":
            if shared.exp.flatten_correction_all:
                for i in range(len(shared.exp.list)):
                    shared.exp.list[i].applied_flatten_order = None
                    shared.exp.list[i].applied_flatten_motor = False
                    shared.exp.list[i].applied_flatten_excluded = ()
                    shared.exp.list[i].flatten_applied = False
            else:
                data.applied_flatten_order = None
                data.applied_flatten_motor = False
                data.applied_flatten_excluded = ()
                data.flatten_applied = False

            self.heights_changed()
            self.BT_remove_flatten.setEnabled(False)
            self.BT_apply_flatten.setEnabled(True)

        elif name == "clear_excluded":
            data.flatten_excluded = []
            self.update_MPL("MPL_canvas_flatten")
            self.update_GUI("buttons_flatten")

        elif name == "radiobutton_flatten_none":
            shared.exp.selected_flatten_order = None
            self.update_MPL("MPL_canvas_flatten")

        elif name == "radiobutton_flatten_order0":
            shared.exp.selected_flatten_order = 0
            self.update_MPL("MPL_canvas_flatten")

        elif name == "radiobutton_flatten_order1":
            shared.exp.selected_flatten_order = 1
            self.update_MPL("MPL_canvas_flatten")

        elif name == "radiobutton_flatten_order2":
            shared.exp.selected_flatten_order = 2
            self.update_MPL("MPL_canvas_flatten")

    def pixel(self, event, inside=True):
        """Pixel (x, y) of a mouse event on the original map, None outside
        of the map; with inside False, the nearest pixel of the map (the
        drag of a rectangle can go beyond its edges)."""
        canvas = self.MPL_canvas_original.canvas
        axes = getattr(canvas, "axes", None)
        if axes is None:
            return None
        if inside:
            if event.inaxes is not axes or event.xdata is None:
                return None
            x, y = event.xdata, event.ydata
        else:
            # Position in the coordinates of the map, also outside of it
            x, y = axes.transData.inverted().transform((event.x, event.y))
        data = shared.exp.current_data
        i = int(numpy.floor(x / canvas.x_size))
        j = int(numpy.floor(y / canvas.y_size))
        return (min(max(i, 0), data.nbr_pixels_x - 1),
                min(max(j, 0), data.nbr_pixels_y - 1))

    def rectangle(self, end):
        """Rectangle (x first, x last, y first, y last) from the start of
        the drag to a pixel."""
        (i0, j0), (i1, j1) = self.drag_start, end
        return (min(i0, i1), max(i0, i1), min(j0, j1), max(j0, j1))

    def mouse_press(self, event):
        """Start of a rectangle (left button)."""
        if event.button == 1:
            self.drag_start = self.pixel(event)

    def mouse_move(self, event):
        """Rectangle being drawn."""
        if self.drag_start is None:
            return
        end = self.pixel(event, inside=False)
        if end is not None:
            self.MPL_canvas_original.canvas.draw_excluded(
                self.rectangle(end))

    def mouse_release(self, event):
        """Rectangle added: the flattened map without its pixels in the
        fits."""
        if self.drag_start is None:
            return
        end = self.pixel(event, inside=False)
        if end is not None:
            shared.exp.current_data.flatten_excluded.append(
                self.rectangle(end))
        self.drag_start = None
        self.update_MPL("MPL_canvas_flatten")
        self.update_GUI("buttons_flatten")

    @staticmethod
    def heights_changed():
        """Topography, maps and values after a change of the flatten (of
        all the files with Apply on all)."""
        if shared.exp.flatten_correction_all:
            datas = list(shared.exp.list)
        else:
            datas = [shared.exp.current_data]
        roi_manager.heights_changed(datas)

    def checkbox_clicked(self, name):
        """Called when a checkbox is clicked."""
        if name == "apply_flatten_on_all":
            state = self.TG_apply_flatten_on_all.isChecked()
            if state:
                shared.exp.flatten_correction_all = True
            else:
                shared.exp.flatten_correction_all = False

        elif name == "flatten_motor":
            shared.exp.selected_flatten_motor = \
                self.TG_flatten_motor.isChecked()
            self.update_MPL("MPL_canvas_flatten")

    def update_MPL(self, name):
        """Updates the plot."""
        if name == "MPL_canvas_flatten":
            # Original piezo_image
            self.MPL_canvas_original.update_plot()

            # Flattened image
            self.MPL_canvas_flatten.update_plot()

    def update_widget(self):
        """Update the widget (GUI and plots)."""
        self.update_MPL("MPL_canvas_flatten")
        self.update_GUI("radiobuttons_flatten_order")
        self.update_GUI("checkbox_flatten_motor")
        self.update_GUI("buttons_flatten")

    def update_GUI(self, what):
        """Update the GUI."""
        if shared.exp.current_data.flatten_applied:
            order = shared.exp.current_data.applied_flatten_order
        else:
            order = shared.exp.selected_flatten_order

        if what == "radiobuttons_flatten_order":
            if order is None:
                self.RV_flatten_none.setChecked(True)
            elif order == 0:
                self.RV_flatten_order0.setChecked(True)
            elif order == 2:
                self.RV_flatten_order2.setChecked(True)
            else:
                self.RV_flatten_order1.setChecked(True)

        elif what == "checkbox_flatten_motor":
            data = shared.exp.current_data
            motor = shared.exp.selected_flatten_motor
            if data.flatten_applied:
                motor = data.applied_flatten_motor
            self.TG_flatten_motor.setEnabled(data.motor_parts is not None)
            self.TG_flatten_motor.setChecked(
                motor and data.motor_parts is not None)

        elif what == "buttons_flatten":
            data = shared.exp.current_data
            self.BT_clear_excluded.setEnabled(bool(data.flatten_excluded))
            # A change of the rectangles can be applied again
            changed = data.flatten_applied and \
                tuple(data.flatten_excluded) != \
                tuple(data.applied_flatten_excluded)
            if changed:
                self.BT_apply_flatten.setEnabled(True)
                self.BT_remove_flatten.setEnabled(True)
                return
            if shared.exp.current_data.flatten_applied is False:
                self.BT_remove_flatten.setEnabled(False)
                self.BT_apply_flatten.setEnabled(True)
            else:
                self.BT_apply_flatten.setEnabled(False)
                self.BT_remove_flatten.setEnabled(True)
