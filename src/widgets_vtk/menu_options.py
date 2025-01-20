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

"""Menu for the 3D widget."""

from PyQt5 import QtCore, QtWidgets
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFButtonGroup
from ..tools import misc_tools
from .slicing import SlicingWidget
from ..widgets.plot_options import PlotOptionsWidget
from ..widgets_small.indentation import WidgetIndentation
from .utils.binary_threshold_widget import ThresholdOptionsWidget
from .. import widgets_list
from .. import shared
from .. import consts


class MenuOptionsWidget(QtWidgets.QWidget):
    """This widget contains the menu with to options for a particular layer."""

    def __init__(self):
        super().__init__()

        VL = QtWidgets.QVBoxLayout()
        VL.setContentsMargins(0, 0, 0, 0)

        # Tabs
        self.tabs = QtWidgets.QTabWidget()

        self.tab_options = WidgetOptions(self)
        self.tab_afm = WidgetAFM(self)
        self.tab_afm2 = WidgetAFM2(self)

        self.tab_HR = None
        self.tab_Optics = None
        self.tab_glass = None
        self.tab_isosurface = None

        self.tabs.addTab(self.tab_options, "Options")
        self.tabs.addTab(self.tab_afm, "AFM")
        self.tabs.addTab(self.tab_afm2, "AFM 2")

        VL.addWidget(self.tabs)

        self.setLayout(VL)

    def update_widget(self):
        """Update the GUI."""
        row = widgets_list.widget_vtk.get_current_layer()

        # Get current tab index to reselect it afterwards
        index = self.tabs.currentIndex()

        # Remove tabs
        if self.tabs.count() == 3:
            self.tabs.removeTab(2)
        self.tabs.removeTab(1)

        # Get the type of tab and display it
        if shared.layer_list[row].type == "afm":
            self.tabs.addTab(self.tab_afm, "AFM 1")
            self.tabs.addTab(self.tab_afm2, "AFM 2")
            self.tab_afm.update_widget()
            self.tab_afm2.update_widget()
        elif shared.layer_list[row].type == "single_image" or \
                shared.layer_list[row].type == "stack":
            self.tab_Optics = WidgetOptics(self)
            self.tabs.addTab(self.tab_Optics, "Image")
            self.tab_Optics.update_widget()
        elif shared.layer_list[row].type == "HR":
            self.tab_HR = WidgetOptics(self)
            self.tabs.addTab(self.tab_HR, "HR")
            self.tab_HR.update_widget()
        elif shared.layer_list[row].type == "glass":
            self.tab_glass = WidgetGlass(self)
            self.tabs.addTab(self.tab_glass, "Glass")
            self.tab_glass.update_widget()
        elif shared.layer_list[row].type == "isosurface":
            self.tab_isosurface = WidgetIsoSurface(self)
            self.tabs.addTab(self.tab_isosurface, "Isosurface")
            self.tab_isosurface.update_widget()

        self.tab_options.update_widget()

        # Select tab
        self.tabs.setCurrentIndex(index)


class WidgetOptions(QtWidgets.QWidget):
    """Tab containing the main options for a layer.

    This tab is common for all layers. It allows to change the position of the
    layer and it's opacity.
    """

    def __init__(self, parent):
        super().__init__()

        self.parent = parent

        self.VL = QtWidgets.QVBoxLayout()

        # Camera actions
        self.VL_actions = QtWidgets.QVBoxLayout()
        self.BG_action = PYAFButtonGroup(self, "RBT_action")
        self.RBT_camera = QtWidgets.QRadioButton("Camera")
        self.RBT_translate = QtWidgets.QRadioButton("Translate")
        self.RBT_rotate = QtWidgets.QRadioButton("Rotate")
        self.RBT_camera.setChecked(True)
        self.BG_action.addButton(self.RBT_camera, 0)
        self.BG_action.addButton(self.RBT_translate, 1)
        self.BG_action.addButton(self.RBT_rotate, 2)
        self.VL_actions.addWidget(self.RBT_camera)
        self.VL_actions.addWidget(self.RBT_translate)
        self.VL_actions.addWidget(self.RBT_rotate)

        # Position
        self.box = QtWidgets.QGroupBox("Position")
        self.vl_box = QtWidgets.QVBoxLayout()

        self.IN_pos_x = PYAFInput(self, "pos_x", "X", width=200)
        self.IN_pos_x.input.setValidator(misc_tools.validator("F"))
        self.IN_pos_y = PYAFInput(self, "pos_y", "Y", width=200)
        self.IN_pos_y.input.setValidator(misc_tools.validator("F"))
        self.IN_pos_z = PYAFInput(self, "pos_z", "Z", width=200)
        self.IN_pos_z.input.setValidator(misc_tools.validator("F"))
        self.IN_angle = PYAFInput(self, "angle", "Angle", width=200)
        self.IN_angle.input.setValidator(misc_tools.validator("F"))

        layer = shared.layer_list[0]
        self.IN_pos_x.changeValue(str(layer.pos_x / 1000.0))
        self.IN_pos_y.changeValue(str(layer.pos_y / 1000.0))
        self.IN_pos_z.changeValue(str(layer.pos_z / 1000.0))
        self.IN_angle.changeValue(str(layer.angle))
        self.vl_box.addWidget(self.IN_pos_x)
        self.vl_box.addWidget(self.IN_pos_y)
        self.vl_box.addWidget(self.IN_pos_z)
        self.vl_box.addWidget(self.IN_angle)
        self.box.setLayout(self.vl_box)

        text = "Opacity (%)"
        self.IN_opacity = PYAFInput(self, "opacity", text, width=80)
        self.IN_opacity.input.setValidator(misc_tools.validator("RI", bottom=0, top=100))

        label = "Display Z scale"
        self.CB_display_zscale = PYAFCheckBox(self, "disp_zscale", label)

        self.VL.addLayout(self.VL_actions)
        self.VL.addWidget(self.box)
        self.VL.addWidget(self.IN_opacity)
        self.VL.addWidget(self.CB_display_zscale)
        self.VL.addStretch(1)

        self.setLayout(self.VL)

    def button_clicked(self, button):
        """Called when a button is clicked."""
        if button == "RBT_action":
            wl = widgets_list.widget_vtk
            index = self.BG_action.checkedId()

            if index == 0:
                wl.action_mode = "camera"
            elif index == 1:
                wl.action_mode = "translate"
                wl.menu_layer_list.button_clicked("lock_top")
            elif index == 2:
                wl.action_mode = "rotate"

    def update_widget(self):
        """Update the GUI."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        # Update the inputs
        self.IN_pos_x.changeValue(str(layer.pos_x / 1000.0))
        self.IN_pos_y.changeValue(str(layer.pos_y / 1000.0))
        self.IN_pos_z.changeValue(str(layer.pos_z / 1000.0))
        self.IN_angle.changeValue(str(layer.angle))
        self.IN_opacity.changeValue(str(layer.opacity))

        # Display or hide the z scale
        self.CB_display_zscale.setChecked(layer.display_z_scale)

    def input_updated(self, name):
        """Called when an input field is updated."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if name == "pos_x" or name == "pos_y" or \
                name == "pos_z" or name == "angle":
            newpos_x = self.IN_pos_x.get_float_value() * 1000.0
            newpos_y = self.IN_pos_y.get_float_value() * 1000.0
            newpos_z = self.IN_pos_z.get_float_value() * 1000.0
            angle = self.IN_angle.get_float_value()

            if angle == layer.angle:
                dx, dy, dz = 0, 0, 0
                if newpos_x != layer.pos_x:
                    dx = newpos_x - layer.pos_x
                if newpos_y != layer.pos_y:
                    dy = newpos_y - layer.pos_y
                if newpos_z != layer.pos_z:
                    dz = newpos_z - layer.pos_z

                layer.translate(dx, dy, dz)

            else:
                layer.rotate(None, None, increment=(angle - layer.angle))

        elif name == "opacity":
            value = self.IN_opacity.get_int_value()
            shared.layer_list[index].set_opacity(value)

    def checkbox_clicked(self, name):
        """Called whenever a checkbox is clicked."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if name == "disp_zscale":
            layer.display_z_scale = self.CB_display_zscale.isChecked()
            # Reinitialize the axes
            layer.init_axes()


class WidgetOptics(QtWidgets.QWidget):
    """Widget with options for the optics layer."""

    def __init__(self, parent):
        super().__init__()

        self.parent = parent

        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]
        self.type = layer.type

        # Load the settings
        self.settings = QtCore.QSettings()

        self.VL = QtWidgets.QVBoxLayout()

        if self.type == "single_image" or self.type == "stack":
            self.box = QtWidgets.QGroupBox()
            self.vl_box = QtWidgets.QVBoxLayout()

            label = "X size (um)"
            self.IN_x_size = PYAFInput(self, "size_x", label, width=100)
            self.IN_x_size.input.setValidator(misc_tools.validator("UF"))
            label = "Y size (um)"
            self.IN_y_size = PYAFInput(self, "size_y", label, width=100)
            self.IN_y_size.input.setValidator(misc_tools.validator("UF"))
            label = "Z size (um)"
            self.IN_z_size = PYAFInput(self, "size_z", label, width=100)
            self.IN_z_size.input.setValidator(misc_tools.validator("UF"))
            label = "Pixel width (um)"
            self.IN_x_res = PYAFInput(self, "res_x", label, width=100)
            self.IN_x_res.input.setValidator(misc_tools.validator("UF"))
            label = "Pixel height (um)"
            self.IN_y_res = PYAFInput(self, "res_y", label, width=100)
            self.IN_y_res.input.setValidator(misc_tools.validator("UF"))
            label = "Pixel height (um)"
            self.IN_z_res = PYAFInput(self, "res_z", label, width=100)
            self.IN_z_res.input.setValidator(misc_tools.validator("UF"))
            self.CB_bind = PYAFCheckBox(self, "bind", "Maintain Aspect Ratio")

            self.BT_apply = PYAFButton(self, "apply", "Apply changes", 200)
            self.BT_reset = PYAFButton(self, "reset", "Reset to default", 200)
            self.vl_box.addWidget(self.IN_x_size)
            self.vl_box.addWidget(self.IN_y_size)
            self.vl_box.addWidget(self.IN_z_size)
            self.vl_box.addWidget(self.IN_x_res)
            self.vl_box.addWidget(self.IN_y_res)
            self.vl_box.addWidget(self.IN_z_res)
            self.vl_box.addWidget(self.CB_bind)
            self.vl_box.addWidget(self.BT_apply)
            self.vl_box.addWidget(self.BT_reset)
            self.box.setLayout(self.vl_box)

            self.VL.addWidget(self.box)

            self.BT_use_LUT = PYAFButton(self, "use_lut", "Use LUT", 200)
            label = "LUT threshold "
            self.IN_threshold = PYAFInput(self, "input_threshold", label)
            self.IN_threshold.input.setValidator(misc_tools.validator("UF"))
            label = "Invert LUT"
            self.BT_invert_LUT = PYAFButton(self, "invert_lut", label, 200)
            label = "LUT max color"
            self.BT_LUT_max_col = PYAFButton(self, "lut_max_color", label, 200)
            self.VL.addWidget(self.BT_use_LUT)
            self.VL.addWidget(self.BT_invert_LUT)
            self.VL.addWidget(self.BT_LUT_max_col)
            self.VL.addWidget(self.IN_threshold)

        elif self.type == "HR":
            self.BT_HR_color = PYAFButton(self, "HR_color", "HR Color", 200)
            self.VL.addWidget(self.BT_HR_color)

        if self.type == "stack":
            self.IN_slide = PYAFInput(self, "input_slide", "Slide", width=100)
            self.IN_slide.input.setValidator(misc_tools.validator("UF"))

            label = "Single opacity"
            name = "slide_op_single"
            self.IN_slide_op_single = PYAFInput(self, name, label, width=100)
            self.IN_slide_op_single.input.setValidator(
                misc_tools.validator("RI", bottom=0, top=100))
            label = "Others opacity"
            name = "slide_op_others"
            self.IN_slide_op_others = PYAFInput(self, name, label, width=100)
            self.IN_slide_op_others.input.setValidator(
                misc_tools.validator("RI", bottom=0, top=100))

            self.slider = QtWidgets.QSlider()
            self.slider.setRange(0, len(layer.actor_list) - 1)
            self.slider.setSingleStep(1)
            self.slider.setValue(layer.slide_pos)
            self.slider.setOrientation(QtCore.Qt.Horizontal)
            self.slider.valueChanged.connect(self.slider_moved)

            self.VL.addWidget(self.slider)
            self.VL.addWidget(self.IN_slide)
            self.VL.addWidget(self.IN_slide_op_single)
            self.VL.addWidget(self.IN_slide_op_others)

            label = "Binary Thresholding"
            self.BT_binary_tresh = PYAFButton(self, "binary_tresh", label, 200)

            self.VL.addWidget(self.BT_binary_tresh)

        # Options to flip the images (X or Y, mirrored flip)
        VL_flip = QtWidgets.QVBoxLayout()

        self.CB_flip_x = PYAFCheckBox(self, "flip_x", "Flip around X axis")
        self.CB_flip_y = PYAFCheckBox(self, "flip_y", "Flip around Y axis")

        VL_flip.addWidget(self.CB_flip_x)
        VL_flip.addWidget(self.CB_flip_y)

        if consts.ADVANCED:
            self.VL.addLayout(VL_flip)

        # Option to smooth the image
        VL_smooth = QtWidgets.QVBoxLayout()

        self.IN_smooth = PYAFInput(self, "smooth", "Smooth (Gaussian)")
        self.IN_smooth.input.setValidator(misc_tools.validator("UF"))

        VL_smooth.addWidget(self.IN_smooth)

        if consts.ADVANCED:
            self.VL.addLayout(VL_smooth)

        self.VL.addStretch(1)

        self.setLayout(self.VL)

    def button_clicked(self, button):
        """Called when a button is clicked."""
        self.setFocus()

        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if button == "reset":
            layer.nbr_pixels_x = layer.original_nbr_pixels_x
            layer.nbr_pixels_y = layer.original_nbr_pixels_y
            layer.nbr_pixels_z = layer.original_nbr_pixels_z
            layer.res_x = layer.original_res_x
            layer.res_y = layer.original_res_y
            layer.res_z = layer.original_res_z
            layer.size_x = layer.original_size_x
            layer.size_y = layer.original_size_y
            layer.size_z = layer.original_size_z

            layer.update_actor()
            self.update_widget()

        elif button == "use_lut":
            if layer.use_lut:
                layer.use_lut = False
            else:
                layer.use_lut = True

            layer.update_actor()

            self.update_widget()

        elif button == "invert_lut":
            if layer.lut_invert:
                layer.lut_invert = False
            else:
                layer.lut_invert = True

            layer.update_actor()

        elif button == "HR_color":
            # self.colorwidget_HR = MPQColorDialog(self, "HR_color")
            pass

        elif button == "apply":
            layer.update_actor()
            self.update_widget()

        elif button == "lut_max_color":
            color = misc_tools.ask_user_for_color(layer.lut_max_color)
            if color is not False:
                layer.lut_max_color = color
                layer.update_actor()

        elif button == "binary_tresh":
            if widgets_list.widget_vtk_binary_threshold is None:
                ThresholdOptionsWidget(self)
                widgets_list.widget_vtk_binary_threshold.show()
            else:
                widgets_list.widget_vtk_binary_threshold.activateWindow()
                widgets_list.widget_vtk_binary_threshold.raise_()

    def slider_moved(self):
        """Called when the slider for the stacks is moved."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        layer.slide_pos = self.slider.value()
        self.IN_slide.changeValue(layer.slide_pos)
        layer.set_opacity()

    def checkbox_clicked(self, name):
        """Called when a checkbox is clicked."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if name == "bind":
            layer.bind = self.CB_bind.isChecked()

        elif name == "flip_x":
            layer.flip_x = self.CB_flip_x.isChecked()
            layer.update_actor()

        elif name == "flip_y":
            layer.flip_y = self.CB_flip_y.isChecked()
            layer.update_actor()

    def input_updated(self, name):
        """Called when a value is changed in an input field."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if name == "size_x":
            new_val = self.IN_x_size.get_float_value() * 1000.0
            layer.size_x = new_val
            layer.res_x = layer.size_x / layer.nbr_pixels_x
            if layer.bind:
                layer.size_y = layer.size_x
                layer.res_y = layer.res_x

        elif name == "size_y":
            new_val = self.IN_y_size.get_float_value() * 1000.0
            layer.size_y = new_val
            layer.res_y = layer.size_y / layer.nbr_pixels_y
            if layer.bind:
                layer.size_x = layer.size_y
                layer.res_x = layer.res_y

        elif name == "size_z":
            new_val = self.IN_z_size.get_float_value() * 1000.0
            layer.size_z = new_val
            layer.res_z = layer.size_z / layer.nbr_pixels_z

        elif name == "res_x":
            new_val = self.IN_x_res.get_float_value() * 1000.0
            layer.res_x = new_val
            layer.size_x = layer.res_x * layer.nbr_pixels_x
            if layer.bind:
                layer.size_y = layer.size_x
                layer.res_y = layer.res_x

        elif name == "res_y":
            new_val = self.IN_y_res.get_float_value() * 1000.0
            layer.res_y = new_val
            layer.size_y = layer.res_y * layer.nbr_pixels_y
            if layer.bind:
                layer.size_x = layer.size_y
                layer.res_x = layer.res_y

        elif name == "res_z":
            new_val = self.IN_z_res.get_float_value() * 1000.0
            layer.res_z = new_val
            layer.size_z = layer.res_z * layer.nbr_pixels_z

        elif name == "input_threshold":
            layer.lut_threshold = self.IN_threshold.get_int_value()

            layer.update_actor()

        elif name == "input_slide":
            layer.slide_pos = self.IN_slide.get_int_value()

            layer.set_opacity()

        elif name == "slide_op_single":
            val = self.IN_slide_op_single.get_int_value()
            layer.single_opacity = val

            layer.set_opacity()

        elif name == "slide_op_others":
            val = self.IN_slide_op_others.get_int_value()
            layer.others_opacity = val

            layer.set_opacity()

        elif name == "smooth":
            layer.smooth_gaussian = self.IN_smooth.get_float_value()
            layer.update_actor()

        self.update_widget()

    def update_widget(self):
        """Updates the GUI."""
        layer_id = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[layer_id]

        if layer.type == "single_image" or layer.type == "stack":
            self.IN_x_size.changeValue(layer.size_x / 1000.0)
            self.IN_y_size.changeValue(layer.size_y / 1000.0)
            self.IN_x_res.changeValue(layer.res_x / 1000.0)
            self.IN_y_res.changeValue(layer.res_y / 1000.0)
            self.IN_threshold.changeValue(layer.lut_threshold)
            if layer.use_lut:
                self.BT_use_LUT.setText("No LUT")
            else:
                self.BT_use_LUT.setText("Use LUT")

            self.CB_bind.setChecked(layer.bind)

            if layer.type == "single_image":
                self.IN_z_size.changeValue("None")
                self.IN_z_size.setEnabled(False)
                self.IN_z_res.changeValue("None")
                self.IN_z_res.setEnabled(False)
            else:
                self.IN_z_size.changeValue(layer.size_z / 1000.0)
                self.IN_z_size.setEnabled(True)
                self.IN_z_res.changeValue(layer.res_z / 1000.0)
                self.IN_z_res.setEnabled(True)

            self.IN_smooth.changeValue(layer.smooth_gaussian)

        if layer.type == "stack":
            self.IN_slide.changeValue(layer.slide_pos)
            self.IN_slide_op_single.changeValue(layer.single_opacity)
            self.IN_slide_op_others.changeValue(layer.others_opacity)
            self.slider.setValue(layer.slide_pos)

        self.CB_flip_x.setChecked(layer.flip_x)
        self.CB_flip_y.setChecked(layer.flip_y)


class WidgetGlass(QtWidgets.QWidget):
    """Widget with options for the glass layers."""

    def __init__(self, parent):
        super().__init__()

        self.parent = parent
        self.type = type

        self.VL = QtWidgets.QVBoxLayout()

        self.box = QtWidgets.QGroupBox()
        self.vl_box = QtWidgets.QVBoxLayout()

        label = "X size (um) : "
        self.IN_x_size = PYAFInput(self, "size_x", label, width=100)
        self.IN_x_size.input.setValidator(misc_tools.validator("UF"))
        label = "Y size (um) : "
        self.IN_y_size = PYAFInput(self, "size_y", label, width=100)
        self.IN_y_size.input.setValidator(misc_tools.validator("UF"))
        self.CB_bind = PYAFCheckBox(self, "bind", "Maintain Aspect Ratio")
        self.BT_reset = PYAFButton(self, "reset", "Reset to default", 200)
        self.vl_box.addWidget(self.IN_x_size)
        self.vl_box.addWidget(self.IN_y_size)
        self.vl_box.addWidget(self.CB_bind)
        self.vl_box.addWidget(self.BT_reset)
        self.box.setLayout(self.vl_box)

        self.VL.addWidget(self.box)

        self.VL.addStretch(1)

        self.setLayout(self.VL)

    def checkbox_clicked(self):
        """Called when a checkbox is clicked."""
        layer_id = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[layer_id]

        layer.bind = self.CB_bind.isChecked()

    def button_clicked(self, button):
        """Called when a button is clicked."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if button == "reset":
            layer.size_x = layer.original_size_x
            layer.size_y = layer.original_size_y

            layer.update_actor()
            self.update_widget()

    def input_updated(self, name):
        """Called when an input is updated."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if name == "size_x":
            new_val = self.IN_x_size.get_float_value() * 1000.0
            if new_val != layer.size_x:
                layer.size_x = new_val
                if layer.bind:
                    layer.size_y = layer.size_x
        elif name == "size_y":
            new_val = self.IN_y_size.get_float_value() * 1000.0
            if new_val != layer.size_y:
                layer.size_y = new_val
                if layer.bind:
                    layer.size_x = layer.size_y

        layer.update_actor()
        self.update_widget()

    def update_widget(self):
        """Update the widget."""
        layer_id = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[layer_id]

        self.IN_x_size.changeValue(layer.size_x / 1000.0)
        self.IN_y_size.changeValue(layer.size_y / 1000.0)

        self.CB_bind.setChecked(layer.bind)


class WidgetAFM(QtWidgets.QWidget):
    """Widget with options for the AFM layer."""

    def __init__(self, parent):
        super().__init__()

        self.parent = parent
        self.open = True

        # Load the settings
        self.settings = QtCore.QSettings()

        self.VL = QtWidgets.QVBoxLayout()

        # parent.parent.parent = MenuWidget.OpenglWidget.DataWidget
        self.afm_colors_widget = PlotOptionsWidget(self, "vtk")
        self.afm_indentation_widget = WidgetIndentation(self)
        self.VL.addWidget(self.afm_colors_widget)
        self.VL.addWidget(self.afm_indentation_widget)

        label = "Refresh background (slower)"
        self.CB_refresh = PYAFCheckBox(self, "refresh_background", label)
        checked = self.settings.value("refreshBackgroundInOpenGL", True)
        checked = checked
        if checked:
            self.CB_refresh.setChecked(True)
        else:
            self.CB_refresh.setChecked(False)

        self.BT_slicing = PYAFButton(self, "slicing", "Slicing")

        self.VL.addWidget(self.CB_refresh)
        self.VL.addWidget(self.BT_slicing)
        self.VL.addStretch(1)

        self.setLayout(self.VL)

    def checkbox_clicked(self, name):
        """Called when a checkbox is clicked."""
        if name == "refresh_background":
            value = self.CB_refresh.isChecked()
            self.settings.setValue("refreshBackgroundInOpenGL", value)

    def button_clicked(self, button):
        """Called when a button is clicked."""
        if button == "slicing":
            if widgets_list.widget_slicing is None:
                # Create new widget
                SlicingWidget(self)
                widgets_list.widget_slicing.setWindowTitle("3D Slicing")
                widgets_list.widget_slicing.show()
            else:
                # Bring to front
                widgets_list.widget_slicing.activateWindow()
                widgets_list.widget_slicing.raise_()

    def update_widget(self):
        """Updates the widget."""
        layer_id = widgets_list.widget_vtk.get_current_layer()
        afm_id = shared.layer_list[layer_id].afm_id
        self.afm_colors_widget.update_widget(data_id=afm_id)
        self.afm_indentation_widget.update_widget(data_id=afm_id)


class WidgetAFM2(QtWidgets.QWidget):
    """Widget with options for the AFM layer."""

    def __init__(self, parent):
        super().__init__()

        self.parent = parent
        self.open = True

        self.VL = QtWidgets.QVBoxLayout()

        # Option to chose between 3D or 2D
        HL_mode = QtWidgets.QHBoxLayout()
        self.GRP_mode = PYAFButtonGroup(self, "RBT_mode")
        self.RBT_3D = QtWidgets.QRadioButton("3D")
        self.RBT_2D = QtWidgets.QRadioButton("2D")
        self.GRP_mode.addButton(self.RBT_3D, 0)
        self.GRP_mode.addButton(self.RBT_2D, 1)
        HL_mode.addWidget(self.RBT_3D)
        HL_mode.addWidget(self.RBT_2D)

        self.VL_radio_buttons = QtWidgets.QVBoxLayout()
        self.GRP_surf_type = PYAFButtonGroup(self, "RBT_action")
        self.RBT_stiffness = QtWidgets.QRadioButton("Elasticity")
        self.RBT_stiffness_sq = QtWidgets.QRadioButton("Elasticity (squares)")
        self.RBT_topo_dots = QtWidgets.QRadioButton("Topo (dots)")
        self.RBT_topo_surf = QtWidgets.QRadioButton("Topo (surf)")
        self.GRP_surf_type.addButton(self.RBT_stiffness, 0)
        self.GRP_surf_type.addButton(self.RBT_stiffness_sq, 1)
        self.GRP_surf_type.addButton(self.RBT_topo_dots, 2)
        self.GRP_surf_type.addButton(self.RBT_topo_surf, 3)
        self.VL_radio_buttons.addWidget(self.RBT_stiffness)
        self.VL_radio_buttons.addWidget(self.RBT_stiffness_sq)
        # self.VL_radio_buttons.addWidget(self.RBT_topo_dots)
        self.VL_radio_buttons.addWidget(self.RBT_topo_surf)

        # Add 3 bouttons to chose the type of smoothing
        VL_smoothing = QtWidgets.QVBoxLayout()
        self.GRP_smoothing = PYAFButtonGroup(self, "RBT_smooth")
        self.RBT_none = QtWidgets.QRadioButton("No smoothing")
        self.RBT_loop = QtWidgets.QRadioButton("LoopSubdivision")
        self.RBT_butterfly = QtWidgets.QRadioButton("ButterflySubdivision")
        self.GRP_smoothing.addButton(self.RBT_none, 0)
        self.GRP_smoothing.addButton(self.RBT_loop, 1)
        self.GRP_smoothing.addButton(self.RBT_butterfly, 2)

        # Add inputs with smoothing options
        name = "Smoothing iterations"
        self.IN_smooth_iter = PYAFInput(self, "smooth_iter", name, width=50)
        self.IN_smooth_iter.input.setValidator(misc_tools.validator("UF"))
        name = "Target decimation"
        self.IN_targ_decim = PYAFInput(self, "targ_decimate", name, width=50)
        self.IN_targ_decim.input.setValidator(misc_tools.validator("UF"))

        BT_apply_smoothing = PYAFButton(self, "apply_smoothing", "Apply")

        VL_smoothing.addWidget(self.RBT_none)
        VL_smoothing.addWidget(self.RBT_loop)
        VL_smoothing.addWidget(self.RBT_butterfly)
        VL_smoothing.addWidget(self.IN_smooth_iter)
        VL_smoothing.addWidget(self.IN_targ_decim)
        VL_smoothing.addWidget(BT_apply_smoothing)

        self.VL.addLayout(HL_mode)
        self.VL.addLayout(self.VL_radio_buttons)
        self.VL.addLayout(VL_smoothing)
        self.VL.addStretch(1)

        self.setLayout(self.VL)

        self.update_widget()

    def _input_updated(self, val):
        """Called when an input is updated."""
        # Can be called only if the widget still exists
        if self.open:
            self.input_updated(val)

    def closeEvent(self, event):
        """Called when the widget is closed."""
        self.open = False
        # Close the widget
        event.accept()

    def button_clicked(self, button):
        """Called when a button is clicked."""
        layer_id = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[layer_id]
        data = shared.exp.list[layer.afm_id]

        if button == "RBT_mode":
            val = self.GRP_mode.checkedId()
            if val == 0:
                data.opengl_afm_flat_view = False
            elif val == 1:
                data.opengl_afm_flat_view = True

            self.update_widget()

            layer.update_tomography()

        elif button == "RBT_action":
            val = self.GRP_surf_type.checkedId()
            widg = widgets_list.widget_results
            if val == 0:
                data.opengl_surf_type = "Stiffness"
                widg.RBT_mesh_stiffness.setChecked(True)
            elif val == 1:
                data.opengl_surf_type = "Stiffness_squares"
                widg.RBT_mesh_stiffness.setChecked(True)
            elif val == 2:
                data.opengl_surf_type = "Topo_dots"
                widg.RBT_mesh_topography.setChecked(True)
            elif val == 3:
                data.opengl_surf_type = "Topo_surf"
                widg.RBT_mesh_topography.setChecked(True)
            widgets_list.widget_results.button_clicked("BTG_meshgrid_type")

            layer.update_tomography()

        elif button == "RBT_smooth":
            # Change the type of smoothing
            val = self.GRP_smoothing.checkedId()
            if val == 0:
                data.vtk_smoothing_type = None
            elif val == 1:
                data.vtk_smoothing_type = "loop"
            elif val == 2:
                data.vtk_smoothing_type = "butterfly"

        elif button == "apply_smoothing":
            # Set focus so that input_updated is called if needed
            # Will update the parameters for the smooting
            self.setFocus()
            layer.update_smoothing()

    def input_updated(self, name):
        """Called when an input is updated."""
        layer_id = widgets_list.widget_vtk.get_current_layer()
        afm_id = shared.layer_list[layer_id].afm_id
        data = shared.exp.list[afm_id]

        if name == "smooth_iter":
            data.vtk_smoothing_iterations = self.IN_smooth_iter.get_int_value()

        elif name == "targ_decimate":
            val = self.IN_targ_decim.get_float_value()
            data.vtk_decimate_target_reduction = val

    def update_widget(self):
        """Updates the widget."""
        layer_id = widgets_list.widget_vtk.get_current_layer()
        afm_id = shared.layer_list[layer_id].afm_id
        data = shared.exp.list[afm_id]

        if data.opengl_afm_flat_view:
            self.RBT_2D.setChecked(True)
            self.RBT_stiffness_sq.setEnabled(True)
        else:
            self.RBT_3D.setChecked(True)
            self.RBT_stiffness_sq.setEnabled(False)
            # Do not allow stiffness squares in 3D
            if data.opengl_surf_type == "Stiffness_squares":
                data.opengl_surf_type = "Stiffness"
                self.RBT_stiffness.setChecked(True)

        # Type of surface to display (set different meshgrid type for color)
        # Changing the meshgrid type for this is not very clean ...
        if data.opengl_surf_type == "Stiffness":
            self.RBT_stiffness.setChecked(True)
            widgets_list.widget_results.RBT_mesh_stiffness.setChecked(True)
        elif data.opengl_surf_type == "Stiffness_squares":
            self.RBT_stiffness_sq.setChecked(True)
            widgets_list.widget_results.RBT_mesh_stiffness.setChecked(True)
        elif data.opengl_surf_type == "Topo_dots":
            self.RBT_topo_dots.setChecked(True)
            widgets_list.widget_results.RBT_mesh_topography.setChecked(True)
        elif data.opengl_surf_type == "Topo_surf":
            self.RBT_topo_surf.setChecked(True)
            widgets_list.widget_results.RBT_mesh_topography.setChecked(True)
        widgets_list.widget_results.button_clicked("BTG_meshgrid_type")

        # Smoothing options
        if data.vtk_smoothing_type is None:
            self.RBT_none.setChecked(True)
        elif data.vtk_smoothing_type == "loop":
            self.RBT_loop.setChecked(True)
        elif data.vtk_smoothing_type == "butterfly":
            self.RBT_butterfly.setChecked(True)

        # Other smoothing options
        self.IN_smooth_iter.changeValue(data.vtk_smoothing_iterations)
        self.IN_targ_decim.changeValue(data.vtk_decimate_target_reduction)


class WidgetIsoSurface(QtWidgets.QWidget):
    """Widget with options for the isosurface layer"""

    def __init__(self, parent):
        super().__init__()

        self.parent = parent

        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]
        self.type = layer.type

        self.VL = QtWidgets.QVBoxLayout()

        self.BT_color = PYAFButton(self, "color", "Isosurface color")
        self.VL.addWidget(self.BT_color)

        self.VL.addStretch(1)

        self.setLayout(self.VL)

    def button_clicked(self, name):
        """Called when a button is clicked."""
        self.setFocus()

        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]

        if name == "color":
            color = misc_tools.ask_user_for_color(layer.isosurface_color)
            if color is not False:
                layer.isosurface_color = color
                layer.update_actor()

    def update_widget(self):
        """Updates the GUI. (Do nothing here for the moment)"""
        pass
