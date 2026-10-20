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

from qtpy import QtCore, QtGui, QtWidgets
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFButtonGroup
from ..tools.gui_tools import PYAFComboBox
from ..tools import map_values
from .utils import vtk_afm
from .utils.vtk_image import COLORMAPS
from ..tools import misc_tools
from .slicing import SlicingWidget
from ..widgets.plot_options import PlotOptionsWidget
from ..widgets_small.indentation import WidgetIndentation
from .utils.binary_threshold_widget import ThresholdOptionsWidget
from .. import widgets_list
from .. import shared
from .. import consts


class MenuOptionsWidget(QtWidgets.QWidget):
    """Options of the selected layer: its position, and the boxes of its
    type (AFM, image, substrate, isosurface)."""

    def __init__(self):
        super().__init__()

        self.VL = QtWidgets.QVBoxLayout()
        self.VL.setContentsMargins(0, 0, 0, 0)

        self.tab_options = WidgetOptions(self)
        self.tab_afm = WidgetAFM(self)
        self.tab_Optics = None
        self.tab_glass = None
        self.tab_isosurface = None

        # Boxes of the type of the selected layer
        self.VL_type = QtWidgets.QVBoxLayout()
        self.VL_type.setContentsMargins(0, 0, 0, 0)
        self.VL_type.addWidget(self.tab_afm)
        self.type_box = None

        # The position box goes in the left column, below the list of the
        # layers (see MenuLayerListWidget.add_position_box)
        self.VL.addLayout(self.VL_type)
        self.VL.addStretch(1)
        self.setLayout(self.VL)
        self.setFixedWidth(300)

    def update_widget(self):
        """Update the GUI: boxes of the type of the selected layer."""
        row = widgets_list.widget_vtk.get_current_layer()
        layer_type = shared.layer_list[row].type

        # Remove the box of the previous layer type
        if self.type_box is not None:
            self.type_box.setParent(None)
            self.type_box.deleteLater()
            self.type_box = None

        # The AFM boxes stay (for the AFM layer of the selected layer, see
        # WidgetAFM.current), the box of the other layer types below
        if any(layer.type == "afm" for layer in shared.layer_list):
            self.tab_afm.update_widget()
        if layer_type != "afm":
            if layer_type in ("single_image", "stack"):
                self.tab_Optics = WidgetOptics(self)
                widget, title = self.tab_Optics, "Image"
            elif layer_type == "glass":
                self.tab_glass = WidgetGlass(self)
                widget, title = self.tab_glass, "Substrate"
            else:
                self.tab_isosurface = WidgetIsoSurface(self)
                widget, title = self.tab_isosurface, "Isosurface"
            self.type_box = QtWidgets.QGroupBox(title)
            layout = QtWidgets.QVBoxLayout()
            layout.setContentsMargins(4, 4, 4, 4)
            layout.addWidget(widget)
            self.type_box.setLayout(layout)
            self.VL_type.addWidget(self.type_box)
            widget.update_widget()

        self.tab_options.update_widget()


class WidgetOptions(QtWidgets.QGroupBox):
    """Position of the selected layer (common to all the layers): mouse
    action, position, angle, opacity and z scale."""

    def __init__(self, parent):
        super().__init__("Position")

        self.parent = parent

        GL = QtWidgets.QGridLayout()
        GL.setContentsMargins(6, 6, 6, 6)
        GL.setVerticalSpacing(4)

        # Mouse action
        HL_actions = QtWidgets.QHBoxLayout()
        self.BG_action = PYAFButtonGroup(self, "RBT_action")
        self.RBT_camera = QtWidgets.QRadioButton("Camera")
        self.RBT_translate = QtWidgets.QRadioButton("Translate")
        self.RBT_rotate = QtWidgets.QRadioButton("Rotate")
        self.RBT_camera.setChecked(True)
        for index, button in enumerate((self.RBT_camera, self.RBT_translate,
                                        self.RBT_rotate)):
            self.BG_action.addButton(button, index)
            HL_actions.addWidget(button)
        HL_actions.addStretch(1)
        self.RBT_camera.setToolTip("The mouse moves the camera")
        self.RBT_translate.setToolTip(
            "The mouse moves the selected layer (view from the top)")
        self.RBT_rotate.setToolTip("The mouse rotates the selected layer")

        # Position [um] and angle [degrees], two per row
        self.IN_pos_x = PYAFInput(self, "pos_x", "X [\u03bcm]", width=60)
        self.IN_pos_y = PYAFInput(self, "pos_y", "Y [\u03bcm]", width=60)
        self.IN_pos_z = PYAFInput(self, "pos_z", "Z [\u03bcm]", width=60)
        self.IN_angle = PYAFInput(self, "angle", "Angle [\u00b0]", width=60)
        for field in (self.IN_pos_x, self.IN_pos_y, self.IN_pos_z,
                      self.IN_angle):
            field.input.setValidator(misc_tools.validator("F"))

        self.IN_opacity = PYAFInput(self, "opacity", "Opacity [%]", width=60)
        self.IN_opacity.input.setValidator(
            misc_tools.validator("RI", bottom=0, top=100))
        self.CB_display_zscale = PYAFCheckBox(self, "disp_zscale",
                                              "Z scale")

        GL.addLayout(HL_actions, 0, 0, 1, 2)
        GL.addWidget(self.IN_pos_x, 1, 0)
        GL.addWidget(self.IN_pos_y, 1, 1)
        GL.addWidget(self.IN_pos_z, 2, 0)
        GL.addWidget(self.IN_angle, 2, 1)
        GL.addWidget(self.IN_opacity, 3, 0)
        GL.addWidget(self.CB_display_zscale, 3, 1)
        self.setLayout(GL)

        self.show_position(shared.layer_list[0])

    def show_position(self, layer):
        """Position of a layer in the inputs (2 decimals, angle 1)."""
        self.IN_pos_x.changeValue("%.2f" % (layer.pos_x / 1000.0))
        self.IN_pos_y.changeValue("%.2f" % (layer.pos_y / 1000.0))
        self.IN_pos_z.changeValue("%.2f" % (layer.pos_z / 1000.0))
        self.IN_angle.changeValue("%.1f" % layer.angle)

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

        self.show_position(layer)
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

            if round(angle, 1) == round(layer.angle, 1):
                # Moves smaller than the displayed precision are ignored
                dx, dy, dz = [new - old if abs(new - old) >= 10 else 0
                              for new, old in ((newpos_x, layer.pos_x),
                                               (newpos_y, layer.pos_y),
                                               (newpos_z, layer.pos_z))]
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


# LUTs of the images from black to a color (RGB); then the colormaps
# (see vtk_image) and the custom colors
LUT_COLORS = {"Gray": (1.0, 1.0, 1.0), "Red": (1.0, 0.0, 0.0),
              "Green": (0.0, 1.0, 0.0), "Blue": (0.0, 0.0, 1.0),
              "Cyan": (0.0, 1.0, 1.0), "Magenta": (1.0, 0.0, 1.0),
              "Yellow": (1.0, 1.0, 0.0)}
LUTS = ["None"] + list(LUT_COLORS) + list(COLORMAPS) + ["Custom"]


def color_icon(color):
    """Square icon of a color (RGB(A) from 0 to 1)."""
    pixmap = QtGui.QPixmap(14, 14)
    pixmap.fill(QtGui.QColor.fromRgbF(*color[:3]))
    return QtGui.QIcon(pixmap)


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
            # Sizes of the image and of its pixels [um], x and y (and z for
            # the stacks) on a line
            GL = QtWidgets.QGridLayout()
            GL.setContentsMargins(0, 0, 0, 0)
            GL.setVerticalSpacing(4)
            GL.addWidget(QtWidgets.QLabel("Size [\u03bcm]"), 0, 0)
            GL.addWidget(QtWidgets.QLabel("Pixel [\u03bcm]"), 1, 0)
            self.IN_x_size = PYAFInput(self, "size_x", "X", width=55)
            self.IN_y_size = PYAFInput(self, "size_y", "Y", width=55)
            self.IN_z_size = PYAFInput(self, "size_z", "Z", width=55)
            self.IN_x_res = PYAFInput(self, "res_x", "X", width=55)
            self.IN_y_res = PYAFInput(self, "res_y", "Y", width=55)
            self.IN_z_res = PYAFInput(self, "res_z", "Z", width=55)
            fields = [[self.IN_x_size, self.IN_y_size],
                      [self.IN_x_res, self.IN_y_res]]
            if self.type == "stack":
                fields[0].append(self.IN_z_size)
                fields[1].append(self.IN_z_res)
            for row, line in enumerate(fields):
                for col, field in enumerate(line):
                    field.input.setValidator(misc_tools.validator("UF"))
                    GL.addWidget(field, row, col + 1)
            self.IN_x_res.setToolTip("Pixel width")
            self.IN_y_res.setToolTip("Pixel height")
            GL.setColumnStretch(len(fields[0]) + 1, 1)
            self.VL.addLayout(GL)

            self.CB_bind = PYAFCheckBox(self, "bind", "Maintain Aspect Ratio")
            self.VL.addWidget(self.CB_bind)
            HL = QtWidgets.QHBoxLayout()
            self.BT_apply = PYAFButton(self, "apply", "Apply changes", 125)
            self.BT_reset = PYAFButton(self, "reset", "Reset", 125)
            self.BT_reset.setToolTip("Sizes of the file")
            HL.addWidget(self.BT_apply)
            HL.addWidget(self.BT_reset)
            HL.addStretch(1)
            self.VL.addLayout(HL)

            # Lookup table (LUT): colors of the gray levels, the background
            # (dark, or bright if inverted) transparent up to the threshold
            self.list_LUT = PYAFComboBox(self, "lut")
            for name in LUTS:
                self.list_LUT.addItem(name)
            self.list_LUT.fit_to_contents()
            self.list_LUT.setToolTip(
                "Colors of the gray levels of the image (None: original "
                "image). Custom: from the min to the max color.")
            self.CB_invert_LUT = PYAFCheckBox(self, "invert_lut", "Invert")
            self.CB_invert_LUT.setToolTip(
                "Colors reversed (bright pixels in the min color); the "
                "threshold makes the bright pixels transparent (e.g. bright "
                "background)")
            HL = QtWidgets.QHBoxLayout()
            HL.addWidget(QtWidgets.QLabel("LUT"))
            HL.addWidget(self.list_LUT)
            HL.addWidget(self.CB_invert_LUT)
            HL.addStretch(1)
            self.VL.addLayout(HL)
            HL = QtWidgets.QHBoxLayout()
            self.BT_LUT_min_col = PYAFButton(self, "lut_min_color",
                                             "Min color", 125)
            self.BT_LUT_max_col = PYAFButton(self, "lut_max_color",
                                             "Max color", 125)
            HL.addWidget(self.BT_LUT_min_col)
            HL.addWidget(self.BT_LUT_max_col)
            HL.addStretch(1)
            self.VL.addLayout(HL)
            label = "Threshold (0-255)"
            self.IN_threshold = PYAFInput(self, "input_threshold", label,
                                          width=55)
            self.IN_threshold.input.setValidator(
                misc_tools.validator("RI", bottom=0, top=255))
            self.IN_threshold.setToolTip(
                "Gray levels of the background made transparent (dark "
                "background, bright if inverted)")
            self.VL.addWidget(self.IN_threshold)

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

        label = "Gaussian smoothing \u03c3 [px]"
        self.IN_smooth = PYAFInput(self, "smooth", label, width=55)
        self.IN_smooth.input.setValidator(misc_tools.validator("UF"))
        self.IN_smooth.setToolTip(
            "Smooths the image (Gaussian filter of standard deviation "
            "\u03c3, in pixels; 0: no smoothing)")

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

        elif button == "apply":
            layer.update_actor()
            self.update_widget()

        elif button in ("lut_min_color", "lut_max_color"):
            # Attributes lut_min_color and lut_max_color of the layer
            color = misc_tools.ask_user_for_color(getattr(layer, button))
            if color is not False:
                setattr(layer, button, color)
                layer.update_actor()
                self.update_widget()

        elif button == "binary_tresh":
            if widgets_list.widget_vtk_binary_threshold is None:
                ThresholdOptionsWidget(self)
                widgets_list.widget_vtk_binary_threshold.show()
            else:
                widgets_list.widget_vtk_binary_threshold.activateWindow()
                widgets_list.widget_vtk_binary_threshold.raise_()

    def list_updated(self, name):
        """A LUT was chosen: colors from black to a color, colormap, or
        custom colors."""
        index = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[index]
        lut = self.list_LUT.currentText()
        layer.use_lut = lut != "None"
        if layer.use_lut:
            layer.lut_name = lut
            if lut in LUT_COLORS:
                layer.lut_min_color = [0.0, 0.0, 0.0, 1.0]
                layer.lut_max_color = list(LUT_COLORS[lut]) + [1.0]
        layer.update_actor()
        self.update_widget()

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

        elif name == "invert_lut":
            layer.lut_invert = self.CB_invert_LUT.isChecked()
            layer.update_actor()

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
            # Sizes with 2 decimals, pixels with 3 significant digits [um]
            self.IN_x_size.changeValue("%.2f" % (layer.size_x / 1000.0))
            self.IN_y_size.changeValue("%.2f" % (layer.size_y / 1000.0))
            self.IN_x_res.changeValue("%.3g" % (layer.res_x / 1000.0))
            self.IN_y_res.changeValue("%.3g" % (layer.res_y / 1000.0))
            if layer.type == "stack":
                self.IN_z_size.changeValue("%.2f" % (layer.size_z / 1000.0))
                self.IN_z_res.changeValue("%.3g" % (layer.res_z / 1000.0))
            self.CB_bind.setChecked(layer.bind)

            self.list_LUT.setCurrentText(
                layer.lut_name if layer.use_lut else "None")
            self.CB_invert_LUT.setChecked(layer.lut_invert)
            self.IN_threshold.changeValue(layer.lut_threshold)
            for button, color in ((self.BT_LUT_min_col, layer.lut_min_color),
                                  (self.BT_LUT_max_col, layer.lut_max_color)):
                button.setIcon(color_icon(color))
            for widget in (self.CB_invert_LUT, self.IN_threshold):
                widget.setEnabled(layer.use_lut)
            # Colors chosen only for the custom LUT
            custom = layer.use_lut and layer.lut_name == "Custom"
            self.BT_LUT_min_col.setEnabled(custom)
            self.BT_LUT_max_col.setEnabled(custom)

            self.IN_smooth.changeValue(layer.smooth_gaussian)

        if layer.type == "stack":
            self.IN_slide.changeValue(layer.slide_pos)
            self.IN_slide_op_single.changeValue(layer.single_opacity)
            self.IN_slide_op_others.changeValue(layer.others_opacity)
            self.slider.setValue(layer.slide_pos)

        self.CB_flip_x.setChecked(layer.flip_x)
        self.CB_flip_y.setChecked(layer.flip_y)


class WidgetGlass(QtWidgets.QWidget):
    """Widget with options for the substrate layers."""

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

    def checkbox_clicked(self, _name=None):
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
    """Boxes of the AFM layers: Display (map, height, 3D or flat view,
    smoothing, slicing) and Colors (color options, apply to all,
    indentation segment)."""

    def __init__(self, parent):
        super().__init__()

        self.parent = parent

        VL = QtWidgets.QVBoxLayout()
        VL.setContentsMargins(0, 0, 0, 0)

        # --- Display --------------------------------------------------------
        box_display = QtWidgets.QGroupBox("Display")
        GL = QtWidgets.QGridLayout()
        GL.setContentsMargins(6, 6, 6, 6)
        GL.setVerticalSpacing(4)

        self.list_maps = PYAFComboBox(self, "map")
        self.list_maps.setToolTip(
            "Map drawn on the surface (the computed maps of the results "
            "tab).\nThe elasticity computed with the Tomography is drawn as "
            "a volume (one value per indentation segment).")
        for combo in ("list_maps",):
            getattr(self, combo).setSizeAdjustPolicy(
                QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
            getattr(self, combo).setMinimumContentsLength(14)
        self.list_height = PYAFComboBox(self, "height")
        self.list_height.setToolTip(
            "Height of the surface: zero-force topography (computed with "
            "the elasticity) or piezo height")

        HL_mode = QtWidgets.QHBoxLayout()
        self.GRP_mode = PYAFButtonGroup(self, "RBT_mode")
        self.RBT_3D = QtWidgets.QRadioButton("3D")
        self.RBT_2D = QtWidgets.QRadioButton("Flat")
        self.RBT_3D.setToolTip("Surface at its height")
        self.RBT_2D.setToolTip("Flat surface (at the Z position)")
        self.GRP_mode.addButton(self.RBT_3D, 0)
        self.GRP_mode.addButton(self.RBT_2D, 1)
        HL_mode.addWidget(self.RBT_3D)
        HL_mode.addWidget(self.RBT_2D)

        self.CB_smoothing = PYAFCheckBox(self, "smoothing", "Smoothing")
        self.CB_smoothing.setToolTip(
            "Colors interpolated between the pixels and smoother surface "
            "(subdivision of the triangles), slower.\nOff: one color per "
            "pixel, as in the maps.")
        # Smoothing, then 3D or flat, on a line
        HL_mode.insertWidget(0, self.CB_smoothing)
        HL_mode.addStretch(1)
        self.BT_slicing = PYAFButton(self, "slicing", "Slicing")
        self.BT_slicing.setToolTip("Part of the map displayed (in x and y)")

        GL.addWidget(QtWidgets.QLabel("Map"), 0, 0)
        GL.addWidget(self.list_maps, 0, 1)
        GL.addWidget(QtWidgets.QLabel("Height"), 1, 0)
        GL.addWidget(self.list_height, 1, 1)
        GL.addLayout(HL_mode, 2, 0, 1, 2)
        GL.addWidget(self.BT_slicing, 3, 0, 1, 2)
        GL.setColumnStretch(1, 1)
        box_display.setLayout(GL)

        # --- Colors ---------------------------------------------------------
        box_colors = QtWidgets.QGroupBox("Colors")
        VL_colors = QtWidgets.QVBoxLayout()
        VL_colors.setContentsMargins(6, 6, 6, 6)
        self.BT_colors = PYAFButton(self, "colors", "Color options")
        self.BT_colors.setToolTip(
            "Color scale of the map (same options as the map of the "
            "results tab)")
        self.CB_apply_to_all = PYAFCheckBox(self, "apply_to_all",
                                            "Apply to all")
        self.CB_apply_to_all.setToolTip(
            "The map, the colors and the indentation segment are changed "
            "for all the files of the 3D view")
        # Indentation segment of the elasticity (Tomography)
        self.afm_indentation_widget = WidgetIndentation(self)
        HL_colors = QtWidgets.QHBoxLayout()
        HL_colors.addWidget(self.BT_colors)
        HL_colors.addWidget(self.CB_apply_to_all)
        HL_colors.addStretch(1)
        VL_colors.addLayout(HL_colors)
        VL_colors.addWidget(self.afm_indentation_widget)
        box_colors.setLayout(VL_colors)

        VL.addWidget(box_display)
        VL.addWidget(box_colors)
        self.setLayout(VL)

    def current(self):
        """AFM layer of the options and its file: the selected layer, or the
        AFM layer of the file of the selected layer (substrate), or the first
        AFM layer (images)."""
        layer_id = widgets_list.widget_vtk.get_current_layer()
        layer = shared.layer_list[layer_id]
        if layer.type != "afm":
            afm = [lay for lay in shared.layer_list if lay.type == "afm"]
            same = [lay for lay in afm
                    if lay.afm_id == getattr(layer, "afm_id", None)]
            layer = (same or afm)[0]
        return layer, shared.exp.list[layer.afm_id]

    def afm_layers(self, layer):
        """Layers changed by an option: all the AFM layers with apply to all,
        else the selected one."""
        if shared.exp.apply_to_all_data:
            return [lay for lay in shared.layer_list if lay.type == "afm"]
        return [layer]

    def show_map(self, data):
        """Displays the map of the 3D view in the results tab (the color
        options are the ones of the displayed map)."""
        data.meshgrid_type = data.opengl_map_type
        if data is shared.exp.current_data:
            results = widgets_list.widget_results
            button = results.GRP_mesh_type.button(
                misc_tools.get_meshgrid_type_by_id(data.opengl_map_type))
            button.setChecked(True)
            results.button_clicked("BTG_meshgrid_type")
        if widgets_list.widget_meshgrid_options is not None:
            widgets_list.widget_meshgrid_options.update_widget()

    def list_updated(self, name):
        """A list was changed: map or height."""
        layer, data = self.current()
        if name == "map":
            map_type = self.list_maps.currentData()
            for lay in self.afm_layers(layer):
                dt = shared.exp.list[lay.afm_id]
                if map_type in map_values.available_maps(dt):
                    dt.opengl_map_type = map_type
                    self.show_map(dt)
                    lay.update_tomography()
        elif name == "height":
            height = self.list_height.currentData()
            for lay in self.afm_layers(layer):
                dt = shared.exp.list[lay.afm_id]
                dt.opengl_height = height
                lay.update_tomography()
        self.update_widget()

    def button_clicked(self, button):
        """Called when a button is clicked."""
        layer, data = self.current()

        if button == "RBT_mode":
            flat = self.GRP_mode.checkedId() == 1
            for lay in self.afm_layers(layer):
                shared.exp.list[lay.afm_id].opengl_afm_flat_view = flat
                lay.update_tomography()
            self.update_widget()

        elif button == "colors":
            # Color options window (as for the maps), for the file of the
            # selected layer and its map
            self.show_map(data)
            if widgets_list.widget_meshgrid_options is None:
                PlotOptionsWidget(self, "vtk")
                widgets_list.widget_meshgrid_options.setWindowTitle(
                    "Color options")
                widgets_list.widget_meshgrid_options.resize(250, 250)
                widgets_list.widget_meshgrid_options.show()
            widgets_list.widget_meshgrid_options.update_widget(
                data_id=layer.afm_id)
            widgets_list.widget_meshgrid_options.activateWindow()
            widgets_list.widget_meshgrid_options.raise_()

        elif button == "slicing":
            if widgets_list.widget_slicing is None:
                # Create new widget
                SlicingWidget(self)
                widgets_list.widget_slicing.setWindowTitle("3D Slicing")
                widgets_list.widget_slicing.show()
            else:
                # Bring to front
                widgets_list.widget_slicing.activateWindow()
                widgets_list.widget_slicing.raise_()

    def checkbox_clicked(self, name):
        """Called when a checkbox is clicked."""
        layer, data = self.current()

        if name == "smoothing":
            for lay in self.afm_layers(layer):
                dt = shared.exp.list[lay.afm_id]
                dt.vtk_smoothing_type = "loop" \
                    if self.CB_smoothing.isChecked() else None
                # Colors per pixel or interpolated: new surface
                lay.update_tomography()

        elif name == "apply_to_all":
            shared.exp.apply_to_all_data = self.CB_apply_to_all.isChecked()
            widgets_list.widget_results.update_GUI("apply_to_all")
            self.afm_indentation_widget.update_widget(data_id=layer.afm_id)

    def update_widget(self):
        """Updates the boxes for the selected AFM layer."""
        layer, data = self.current()
        if data.opengl_map_type not in map_values.available_maps(data):
            data.opengl_map_type = "stiffness"

        # Maps computed for the file
        self.list_maps.blockSignals(True)
        self.list_maps.clear()
        for map_type in map_values.available_maps(data):
            self.list_maps.addItem(map_values.MAPS[map_type][0], map_type)
        self.list_maps.setCurrentIndex(
            self.list_maps.findData(data.opengl_map_type))
        self.list_maps.blockSignals(False)

        # Heights: topography if computed, piezo height
        self.list_height.blockSignals(True)
        self.list_height.clear()
        if data.stiffness_calculated:
            self.list_height.addItem("Topography (zero force)", "topo")
        self.list_height.addItem("Piezo height", "piezo")
        index = self.list_height.findData(data.opengl_height)
        self.list_height.setCurrentIndex(max(index, 0))
        self.list_height.blockSignals(False)

        flat = data.opengl_afm_flat_view
        (self.RBT_2D if flat else self.RBT_3D).setChecked(True)
        self.CB_smoothing.setChecked(data.vtk_smoothing_type is not None)
        self.CB_apply_to_all.setChecked(shared.exp.apply_to_all_data)

        # Indentation segment: only for the elasticity volume (Tomography)
        self.afm_indentation_widget.update_widget(data_id=layer.afm_id)
        self.afm_indentation_widget.setEnabled(vtk_afm.is_volume(data))


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
