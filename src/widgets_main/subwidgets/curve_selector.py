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

"""Widget used to go through the curves."""

from PyQt5 import QtGui, QtWidgets

from ... import shared
from ... import widgets_list
from ...tools.gui_tools import PYAFButton, PYAFInput
from ...tools.PYAFWidget import PYAFWidget


class ResultsSelectorWidget(PYAFWidget):
    """The curve selector allows to go through the curves.

    There is a random button to move randomly through the dataset, two input
    fields to go to a specific curve, and buttons to go back and forth in the
    last seen curves.
    """

    def __init__(self, parent, name):
        super().__init__(parent, name)
        self.IN_curve_x = PYAFInput(
            self, "input_chose_curve_x", ",")
        self.IN_curve_y = PYAFInput(self, "input_chose_curve_y", "")
        self.IN_curve_x.input.setValidator(validator("x"))
        self.IN_curve_y.input.setValidator(validator("y"))
        self.IN_curve_x.changeValue(shared.exp.meshgrid_click_xpos + 1)
        self.IN_curve_y.changeValue(shared.exp.meshgrid_click_ypos + 1)

        self.button_rand = PYAFButton(self, "button_rand", "Random curve")

        self.BT_back = PYAFButton(self, "back", "<", size=60)
        self.BT_forward = PYAFButton(self, "forward", ">", size=60)

        self.label_info = QtWidgets.QLabel()

        self.HL_select_curve = QtWidgets.QHBoxLayout()
        self.HL_select_curve.addStretch(1)
        self.HL_select_curve.addWidget(self.IN_curve_x)
        self.HL_select_curve.addWidget(self.IN_curve_y)
        self.HL_select_curve.addWidget(self.button_rand)
        self.HL_select_curve.addWidget(self.BT_back)
        self.HL_select_curve.addWidget(self.BT_forward)
        self.HL_select_curve.addStretch(1)

        self.setLayout(self.HL_select_curve)

        self.update_widget()

    def update_widget(self):
        """Update the widget."""
        data = shared.exp.list[shared.exp.id_selected]

        # Get the position in user coordinates
        xpos = shared.exp.meshgrid_click_xpos + 1
        ypos = shared.exp.meshgrid_click_ypos + 1

        # Update the inputs
        self.IN_curve_x.changeValue(xpos)
        self.IN_curve_y.changeValue(ypos)
        self.IN_curve_x.changeValue(xpos)
        self.IN_curve_y.changeValue(ypos)

        length = len(shared.exp.last_ten_curves)
        pos = shared.exp.pos_in_last_ten_curves

        # Update the back and forward buttons
        if pos < length + 1 and pos != 0:
            self.BT_back.setEnabled(True)
        else:
            self.BT_back.setEnabled(False)

        if shared.exp.pos_in_last_ten_curves + 1 == length:
            self.BT_forward.setEnabled(False)
        else:
            self.BT_forward.setEnabled(True)

        # Display a message for Bruker files where there are more pixels
        # than force curves. The reduction of the number of pixels is done
        # during loading.
        tp = (data.file_type == "Nanoscope (Force Volume)" or
              data.file_type == "Nanoscope (Peak Force)")
        if tp:
            pixels = data.piezo_image_samps_per_line
            force_curves = data.nbr_pixels_x
            div = pixels / force_curves
            if div != 1:
                text = ("Original number of pixels : "
                        + str(pixels) + "x" + str(pixels) + ". Reduced to : "
                        + str(force_curves) + "x" + str(force_curves) + ".")

                # Display the text
                self.label_info.setText(text)

            else:
                self.label_info.setText("")

        else:
            self.label_info.setText("")

    def input_updated(self, field):
        """Called when an input field is updated."""
        if field == "input_chose_curve_x" or field == "input_chose_curve_y":
            x = self.IN_curve_x.get_int_value()
            y = self.IN_curve_y.get_int_value()
            widgets_list.widget_main.change_curve(x, y)

    def button_clicked(self, button):
        """Called when a button is clicked."""
        if button == "button_rand":
            widgets_list.widget_main.change_curve(0, 0, "rand")

        elif button == "back":
            ls = shared.exp.last_ten_curves
            shared.exp.pos_in_last_ten_curves -= 1
            posx = ls[shared.exp.pos_in_last_ten_curves][0]
            posy = ls[shared.exp.pos_in_last_ten_curves][1]

            widgets_list.widget_main.change_curve(posx, posy,
                                                  save_in_last_ten=False)
            self.update_widget()

        elif button == "forward":
            ls = shared.exp.last_ten_curves
            shared.exp.pos_in_last_ten_curves += 1
            posx = ls[shared.exp.pos_in_last_ten_curves][0]
            posy = ls[shared.exp.pos_in_last_ten_curves][1]
            widgets_list.widget_main.change_curve(posx, posy,
                                                  save_in_last_ten=False)
            self.update_widget()


def validator(nb):
    """Validate the position.

    Check if the inputed value is an int and lies in the right boundaries.
    """
    data = shared.exp.current_data

    if nb == "x":
        x = int(data.nbr_pixels_x)
        val = QtGui.QIntValidator(1, x)
    if nb == "y":
        y = int(data.nbr_pixels_y)
        val = QtGui.QIntValidator(1, y)

    return val
