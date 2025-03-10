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

"""See the values of the results for a single curve."""

from PyQt5 import QtGui, QtWidgets
from ..tools.gui_tools import PYAFButton as PYAFButton
from ..tools.gui_tools import ClearWidgetsFromLayout as ClearWidgetsFromLayout
from ..tools.PYAFWidget import PYAFWidget
from ..tools import curve_tools
from .. import shared


class DetailsWidget(PYAFWidget):
    """Widget with results."""

    def __init__(self, parent):
        super().__init__(parent, "widget_details")

        self.data = shared.exp.current_data

        self.HL = QtWidgets.QHBoxLayout()

        # Stiffness
        self.box_stiffness = QtWidgets.QGroupBox()
        self.VL_box_stiffness = QtWidgets.QVBoxLayout()
        self.VL_stiffness_labels = QtWidgets.QVBoxLayout()
        self.button_save = PYAFButton(self, "button_save", "Save as ... ")
        self.VL_box_stiffness.addLayout(self.VL_stiffness_labels)
        self.VL_box_stiffness.addWidget(self.button_save)
        self.box_stiffness.setLayout(self.VL_box_stiffness)

        # Work and rupture_force1
        self.box_WR1 = QtWidgets.QGroupBox("Work and rupture force")
        self.VL_box_WR1 = QtWidgets.QVBoxLayout()
        self.VL_WR1_labels = QtWidgets.QVBoxLayout()
        self.VL_box_WR1.addLayout(self.VL_WR1_labels)
        self.box_WR1.setLayout(self.VL_box_WR1)

        # Events
        self.box_events = QtWidgets.QGroupBox("Events")
        self.VL_box_events = QtWidgets.QVBoxLayout()
        self.VL_events_labels = QtWidgets.QVBoxLayout()
        self.VL_box_events.addLayout(self.VL_events_labels)
        self.box_events.setLayout(self.VL_box_events)

        # Finish
        self.HL.addWidget(self.box_stiffness)
        self.HL.addWidget(self.box_WR1)
        self.HL.addWidget(self.box_events)
        self.setLayout(self.HL)

        self.update_widget()

    def update_widget(self):
        """Update the GUI."""
        data = shared.exp.current_data
        i = shared.exp.meshgrid_click_xpos
        j = shared.exp.meshgrid_click_ypos

        if data.used_stiffness_model_selected == 3:
            # Slope
            self.box_stiffness.setTitle("Apparent stiffness")
        else:
            # Stiffness
            self.box_stiffness.setTitle("Elasticity")

        # --- Stiffness ------------------------------------------------
        # Clear the layout with the labels
        ClearWidgetsFromLayout(self.VL_stiffness_labels)
        if data.stiffness_calculated:
            values = shared.exp.current_data.stiffness_array[i][j]

            if len(values) > 20:
                steps = 20
            else:
                steps = len(values)

            for ind in range(steps):
                label_value = QtWidgets.QLabel()
                label_value.setFont(QtGui.QFont("defaultFamily", 12))
                label_value.setText(get_one_line(i, j, ind, values[ind]))
                self.VL_stiffness_labels.addWidget(label_value)

            if len(values) > 20:
                label_value = QtWidgets.QLabel()
                label_value.setFont(QtGui.QFont("defaultFamily", 12))
                label_value.setText("See more in save ... ")
                self.VL_stiffness_labels.addWidget(label_value)

            self.button_save.setEnabled(True)
        else:
            self.button_save.setEnabled(False)
        self.VL_stiffness_labels.addStretch(1)

        # --- Work and rupture_force1 ----------------------------------
        ClearWidgetsFromLayout(self.VL_box_WR1)
        if data.work_and_rupture_force1_calculated:
            label_work = QtWidgets.QLabel()
            label_rf1 = QtWidgets.QLabel()
            label_work.setFont(QtGui.QFont("defaultFamily", 12))
            label_rf1.setFont(QtGui.QFont("defaultFamily", 12))
            label_work.setText(str(data.work[i][j] * 1e15) + " fJ")
            value = str(data.rupture_force1[i][j] * 1e12)
            label_rf1.setText(value + " pN")
            self.VL_box_WR1.addWidget(label_work)
            self.VL_box_WR1.addWidget(label_rf1)
            self.VL_box_WR1.addStretch(1)

        # --- Events ---------------------------------------------------
        ClearWidgetsFromLayout(self.VL_events_labels)
        if data.events_calculated:
            for a in range(len(data.events_forces[i][j])):
                label_event = QtWidgets.QLabel()
                label_event.setFont(QtGui.QFont("defaultFamily", 12))
                val = str(data.events_forces[i][j][a] * 1e12)
                txt = "Force : " + val + " pN"
                label_event.setText(txt)
                self.VL_events_labels.addWidget(label_event)
            self.VL_events_labels.addStretch(1)

    def button_clicked(self, what):
        """Called when a button is clicked."""
        if what == "button_save":
            filename, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Save data", "", "Text files (*.txt)")

            if filename:
                newfile = open(str(filename), "w")

                i = shared.exp.meshgrid_click_xpos
                j = shared.exp.meshgrid_click_ypos
                values = shared.exp.current_data.stiffness_array[i][j]

                for ind in range(len(values) - 1):
                    newfile.write(get_one_line(i, j, ind, values[ind], True))

                newfile.close()


def get_one_line(i, j, ind, value, export=False):
    """Return a formatted line with values to be written."""
    data = shared.exp.current_data
    step = data.used_indentation_step
    start = data.used_indentation_start

    if export:
        sep = "\t"
        endl = "\r\n"
    else:
        sep = " "
        endl = ""

    # End value
    if data.used_indentation_stop == 0 and step == 0:
        # Get the max value (last value of the force curve)
        approach, _, _, _, _, _ = curve_tools.get_curve(
            data, [i, j], mode="details_widget")

        tr = curve_tools.get_force_curves(
            approach[0],
            approach[1],
            data.pocs_real[i][j],
            data.spring_constant)
        end = str(tr[0][-1])

    else:
        end = str(ind * step + step + start)

    # Used tomography mode or not
    if not data.used_tomography:
        start_z = str(start)
    else:
        start_z = str(ind*step + start)

    txt = str(ind) + sep
    txt += " (" + start_z + " - " + end + " nm)" + sep

    if data.used_stiffness_model_selected == 3:
        # Slope
        txt += str(value).replace('.', ',') + sep + "N/m" + endl
    else:
        # Stiffness
        txt += str(value / 1000.0).replace('.', ',') + sep + "kPa" + endl

    return txt
