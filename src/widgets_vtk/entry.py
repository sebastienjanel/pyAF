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

"""Small widget asking the user which AFM tomographies he wants to display in VTK."""

from PyQt5 import QtWidgets
from ..tools.gui_tools import PYAFButton as PYAFButton
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFCheckBox
from .main import OpenglWidget
from .. import shared
from .. import widgets_list
from .. import consts


class VTKEntryWidget(PYAFWidget):
    """VTKEntryWidget asks the user which tomography to display."""

    def __init__(self, parent, mode=None):
        super().__init__(parent, "widget_vtk_entry")

        self.parent = parent
        self.CB_list = []
        self.mode = mode

        VL = QtWidgets.QVBoxLayout()

        # Label
        label = QtWidgets.QLabel()
        label.setText("Choose the files you want to display in 3D :")

        # Display a list of checkboxes with the names in a scroll area
        scroll_area = QtWidgets.QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setMinimumSize(150, 250)
        widget = QtWidgets.QWidget()
        VL_list = QtWidgets.QVBoxLayout()

        for i in range(len(shared.exp.list)):
            data = shared.exp.list[i]
            if data.stiffness_calculated:
                CB = PYAFCheckBox(self, str(data.unique_id), data.filename)
                CB.setChecked(False)
                self.CB_list.append(CB)
                VL_list.addWidget(CB)

        VL_list.addStretch(1)
        widget.setLayout(VL_list)

        scroll_area.setWidget(widget)

        HL_buttons = QtWidgets.QHBoxLayout()
        BT_continue = PYAFButton(self, "continue", "OK")
        BT_cancel = PYAFButton(self, "cancel", "Cancel")
        HL_buttons.addWidget(BT_cancel)
        HL_buttons.addStretch(1)
        HL_buttons.addWidget(BT_continue)

        VL.addWidget(label)
        VL.addWidget(scroll_area)
        VL.addLayout(HL_buttons)

        self.setLayout(VL)

        if consts.AUTO3D or consts.UNIT_TESTING:
            self.CB_list[0].setChecked(True)
            self.button_clicked("continue")

    def button_clicked(self, what):
        """Called when a button is clicked."""
        if what == "continue":
            # Make a list of result id's which will be loaded.
            load_list = []
            for CB in self.CB_list:
                if CB.isChecked():
                    load_list.append(int(CB.name))

            self.close()

            if load_list == []:
                # Close and do nothing if no file was selected.
                return False

            if self.mode is None:
                # When we come from the results tab and the user has clicked
                # on the 3D button, open the 3D widget
                OpenglWidget(self.parent, load_list)

                # The show() call is inside the VTK widget

            elif self.mode == "add_more":
                # In case the user used the menu to add more tomographies to
                # the 3D widget.

                widgets_list.widget_vtk.canvas.add_more_afm_layers(load_list)

        elif what == "cancel":
            self.close()

    def checkbox_clicked(self, what):
        """Called when a checkbox is clicked.

        Is not needed here as the values are only saved when the OK button is
        pushed.
        """
        pass
