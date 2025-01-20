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

"""Widget used to rename the groups."""

from PyQt5 import QtWidgets
from ..tools.gui_tools import PYAFButton
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFInput
from .. import shared
from .. import widgets_list


class RenameConditionsWidget(PYAFWidget):
    """Shows a list of groups and allows them to be renamed."""

    def __init__(self, parent):
        name = "widget_rename_conditions"
        super().__init__(parent, name)

        self.input_list = []

        VL = QtWidgets.QVBoxLayout()

        for i in range(1, len(shared.exp.conditions_list), 1):
            field = PYAFInput(self, "input", "Name", width=200)
            field.changeValue(shared.exp.conditions_list[i].name)
            self.input_list.append(field)
            VL.addWidget(field)

        self.button_ok = PYAFButton(self, "button_ok", "OK")

        VL.addWidget(self.button_ok)
        VL.addStretch(1)
        self.setLayout(VL)

    def input_updated(self, _):
        """Called when an input field is updated.

        For this widget, do nothing in this case.
        """
        pass

    def button_clicked(self, button):
        """Called when a button is clicked."""
        if button == "button_ok":
            # Save the new names
            for i in range(1, len(shared.exp.conditions_list), 1):
                name = str(self.input_list[i - 1].get_str_value())
                shared.exp.conditions_list[i].name = name

            # Update the names in the histogram's list
            widgets_list.widget_results_single.reset_table()
            widgets_list.widget_results_single.update_widget()

            tw = widgets_list.widget_results_experiment.tableWidget

            # Update the names in the histogram's group's list
            for i in range(1, len(shared.exp.conditions_list), 1):
                name = shared.exp.conditions_list[i].name
                tw.item(i - 1, 1).setText(name)

            # Update the groups plot
            widgets_list.widget_results_experiment.update_MPL()

            # Close widget
            self.close()
