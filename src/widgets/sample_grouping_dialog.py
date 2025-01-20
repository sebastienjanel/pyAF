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

"""Advanced options for statistical analysis, called from Experiment Plots tab."""

from .. import shared
from PyQt5 import QtCore, QtWidgets
from ..tools.gui_tools import PYAFButton, PYAFComboBox
from ..tools.PYAFWidget import PYAFWidget
from .. import widgets_list


class SampleGroupingDialog(PYAFWidget):
    """Widget used to add comparison rules used in the statistical analysis manually."""

    def __init__(self, parent):
        name = "widget_sample_grouping_dialog"
        super().__init__(parent, name)

        BOX_add_rule = QtWidgets.QGroupBox("Add rule")

        label_condition_a = QtWidgets.QLabel("A", self)     # Implement dynamic management based on selected test.
        label_condition_b = QtWidgets.QLabel("B", self)

        self.list_condition_a = PYAFComboBox(self, "list_condition_a")
        self.list_condition_b = PYAFComboBox(self, "list_condition_b")

        FL = QtWidgets.QFormLayout()

        for condition in shared.exp.conditions_list[1:]:
            self.list_condition_a.addItem(condition.name)
            self.list_condition_b.addItem(condition.name)

        self.button_add_another_rule = PYAFButton(self, "button_add_another_rule", "Add another rule")
        self.button_done = PYAFButton(self, "button_done", "Done")

        FL.addRow(label_condition_a, self.list_condition_a)
        FL.addRow(label_condition_b, self.list_condition_b)
        FL.addRow(self.button_add_another_rule, self.button_done)
        FL.setFieldGrowthPolicy(2)  # Let the fields expand
        BOX_add_rule.setLayout(FL)

        buttonBox = QtWidgets.QHBoxLayout()
        buttonBox.addWidget(self.button_add_another_rule)
        buttonBox.addWidget(self.button_done)

        mainLayout = QtWidgets.QVBoxLayout()
        mainLayout.addWidget(BOX_add_rule)
        mainLayout.addLayout(buttonBox)
        self.setLayout(mainLayout)

    def list_updated(self, _):
        """Called when an list field is updated.

        For this widget, do nothing in this case.
        """
        pass

    def button_clicked(self, button):
        """Called when a button is clicked."""
        if button == "button_add_another_rule":
            """Adds rule and keeps the widget open, allowing the user to add another rule"""
            self.add_rule()

        if button == "button_done":
            """Adds rule and closes the widget"""
            self.add_rule()

            # Close widget
            self.close()

    def add_rule(self):
        """Add condition comparison based on the combo boxes input"""
        index_a = self.list_condition_a.currentIndex()
        index_b = self.list_condition_b.currentIndex()
        if index_a != index_b:
            sample_pair = self.order_sample_pair(
                (index_a + 1, index_b + 1))
            widgets_list.widget_results_experiment.add_sample_group(sample_pair)

    @staticmethod
    def order_sample_pair(sample_pair):
        """Order condition pairs"""
        return tuple(
            (sample_pair[0], sample_pair[1]) if sample_pair[0] <= sample_pair[1] else (sample_pair[1], sample_pair[0]))
