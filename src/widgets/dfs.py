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

"""Display a list with results of DFS experiments."""

from PyQt5 import QtCore, QtWidgets
from ..tools import gui_tools
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFButtonGroup
from ..tools import math_tools
import math
import numpy
from ..tools.PYAFWidget import PYAFWidget
from ..widgets.progressbar import Progressbar
from .. import widgets_list
from .. import shared
from ..tools import results_sorting


class DFSWidget(PYAFWidget):
    """Class displaying a widget with DFS results and compute button.

    The widget is called from the parent widget (results.py). It displays a
    list with the results from previous DFS computations (size of the potential
    barrier (x) in nm and dissociation constant ku (s^-1).

    The compute button lets you start a new computation with the currently
    selected rows in the result table from the parent widget.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_dfs")

        self.VL = QtWidgets.QVBoxLayout()

        self.HL = QtWidgets.QHBoxLayout()
        self.GRP_buttons = PYAFButtonGroup(self, "GRP_buttons")
        self.RBT_single = QtWidgets.QRadioButton("Single results")
        self.RBT_groups = QtWidgets.QRadioButton("Grouped results")
        self.GRP_buttons.addButton(self.RBT_single, 0)
        self.GRP_buttons.addButton(self.RBT_groups, 1)
        self.HL.addStretch(1)
        self.HL.addWidget(self.RBT_single)
        self.HL.addWidget(self.RBT_groups)
        self.HL.addStretch(1)

        self.tableWidget = QtWidgets.QTableWidget(0, 3)
        policy = QtCore.Qt.ScrollBarAlwaysOff
        self.tableWidget.setHorizontalScrollBarPolicy(policy)
        behavior = QtWidgets.QAbstractItemView.SelectRows
        self.tableWidget.setSelectionBehavior(behavior)
        self.tableWidget.horizontalHeader().setHighlightSections(False)
        self.tableWidget.setHorizontalHeaderLabels(["", "x [nm]", "ku [s^-1]"])
        self.tableWidget.setFixedWidth(300)
        self.tableWidget.setColumnWidth(0, 40)
        self.tableWidget.setColumnWidth(1, 130)
        self.tableWidget.setColumnWidth(2, 130)
        self.tableWidget.verticalHeader().hide()

        self.BT_get_dfs = PYAFButton(self, "get_dfs", "Get values")
        self.BT_remove = PYAFButton(self, "remove", "Remove selected")
        self.HL_buttons = QtWidgets.QHBoxLayout()
        self.HL_buttons.addWidget(self.BT_get_dfs)
        self.HL_buttons.addWidget(self.BT_remove)

        self.VL.addLayout(self.HL)
        self.VL.addWidget(self.tableWidget)
        self.VL.addLayout(self.HL_buttons)
        self.setLayout(self.VL)

        self.update_widget()

    def input_updated(self, val):
        """Do nothing here."""
        pass

    def button_clicked(self, button):
        """Called when the Get values button is clicked.

        Will compute the two values (ku and x) from the DFS experiment.
        """
        # Boltzmann constant
        KB = 1.3806503 * 1e-23

        if button == "get_dfs":
            Progressbar()
            widgets_list.widget_progressbar.set_label("Fetching data")

            medians_force = []
            medians_lr = []

            if shared.exp.dfs_mode == "single":
                wl = widgets_list.widget_results_single
                table = wl.tableWidget.selectionModel().selectedRows()
            elif shared.exp.dfs_mode == "groups":
                wl = widgets_list.widget_results_groups
                table = wl.tableWidget.selectionModel().selectedRows()

            # Get the range of the progressbar
            progressbar_range = 0
            for item in table:
                row_index = item.row()
                if shared.exp.dfs_mode == "single":
                    data_id = shared.exp.results_list[row_index].data_id
                    lr = shared.exp.list[data_id].loading_rates
                    progressbar_range += len(lr) * 2
                elif shared.exp.dfs_mode == "groups":
                    group_id = row_index + 1
                    for i in range(len(shared.exp.results_list)):
                        if group_id == shared.exp.results_list[i].group:
                            data_id = shared.exp.results_list[i].data_id
                            lr = shared.exp.list[data_id].loading_rates
                            progressbar_range += len(lr) * 2

            widgets_list.widget_progressbar.set_range(0, progressbar_range)

            results_util = results_sorting.GetResults()

            if shared.exp.dfs_mode == "single":
                for item in table:
                    row_index = item.row()
                    result_id = shared.exp.results_list[row_index].result_id
                    res = results_util.get_results_single(
                        result_id=result_id,
                        force_type="loading_rates")
                    array, values = res[0], res[1]
                    median_lr = values[1]
                    median_force = numpy.median(array[0])
                    medians_force.append(median_force * 1e-12)
                    medians_lr.append(median_lr * 1e-12)

            if shared.exp.dfs_mode == "groups":
                res = results_util.get_results_groups(
                    force_type="loading_rates")
                groups, values = res[0], res[1]

                for item in table:
                    row_index = item.row()
                    group_id = row_index + 1
                    if groups[group_id] != [[], []]:
                        val = groups[group_id][0]
                        medians_force.append(numpy.median(val) * 1e-12)
                        medians_lr.append(values[group_id][1] * 1e-12)

            # F = KBT/x * ln ( r * x / KBTku )

            # F = KBT/x * [ ln (r) + ln ( x/KBTku) ]
            # F = KBT/x * ln (r) + KBT / x [ ln ( x/KBTku) ]
            # F = a * ln(r) + b
            # a = KBT/x and b = KBT / x [ ln ( x/KBTku) ]

            # a = KBT/x --> x = KBT/a
            #
            # b = KBT / x [ ln ( x/KBTku) ]
            # exp(bx/KBT) =  x/KBTku
            # ku = x/KBT*exp(-bx/KBT)

            T = shared.exp.list[0].temperature

            medians_force = numpy.sort(medians_force, -1)
            medians_lr = numpy.sort(medians_lr, -1)

            if medians_lr:
                # Get a fit
                coeffs, _ = \
                    math_tools.fit_linear(numpy.log(medians_lr), medians_force)

                x = (T * KB) / coeffs[0]
                ku = (x / (KB * T)) * math.exp((-coeffs[1] * x) / (KB * T))

                # a, b, minx, maxx, displayed = True
                coeffs, _ = math_tools.fit_linear(medians_lr, medians_force)
                params = [coeffs[0], coeffs[1] * 1e12,
                          numpy.amin(medians_lr) * 1e12,
                          numpy.amax(medians_lr) * 1e12,
                          True]
                result = [x, ku]

                if shared.exp.dfs_mode == "single":
                    shared.exp.list_dfs_single_fits.append(params)
                    shared.exp.list_dfs_single_results.append(result)
                elif shared.exp.dfs_mode == "groups":
                    shared.exp.list_dfs_groups_fits.append(params)
                    shared.exp.list_dfs_groups_results.append(result)

                self.update_widget()

                # Update plot (display fits)
                if shared.exp.dfs_mode == "single":
                    widgets_list.widget_results_single.update_MPL()
                elif shared.exp.dfs_mode == "groups":
                    widgets_list.widget_results_groups.update_MPL()

            # Close progressbar
            widgets_list.widget_progressbar.close()

        elif button == "GRP_buttons":
            if self.GRP_buttons.checkedId() == 0:
                shared.exp.dfs_mode = "single"
            elif self.GRP_buttons.checkedId() == 1:
                shared.exp.dfs_mode = "groups"

            self.update_widget()

        elif button == "remove":
            table = self.tableWidget.selectionModel().selectedRows()
            rows_to_delete = []
            for item in table:
                rows_to_delete.append(item.row())

            # It depends which way you have selected
            # so I need to sort the index array :
            rows_to_delete = \
                numpy.sort(numpy.unique(rows_to_delete), -1)[::-1].tolist()

            for row_to_delete in rows_to_delete:
                self.tableWidget.removeRow(row_to_delete)
                if shared.exp.dfs_mode == "single":
                    del shared.exp.list_dfs_single_fits[row_to_delete]
                    del shared.exp.list_dfs_single_results[row_to_delete]
                    widgets_list.widget_results_single.update_MPL()
                elif shared.exp.dfs_mode == "groups":
                    del shared.exp.list_dfs_groups_fits[row_to_delete]
                    del shared.exp.list_dfs_groups_results[row_to_delete]
                    widgets_list.widget_results_groups.update_MPL()

    def update_widget(self):
        """update the GUI."""
        # Check the right radio button
        if shared.exp.dfs_mode == "single":
            self.RBT_single.setChecked(True)
            length = len(shared.exp.list_dfs_single_fits)
        elif shared.exp.dfs_mode == "groups":
            self.RBT_groups.setChecked(True)
            length = len(shared.exp.list_dfs_groups_fits)

        # Remove old rows from current table
        for i in range(self.tableWidget.rowCount(), -1, -1):
            self.tableWidget.removeRow(i)

        for i in range(length):
            self.tableWidget.insertRow(i)
            self.tableWidget.setRowHeight(i, 25)

            # Checkbox
            if shared.exp.dfs_mode == "single":
                display = shared.exp.list_dfs_single_fits[i][4]
                xvalue = shared.exp.list_dfs_single_results[i][0] * 1e9
                kuvalue = shared.exp.list_dfs_single_results[i][1]
            elif shared.exp.dfs_mode == "groups":
                display = shared.exp.list_dfs_groups_fits[i][4]
                xvalue = shared.exp.list_dfs_groups_results[i][0] * 1e9
                kuvalue = shared.exp.list_dfs_groups_results[i][1]

            cb = gui_tools.CenteredCellCheckbox(self, "display_fit", i)
            self.tableWidget.setCellWidget(i, 0, cb)
            self.tableWidget.cellWidget(i, 0).checkbox.setChecked(display)

            # x
            label = QtWidgets.QLabel()
            label.setFont(self.parent.parent.smallfont)
            label.setText(str(xvalue))
            centered_label = gui_tools.CenteredCellwidget(label)
            self.tableWidget.setCellWidget(i, 1, centered_label)

            # ku
            label = QtWidgets.QLabel()
            label.setFont(self.parent.parent.smallfont)
            label.setText(str(kuvalue))
            centered_label = gui_tools.CenteredCellwidget(label)
            self.tableWidget.setCellWidget(i, 2, centered_label)

    def checkbox_clicked(self, checkbox, row):
        """Called when a checkbox is clicked.

        Hides or displays a DFS fit.
        """
        if checkbox == "display_fit":
            value = self.tableWidget.cellWidget(row, 0).checkbox.isChecked()

            if shared.exp.dfs_mode == "single":
                shared.exp.list_dfs_single_fits[row][4] = value
                widgets_list.widget_results_single.update_MPL()
            elif shared.exp.dfs_mode == "groups":
                shared.exp.list_dfs_groups_fits[row][4] = value
                widgets_list.widget_results_groups.update_MPL()
