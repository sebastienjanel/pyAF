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

"""Functions used in the results_single tab.

You can duplicate, rename or remove a result in the results table.
You can export your results to a folder.
"""

import copy
import os
from PyQt5 import QtWidgets

from .. import consts
from .. import shared
from .. import widgets_list
from ..tools import export, misc_tools, stat_tools


def duplicate_results():
    """Duplicates the selected rows in the results table."""
    widget = widgets_list.widget_results_single

    # Save the currently chosen parameters and options
    widget.save_hist_prefs()

    # Go through the selected rows of the tableWidget
    for item in widget.tableWidget.selectionModel().selectedRows():
        list_id = item.row()
        result = shared.exp.results_list[list_id]

        # Create a new result (copy the old parameters)
        newresult = copy.deepcopy(result)
        newresult.result_id = len(shared.exp.results_list)

        # Change the name
        newresult.name = result.name + " (Copy)"

        # Append to the results list
        shared.exp.results_list.append(newresult)

    # Update the table and the values
    widget.reset_table()
    stat_tools.get_values()
    widget.update_widget()


def rename_result():
    """Renames the selected row in the results table.

    You can only rename one row at a time.
    """
    widget = widgets_list.widget_results_single

    isok = False
    if (len(widget.tableWidget.selectionModel().selectedRows())) > 1:
        # Tell the user to select only one row
        QtWidgets.QMessageBox.warning(
            widget,
            "Error",
            "Select only one row please")
    else:
        for item in widget.tableWidget.selectionModel().selectedRows():
            result_id = item.row()
            oldname = shared.exp.results_list[result_id].name
            text = "Rename as : "
            mode = QtWidgets.QLineEdit.Normal
            newname, isok = QtWidgets.QInputDialog.getText(widget, "Input", text,
                                                       mode, oldname)
    if isok:
        # Update the value
        newname = str(newname)
        shared.exp.results_list[result_id].name = newname

        # Update the GUI
        widget.tableWidget.cellWidget(result_id, 0).item.setText(newname)
        widget.update_MPL()


def remove_result():
    """Removes the selected row in the results table.

    You can only remove one row at a time.
    """
    widget = widgets_list.widget_results_single

    nbr_found = 0
    if (len(widget.tableWidget.selectionModel().selectedRows())) > 1:
        # Tell the user to select only one row
        QtWidgets.QMessageBox.warning(
            widget,
            "Error",
            "Select only one row please")
    else:
        for item in widget.tableWidget.selectionModel().selectedRows():
            old_id = item.row()

            data_id = shared.exp.results_list[old_id].data_id

            # Check if it's the last result for this data_id
            for result in shared.exp.results_list:
                if result.data_id == data_id:
                    nbr_found += 1
        # Only delete if there is still a file left for this result
        if nbr_found > 1:
            # Delete result
            del shared.exp.results_list[old_id]

            # Shift indices
            for i in range(len(shared.exp.results_list)):
                shared.exp.results_list[i].result_id = i

            # Update GUI
            widget.reset_table()
            stat_tools.get_values()
            widget.update_widget()
        else:
            text = "You can't remove this line, "\
                + "it's the last one linked to the "\
                + shared.exp.list[data_id].filename + " file."
            QtWidgets.QMessageBox.warning(widget, "Error", text)


def export_results(mode):
    """Exports the data as text files.

    Is called from the menu. Will ask the user where to save the files.

    Warn the user that if he checked the log scale mode in the
    histograms, the data will not be exported as log(E) but as E. He
    can do the log computation himself on his data.
    """
    widget = widgets_list.widget_results_single

    # Get the save path from the settings
    save_path = misc_tools.get_user_path(widget.settings)

    # Check for log mode
    if shared.exp.hist_log:
        text = "You have checked the log scale checkbox for the stiffness "\
               "results. Your data will not be exported as log(E), but you "\
               "can still do this computation yourself afterwards."
        QtWidgets.QMessageBox.warning(widget, "Exporting with log scale", text)

    # Open dialog the save the file
    text = "Export checked"
    folder, _ = QtWidgets.QFileDialog.getSaveFileName(widget, text, save_path)

    if folder:
        widgets_list.widget_results_single.save_hist_prefs()
        folder = str(folder)

        # Save current path
        misc_tools.set_user_path(widget.settings, os.path.dirname(folder))

        # Export the data
        if mode == "classic_format":
            export.export_results(folder)

        elif mode == "r_format":
            export.export_results_r(folder)

        QtWidgets.QApplication.beep()
        QtWidgets.QMessageBox.warning(widget, "Information", "Exporting done.")


def change_color():
    """Change the color of the currently selected histogram."""
    # On which tab are we ?
    tab = widgets_list.widget_results_single.parent.tabs.currentIndex()
    if tab == 3:
        tw = widgets_list.widget_results_single.tableWidget
    elif tab == 4:
        tw = widgets_list.widget_results_groups.tableWidget

    # Allow only one row
    if (len(tw.selectionModel().selectedRows())) > 1:
        # Tell the user to select only one row
        text = "Select only one row please"
        widget = widgets_list.widget_results_single
        QtWidgets.QMessageBox.warning(widget, "Error", text)
    elif (len(tw.selectionModel().selectedRows())) == 0:
        # Tell the user to select at least one row
        text = "Select at least on row please"
        widget = widgets_list.widget_results_groups
        QtWidgets.QMessageBox.warning(widget, "Error", text)

    else:
        # Get the row
        for item in tw.selectionModel().selectedRows():
            row_id = item.row()
        if tab == 3:
            current_color = shared.exp.results_list[row_id].color
        elif tab == 4:
            # Skip group 0 == None group
            current_color = shared.exp.groups_list[row_id + 1].color

        if consts.UNIT_TESTING is False:
            color = misc_tools.ask_user_for_color(current_color)
            if color is not False:
                if tab == 3:
                    shared.exp.results_list[row_id].color = color
                elif tab == 4:
                    # Skip group 0 == None group
                    shared.exp.groups_list[row_id + 1].color = color

        # Update the plot
        if tab == 3:
            widgets_list.widget_results_single.update_MPL()
        elif tab == 4:
            widgets_list.widget_results_groups.update_MPL()
