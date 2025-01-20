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

"""Tools used to discard curves."""


from .. import shared
from .. import widgets_list
from PyQt5 import QtWidgets


def add_discard_actions_to_menu(parent, menu):
    """Called by utils/meshgrid_menu. Adds entries about discarding in the menu."""
    data = shared.exp.current_data
    xpos = shared.exp.meshgrid_click_xpos
    ypos = shared.exp.meshgrid_click_ypos

    if data.discarded_curves[xpos][ypos] == 1:
        single_curve = "Keep curve"
    else:
        single_curve = "Discard curve"

    if data.discard_corrupted_curves_applied:
        single_file = "Keep corrupted (this file)"
    else:
        single_file = "Discard corrupted (this file)"

    if shared.exp.discard_corrupted_curves_applied_all:
        all_files = "Keep corrupted (all files)"
    else:
        all_files = "Discard corrupted (all files)"

    # Change curve's status (discarded or not)
    action = QtWidgets.QAction(single_curve, parent)
    action.triggered.connect(change_single_curve_status)
    menu.addAction(action)
    action = QtWidgets.QAction(single_file, parent)
    action.triggered.connect(change_single_file_status)
    menu.addAction(action)
    action = QtWidgets.QAction(all_files, parent)
    action.triggered.connect(change_all_files_status)
    menu.addAction(action)
    menu.addSeparator()


def change_single_curve_status():
    """Discard a single curve. Takes the currently displayed curve position."""
    data = shared.exp.list[shared.exp.id_selected]
    xpos = shared.exp.meshgrid_click_xpos
    ypos = shared.exp.meshgrid_click_ypos
    if data.discarded_curves[xpos][ypos] == 1:
        data.discarded_curves[xpos][ypos] = 0
    else:
        data.discarded_curves[xpos][ypos] = 1

    refresh_after_status_change()


def change_single_file_status():
    """Remember if all the corrupted curves have been discared or not in a file."""
    data = shared.exp.list[shared.exp.id_selected]
    if data.discard_corrupted_curves_applied:
        data.discard_corrupted_curves_applied = False
        if data.discard_corrupted_curves_applied != \
                shared.exp.discard_corrupted_curves_applied_all:
            shared.exp.discard_corrupted_curves_applied_all = True
    else:
        data.discard_corrupted_curves_applied = True
        if data.discard_corrupted_curves_applied != \
                shared.exp.discard_corrupted_curves_applied_all:
            shared.exp.discard_corrupted_curves_applied_all = False

    discard_corrupted_curves(shared.exp.id_selected)
    refresh_after_status_change()


def change_all_files_status():
    """Discard all the corrupted curves from a file. (Or re-enable them)."""
    orig_value = shared.exp.discard_corrupted_curves_applied_all
    if shared.exp.discard_corrupted_curves_applied_all:
        shared.exp.discard_corrupted_curves_applied_all = False
    else:
        shared.exp.discard_corrupted_curves_applied_all = True
    for data_id in range(len(shared.exp.list)):
        data = shared.exp.list[data_id]
        if orig_value:
            data.discard_corrupted_curves_applied = False
        else:
            data.discard_corrupted_curves_applied = True
        discard_corrupted_curves(data_id)

    refresh_after_status_change()


def refresh_after_status_change():
    """Refresh the plots."""
    widgets_list.widget_data.meshgrid.canvas.update_blit("discarded")
    widgets_list.widget_data.curve.update_plot()

    # Refresh only the results plots if something has been computed
    widgets_list.widget_results.MPL_canvas1.canvas.update_blit("discarded")
    widgets_list.widget_results.MPL_canvas2.update_plot()

    if widgets_list.widget_curve_mod is not None:
        widgets_list.widget_curve_mod.update_widget()


def discard_corrupted_curves(data_id):
    """Do the real discarding of all corrupted curves.

    Is called by the change_all_files_status function.
    """
    data = shared.exp.list[data_id]
    value = data.discard_corrupted_curves_applied
    # Make sure the value is an integer
    if value:
        value = 1
    else:
        value = 0
    if data.corrupted_curves is not None:
        for curve in data.corrupted_curves:
            data.discarded_curves[curve[0]][curve[1]] = value


def display_discarded():
    """Display or hide the discarded curves from the meshgrid."""
    if shared.exp.meshgrid_display_discarded:
        shared.exp.meshgrid_display_discarded = False
    else:
        shared.exp.meshgrid_display_discarded = True

    refresh_after_status_change()
