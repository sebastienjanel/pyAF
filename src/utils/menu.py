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

"""Main menu of PYAF."""

import sys
import webbrowser
from PyQt5 import QtWidgets
from ..load_and_save.save import Save
from . import tablewidget
from .export_figure import export_figure
from ..widgets_vtk.utils.save import save_layers
from ..widgets_vtk.utils.load import load_layers
from ..widgets_vtk.utils import menu_actions
from .. import widgets_list


def create_menu(parent):
    """Creates the menu of PYAF with different actions."""
    if sys.platform == "darwin":
        # No parent for the OS X's menu bar, to have it on every window
        parent.menu = QtWidgets.QMenuBar()
    else:
        parent.menu = QtWidgets.QMenuBar(parent)

    # The list actions is returned to the main pyaf.py module, so that
    # depending on which tab you are on the different menu elements are
    # enabled or disabled
    list_actions = {}

    # File menu
    file_menu = parent.menu.addMenu("File")

    action = QtWidgets.QAction("Save", parent)
    action.setShortcut("Ctrl+S")
    action.triggered.connect(lambda: Save(parent))
    file_menu.addAction(action)

    action = QtWidgets.QAction("Load more ...", parent)
    action.setShortcut("Ctrl+O")
    action.triggered.connect(parent.load_more)
    file_menu.addAction(action)

    # Add prefs (on OSX will be added to the Python menu, not to File)
    action = QtWidgets.QAction("Preferences", parent)
    action.triggered.connect(parent.open_pref_widget)
    file_menu.addAction(action)

    # Misc menu
    misc_menu = parent.menu.addMenu("Misc")

    action = QtWidgets.QAction("Random Curve", parent)
    action.setShortcut("R")
    w_main = widgets_list.widget_main
    action.triggered.connect(lambda: w_main.change_curve(0, 0, "rand"))
    misc_menu.addAction(action)

    action = QtWidgets.QAction("Load plugin", parent)
    action.triggered.connect(parent.open_plugin_window)
    misc_menu.addAction(action)

    # Results menu
    results_menu = parent.menu.addMenu("Plots")

    action = QtWidgets.QAction("Export selected", parent)
    action.triggered.connect(lambda: tablewidget.export_results("classic_format"))
    results_menu.addAction(action)
    list_actions["export_selected"] = action

    action = QtWidgets.QAction("Export selected for R", parent)
    action.triggered.connect(lambda: tablewidget.export_results("r_format"))
    results_menu.addAction(action)
    list_actions["export_selected_r"] = action

    action = QtWidgets.QAction("Rename row", parent)
    action.triggered.connect(tablewidget.rename_result)
    results_menu.addAction(action)
    list_actions["rename_row"] = action

    action = QtWidgets.QAction("Duplicate row", parent)
    action.triggered.connect(tablewidget.duplicate_results)
    results_menu.addAction(action)
    list_actions["duplicate_row"] = action

    action = QtWidgets.QAction("Remove row", parent)
    action.triggered.connect(tablewidget.remove_result)
    results_menu.addAction(action)
    list_actions["remove_row"] = action

    action = QtWidgets.QAction("Change color of row", parent)
    action.triggered.connect(tablewidget.change_color)
    results_menu.addAction(action)
    list_actions["row_color"] = action

    # Figures menu
    figure_menu = parent.menu.addMenu("Figures")

    action = QtWidgets.QAction("Open meshgrid", parent)
    action.triggered.connect(lambda: export_figure("meshgrid"))
    figure_menu.addAction(action)
    list_actions["open_meshgrid"] = action

    action = QtWidgets.QAction("Open curve", parent)
    action.triggered.connect(lambda: export_figure("curve"))
    figure_menu.addAction(action)
    list_actions["open_curve"] = action

    action = QtWidgets.QAction("Open results", parent)
    action.triggered.connect(lambda: export_figure("results"))
    figure_menu.addAction(action)
    list_actions["open_results"] = action

    # VTK (3D) menu
    vtk_menu = parent.menu.addMenu("3D")

    action = QtWidgets.QAction("Save layers", parent)
    action.triggered.connect(lambda: save_layers(parent))
    vtk_menu.addAction(action)

    action = QtWidgets.QAction("Load layers", parent)
    action.triggered.connect(lambda: load_layers(parent))
    vtk_menu.addAction(action)

    action = QtWidgets.QAction("Remove layer", parent)
    action.triggered.connect(menu_actions.remove_layer)
    vtk_menu.addAction(action)

    action = QtWidgets.QAction("Add AFM layer", parent)
    action.triggered.connect(menu_actions.ask_user_for_afm_layers)
    vtk_menu.addAction(action)
    list_actions["add_afm_layer"] = action

    action = QtWidgets.QAction("Copy", parent)
    action.setShortcut("Ctrl+C")
    action.triggered.connect(menu_actions.copy_pos)
    vtk_menu.addAction(action)

    action = QtWidgets.QAction("Paste", parent)
    action.setShortcut("Ctrl+V")
    action.triggered.connect(menu_actions.paste_pos)
    vtk_menu.addAction(action)

    help_menu = parent.menu.addMenu("Help")

    action = QtWidgets.QAction("About pyAF", parent)
    action.triggered.connect(parent.open_about)
    help_menu.addAction(action)

    action = QtWidgets.QAction("Documentation", parent)
    action.triggered.connect(lambda: open_link("docs"))
    help_menu.addAction(action)

    action = QtWidgets.QAction("Report a bug", parent)
    action.triggered.connect(lambda: open_link("bugs"))
    help_menu.addAction(action)

    if sys.platform != "darwin":
        # Do not set the menubar on OS X to have it on every window
        parent.setMenuBar(parent.menu)

    return list_actions


def open_link(link):
    """Opens a link in the browser."""
    if link == "bugs":
        webbrowser.open("https://bitbucket.org/adujardin/pyaf/issues/new")
    elif link == "docs":
        # Not defined for the moment
        pass
