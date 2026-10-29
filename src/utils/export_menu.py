# Copyright Sébastien Janel (2024- ) sebastien.janel@cnrs.fr
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

"""Export menu: dialogs of the exports of tools.export (the parameters of
each data export are written next to it, see tools.params_export)."""

import os

from qtpy import QtCore, QtWidgets

from .. import shared
from .. import widgets_list
from ..tools import export, gui_tools, map_values, misc_tools


# Items of the Export menu: (key, text, tooltip)
ITEMS = (
    ("data_table", "Data table (CSV)...",
     "Values of the checked rows of the Plots tab, one row per curve (or "
     "per event), with their file, replicate and condition (SI units)"),
    ("summary", "Summary per row (CSV)...",
     "Descriptive statistics of each checked row of the Plots tab"),
    ("superplot", "SuperPlot values and tests (CSV)...",
     "Value of each file of the SuperPlot and the tests of the conditions"),
    ("dfs", "DFS fits (CSV)...",
     "Fits of the DFS window (files or speeds) with their points"),
    ("map", "Map values (32-bit TIFF, for analysis)...",
     "Values of the map displayed in the Results tab, for its file: "
     "32-bit TIFF for ImageJ / Fiji (calibrated in um, with the colour "
     "table and range of pyAF; NaN without value), to measure"),
    ("map_figure", "Map figure (PNG, PDF or SVG)...",
     "The map displayed in the Results tab, with its axes and colour bar "
     "(300 dpi for PNG)"),
    (None, None, None),
    ("parameters", "Parameters (JSON)...",
     "Parameters of all the files: calibration, acquisition and the "
     "parameters used for the computations"),
    ("all", "Everything to a folder...",
     "Everything computed, whatever is displayed: data table and summary "
     "of all the rows of the Plots tab, SuperPlot of each result type, DFS "
     "fits, all the maps of the files (values and figures) and the "
     "parameters, in a new folder"),
)


def add_menu(parent, menu_bar):
    """Adds the Export menu to the menu bar."""
    menu = menu_bar.addMenu("Export")
    menu.setToolTipsVisible(True)
    for key, text, tooltip in ITEMS:
        if key is None:
            menu.addSeparator()
            continue
        action = QtWidgets.QAction(text, parent)
        action.setToolTip(tooltip + ".\nThe parameters are saved next to "
                          "it (.params.json)." if key != "parameters"
                          else tooltip + ".")
        action.triggered.connect(
            lambda _=False, k=key: export_clicked(parent, k))
        menu.addAction(action)
    return menu


def ask_path(parent, title, name, file_filter):
    """Asks where to save, in the last folder used. None if cancelled."""
    settings = QtCore.QSettings()
    folder = misc_tools.get_user_path(settings)
    path, _ = QtWidgets.QFileDialog.getSaveFileName(
        parent, title, os.path.join(folder, name), file_filter)
    if not path:
        return None
    misc_tools.set_user_path(settings, os.path.dirname(path))
    return path


def base_name():
    """Name of the exports: results type of the plots."""
    return str(shared.exp.results_type)


def plots_ready(parent):
    """True if the Plots tab has values and checked rows (else a message)."""
    # The rows checked in the table
    if shared.single_data is not None:
        widgets_list.widget_results_single.save_hist_prefs()
    if shared.single_data is None or not export.checked_rows():
        QtWidgets.QMessageBox.information(
            parent, "Export", "Check the rows to export in the Plots tab "
            "(the values of the displayed results are exported).")
        return False
    return True


def export_clicked(parent, key):
    """An item of the Export menu was chosen."""
    csv_filter = "CSV (*.csv)"
    files = None
    if key == "data_table":
        if plots_ready(parent):
            path = ask_path(parent, "Export the data table",
                            base_name() + "_data.csv", csv_filter)
            files = path and export.export_data_table(path)
    elif key == "summary":
        if plots_ready(parent):
            path = ask_path(parent, "Export the summary",
                            base_name() + "_summary.csv", csv_filter)
            files = path and export.export_summary(path)
    elif key == "superplot":
        files = export_superplot(parent)
    elif key == "dfs":
        if not [f for f in export.dfs_tools.get_fits(shared.exp.dfs_mode)
                if f.get("method") != "old"]:
            QtWidgets.QMessageBox.information(
                parent, "Export", "No fit in the DFS window (files or "
                "speeds, as selected in it).")
            return
        path = ask_path(parent, "Export the DFS fits", "dfs_fits.csv",
                        csv_filter)
        files = path and export.export_dfs_fits(path)
    elif key == "map":
        data = shared.exp.current_data
        map_type = data.meshgrid_type
        if map_type not in map_values.available_maps(data):
            map_type = "piezo"
        name = os.path.splitext(os.path.basename(str(data.filename)))[0]
        path = ask_path(parent, "Export the map values",
                        name + "_" + export.map_type_name(map_type) + ".tif",
                        "TIFF (*.tif *.tiff)")
        if path:
            if not path.lower().endswith((".tif", ".tiff")):
                path += ".tif"
            files = export.export_map(path, shared.exp.id_selected,
                                      map_type)
    elif key == "map_figure":
        data = shared.exp.current_data
        name = os.path.splitext(os.path.basename(str(data.filename)))[0]
        path = ask_path(parent, "Export the map figure",
                        name + "_" + export.map_type_name(
                            data.meshgrid_type) + ".png",
                        "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if path:
            if not path.lower().endswith((".png", ".pdf", ".svg")):
                path += ".png"
            files = export.export_map_figure(path, shared.exp.id_selected)
    elif key == "parameters":
        path = ask_path(parent, "Export the parameters", "parameters.json",
                        "JSON (*.json)")
        files = path and export.export_parameters(path)
    elif key == "all":
        files = export_all(parent)
    if files:
        done(parent, files)


def export_superplot(parent):
    """SuperPlot values and tests (also the Export CSV button of the
    SuperPlots tab). Returns the written files."""
    if not export.superplot_tools.get_data():
        QtWidgets.QMessageBox.information(
            parent, "Export", "No SuperPlot: put the files in replicates "
            "(Plots tab) and the replicates in conditions (Replicates tab).")
        return None
    path = ask_path(parent, "Export the SuperPlot values",
                    "superplot_" + base_name() + ".csv", "CSV (*.csv)")
    return path and export.export_superplot(path)


def export_all(parent):
    """Everything in a new folder (or an existing one, after asking).
    Returns the written files."""
    path = ask_path(parent, "Name of the folder of the export",
                    "pyAF_export_" + base_name(), "Folder (*)")
    if not path:
        return None
    if os.path.isfile(path):
        QtWidgets.QMessageBox.warning(
            parent, "Export", "A file has this name: choose a folder name.")
        return None
    if os.path.isdir(path) and os.listdir(path):
        answer = QtWidgets.QMessageBox.question(
            parent, "Export", "The folder " + os.path.basename(path) +
            " is not empty: write the files in it (the files of the same "
            "names are replaced, the others are kept)?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if answer != QtWidgets.QMessageBox.Yes:
            return None
    if shared.single_data is not None:
        widgets_list.widget_results_single.save_hist_prefs()
    return export.export_all(path)


def done(parent, files):
    """Message with the written files."""
    folder = os.path.dirname(files[0])
    names = [os.path.relpath(f, folder) for f in files]
    gui_tools.show_success(
        parent, "Export complete", str(len(files)) + " file" +
        ("s" if len(files) > 1 else "") + " written in " + folder + ":\n" +
        "\n".join(names))
