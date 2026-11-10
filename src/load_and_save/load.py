# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
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

"""Module containing the class to load AFM files in pyAF."""

import logging
import copy
import os
from qtpy import QtWidgets
from ..tools.temp_file import TempFile
from ..widgets.progressbar import Progressbar
from ..tools import gui_tools
from ..tools import misc_tools
from ..tools import stat_tools
from ..load_and_save.file_loader.pyaf import LoaderPYAF
from ..load_and_save.file_loader.nanoscope import LoaderNanoscope
from ..load_and_save.file_loader.jpk import LoaderJPK
from ..load_and_save.file_loader.jpk_smart import LoaderJPKSmart
from ..tools import events_refresher
from .. import widgets_list
from .. import shared
from .. import consts

# Maps of the viscoelasticity with the same color scale in all the files at
# load (the moduli, as the elasticity, keep the scale of each file)
VISCO_MAPS = ("visco_fluidity", "visco_tan_delta")
# Maximum number of files listed in the message about the short curves
MAX_LISTED_FILES = 10


class Load:
    """Class to load AFM files or .pyaf files in pyAF.

    Gets as input the files_info list (which was defined in pre_check).
    This list contains dictionnarys with informations about the files to
    load. (Their type, path, postion, and if they should be loaded or
    not).
    """

    def __init__(self, parent, files_info, load_more=False):
        self.logger = logging.getLogger()
        self.logger.debug("Loading files")

        self.parent = parent
        self.load_more = load_more

        # Create a progressbar
        Progressbar()

        # Check for non selected files. These need to be removed from
        # the list
        to_delete = []
        for i in range(len(files_info)):
            if not files_info[i]["checked"]:
                to_delete.append(i)

        if to_delete:
            for i in range(len(to_delete) - 1, -1, -1):
                del files_info[to_delete[i]]

        # The first time we create the temp file
        if not self.load_more:
            # New temporary file object
            shared.exp.temp_file = TempFile()
            # Delete old temporary files from tmp folder
            shared.exp.temp_file.clean_temp_folder()
            # Create a temporary file
            # (only if the first file is not a pyAF file)
            if files_info[0]["file_type"] != "pyAF":
                shared.exp.temp_file.create_new_file()

        # Store the old length
        offset = len(shared.exp.list)

        # Files loaded here with curves too short to be fitted
        self.files_short_curves = []
        # Files loaded here with curves whose position is not in their pixel
        # ([file name, number of curves])
        self.files_misplaced = []

        # Go through the list of files
        for i in range(len(files_info)):
            self.load_single_file(files_info[i])

        if not self.load_more:
            # Save multiprocessing options at first load only
            cores = misc_tools.get_cores()
            shared.exp.cores = cores

        if self.load_more:
            # Update the color scales
            if shared.exp.apply_to_all_data:
                update_colors_at_load()

            # Recreate data list
            widgets_list.widget_data.create_new_list()

            # Select the last file in the list and update GUI
            shared.exp.id_selected = len(shared.exp.list) - 1
            self.parent.file_changed(option="load_more")

            # Enable remove button (we have more than 2 files now)
            widgets_list.widget_data.BT_remove_file.setEnabled(True)

            # Refetch results to fill all the shared data for the tables
            widgets_list.widget_results_single.reset_table()
            stat_tools.get_values()

            # Refetch results for the groups
            stat_tools.fetch_group_data()
            widgets_list.widget_results_groups.update_widget()

            # Refetch results for the conditions
            widgets_list.widget_results_experiment.update_widget()

            # Update the widgets
            widgets_list.widget_results_single.update_widget()
            widgets_list.widget_results_groups.update_widget()

        # If more .pyaf files are loaded, it can happen that the positions in
        # the first file are outside of the bounds of the meshgrids of the
        # second file. So it's safer to reset the current position to 0,0
        shared.exp.meshgrid_click_xpos = 0
        shared.exp.meshgrid_click_ypos = 0
        # (in the region of the map of the selected file)
        if shared.exp.list:
            shared.exp.meshgrid_click_xpos, shared.exp.meshgrid_click_ypos = \
                shared.exp.current_data.nearest_inside(0, 0)

        # Display a message if some curves are too short to be fitted
        lines = ""
        for data in self.files_short_curves[:MAX_LISTED_FILES]:
            # Short name (as in the plots)
            from ..experiment import auto_short_name
            name = auto_short_name(os.path.basename(str(data.filename)))
            total = data.nbr_pixels_x * data.nbr_pixels_y
            lines += "<li>" + name + " : " + \
                str(len(data.corrupted_curves)) + " of " + str(total) + \
                " curves</li>"
        others = len(self.files_short_curves) - MAX_LISTED_FILES
        if others > 0:
            lines += "<li>... and " + str(others) + " other file" + \
                ("s" if others > 1 else "") + "</li>"

        if lines != "" and consts.UNIT_TESTING is False:
            text = "Some curves are shorter than the baseline window of " \
                "the point of contact, or have no data, and may not be " \
                "fitted. The Clean up widget lists these curves: <ul>" + \
                lines + "</ul>"

            # Warning (orange sign, as the other messages)
            gui_tools.show_success(None, "Short curves", details=text,
                                   warning=True)

        # Display a message if the position of some curves is not in their
        # pixel of the map
        if self.files_misplaced and consts.UNIT_TESTING is False:
            from ..experiment import auto_short_name
            lines = "".join(
                "<li>" + auto_short_name(name) + " : " + str(nbr) +
                " curves</li>"
                for name, nbr in self.files_misplaced[:MAX_LISTED_FILES])
            text = "The position of some curves recorded by the " \
                "microscope is not in their pixel of the map: <ul>" + \
                lines + "</ul>The curves are placed with their index on " \
                "the map. Check the maps of these files."
            gui_tools.show_success(None, "Positions of the curves",
                                   details=text, warning=True)

        # Check files with missmatch in the number of segments described in the header and the segments found
        # while loading the file.

        str_missmatch = ""
        # print(shared.exp.segment_handling)
        for file in shared.exp.segment_handling:
            str_missmatch = str_missmatch + "<li> " + f"{file[0]} <ul>" \
                                          + f"<li>{file[1]} segments are described. </li>" \
                                          + f"<li>{file[2]} segments where found. </li> </ul>" + "</li>"

        if str_missmatch != "" and consts.UNIT_TESTING is False:
            text = "On the header of the following files: <ul>"\
                    + str_missmatch + "</ul>" + "<br>"\
                    + "[!] Only the found segments will be loaded."

            msg_corr = QtWidgets.QMessageBox()
            msg_corr.setText("Error (Number of segments described on header does not match number of segments found)")
            msg_corr.setInformativeText(text)
            msg_corr.setIcon(QtWidgets.QMessageBox.Critical)
            msg_corr.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
            msg_corr.exec_()

        # Reset focus to main window
        if self.parent.isVisible() and not self.parent.isActiveWindow():
            self.parent.activateWindow()

        self.logger.debug("Loading successful")

    def load_single_file(self, file_info):
        """Loads a single file, depending on it's type.

        Will call the classes in the file_loader folder and load a file.
        Some update and cleanup actions are done after the loading.
        """
        offset = len(shared.exp.list)

        file_type = file_info["file_type"]
        filename = file_info["filename"]
        self.logger.debug("Loading : %s, %s", filename, file_type)

        if widgets_list.widget_progressbar is None:
            widgets_list.widget_progressbar = Progressbar()
        # Define progressbar's range and label and reset it to 0
        widgets_list.widget_progressbar.set_label("Loading " + filename)
        widgets_list.widget_progressbar.reset()

        if file_type != "pyAF":
            # Set the range. For pyAF files this is done in the LoaderPYAF
            # class.
            nbr = file_info["nbrcurves"]
            widgets_list.widget_progressbar.set_range(0, nbr)

        # Add group headers to temp_file with an unique id
        if file_type != "pyAF":
            shared.exp.temp_file.create_storage_group(offset)
            shared.exp.temp_file.flush_file()

        # Add empty dataset to the experiment
        # For PYAF files this is done inside the pyaf.py loader, because there
        # is more stuff to do than just create an empty data objet in
        # shared.exp
        if file_type != "pyAF":
            shared.exp.addData(filename, offset)
            shared.exp.list[offset].file_type = file_type

        if file_type == "pyAF":
            if offset != 0:
                load_more = True
            else:
                load_more = False

            LoaderPYAF(file_info, offset, load_more)

        elif file_type == "JPK (Force Map)" or\
                file_type == "JPK (QI)" or\
                file_type == "JPK (Single File)":
            loader = LoaderJPK(file_info, offset)
            if loader.misplaced:
                self.files_misplaced.append([filename, loader.misplaced])

        elif file_type == "JPK (Smart Mapping)":
            loader = LoaderJPKSmart(file_info, offset)
            if loader.misplaced:
                self.files_misplaced.append([filename, loader.misplaced])

        elif file_type == "Nanoscope (Peak Force)" or\
                file_type == "Nanoscope (Force Volume)" or\
                file_type == "Nanoscope (Single File)":
            LoaderNanoscope(file_info, offset)

            data = shared.exp.list[offset]

            if data.version == 5120000 and data.nbr_pixels_x >= 128:
                data.version = file_info["version"]

        # Update the data (reset discarded curves array if needed)
        if file_type != "pyAF":
            data = shared.exp.list[offset]
            data.reset_discarded()
            data.update(color=False)
        else:
            # For pyaf files update all the data
            # First, update the events arrays, this is slow so this is done
            # while displaying a progressbar
            events_refresher.update_events("all")
            for data in shared.exp.list:
                data.update(color=False)

        # Add to parent's lists (only if load more, for the first load the
        # lists are populated during UI initialisation in main.py). If it's a
        # .pyaf file this is done inside the pyaf.py loader.
        if self.load_more and file_type != "pyAF":
            # Loading single file
            widgets_list.widget_results.list_exp.addItem(filename)
            widgets_list.widget_compute.list_exp2.addItem(filename)

        # Update some values
        shared.exp.update_global_values()

        if file_type != "pyAF":
            # Curves too short to be fitted, listed to be discarded with the
            # clean up tool: shorter than the baseline window of the point of
            # contact (height range of the valid points), or without data. The number of points of a curve is
            # not a criterion: JPK curves can have far fewer points than the
            # longest one and be complete.

            data_set = shared.exp.list[offset]

            # Baseline window, at most half of the ramp: with a window longer
            # than the ramp (default values, short ramp), all the curves would
            # be listed, and the window has to be changed for the file.
            limit = (data_set.fitparam_poc_skip_start +
                     data_set.fitparam_poc_fit_length)
            if data_set.ramp_size:
                limit = min(limit, abs(data_set.ramp_size) / 2.0)

            # List of curves too short
            corrupted_curves = []

            # Positions are determined in the respective loading class. These
            # are the positions of the first and last valid values in the
            # curves.
            positions = data_set.approach_positions
            for i in range(data_set.nbr_pixels_x):
                heights = data_set.curves_approach[i, :, 0, :]
                for j in range(data_set.nbr_pixels_y):
                    # No curve outside of the region of the map
                    if data_set.is_outside(i, j):
                        continue
                    app_start = int(positions[i][j][0])
                    app_end = int(positions[i][j][1])
                    if app_end - app_start < 2 or \
                            abs(heights[j][app_end - 1] -
                                heights[j][app_start]) < limit:
                        corrupted_curves.append([i, j])

            # Store the curves in a list
            if corrupted_curves == []:
                data_set.corrupted_curves = None
            else:
                data_set.corrupted_curves = corrupted_curves
                self.files_short_curves.append(data_set)

            # Create tables to store results (Already done for .pyaf files)
            shared.exp.temp_file.create_tables_for_results(offset, "all")

            # Piezo height: scale of each file (autoscale), the heights of
            # different files are not comparable
            data_set.color_opts_piezo[1] = 0
            data_set.color_opts_piezo[2] = data_set.max_piezo
            data_set.color_opts_piezo[3] = data_set.max_piezo / 2.0


def update_colors_at_load():
    """Update colorscales during the loading procedure.

    Lets say you load a first file, without any stiffness computation.
    The colorscale is not set for this dataset. Now you load a .pyaf file
    with the stiffness computed. The apply_to_all (called during the GUI
    updated) function will then apply the empty colorscale to every file,
    which will make the app bug ...
    The best solution is to take the first non-empty value found in the
    dataset and apply it everywhere (before the apply_to_all part).

    This solution is not perfect and can lead to not adapted values for the
    colorscales, but it is still better than no values at all.

    I use deepcopy because copying a list will keep the values linked to the
    same memory adress, messing up the colorscales afterwards.
    """
    found_stiffness = False
    found_work_and_rupture_force1_calculated = False
    found_events_calculated = False
    visco = None

    for i in range(len(shared.exp.list)):
        data = shared.exp.list[i]

        if data.stiffness_calculated and not found_stiffness:
            found_stiffness = True
            topo = copy.deepcopy(data.color_opts_topo)
            stiffness = copy.deepcopy(data.color_opts_stiffness)

        cond1 = data.work_and_rupture_force1_calculated
        cond2 = found_work_and_rupture_force1_calculated

        if cond1 and cond2 is False:
            found_work_and_rupture_force1_calculated = True
            work = copy.deepcopy(data.color_opts_work)
            force1 = copy.deepcopy(data.color_opts_rupture_force1)

        if data.visco_calculated and visco is None:
            visco = [copy.deepcopy(data.get_color_opts(mt))
                     for mt in VISCO_MAPS]

        if data.events_calculated and found_events_calculated is False:
            found_events_calculated = True
            events = copy.deepcopy(data.color_opts_events_per_curve)
            force2 = copy.deepcopy(data.color_opts_rupture_force2)

    # Set the values
    for i in range(len(shared.exp.list)):
        dt = shared.exp.list[i]

        if found_stiffness:
            # The elasticity keeps the scale of each file (computed or
            # saved), see apply_to_all
            dt.color_opts_topo = copy.deepcopy(topo)
            if not dt.stiffness_calculated:
                dt.color_opts_stiffness = copy.deepcopy(stiffness)

        if found_work_and_rupture_force1_calculated:
            dt.color_opts_work = copy.deepcopy(work)
            dt.color_opts_rupture_force1 = copy.deepcopy(force1)

        if found_events_calculated:
            dt.color_opts_events_per_curve = copy.deepcopy(events)
            dt.color_opts_rupture_force2 = copy.deepcopy(force2)

        if visco is not None:
            for mt, opts in zip(VISCO_MAPS, visco):
                dt.get_color_opts(mt)[:] = copy.deepcopy(opts)