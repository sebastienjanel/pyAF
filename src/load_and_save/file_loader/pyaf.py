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

"""Loads .pyaf files."""

import os
import tables

from ... import consts
import logging
try:
    import pickle as pickle
except ImportError:
    import pickle
from ... import widgets_list
from ... import shared
from ...tools import misc_tools

class LoaderPYAF:
    """Class used to load pyAF files."""

    def __init__(self, file_info, more_offset, load_more):
        logger = logging.getLogger()

        wl = widgets_list

        self.path = file_info["path"]
        self.load_more = load_more

        # Get the number of files to load, and the bytesoffsets
        inputfile = open(self.path, "rb")
        line = inputfile.readline().decode()

        while line.find("Version=") != 0:
            line = inputfile.readline().decode()
        f_version = line.split("=")[1].splitlines()[0]
        ver = consts.VERSION
        logger.debug("pyAF file version %s, pyAF version %s", f_version, ver)
        while line.find("Nbr_of_files=") != 0:
            line = inputfile.readline().decode()
        nbr_files = int(line.split("=")[1].splitlines()[0])
        while line.find("Global_Prefs_Offset=") != 0:
            line = inputfile.readline().decode()
        globaloffset = int(line.split("=")[1].splitlines()[0])

        # Get the offsets for each prefs
        single_offsets = []
        singleoffset_total = 0
        for i in range(nbr_files):
            line = inputfile.readline().decode()
            offset = int(line.split("=")[1].splitlines()[0])
            single_offsets.append(offset)
            singleoffset_total += offset
        offset_prefs = inputfile.tell()  # Get the current position

        # Skip the prefs files to directly go to the hdf5 part
        inputfile.seek(globaloffset + singleoffset_total, 1)

        # Get size of the file to load it piecewise
        size = os.path.getsize(self.path) - inputfile.tell()
        blocksize = 100000.0

        tf = shared.exp.temp_file

        if not self.load_more and len(shared.exp.list) == 0:
            more_offset = 0
            wl.widget_progressbar.set_range(0, int(size / blocksize))
            tf.load_old_file(inputfile, size, blocksize)
            # Delete old temporary files from tmp folder
            shared.exp.temp_file.clean_temp_folder()
        else:
            # Copy data (can be slow)
            wl.widget_progressbar.set_range(0, int(size / blocksize))
            tf.load_more_old_file(inputfile, shared.exp,
                                  nbr_files, more_offset, size, blocksize)

        # Flush file to be sure everything is written on the disc
        shared.exp.temp_file.flush_file()

        # Go back to the prefs part
        inputfile.seek(offset_prefs)

        # Load main prefs (experiment class)
        general_prefs = inputfile.read(globaloffset)

        try:
            general_prefs = pickle.loads(general_prefs, encoding=consts.ENCODING)

        except ImportError:
            # In previous versions of pyAF, when you pickle dump the general
            # preferences, you dump also some objects linked to classes in
            # the experiment module. (for example groups or results).
            # The imports were then changed during 1.4 dev to use relative
            # paths, meaning that the module name is now src.experiment.
            # In this case we get an import error as pickle is looking for the
            # module name "experiment". We can "fake" the module name by adding
            # an experiment module to the sys.modules keys. This will allow the
            # loading of old files. See error #343
            # Could be removed for pyAF version 1.6 as this is a dirty hack
            # to get this working.
            from ... import experiment
            import sys
            sys.modules["experiment"] = experiment
            general_prefs = pickle.loads(general_prefs, encoding=consts.ENCODING)

        if not self.load_more:
            # Load general prefs (only for the first load)
            # In load more, the preferences from the first .pyAF file are kept.
            do_not_copy = ["filenames_list", "results_list"]
            # Some times while decoding general_pref they can be decoded into bytes.
            # This will cause an error.
            for key in list(general_prefs.keys()):
                # Do not copy some stuff, those are updated separatly because
                # we don't want to overwrite those lists.
                if key not in do_not_copy:
                    value = general_prefs[key]
                    misc_tools.setattr_special(shared.exp, key, value)

        # Update the results_list, the id's have to be shifted !
        new_list = general_prefs["results_list"]
        old_list = shared.exp.results_list
        max_id = len(old_list)

        added_one_result = False

        old_list = shared.exp.groups_list
        group_max_id = len(old_list)

        for i in old_list:
            misc_tools.add_attributes(i)


        # Add to results
        for i in range(len(new_list)):

            misc_tools.add_attributes(new_list[i])

            # Append the result
            shared.exp.results_list.append(new_list[i])
            # Update the id's
            result = shared.exp.results_list[max_id + i]
            result.result_id = max_id + i
            result.data_id = more_offset + new_list[i].data_id
            added_one_result = True
            if self.load_more:
                # Update the group id for the new results. The corresponding
                # groups will be added below.
                result.group += group_max_id - 1

        # At least one result was added, if there was none in the previous
        # file, we need to update the result's type to have one. Else it
        # will remain set to None ...
        if added_one_result:
            shared.exp.results_type = general_prefs["results_type"]

        if self.load_more:
            fname = os.path.basename(self.path)

            # Update the groups list. Create new groups from the new file.
            # The groups names are updated with the name of the file they are
            # comming from.

            # The groups ids in the results_list have already been updated
            # above.

            # Append new groups (skip first one, is None group)
            for i in range(1, len(general_prefs["groups_list"]), 1):

                value = general_prefs["groups_list"][i]

                misc_tools.add_attributes(value)

                shared.exp.groups_list.append(value)

                # Update the name before adding the row to the tablewidget
                nm = " (from " + fname + ")"
                shared.exp.groups_list[group_max_id + i - 1].name += nm

                # - 2 because : -1 is because we start at 1 (skip None group)
                # The other -1 is because insertRow needs the position of the
                # last row.
                if wl.widget_results_groups is not None:
                    pos = group_max_id + i - 2
                    wl.widget_results_groups.tableWidget.insertRow(pos)
                    wl.widget_results_groups.add_row_in_tablewidget_groups(pos)

        # Create the datasets
        for i in range(more_offset, more_offset + nbr_files, 1):
            filename = general_prefs["filenames_list"][i - more_offset]
            shared.exp.addData(filename, i)

            if self.load_more and wl.widget_results is not None:
                # Add the names to the list (only when loading more, on first
                # load this is done during UI initialisation)
                wl.widget_results.list_exp.addItem(filename)
                wl.widget_compute.list_exp2.addItem(filename)

            # Load single prefs in the datasets
            value = inputfile.read(single_offsets[i - more_offset])

            single_prefs = pickle.loads(value, encoding=consts.ENCODING)

            skip = False
            for key in single_prefs:
                value = single_prefs[key]
                # COMPATIBILITY MODE 11/06/14
                if key == "nbr_points_x":
                    key = "nbr_pixels_x"
                if key == "nbr_points_y":
                    key = "nbr_pixels_y"
                # COMPATIBILITY MODE 10/07/14
                if key == "force_samples_per_pause":
                    key = "nbr_points_per_pause_curve"
                if key == "force_samples_per_curve":
                    key = "nbr_points_per_curve_approach"
                    # Add also the retraction with the same value
                    misc_tools.setattr_special(
                        shared.exp.list[i],
                        "nbr_points_per_curve_retraction",
                        value)
                if key == "force_samples_per_curve_real":
                    key = "nbr_points_per_curve_approach_real"
                    # Add also the retraction with the same value
                    misc_tools.setattr_special(
                        shared.exp.list[i],
                        "nbr_points_per_curve_retraction_real",
                        value)
                # COMPATIBILITY MODE 21/08/14
                if key == "used_spring_constant":
                    misc_tools.setattr_special(
                        shared.exp.list[i],
                        "_used_spring_constant",
                        value)
                    skip = True
                # ---------------------------
                if key != "dontsave" and skip is False:
                    # dontsave key was added to the dontsave list during
                    # 1.5 development. For all the files before, the dontsave
                    # list was saved and would overwrite the new one.
                    # So don't overwrite it. This if statement can be removed
                    # for the 1.6 release.
                    misc_tools.setattr_special(shared.exp.list[i], key, value)
                # Reset skip
                skip = False

            # Update the id
            data = shared.exp.list[i]
            data.unique_id = i

            # COMPATIBILITY MODE 17/01/14 add new vtk values
            try:
                _ = data.vtk_smoothing_type
            except AttributeError:
                data.vtk_smoothing_type = None
                data.vtk_smoothing_iterations = 1
                data.vtk_decimate_target_reduction = 0.2

        # COMPATIBILITY MODE 09/04/14 change pocs_corrected to topography
        for i in range(more_offset, more_offset + nbr_files, 1):
            # Check if node exists
            fi = shared.exp.temp_file
            st = "/data/_" + str(i) + "/results/"
            try:
                fi.file.get_node(st + "topography")
            except tables.exceptions.NoSuchNodeError:
                # Rename node
                fi.file.rename_node(st + "pocs_corrected", "topography")

        # COMPATIBILITY MODE 12/08/14 add events joc fits
        for i in range(more_offset, more_offset + nbr_files, 1):
            # Check if node exists
            fi = shared.exp.temp_file
            st = "/data/_" + str(i) + "/results/"
            try:
                fi.file.get_node(st + "events_fits_joc")
            except tables.exceptions.NoSuchNodeError:
                # Rename node
                fi.file.create_carray(
                    st, "events_fits_joc",
                    tables.Float32Atom(),
                    shape=(data.nbr_pixels_x, data.nbr_pixels_y, 3),
                    filters = tables.Filters(complevel=9, complib="blosc"))
                values = fi.file.get_node(st + "events_fits_joc")
                # Fill with empty values, the user will have to recompute
                for i in range(data.nbr_pixels_x):
                    for j in range(data.nbr_pixels_y):
                        values[i, j] = [0, 0, 0]

        # COMPATIBILITY MODE 22/02/14 add color to results and groups
        i = 0
        for res in shared.exp.results_list:
            try:
                _ = res.color
            except AttributeError:
                res.color = shared.colors_list[i]
            if i < len(shared.colors_list) - 1:
                i = i + 1
            else:
                i = 0
        # First "None" group
        try:
            _ = shared.exp.groups_list[0].color
        except AttributeError:
            shared.exp.groups_list[0].color = None
        # The other groups
        i = 0
        for group_id in range(1, len(shared.exp.groups_list)):
            group = shared.exp.groups_list[group_id]
            try:
                _ = group.color
            except AttributeError:
                group.color = shared.colors_list[i]
            if i < len(shared.colors_list) - 1:
                i = i + 1
            else:
                i = 0

        # Close the inputfile
        inputfile.close()
