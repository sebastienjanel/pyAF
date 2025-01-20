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

"""
The temporary file is a hdf5 file containing all the data and the results.

It is saved in the default tmp folder of the OS, or if this folder is not
available, in a folder specified by the user.

The advantage of using a hdf5 file is that the data does not need to be loaded
in the RAM, and the saving of the file is fast because it only implies copying
the tmp file to a new destination.

"""

import os
import sys
import time
import tables
from PyQt5 import QtCore
import tempfile
from .. import widgets_list
from .. import shared
from packaging.version import Version


class TempFile:
    """Class with all the methods for the usage of the temp file."""

    def __init__(self, path=None):
        if path is None:
            # Get the settings
            self.settings = QtCore.QSettings()
            default_path = tempfile.gettempdir() + os.sep
            temp_path = self.settings.value("tempPath", default_path)
            # Add a slash if it's missing
            if temp_path[-1] != os.sep:
                temp_path += os.sep
        else:
            temp_path = path

        # Path where the tmp file will be stored
        self.dirpath = str(temp_path)

        self.date = str(time.time()).split(".")[0]
        self.filename = "tmp_" + self.date
        self.filepath = None
        self.file = None

    def create_new_file(self):
        """Create a new (empty) hdf5 file.

        The file will have as name : tmp_date.hdf5. The date is the number of
        seconds passed since the epoch (1970 on Unix).
        See : http://docs.python.org/2/library/time.html

        This allows to know when the file was createad.
        """
        self.filepath = self.dirpath + self.filename + ".hdf5"
        self.file = tables.open_file(self.filepath, "w")
        self.file.create_group("/", "data", "Data")

    def load_old_file(self, inputfile, size, blocksize):
        """Load an saved PYAF file."""
        self.filepath = self.dirpath + self.filename + ".hdf5"
        tosave = open(self.filepath, "wb")

        # Copy blockwise because reading the whole file can result in mmap
        # errors for huge files. Besides, this allows the displaying of a
        # progressbar.
        for _ in range(int(size / blocksize) + 1):
            tosave.write(inputfile.read(int(blocksize)))

            # Update progressbar
            widgets_list.widget_progressbar.update()

        tosave.close()
        self.open_file()

    def load_more_old_file(
            self, inputfile, exp, nbr_files, more_offset, size, blocksize):
        r"""
        Loadind a second PYAF file.

        The difference here is that the nodes need to be copied from the second
        file to the first file.

        The temporary file created to store the data during the copying should
        have a different name than the other tmp files, thats why it is called
        tmp2\_. In some very rare cases where the file loaded just before this
        one is a single curve, time.time returns the same value as the one
        before, and this part of code would then overwrite the first file ...

        """

        # Slower than normal loadOldFile because there is a double file copy

        # Copy hdf5 file to new tmp file
        date = str(time.time()).split(".")[0]
        filepath = self.dirpath + "tmp2_" + date + ".hdf5"
        tosave = open(filepath, "wb")

        # Copy blockwise because reading the whole file can result in mmap
        # errors for huge files
        for _ in range(int(size / blocksize) + 1):
            tosave.write(inputfile.read(int(blocksize)))

            # Update progressbar
            widgets_list.widget_progressbar.update()

        tosave.close()

        # Set a moving progressbar without range
        widgets_list.widget_progressbar.reset()
        widgets_list.widget_progressbar.set_range(0, 0)
        widgets_list.widget_progressbar.update()

        # Import the data in actual tmp file
        oldfile = tables.open_file(filepath, "r+")

        oldfile_id = 0
        for unique_id in range(more_offset, more_offset + nbr_files, 1):
            # Create group
            exp.temp_file.create_storage_group(unique_id)
            # Copy everything
            dstgroup = exp.temp_file.file.get_node("/data/_" + str(unique_id))
            srcgroup = oldfile.get_node("/data/_" + str(oldfile_id))
            exp.temp_file.file.copy_children(srcgroup, dstgroup,
                                             overwrite=True, recursive=True)
            exp.temp_file.flush_file()
            oldfile_id += 1

        # Close old file
        oldfile.close()

        # Delete tmp file
        os.remove(filepath)

    def open_file(self):
        """Open the hdf5 file.

        The file can then be accessed throught the self.file variable.
        """
        self.file = tables.open_file(self.filepath, "r+")

    def flush_file(self):
        """Flush the data to the file."""
        self.file.flush()

    def close_file(self):
        """Close the file."""
        self.file.close()

    def delete_temp_file(self):
        """Delete the file."""
        os.remove(self.filepath)

    def clean_temp_folder(self):
        """Remove old files from the tmp folder.

        It can happen that the software crashed. In this case, the temporary
        file is not removed. This method checks at each launch of the softzare
        if the tmp folder contains old files to remove. We don't want to
        cripple the user's hard drives with old broken hdf5 files.

        The deletion takes places after 12 hours.
        """
        for path in os.listdir(self.dirpath):
            if path.split(".")[-1] == "hdf5":
                # 12 h
                if time.time() - int(path.split(".")[0].split("_")[1]) > 43200:
                    try:
                        os.remove(self.dirpath + path)
                    except OSError:
                        # In some rare cases the rights of the folder could
                        # have changed ... In this case, do nothing.
                        pass

    def create_storage_group(self, unique_id):
        """Create groups to store data.

        The positions group will contain the approach and retraction positions.
        """
        unique_id = str(unique_id)
        pre = "/data/_" + unique_id + "/"
        fi = self.file
        fi.create_group("/data/", "_" + unique_id, "Single Data")
        fi.create_group(pre, "curves", "Curves")
        fi.create_group(pre, "piezo_image", "Piezo Image")
        fi.create_group(pre, "positions", "Positions")
        fi.create_group(pre, "results", "Results")

    def copy_data(self, data, old_id, new_id):
        """Copy the data from one dataset to another."""
        st1 = "/data/_" + str(new_id)
        st2 = "/data/_" + str(old_id)

        self.cp_node(st2 + "/curves/approach", st1 + "/curves")
        self.cp_node(st2 + "/curves/retraction", st1 + "/curves")

        if data.nbr_points_per_pause_curve is not None:
            self.cp_node(st2 + "/curves/pause", st1 + "/curves")

        if data.nbr_points_per_modulation_curve is not None:
            self.cp_node(st2 + "/curves/modulation", st1 + "/curves")

        self.cp_node(st2 + "/piezo_image/piezo_image", st1 + "/piezo_image")
        node = st2 + "/positions/approach_positions"
        self.cp_node(node, st1 + "/positions")
        node = st2 + "/positions/retraction_positions"
        self.cp_node(node, st1 + "/positions")

        st2 = st2 + "/results/"
        st1 = st1 + "/results"

        if data.stiffness_calculated:
            self.cp_node(st2 + "pocs_indices", st1)
            self.cp_node(st2 + "pocs_real", st1)
            self.cp_node(st2 + "topography", st1)
            self.cp_node(st2 + "fits_poc", st1)
            self.cp_node(st2 + "indentation", st1)
            self.cp_node(st2 + "stiffness_array", st1)
            # self.cp_node(st2 + "yf", st1)

        if data.work_and_rupture_force1_calculated:
            self.cp_node(st2 + "jocs1_indices", st1)
            self.cp_node(st2 + "jocs2_indices", st1)
            self.cp_node(st2 + "jocs1_real", st1)
            self.cp_node(st2 + "jocs2_real", st1)
            self.cp_node(st2 + "work", st1)
            self.cp_node(st2 + "rupture_force1", st1)
            self.cp_node(st2 + "fits_joc", st1)

        if data.events_calculated:
            self.cp_node(st2 + "events_positions_middle", st1)
            self.cp_node(st2 + "events_positions_start", st1)
            self.cp_node(st2 + "events_positions_stop", st1)
            self.cp_node(st2 + "events_forces", st1)
            self.cp_node(st2 + "events_slopes", st1)
            self.cp_node(st2 + "events_jocs2_indices", st1)
            self.cp_node(st2 + "events_jocs2_real", st1)
            self.cp_node(st2 + "events_forces_distance", st1)
            self.cp_node(st2 + "events_fits_joc", st1)

        if data.stiffness_corrected:
            self.cp_node(st2 + "stiffness_corrected", st1)

        if data.loading_rates_calculated:
            self.cp_node(st2 + "loading_rates", st1)

    def cp_node(self, st2, st1):
        """Copy single node."""
        fi = self.file
        fi.get_node(st2)._f_copy(st1, overwrite=True)

    def create_data_tables(self, unique_id):
        """Create table for the data.

        Data is stored as 32 bit floats. I tested with 64, it's twice as slow,
        and such a precision is not needed for AFM data.
        """
        dt = shared.exp.list[int(unique_id)]
        st = "/data/_" + str(unique_id)

        if dt.file_type in ("JPK (Single File)", "JPK (Force Map)", "JPK (QI)"):
            n = 3

        else:
            n = 2

        # Create tables to store the data
        # Use compression to reduce the file size
        comp = tables.Filters(complevel=9, complib="blosc")
        fi = self.file

        # Curves (approach)
        shape = (dt.nbr_pixels_x,
                 dt.nbr_pixels_y,
                 n,
                 dt.nbr_points_per_curve_approach)
        fi.create_carray(st + "/curves", "approach", tables.Float32Atom(),
                         shape=shape, filters=comp)

        # Curves (retraction)
        shape = (dt.nbr_pixels_x,
                 dt.nbr_pixels_y,
                 n,
                 dt.nbr_points_per_curve_retraction)
        fi.create_carray(st + "/curves", "retraction", tables.Float32Atom(),
                         shape=shape, filters=comp)

        # Curves (pause)
        if dt.nbr_points_per_pause_curve is not None:
            shape = (
                dt.nbr_pixels_x,
                dt.nbr_pixels_y,
                n,
                dt.nbr_points_per_pause_curve)
            fi.create_carray(st + "/curves", "pause", tables.Float32Atom(),
                             shape=shape, filters=comp)

        # Curves (modulation)
        if dt.nbr_points_per_modulation_curve is not None:
            shape = (
                dt.nbr_pixels_x,
                dt.nbr_pixels_y,
                n,
                dt.nbr_points_per_modulation_curve)
            fi.create_carray(st + "/curves", "modulation", tables.Float32Atom(),
                             shape=shape, filters=comp)

        # Piezo image
        fi.create_carray(st + "/piezo_image", "piezo_image",
                         tables.Float32Atom(),
                         shape=(dt.nbr_pixels_x, dt.nbr_pixels_y),
                         filters = comp)

        # Positions (corrupted parts of the curve)
        fi.create_carray(st + "/positions", "approach_positions",
                         tables.Int32Atom(),
                         shape=(dt.nbr_pixels_x, dt.nbr_pixels_y, 2),
                         filters = comp)
        fi.create_carray(st + "/positions", "retraction_positions",
                         tables.Int32Atom(),
                         shape=(dt.nbr_pixels_x, dt.nbr_pixels_y, 2),
                         filters = comp)

    def delete_tables_for_results(self, unique_id, tabletypes):
        """Delete results tables."""
        st = "/data/_" + str(unique_id)

        if tabletypes == "all" or tabletypes == "stiffness":
            self.file.remove_node(st + "/results/pocs_indices")
            self.file.remove_node(st + "/results/pocs_real")
            self.file.remove_node(st + "/results/topography")
            self.file.remove_node(st + "/results/fits_poc")
            self.file.remove_node(st + "/results/indentation")
            self.file.remove_node(st + "/results/stiffness_array")
        if tabletypes == "all" or tabletypes == "work":
            self.file.remove_node(st + "/results/jocs1_indices")
            self.file.remove_node(st + "/results/jocs2_indices")
            self.file.remove_node(st + "/results/jocs1_real")
            self.file.remove_node(st + "/results/jocs2_real")
            self.file.remove_node(st + "/results/work")
            self.file.remove_node(st + "/results/rupture_force1")
            self.file.remove_node(st + "/results/fits_joc")
        if tabletypes == "all" or tabletypes == "events":
            self.file.remove_node(
                st + "/results/events_positions_middle")
            self.file.remove_node(st + "/results/events_positions_start")
            self.file.remove_node(st + "/results/events_positions_stop")
            self.file.remove_node(st + "/results/events_forces")
            self.file.remove_node(st + "/results/events_slopes")
            self.file.remove_node(st + "/results/events_jocs2_indices")
            self.file.remove_node(st + "/results/events_jocs2_real")
            self.file.remove_node(st + "/results/events_forces_distance")
            self.file.remove_node(st + "/results/events_fits_joc")
        if tabletypes == "all" or tabletypes == "stiffness_corrected":
            self.file.remove_node(st + "/results/stiffness_corrected")
        if tabletypes == "all" or tabletypes == "loading_rates":
            self.file.remove_node(st + "/results/loading_rates")

    def create_tables_for_results(
            self,
            unique_id,
            tabletypes,
            recalc_pocs=True):
        """
        Create tables for the results.

        """

        dt = shared.exp.list[int(unique_id)]
        fi = self.file
        ty_int = tables.Int32Atom()
        ty_fl = tables.Float32Atom()
        st = "/data/_" + str(unique_id) + "/results"

        # Use compression to reduce the file size
        comp = tables.Filters(complevel=9, complib="blosc")

        if tabletypes == "all" or tabletypes == "stiffness":
            if recalc_pocs:
                fi.create_carray(st, "pocs_indices",
                                 ty_int,
                                 shape=(dt.nbr_pixels_x, dt.nbr_pixels_y),
                                 filters = comp)
                fi.create_carray(st, "pocs_real",
                                 ty_fl,
                                 shape=(dt.nbr_pixels_x, dt.nbr_pixels_y, 2),
                                 filters = comp)
                fi.create_carray(st, "topography",
                                 ty_fl,
                                 shape=(dt.nbr_pixels_x, dt.nbr_pixels_y),
                                 filters = comp)
                fi.create_carray(st, "fits_poc",
                                 ty_fl,
                                 shape=(dt.nbr_pixels_x, dt.nbr_pixels_y, 3),
                                 filters = comp)
                fi.create_carray(st, "indentation",
                                 ty_fl,
                                 shape=(dt.nbr_pixels_x, dt.nbr_pixels_y),
                                 filters = comp)

            # Here it is tricky, we need to guess the number of segments that
            # will be stored in the array, to optimize speed.
            # The current value was used as is for the last 2-3 years without
            # problem, but it could of course be optimized by determining
            # on some randoms curves the number of segments, and the using a
            # mean value ...
            fi.create_vlarray(
                st,
                "stiffness_array",
                ty_fl,
                expectedrows=10 *
                dt.nbr_pixels_x *
                dt.nbr_pixels_y,
                filters=comp)

        if tabletypes == "all" or tabletypes == "work":
            fi.create_carray(st, "jocs1_indices",
                             ty_int,
                             shape=(dt.nbr_pixels_x, dt.nbr_pixels_y),
                             filters = comp)
            fi.create_carray(st, "jocs2_indices",
                             ty_int,
                             shape=(dt.nbr_pixels_x, dt.nbr_pixels_y),
                             filters = comp)
            fi.create_carray(st, "jocs1_real",
                             ty_fl,
                             shape=(dt.nbr_pixels_x, dt.nbr_pixels_y, 2),
                             filters = comp)
            fi.create_carray(st, "jocs2_real",
                             ty_fl,
                             shape=(dt.nbr_pixels_x, dt.nbr_pixels_y, 2),
                             filters = comp)
            fi.create_carray(st, "fits_joc",
                             ty_fl,
                             shape=(dt.nbr_pixels_x, dt.nbr_pixels_y, 3),
                             filters = comp)
            fi.create_carray(st, "work",
                             ty_fl,
                             shape=(dt.nbr_pixels_x, dt.nbr_pixels_y),
                             filters = comp)
            fi.create_carray(st, "rupture_force1",
                             ty_fl,
                             shape=(dt.nbr_pixels_x, dt.nbr_pixels_y),
                             filters = comp)

        if tabletypes == "all" or tabletypes == "events":
            fi.create_vlarray(
                st,
                "events_positions_middle",
                ty_int,
                expectedrows=5 *
                dt.nbr_pixels_x *
                dt.nbr_pixels_y,
                filters=comp)
            fi.create_vlarray(
                st,
                "events_positions_start",
                ty_int,
                expectedrows=5 *
                dt.nbr_pixels_x *
                dt.nbr_pixels_y,
                filters=comp)
            fi.create_vlarray(
                st,
                "events_positions_stop",
                ty_int,
                expectedrows=5 *
                dt.nbr_pixels_x *
                dt.nbr_pixels_y,
                filters=comp)
            fi.create_vlarray(
                st,
                "events_forces",
                ty_fl,
                expectedrows=5 *
                dt.nbr_pixels_x *
                dt.nbr_pixels_y,
                filters=comp)
            fi.create_vlarray(
                st,
                "events_slopes",
                ty_fl,
                expectedrows=5 *
                dt.nbr_pixels_x *
                dt.nbr_pixels_y,
                filters=comp)
            fi.create_carray(st, "events_jocs2_indices",
                             ty_int,
                             shape=(dt.nbr_pixels_x, dt.nbr_pixels_y),
                             filters = comp)
            fi.create_carray(st, "events_jocs2_real",
                             ty_fl,
                             shape=(dt.nbr_pixels_x, dt.nbr_pixels_y, 2),
                             filters = comp)
            fi.create_vlarray(
                st,
                "events_forces_distance",
                ty_fl,
                expectedrows=5 *
                dt.nbr_pixels_x *
                dt.nbr_pixels_y,
                filters=comp)
            fi.create_carray(
                st, "events_fits_joc",
                ty_fl,
                shape=(dt.nbr_pixels_x, dt.nbr_pixels_y, 3),
                filters = comp)

        if tabletypes == "all" or tabletypes == "stiffness_corrected":
            fi.create_vlarray(
                st,
                "stiffness_corrected",
                ty_fl,
                expectedrows=10 *
                dt.nbr_pixels_x *
                dt.nbr_pixels_y,
                filters=comp)

        if tabletypes == "all" or tabletypes == "loading_rates":
            fi.create_vlarray(
                st,
                "loading_rates",
                ty_fl,
                expectedrows=5 *
                dt.nbr_pixels_x *
                dt.nbr_pixels_y,
                filters=comp)


def display_close_file_message(verbose=False):
    """Hides or displays the closing message in the terminal.

    See : http://www.pytables.org/moin/UserDocuments/AtexitHooks
    """
    open_files = tables.file._open_files

    are_open_files = len(open_files) > 0

    if verbose and are_open_files:
        sys.stderr.write("Closing remaining open files:")

    if Version(tables.__version__) >= Version("3.1.0"):
        # make a copy of the open_files.handlers container for the iteration
        handlers = list(open_files.handlers)
    else:
        # for older versions of pytables, setup the handlers list from the
        # keys
        keys = list(open_files.keys())
        handlers = []
        for key in keys:
            handlers.append(open_files[key])

    for fileh in handlers:
        if verbose:
            sys.stderr.write("%s..." % fileh.filename)

        fileh.close()

        if verbose:
            sys.stderr.write("done")

    if verbose and are_open_files:
        sys.stderr.write("\n")
