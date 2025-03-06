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

"""Module containing the class to save the data to .pyaf files."""

import os
try:
    import pickle as pickle
except ImportError:
    import pickle
from ..tools import misc_tools
from ..tools import math_tools
from ..tools import temp_file
from ..tools import apply_to_all
from PyQt5 import QtCore, QtWidgets
import shutil
from .. import consts
from ..widgets.progressbar import Progressbar
from .. import widgets_list
from .. import shared


class Save:
    """This class lets you save the data to a single .pyaf file.

    The hdf5 temporary file is simply copied to a new destination, and
    the preferences coming from the experiment are pickled and written
    as a header of the .pyaf file.

    The file is copied piecewise (100000 bytes at once) to be able to
    have a progressbar during the save.

    A .pyaf file starts with at least these 3 lines ::
        Version=xxx (this line doesn't exist for 1.0.x versions)
        Nbr_of_files=[int]
        Global_Prefs_Offset=[int bytes]
        Single_Prefs_Offset=[int bytes]

    This lets the Loader() class know the number of files to load, and
    were to look for the hdf5 file (After the last Single_Prefs_Offset).
    """

    def __init__(self, parent, filename=None):
        # Setting focus on main window, so no editingFinished signals
        # will be send from the inputs before opening save window
        # (IMPORTANT).
        parent.setFocus()

        if shared.error_triggered:
            # In this case, there was an error while the user worked with pyAF
            # There were many cases were the user saved the file, and tryied
            # to reopen this. The re-opening failed and the user got confused
            # as he did not know from were the bug was comming, making it
            # difficult to debug. To prevent this, tell the user there is a
            # problem, but also allow him to save his file (sometimes the
            # errors are minor, or coming from a plugin).
            # See also #476
            text = (
                "An error occured while you were working with pyAF. It is "
                "not recommended to save your file as it could have been "
                "corrupted.")

            msg = QtWidgets.QMessageBox()
            msg.setText("Error (Saving)")
            msg.setInformativeText(text)
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.addButton("Abort", QtWidgets.QMessageBox.RejectRole)
            msg.addButton("Force Save", QtWidgets.QMessageBox.AcceptRole)
            response = msg.exec_()

            if response == 0:
                return None

        # Save the currently defined parameters if apply to all is checked
        apply_to_all.apply_to_all("all")

        # Get the saved save path from the settings
        self.settings = QtCore.QSettings()
        self.settings_path = misc_tools.get_user_path(self.settings)

        # Get the global preferences
        global_prefs_experiment = shared.exp.__dict__
        global_prefs_to_save = {}
        for key in list(global_prefs_experiment.keys()):
            if not (math_tools.in_list(shared.exp.dontsave, key)
                    and key != "dontsave"):
                global_prefs_to_save[key] = global_prefs_experiment[key]

        newlist = []
        for i in range(len(shared.exp.list)):
            newlist.append(shared.exp.list[i].filename)
        global_prefs_to_save["filenames_list"] = newlist

        # Get the single preferences
        prefs_single = []
        id_list = []
        for dataset in shared.exp.list:
            id_list.append(dataset.unique_id)
            datasets_prefs = {}
            for key in list(dataset.__dict__.keys()):
                if not (math_tools.in_list(dataset.dontsave, key)
                        and key != "dontsave"):
                    datasets_prefs[key] = dataset.__dict__[key]
            prefs_single.append(datasets_prefs)

        if shared.files_removed and consts.UNIT_TESTING is False:
            # Tell the user that this action could be slow

            text = ("You removed files. pyAF will repack the data before "
                    "saving. This could be slow and the progressbar will not "
                    "update. Please wait until the end of the saving process.")

            QtWidgets.QMessageBox.warning(None, "Info", text)

        # Open dialog to save the file
        if consts.UNIT_TESTING is False or filename is None:
            settings = self.settings_path
            filename, _ = QtWidgets.QFileDialog.getSaveFileName(
                parent, "Save data", settings, "Text files (*.pyaf)")
        #else: filename stays as is

        # Check before if we can write to the folder. Can happen it you try
        # to write to a folder were you don't have any write rights
        # Happened for example on Fedora while wanting to write to a HFS+
        # disk. There were read rights but no write rights ... See bug #363
        allow = True
        save_folder = os.path.dirname(os.path.normpath(filename))
        if filename and not os.access(save_folder, os.W_OK):
            if not consts.UNIT_TESTING:
                text = ("You can not write to this folder !")
                QtWidgets.QMessageBox.warning(None, "Info", text)
            allow = False

        # Save the file
        if filename and allow:
            # Restart the timer displaying the save warning
            if consts.UNIT_TESTING is False:
                parent.timer.start()

            # Add extension (for Linux, on OS X it will work out of the box)
            if not filename.endswith(".pyaf"):
                filename += ".pyaf"

            # Get filename as string
            filename = str(filename)

            # Save current path
            misc_tools.set_user_path(self.settings, os.path.dirname(filename))

            # Get filepath of hd5f file
            filepath = shared.exp.temp_file.filepath

            # Make a temporay folder to store the header file
            tmp_dir_path = shared.exp.temp_file.dirpath\
                + shared.exp.temp_file.date
            os.mkdir(tmp_dir_path)
            headerfile = open(tmp_dir_path + "/header.txt", "wb")

            # Save version number
            headerfile.write(str.encode("Version=" + consts.VERSION + "\n"))

            # Save the header
            headerfile.write(str.encode(
                "Nbr_of_files=" + str(len(prefs_single)) + "\n"))

            # Save global prefs
            pickle.dump(global_prefs_to_save,
                        open(tmp_dir_path + "/global_prefs.txt", "wb"),
                        pickle.HIGHEST_PROTOCOL)

            # Write offset for global prefs
            size = os.path.getsize(tmp_dir_path + "/global_prefs.txt")
            headerfile.write(str.encode("Global_Prefs_Offset=" +
                str(size) + "\n"))

            # Save single prefs
            for i in range(len(prefs_single)):
                pref_path = tmp_dir_path + "/single_prefs_" + str(i) + ".txt"
                pickle.dump(
                    prefs_single[i],
                    open(pref_path, "wb"),
                    pickle.HIGHEST_PROTOCOL)

                # Write offset for single prefs
                size = os.path.getsize(pref_path)
                headerfile.write(str.encode("Single_Prefs_Offset=" +
                    str(size) + "\n"))
            headerfile.close()

            # Create a progressbar which is displayed during saving
            Progressbar()

            if shared.files_removed:
                # If a file has been removed, create a new temp_file in
                # the temporary folder. This will remove dead nodes.

                # Set the progressbar's label
                text = "Cleaning up ... please wait"
                widgets_list.widget_progressbar.set_label(text)

                newfile = temp_file.TempFile(tmp_dir_path)
                newfile.create_new_file()

                for data_id in id_list:
                    old = shared.exp.temp_file.file
                    new = newfile.file
                    st = "_" + str(data_id)
                    new.create_group("/data/", st, "Single Data")

                    node = new.get_node("/data/")
                    st = "/data/_" + str(data_id)
                    old.copy_node(st, newparent=node,
                                 overwrite=True, recursive=True)

                newpath = newfile.filepath

                # Close new file before copying
                newfile.close_file()

            # Set the progressbar's label
            widgets_list.widget_progressbar.set_label("Saving ...")

            # Close hdf5 file before copying
            shared.exp.temp_file.close_file()

            # Get the hdf5 file to copy. In normal mode just get the
            # file from the temp folder. In case if files have been
            # removed, use the newly created file for the occasion,
            # which does not contain dead nodes
            if not shared.files_removed:
                hdf5file = open(filepath, "rb")
                size = os.path.getsize(filepath)
            else:
                hdf5file = open(newpath, "rb")
                size = os.path.getsize(newpath)
                # Reset the value for the next time
                shared.files_removed = False

            # Open destination file
            destination = open(filename, "wb")
            blocksize = 100000.0

            # Set range of progressbar
            prog = widgets_list.widget_progressbar
            prog.set_range(0, int(size / blocksize) + 1)

            # Copy the header and the global prefs
            destination.write(open(tmp_dir_path + "/header.txt", "rb").read())
            prefs_path = tmp_dir_path + "/global_prefs.txt"
            destination.write(open(prefs_path, "rb").read())

            # Copy the single prefs
            for i in range(len(prefs_single)):
                prefs_path = tmp_dir_path + "/single_prefs_" + str(i) + ".txt"
                destination.write(open(prefs_path, "rb").read())

            # Copy the hdf5 file piecewise
            for i in range(int(size / blocksize) + 1):
                destination.write(hdf5file.read(int(blocksize)))

                # Update progressbar
                widgets_list.widget_progressbar.update()

            # Close files
            destination.close()
            hdf5file.close()

            # Close progressbar
            widgets_list.widget_progressbar.close()

            # Delete temporary directory
            shutil.rmtree(tmp_dir_path)

            # Reopen temp_file and update values in exp
            shared.exp.temp_file.open_file()

        # Reset focus to main window
        parent.activateWindow()
