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

"""Module to get the type of file to load."""

import os
import logging
import zipfile

try:
    import pickle as pickle
except ImportError:
    import pickle
from .. import consts
from .. import widgets_list
from .file_loader.nanoscope import extract_text, extract_num
from .file_loader.jpk import extract_int
from ..widgets.progressbar import Progressbar

encodings = ["latin1", "ASCII", "utf-8"]


class CheckFiles:
    """Class to check the files before loading.

    This class allows to pre-check files before loading.
    The check_single_file method allows to see if a file is a valid.
    It tries to detect if the file is corrupted, incomplete or has any
    other problem which would make pyAF crash. It is of course impossible
    to predict when a file will make pyAF crash, but the most common cases
    are taken into account here.

    (Thanks to the CMIP/BICEL teams in in Lille, which provided most of
    those files during the development of pyAF).

    The self.infos variable is filled with all the needed informations about
    the loaded files and can be accessed from outside the class.
    """

    def __init__(self, files):
        self.logger = logging.getLogger()
        self.logger.debug("Pre-checking files")

        self.infos = []

        # Set first label for the progressbar
        # Get path as string, unicode is there for paths with
        # non ASCII characters
        path = str(files[0])

        Progressbar("Checking files", os.path.basename(path))

        # Check each file
        for i in range(len(files)):
            path = str(files[i])
            self.logger.debug("Checking %s", path)

            widgets_list.widget_progressbar.set_label(os.path.basename(path))

            # Check the file
            result = check_single_file(path)
            self.infos.append(result)
            self.logger.debug("File info : %s", result)

        widgets_list.widget_progressbar.close()


def check_single_file(path):
    """Get the information and type of a single file."""
    filename = os.path.basename(path)

    result = {"enabled": False,
              "checked": False,
              "filename": filename,
              "file_type": "",
              "error": "",
              "version": None,
              "path": path,
              "nbrcurves": 0}

    # Check if the file can be opened
    try:
        with open(path, "rb"):
            pass
    except IOError as err:
        if err.errno == 13:
            result["error"] = "Permission denied"
        else:
            result["error"] = "IOError"
        return result

    file_type = get_file_type(path)
    result["file_type"] = file_type
    nbrcurves = 0

    if file_type == "":
        result["error"] = "Not a valid AFM file"
        return result

    else:
        # If it's a Nanoscope file, check for version and number of points
        if file_type == "Nanoscope (Force Volume)" \
                or file_type == "Nanoscope (Single File)" \
                or file_type == "Nanoscope (Peak Force)":
            nanofile = open(path, mode="r", encoding="latin_1")
            for line in nanofile:
                if line.find("\\Version:") + 1:
                    version = extract_num(line, -1)

                if line.find("\\force/line") + 1:
                    nbr_points_per_line = extract_num(line)
                    break

            nanofile.close()

            if file_type == "Nanoscope (Single File)":
                nbr_points_per_line = 1

            # 128x128 (and more) files have often a wrong version number.
            # This is not a bug, it's only due to bad Veeco/Bruker code
            if version == 5120000 and nbr_points_per_line >= 128:
                result["error"] = "Nanoscope version error"

                if consts.AUTOLOAD:
                    # Disable the error in case of debugging, we want the
                    # file to load anyway.
                    result["error"] = ""

            # The nanoscope files can only be squares
            nbrcurves = nbr_points_per_line * nbr_points_per_line

        if file_type == "JPK (Force Map)" or file_type == "JPK (QI)":
            # Check if the file is a valid zip file
            # It can happen that files get corrupted during transfers
            # between computers or with USB Keys ... often there are some
            # bytes missing in the file.
            try:
                jpkfile = zipfile.ZipFile(path, "r")
            except zipfile.BadZipfile:
                result["error"] = "Corrupted file"

            # Check if there are no curves missing as it is possible to
            # save incomplete files in the JPK software
            if not result["error"]:
                if file_type == "JPK (Force Map)":
                    prefix = "force-scan-map"
                elif file_type == "JPK (QI)":
                    prefix = "quantitative-imaging-map"

                # Get the wanted number of pixels and the real number of pixels
                for line in jpkfile.read("header.properties").splitlines():
                    line = line.decode()
                    theline = prefix + ".position-pattern.grid.ilength"
                    if line.find(theline) != -1:
                        nbr_pixels_x = extract_int(line)
                    theline = prefix + ".position-pattern.grid.jlength"
                    if line.find(theline) != -1:
                        nbr_pixels_y = extract_int(line)
                    theline = prefix + ".indexes.max"
                    if line.find(theline) != -1:
                        real_nbr_pixels = extract_int(line)

                nbrcurves = nbr_pixels_x * nbr_pixels_y

                if real_nbr_pixels != nbrcurves - 1:
                    result["error"] = "Incomplete file"
                    result["enabled"] = True
                    result["checked"] = True

                jpkfile.close()

        if file_type == "pyAF":
            # Read the first line and check if there is a version number
            # If not, then the file could be corrupted (happens sometimes
            # while using USB sticks to transfer files, or due to bad
            # transfers through networks. The filesize can be exactly the
            # same but the content can be messed up. This results in the
            # fact that you can't read the version number from the file.)
            pyaffile = open(path, "rb")
            line = pyaffile.readline().decode()

            if line.split("=")[0] != "Version":
                result["error"] = "Corrupted file"

            # Check if the HDF5 part of the file can be unpickled correctly.
            while line.find("Nbr_of_files=") != 0:
                line = pyaffile.readline().decode()
            nbr_files = int(line.split("=")[1].splitlines()[0])
            while line.find("Global_Prefs_Offset=") != 0:
                line = pyaffile.readline().decode()
            globaloffset = int(line.split("=")[1].splitlines()[0])

            # Get the offsets for each prefs
            single_offsets = []
            singleoffset_total = 0
            for i in range(nbr_files):
                line = pyaffile.readline().decode()
                offset = int(line.split("=")[1].splitlines()[0])
                single_offsets.append(offset)
                singleoffset_total += offset
            offset_prefs = pyaffile.tell()  # Get the current position

            # Go back to the prefs part
            pyaffile.seek(offset_prefs)

            general_prefs = pyaffile.read(globaloffset)

            pyaffile.close()

            for i in range(len(encodings)):
                try:
                    general_prefs = pickle.loads(general_prefs, encoding=encodings[i])
                    consts.ENCODING = encodings[i]
                    break

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
                    from .. import experiment
                    import sys
                    sys.modules["experiment"] = experiment
                    general_prefs = pickle.loads(general_prefs, encoding=encodings[i])
                    consts.ENCODING = encodings[i]
                    break

                except UnicodeDecodeError:
                    if i + 1 < len(encodings):
                        # Go through supported encodings.
                        continue

                    else:
                        # If there are no more encodings, raise the error.
                        result["error"] = "Decode error"


        if result["error"] == "":
            # If no error is found, enable the loading of this file and
            # check the checkbox in the GUI.
            result["enabled"] = True
            result["checked"] = True

        result["nbrcurves"] = nbrcurves
        return result


def get_file_type(path):
    """Tries to guess the file type.

    Different strategies are used until no corresponging file type is found.
    When this happens the function returns and empty string.
    """
    file_type = ""

    # Check if file is not a dir (can happend for .app folder in OS X)
    if os.path.isdir(path):
        return file_type

    # Split the name with a dot to get the extension
    splitted_name = os.path.basename(path).split(".")

    # If there are more dots in the name we have to take the last element of
    # the splitted list to get the extension
    # (JPK has these strange filenames for example)
    extension = splitted_name[-1]

    if extension == "fva" or extension == "pyaf":
        # fva was the old extension, was changed during the 1.3 to 1.4
        # development.
        file_type = "pyAF"
    elif extension == "jpk-force":
        file_type = "JPK (Single File)"
    elif extension == "jpk-force-map":
        file_type = "JPK (Force Map)"
    elif extension == "jpk-qi-data":
        file_type = "JPK (QI)"
    elif extension == "pfc":
        file_type = "Nanoscope (Peak Force)"
    else:
        is_nanoscope = False

        thefile = open(path, "r")
        first_line = thefile.readline()
        if first_line.find("\\*Force file list") + 1:
            is_nanoscope = True
        thefile.close()

        if is_nanoscope:
            thefile = open(path, "r")
            for line in thefile:
                if line.find("\\Start context:") + 1:
                    value = extract_text(line)
                    if value == "FOL":
                        file_type = "Nanoscope (Single File)"
                        break
                    elif value == "FVOL":
                        file_type = "Nanoscope (Force Volume)"
                        break
            thefile.close()

    return file_type
