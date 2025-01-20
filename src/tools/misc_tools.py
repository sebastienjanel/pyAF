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

"""Miscellaneous tools for pyAF."""

import os
import sys
import copy
import time
import numpy
import ctypes
import platform
import logging
import logging.handlers
from PIL import Image
#from PIL.ImageQt import ImageQt
from PyQt5 import QtCore, QtGui, QtWidgets
import matplotlib
import multiprocessing
from .. import shared
from .. import consts
import psutil
from . import git_tools
from . import events_tools
from ..experiment import ResultGroup, Result



def get_app_path():
    """Gets the current path of the app.

    Split the current's file path and go back 2 levels, to be in main folder.
    """
    splitted_path = os.path.realpath(__file__).split(os.sep)
    path = ""
    for i in range(1, len(splitted_path) - 2):
        path += os.sep + splitted_path[i]

    return path


def get_base_window_title():
    """Get the basic window title for pyAF."""
    # Get git info for the title
    branch, short_sha = git_tools.get_git_info()
    if branch is not None:
        title = " (" + branch + ", " + short_sha + ")"
    else:
        title = ""

    return "pyAF" + title


def create_new_logger():
    """Creates a logger for debug messages.

    The messages are written to log files inside the log's folder.
    Users can then create a new issue in the bugtracker.
    (See main.py in the my_excepthook function, where the LogSender is opened)
    """
    if platform.uname()[0] == "Darwin":
        logfolder = os.path.expanduser("~") + "/Library/Logs/"
        path = logfolder + "pyAF/"

    elif platform.uname()[0] == "Linux":
        logfolder = os.path.expanduser("~") + "/.pyAF/"
        path = logfolder + "logs/"

    elif platform.uname()[0] == "Windows":
        logfolder = os.path.expanduser("~") + r"\AppData\Roaming"
        path = logfolder + r"\pyAF\\"

    disable = False

    # Create directory to store the logs
    try:
        os.makedirs(path)

    except (OSError, IOError):
        try:
            # The only way to know if the folder is writable is to try writing
            # to it ...
            afile = open(path + "test.txt", "w")
            afile.close()
            os.remove(path + "test.txt")

        except (OSError, IOError):
            disable = True

    if disable:
        fh = logging.NullHandler()

    else:
        # Check for old log files from 48 hours before
        for spath in os.listdir(path):
            if spath.split(".")[-1] == "log":
                # 48 h
                splitted = spath.split(".")[0].split("_")
                if len(splitted) == 2:
                    value = int(splitted[1])
                    if time.time() - value > 172800:
                        try:
                            os.remove(path + spath)
                        except OSError:
                            # In some rare cases the rights of the folder could
                            # have changed ... In this case, do nothing.
                            pass

        # Define a logger
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

        date = str(time.time()).split(".")[0]

        # Add a file handler
        path = path + "pyAF_" + date + ".log"
        fh = logging.FileHandler(path, mode="w")
        fh.setLevel(logging.DEBUG)

        # Create a formatter and set the formatter for the handler.
        fstring = "%(asctime)s, %(module)s, Line : %(lineno)s, %(message)s"
        frmt = logging.Formatter(fstring)
        fh.setFormatter(frmt)

    # Add the Handler to the logger
    logger.addHandler(fh)

    return logger


def find_logger_filename(logger):
    """Finds the file where the logger is writing to.

    http://stackoverflow.com/a/7787832
    """
    log_file = None
    h = logger.__dict__['handlers'][0]
    if h.__class__.__name__ == 'FileHandler':
        log_file = h.baseFilename

    return log_file


def setattr_special(obj, name, value):
    """Set the value of an argument with the name of the argument as a string.

    Usage in pyAF :
    Used to load the saved data from the prefs in the experiment class.

    Important : deepcopy is used on the value to have a different
    memory adress for the newly created argument !

    Note : This function is a trick found on the following page :
    http://stackoverflow.com/questions/295058/
    convert-string-to-preexisting-variable-names
    """
    parts = name.split(".")
    for attr in parts[:-1]:
        obj = getattr(obj, attr)
    setattr(obj, parts[-1], copy.deepcopy(value))


def add_attributes(obj):
    """Method used for adding new attributes to old specific pickled objects,
     making them compatible with new features.
     """

    # Add condition attribute to old Result and ResultGroup classes.
    if not hasattr(obj, 'condition') and (isinstance(obj, ResultGroup) or isinstance(obj, Result)):
        default = 0
        setattr(obj, 'condition', default)


def get_free_space_mb(folder):
    """Return folder/drive free space (in bytes)

    http://stackoverflow.com/questions/51658/
    cross-platform-space-remaining-on-volume-using-python
    """
    if platform.uname()[0] == "Windows":
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(str(folder, 'utf-8')),
                                                   None, None,
                                                   ctypes.pointer(free_bytes))

        return free_bytes.value / 1024 / 1024

    else:
        st = os.statvfs(folder)

        return st.f_bavail * st.f_frsize / 1024 / 1024


def save_used_params(data, mode):
    """Saved current parameters to shared.exp."""
    if mode == "stiffness" or mode == "all":
        data.used_stiffness_model_selected = data.stiffness_model_selected
        data.used_tip_radius = data.tip_radius
        data.used_tip_angle = data.tip_angle
        data.used_poisson_ratio = data.poisson_ratio
        data.used_indentation_start = data.indentation_start
        data.used_indentation_stop = data.indentation_stop
        data.used_indentation_step = data.indentation_step
        data.used_force_start = data.force_start
        data.used_force_stop = data.force_stop
        data.used_strict_stop = data.strict_stop
        data.used_tomography = data.tomography
        data.used_fitparam_poc_skip_start = data.fitparam_poc_skip_start
        data.used_fitparam_poc_fit_length = data.fitparam_poc_fit_length
        data.used_fitparam_poc_refit_option = data.fitparam_poc_refit_option
        data.used_fitparam_poc_noise_multiplicator = \
            data.fitparam_poc_noise_multiplicator
        data.used_fitparam_poc_refit_times = data.fitparam_poc_refit_times

    if mode == "work_and_rupture_force" or mode == "all":
        data.used_fitparam_joc_skip_start = data.fitparam_joc_skip_start
        data.used_fitparam_joc_fit_length = data.fitparam_joc_fit_length
        data.used_fitparam_joc_refit_option = data.fitparam_joc_refit_option
        data.used_fitparam_joc_noise_multiplicator = \
            data.fitparam_joc_noise_multiplicator
        data.used_fitparam_joc_refit_times = data.fitparam_joc_refit_times

    if mode == "events" or mode == "all":
        data.used_fitparam_events_joc_skip_start = \
            data.fitparam_events_joc_skip_start
        data.used_fitparam_events_joc_fit_length = \
            data.fitparam_events_joc_fit_length
        data.used_fitparam_events_joc_refit_option = \
            data.fitparam_events_joc_refit_option
        data.used_fitparam_events_joc_noise_multiplicator = \
            data.fitparam_events_joc_noise_multiplicator
        data.used_fitparam_events_joc_refit_times = \
            data.fitparam_events_joc_refit_times
        data.used_events_kernel_detection_threshold = \
            data.events_kernel_detection_threshold
        data.used_events_msf_detection_threshold = \
            data.events_msf_detection_threshold
        data.used_events_msf_window_size = \
            data.events_msf_window_size
        data.used_events_algorithm = data.events_algorithm
        data.used_events_fit_size = data.events_fit_size
        data.used_kernel = events_tools.create_kernel(data.kernel_size)
        data.used_kernel_size = data.kernel_size
        data.used_adaptive_threshold_option = \
            data.adaptive_threshold_option
        data.used_adaptive_smoothing_window = \
            data.adaptive_smoothing_window
        data.used_adaptive_smoothing_order = \
            data.adaptive_smoothing_order

        data.used_adaptive_threshold_option = \
            data.adaptive_threshold_option
        data.used_events_kernel_adaptive_threshold = \
            data.events_kernel_adaptive_threshold

    # Values saved for all three types
    # Spring constant is given in N/nm, convert it back to N/m to save it
    data.used_spring_constant = data.spring_constant * 1e9
    data.used_deflection_sensitivity = data.deflection_sensitivity
    data.used_tilt_limit_1 = data.tilt_limit_1
    data.used_tilt_limit_2 = data.tilt_limit_2
    data.used_tilt_applied = data.tilt_applied
    data.used_stretch_app_lim1 = data.stretch_app_lim1
    data.used_stretch_app_lim2 = data.stretch_app_lim2
    data.used_stretch_len_app = data.stretch_len_app
    data.used_stretch_ret_lim1 = data.stretch_ret_lim1
    data.used_stretch_ret_lim2 = data.stretch_ret_lim2
    data.used_stretch_len_ret = data.stretch_len_ret
    data.used_sg_smoothing_enabled = data.sg_smoothing_enabled
    data.used_sg_smoothing_order = data.sg_smoothing_order
    data.used_sg_smoothing_width = data.sg_smoothing_width
    data.used_sg_smoothing_uniform = data.sg_smoothing_uniform


def get_hist_refresh_range(exp, loading=False):
    """Get a range for the progressbar during the results refresh.

    The grouped results preparation is very fast so no need to display a
    progressbar for it.
    """
    # At load, if the results are not refreshed, no need to display a
    # progressbar
    if loading:
        settings = QtCore.QSettings()
        refresh_hist_at_load = settings.value("refreshHistsAtLoad")
        if not refresh_hist_at_load:
            return 0

    # Progressbar for single results
    progressbar_range = 0
    for i in range(len(exp.results_list)):
        result = exp.results_list[i]
        data_id = result.data_id
        if exp.results_type == "stiffness":
            if exp.list[data_id].stiffness_calculated:
                length = len(exp.list[data_id].stiffness_array[0])
                progressbar_range = progressbar_range + length
        elif exp.results_type == "stiffness_corr":
            if exp.list[data_id].stiffness_corrected:
                length = len(exp.list[data_id].stiffness_array[0])
                progressbar_range = progressbar_range
        elif exp.results_type == "work" or exp.results_type == "rupture_force":
            if exp.list[data_id].work_and_rupture_force1_calculated:
                length = len(exp.list[data_id].work[0])
                progressbar_range = progressbar_range + length
        elif exp.results_type == "events_forces" or \
                exp.results_type == "events_per_curve" or \
                exp.results_type == "events_rupture_force":
            if exp.list[data_id].events_calculated:
                length = len(exp.list[data_id].events_forces[0])
                progressbar_range = progressbar_range + length
        elif exp.results_type == "loading_rates":
            if exp.list[data_id].loading_rates_calculated:
                length = len(exp.list[data_id].loading_rates[0])
                progressbar_range = progressbar_range + length
        elif exp.results_type == "events_distance":
            if exp.list[data_id].events_calculated:
                length = len(exp.list[data_id].events_forces_distance[0])
                progressbar_range = progressbar_range + length

    return progressbar_range


def get_user_path(settings, first_load=False):
    """Get the saved path from the settings."""
    mpl_path = matplotlib.rcParams["savefig.directory"]
    if first_load:
        mpl_path = None

    # Get the saved save path from the settings
    settings_path = str(settings.value("userPath",
                                       os.path.expanduser("~")))
    if settings_path != mpl_path and mpl_path is not None:
        settings_path = mpl_path

    if os.path.exists(settings_path):
        path = settings_path
    else:
        path = os.path.expanduser("~")

    return path


def set_user_path(settings, path):
    """Save a path in the settings."""
    settings.setValue("userPath", path)

    # Save path for matplotlib
    matplotlib.rcParams["savefig.directory"] = path


def get_meshgrid_type_as_string(mesh_id):
    """Given a meshgrid id, return the type as a string."""
    if mesh_id == 0:
        return "piezo"
    elif mesh_id == 1:
        return "topo"
    elif mesh_id == 2:
        return "stiffness"
    elif mesh_id == 3:
        return "stiffness_slice"
    elif mesh_id == 4:
        return "stiffness_corr"
    elif mesh_id == 5:
        return "stiffness_corr_slice"
    elif mesh_id == 6:
        return "work"
    elif mesh_id == 7:
        return "rupture_force"
    elif mesh_id == 8:
        return "events_per_curve"
    elif mesh_id == 9:
        return "events_rupture_force"


def get_meshgrid_type_by_id(meshgrid_type):
    """Given a meshgrid type, return the type as an id."""
    if meshgrid_type == "piezo":
        mesh_id = 0
    elif meshgrid_type == "topo":
        mesh_id = 1
    elif meshgrid_type == "stiffness":
        mesh_id = 2
    elif meshgrid_type == "stiffness_slice":
        mesh_id = 3
    elif meshgrid_type == "stiffness_corr":
        mesh_id = 4
    elif meshgrid_type == "stiffness_corr_slice":
        mesh_id = 5
    elif meshgrid_type == "work":
        mesh_id = 6
    elif meshgrid_type == "rupture_force":
        mesh_id = 7
    elif meshgrid_type == "events_per_curve":
        mesh_id = 8
    elif meshgrid_type == "events_rupture_force":
        mesh_id = 9
    return mesh_id


def get_color_opts_by_meshgrid_type(data, meshgrid_type):
    """Return the colorscale depending on the meshgrid type."""
    if meshgrid_type == "piezo":
        return data.color_opts_piezo
    elif meshgrid_type == "topo":
        return data.color_opts_topo
    elif meshgrid_type == "stiffness":
        return data.color_opts_stiffness
    elif meshgrid_type == "work":
        return data.color_opts_work
    elif meshgrid_type == "rupture_force":
        return data.color_opts_rupture_force1
    elif meshgrid_type == "events_per_curve":
        return data.color_opts_events_per_curve
    elif meshgrid_type == "events_rupture_force":
        return data.color_opts_rupture_force2
    elif meshgrid_type == "stiffness_slice":
        return data.color_opts_stiffness
    elif meshgrid_type == "stiffness_corr":
        return data.color_opts_stiffness
    elif meshgrid_type == "stiffness_corr_slice":
        return data.color_opts_stiffness


def get_results_type_from_id(results_type):
    """Return the results type as an id, given a string."""
    if results_type == 0:
        result_id = "stiffness"
    elif results_type == 1:
        result_id = "work"
    elif results_type == 2:
        result_id = "rupture_force"
    elif results_type == 3:
        result_id = "events_forces"
    elif results_type == 4:
        result_id = "events_per_curve"
    elif results_type == 5:
        result_id = "events_rupture_force"
    elif results_type == 6:
        result_id = "loading_rates"
    elif results_type == 7:
        result_id = "events_distance"

    return result_id


def get_cores():
    """Function used to determine the number of cores aviable,"""
    cores = multiprocessing.cpu_count()

    if consts.UNIT_TESTING is False:
        return int(cores)
    else:
        # The unit tests segfault often due to problems
        # with process forks on OS X. Remove the segfaults
        # by setting the cores to 1.
        return 1


def get_processes(nbr_pixels):
    """Get the number of processes the software can spawn for the multiprocessing.

    The data has to be evenly shared through the processes, so the number of
    processes will depend on the number of curves and the numbers of cores.

    On windows the maximum processes is set to 4 for performance reasons.
    (See documentation).
    """
    cores = shared.exp.cores

    # Define the number of processes to spawn.
    if nbr_pixels < cores:
        nbr_processes = nbr_pixels
    else:
        nbr_processes = cores

    # Reduce number of processes if needed
    # (for example using 24 cores on 4096 curves will not work).
    # The data has to be evenly distributed to every process.
    while nbr_pixels % nbr_processes != 0:
        nbr_processes -= 1
        if nbr_processes == 1:
            break

    # On windows allow only the creation of maximum 4 processes, spawning
    # more processes is slow.
    if platform.uname()[0] == "Windows":
        if nbr_processes != 1:
            nbr_processes = 4
            while nbr_processes > cores:
                nbr_processes = nbr_processes / 2

    return int(nbr_processes)


def get_nbr_blocks(cols, rows):
    """Get the number of blocks for the computations.

    When computing the stiffness and the work, the data is read and written
    in blocks. Define how many blocks there are here. Having too many blocks
    slows down the computation.
    """
    nbr_total = cols * rows
    free_mem = get_free_memory()

    if free_mem / 1000.0 > 6000:
        max_curves = 128 * 128
    elif free_mem / 1000.0 > 2000:
        max_curves = 64 * 64
    else:
        max_curves = 64 * 32

    # Max size of a block : 600 Meg (seems reasonable)
    # (600*1e6)/(24*2*2*3000) = 2083.3333333333335 curves
    # 24 bytes = 1 float in memory
    # 2 = approach - retraction
    # 2 0 defl - ext
    # 3000 points per curve
    # Lets say maximum 2048 curves per block (64x32 file)

    nbr_of_blocks = int(nbr_total) // max_curves
    if nbr_of_blocks == 0:
        nbr_of_blocks = 1

    return nbr_of_blocks


def get_free_memory():
    """Return the quantity of free RAM in megabytes.

    Needs the psutil library.
    """
    return int(psutil.virtual_memory()[4]) / 1e6


def get_tiff_calibration(path):
    """Gets the calibration of Tiff files coming from various microscopes.

    These tiffs have headers with a lot of \x00 (null) characters. The
    function will parse the text and try to get the calibration out of this
    mess. It is not the cleanest solution but it works ...

    If no calibration is found, the function will return None, None.
    Else it returns the size of a pixel in meters.

    The supported files are tiff files coming from the TEM of the Bicel
    Platform in Lille, and LSM/CZI files converted in Tiffs with imageJ.
    """
    string = ""

    with open(path, "rb") as f:
        # The information is stored in the first 1000 bytes
        nbr_bytes = f.read(1000)

        val = nbr_bytes.split(b'\x00')
        string = ""
        for char in val:
            string += char.decode("latin-1")

        splitted = repr(string).split("XpixCal=")

        if len(splitted) > 1:
            # It's a TEM file

            splitted = repr(string).split("XpixCal=")[1]

            splitted = splitted.split("Unit")[0]
            xcal = splitted.split("\\r")[0]
            ycal = splitted.split("\\r")[1].split("=")[1]

            xcal = round(1000.0 / float(xcal.replace(",", ".")), 3) * 1e-9
            ycal = round(1000.0 / float(ycal.replace(",", ".")), 3) * 1e-9

            return [xcal, ycal, 0]

        try:
            # Is it an LSM file converted in Tiff ?

            # Go back to 0
            f.seek(0)

            # The information is stored in the first 3000 bytes
            nbr_bytes = f.read(3000)
            val = nbr_bytes.split(b'\x00')
            string = ""
            for char in val:
                string += char.decode("latin-1")
            string = repr(string)

            # \xb5m == um
            string1 = string.split("Voxel_size_X: ")[1]
            xcal = float(string1.split(" \\xb5m")[0]) * 1e-6
            string2 = string.split("Voxel_size_Y: ")[1]
            ycal = float(string2.split(" \\xb5m")[0]) * 1e-6
            string3 = string.split("Voxel_size_Z: ")[1]
            zcal = float(string3.split(" \\xb5m")[0]) * 1e-6

            return [xcal, ycal, zcal]

        except IndexError:
            pass

        try:
            # Is it a CZI file converted in Tiff ?

            # Go back to 0
            f.seek(0)

            header = f.read(22000)
            header = repr(header).replace("\\x00", "")
            header = header.split("Scaling|Distance|Id")[1]
            header = header.split("Scaling|Distance|Value")[1]
            xcal = float(header.split("=")[1].split("\\n")[0].split(" ")[1])
            ycal = xcal

            return [xcal, ycal, 0]

        except IndexError:
            pass

        # Return nothing
        return [None, None, None]


#def get_pixmap(array):
#    """Makes a QPixmap from a numpy array."""
#    image = Image.fromarray(numpy.uint8(array))
#    return QtGui.QPixmap.fromImage(ImageQt.ImageQt(image))

def get_pixmap(array):
    """Makes a QPixmap from a numpy array."""
    # Convert the numpy array to a PIL image
    image = Image.fromarray(numpy.uint8(array))

    # Convert the PIL image to a QImage
    image = image.convert("RGBA")
    data = image.tobytes("raw", "RGBA")
    qimage = QtGui.QImage(data, image.size[0], image.size[1], QtGui.QImage.Format_RGBA8888)

    # Convert the QImage to a QPixmap
    return QtGui.QPixmap.fromImage(qimage)


def make_path(pos_a, pos_b):
    """Get the coordinates of the points for a profile.

    Original function by Charles Roduit (OpenFovea)
    """
    pos_a[0] = int(pos_a[0])
    pos_a[1] = int(pos_a[1])
    pos_b[0] = int(pos_b[0])
    pos_b[1] = int(pos_b[1])

    if pos_b[0] == pos_a[0]:
        if pos_a[1] > pos_b[1]:
            pos_a, pos_b = pos_b, pos_a
        length = pos_b[1] - pos_a[1] + 1
        y = list(range(pos_a[1], pos_b[1] + 1))
        x = [pos_b[0]] * length
    elif pos_b[1] == pos_a[1]:
        if pos_a[0] > pos_b[0]:
            pos_a, pos_b = pos_b, pos_a
        length = pos_b[0] - pos_a[0] + 1
        x = list(range(pos_a[0], pos_b[0] + 1))
        y = [pos_b[1]] * length
    else:
        a = float(pos_b[1] - pos_a[1]) / (pos_b[0] - pos_a[0])
        b = pos_a[1] - a * pos_a[0]

        if abs(pos_b[0] - pos_a[0]) > abs(pos_b[1] - pos_a[1]):
            if pos_a[0] > pos_b[0]:
                pos_a, pos_b = pos_b, pos_a
            # x0 -> x1 is bigger to y0 -> y1
            x = list(range(pos_a[0], pos_b[0] + 1))
            y = [round(a * item + b) for item in x]
        else:
            if pos_a[1] > pos_b[1]:
                pos_a, pos_b = pos_b, pos_a
            # y0 -> y1 is bigger than x0 -> x1
            y = list(range(pos_a[1], pos_b[1] + 1))
            x = [round((item - b) / a) for item in y]

    return numpy.asarray([x, y])


def reset_computation(data_id, computation_type):
    """Resets results arrays."""
    data = shared.exp.list[data_id]

    if computation_type == "stiffness" or "all":
        if data.stiffness_calculated:
            data.stiffness_calculated = False

            shared.exp.temp_file.delete_tables_for_results(str(data_id),
                                                           "stiffness", True)
            shared.exp.temp_file.create_tables_for_results(str(data_id),
                                                           "stiffness", True)

    if computation_type == "stiffness_corrected" or "all":
        if data.stiffness_corrected:
            data.stiffness_corrected = False

            shared.exp.temp_file.delete_tables_for_results(
                str(data_id),
                "stiffness_corrected")
            shared.exp.temp_file.create_tables_for_results(
                str(data_id),
                "stiffness_corrected")

    if computation_type == "work_and_rupture_force1" or "all":
        if data.work_and_rupture_force1_calculated:
            data.work_and_rupture_force1_calculated = False

            shared.exp.temp_file.delete_tables_for_results(
                str(data_id),
                "work")
            shared.exp.temp_file.create_tables_for_results(
                str(data_id),
                "work")

    if computation_type == "events" or "all":
        if data.events_calculated:
            data.events_calculated = False

            shared.exp.temp_file.delete_tables_for_results(
                str(data_id),
                "events")
            shared.exp.temp_file.create_tables_for_results(
                str(data_id),
                "events")

    if computation_type == "loading_rates" or "all":
        if data.loading_rates_calculated:
            data.loading_rates_calculated = False

            shared.exp.temp_file.delete_tables_for_results(
                str(data_id),
                "loading_rates")
            shared.exp.temp_file.create_tables_for_results(
                str(data_id),
                "loading_rates")

    data.display_curve_type = "defl_ext"
    data.update()


def get_position_on_meshgrid(event):
    """Method to get the position of a click on the meshgrid.

    The curve number is computed and 1 is added (because curves are numbered
    starting at 1 and not 0).
    """
    data = shared.exp.current_data

    if data.meshgrid_units == "um":
        meshfactor = 1000.0
    elif data.meshgrid_units == "nm":
        meshfactor = 1.0

    if data.scan_size_x != 0:
        scan_size_x = float(data.scan_size_x / meshfactor)
        scan_size_y = float(data.scan_size_y / meshfactor)

        xsize = scan_size_x / float(data.nbr_pixels_x)
        ysize = scan_size_y / float(data.nbr_pixels_y)

    else:
        xsize = 1
        ysize = 1

    return int(event.xdata / xsize) + 1, int(event.ydata / ysize) + 1


def validator(name, reg=None, bottom=None, top=None, decimal=None):
    """Validate the inputs.

    For unsigned values, just set the minimum to 0, the maxium stays at
    +infinity
    """
    # Positive integer
    if name == "UI":
        val = QtGui.QIntValidator()
        val.setBottom(0)

    # Integer between 'bottom' and 'top'
    elif name == "RI":
        val = QtGui.QIntValidator()
        val.setRange(bottom, top)

    # All integers
    elif name == "I":
        val = QtGui.QIntValidator()

    # Positive float
    elif name == "UF":
        val = QtGui.QDoubleValidator()
        val.setBottom(0)

    # Float between 'bottom' and 'top' and max decimal of 'decimal'
    elif name == "RF":
        val = QtGui.QDoubleValidator()
        val.setRange(bottom, top, decimal)

    # All floats
    elif name == "F":
        val = QtGui.QDoubleValidator()

    # Specific regular expression
    elif name == "R":
        val = QtGui.QRegExpValidator(reg)

    return val


def ask_user_for_color(current_color):
    """Opens a color dialog asking the user to chose  color."""
    current_color = QtGui.QColor(
        int(current_color[0] * 255), int(current_color[1] * 255), int(current_color[2] * 255))

    color = QtWidgets.QColorDialog.getColor(current_color)
    if color.isValid():
        color = color.getRgb()[:-1]
        return [color[0] / 256.0, color[1] / 256.0, color[2] / 256.0, 1.0]
    else:
        return False


def update_indentation_list(unique_id):
    """Update the values in the indentation list."""
    data = shared.exp.list[unique_id]

    # Create an array with text, defining slices of indent
    # for the list in the results tab
    data.indentation_list = []
    step = data.indentation_step
    start = data.indentation_start
    for indentation in range(data.max_indentation_index):
        if data.tomography:
            start_z = str(indentation * step + start)
        else:
            start_z = str(0 + start)
        if step:
            stop = str(indentation * step + step + start) + " nm"
        else:
            stop = "max"
        text = str(start_z) + " - " + stop
        data.indentation_list.append(text)


def update_slice_position(unique_id):
    """Update the position of the slice after a computation.

    Check if the currently selected slice still can be
    selected. It can be that this time we have less.
    In this case go back to the 0 (all) slice
    Check also the special case were we have only one slice
    then we also go to 0.
    """
    data = shared.exp.list[unique_id]

    for res in shared.exp.results_list:
        if res.data_id == unique_id:
            if res.slice > len(data.indentation_list):
                res.slice = 0
            if len(data.indentation_list) == 1:
                res.slice = 0


def get_resource_path(relative_path):

    extension = os.path.splitext(relative_path)[1]

    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)

    elif extension == ".png":
        return os.path.join(os.path.abspath("."), "src/images", relative_path)

    elif extension == ".mplstyle":
        return os.path.join(os.path.abspath("."), "src/gui_styles/plot_styles", relative_path)

    elif extension == ".qss":
        return os.path.join(os.path.abspath("."), "src/gui_styles/stylesheets", relative_path)