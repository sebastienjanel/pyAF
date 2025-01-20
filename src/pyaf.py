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

"""Main pyAF file, entry point to launch the app.

Once main() is called, the entry widget is displayed. When the user has
selected the file to open, the PyQT app is created, the data is loaded, then
the GUI is loaded. All the main widgets are created here. Widgets which will
open in a new window are registered in wl.

The libraries (except multiprocessing and QtGui) are only imported in the main
thread. This is due to the fact that on windows, each computation process is
spawns a new python process, which will reload all the libraries, slowing down
the computation.
QtGui has to be loaded everytime because the MainWindow class subclasses
a QtWidgets.QMainWindow.
"""

from pyAF.src import consts


class FetchStderr:
    """Utilitary class to hide ITK warnings at load.

    Should be removed once ITK warnings are fixed. The warnings can be
    displayed when in DEBUG mode.
    """

    def write(self, err):
        """Writes the error to stdout."""
        # pylint: disable=R0201
        if consts.DEBUG:
            sys.stdout.write(err)

import multiprocessing

from PyQt5 import QtGui, QtWidgets

if multiprocessing.current_process().name == "MainProcess":
    # Import only once (for performance issue with multiprocessing on windows)
    import sys
    import os
    import logging
    from PyQt5 import QtCore

    if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
        QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

    if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
        QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

    # Import matplotlib the first time to define the backend (You won't be able
    # to open plots in new windows in Linux if you don't define
    # this here at the beginning of the script.
    # The backend change is conditional because Sphinx (documentation system)
    # will display an error message if we want to apply the backend two times.
    import matplotlib
    matplotlib.use("Qt5Agg")

    # Check for the versions. See PEP 386 for the check
    # [AD: deprecated under PEP 386 (!)]
    from packaging.version import Version
    import tables
    import platform
    from pyAF.src.tools.utils import module_exists

    """if Version(tables.__version__) < Version("3.0"):
        message = "PyTables version is " + tables.__version__ + "." \
                  " Version 3.0.x is needed."
        sys.exit(message)"""
    """if Version(matplotlib.__version__) < Version("1.4.0"):
        message = "Matplotlib's version is " + matplotlib.__version__ + "." \
            " Version 1.4.0 is needed."
        sys.exit(message)"""
    qt_version = Version(QtCore.QT_VERSION_STR)
    if platform.uname()[0] == "Darwin" and qt_version < Version("4.8.6"):
        # For mac we need 4.8.6
        message = "Qt's version is " + QtCore.QT_VERSION_STR + "." \
            " Version 4.8.6 is needed."
        sys.exit(message)
    else:
        if qt_version < Version("4.8.5"):
            message = "Qt's version is " + QtCore.QT_VERSION_STR + "." \
                " Version 4.8.5 is needed."
            sys.exit(message)
    if module_exists("vtk"):
        import vtk
        if vtk.VTK_MAJOR_VERSION <= 5:
            val = str(vtk.VTK_MAJOR_VERSION)
            message = "VTK's version is " + val + ". Version 6.x is needed."
            print(message)
            del vtk
        else:
            consts.ALLOW_VTK = True
    if module_exists("itk"):
        original_stderr = sys.stderr  # Keep a reference to STDERR
        sys.stderr = FetchStderr()  # Redirect the real STDERR
        import itk
        try:
            if Version(itk.Version.GetITKVersion()) < Version("4.5.0"):
                val = str(itk.Version.GetITKVersion())
                message = "ITK's version is " + \
                    val + ". Version 4.5.x is needed."
                print(message)
                del itk
            else:
                consts.ALLOW_ITK = True
        except ImportError:
            message = "ITK is installed but seems misconfigured."
            consts.ALLOW_ITK = False
        sys.stderr = original_stderr  # Set back STDERR

    from pyAF.src.utils import menu
    from pyAF.src.widgets_main.data import DataWidget
    from pyAF.src.widgets_main.compute import ComputeWidget
    from pyAF.src.widgets_main.results import ResultsWidget
    from pyAF.src.widgets_main.results_single import ResultsSingleWidget
    from pyAF.src.widgets_main.results_groups import ResultsGroupsWidget
    from pyAF.src.widgets_main.results_experiment import ResultsExperimentWidget
    from pyAF.src.widgets.entry import EntryWidget
    from pyAF.src.widgets.log_sender import LogSenderWidget
    from pyAF.src.widgets.about import AboutWidget
    from pyAF.src.tools import apply_to_all
    from pyAF.src.tools import misc_tools
    from pyAF.src.tools import math_tools as mt
    from pyAF.src.tools import temp_file
    from pyAF.src.widgets.progressbar import Progressbar
    from pyAF.src.experiment import Experiment
    from pyAF.src.load_and_save.load import Load
    import random
    random.seed()  # Will use current system time
    from pyAF.src import widgets_list as wl
    from pyAF.src import shared
    from pyAF.src.load_and_save.pre_check import CheckFiles
    from pyAF.src.widgets.plugins import PluginsWidget
    from pyAF.src.widgets.preferences import PreferencesWidget
    from pyAF.src.utils import discarding
    from pyAF.src.gui_styles.theme_handler import theme_handler


class MainWindow(QtWidgets.QMainWindow):
    """The MainWindow class is the main window of pyAF.

    It will hold the shared.exp instance, containing all the force map
    data, the results and all the options and preferences defined by the
    user.

    The MainWindow is divided in 5 tabs :
    Data, Compute, Results, Results Single, Results Groups.
    """

    def __init__(self, files_info):
        super().__init__()
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Minimum,
            QtWidgets.QSizePolicy.Minimum)

        # Register the window in widgets_list
        wl.widget_main = self

        # Some operations are done only on first GUI load
        self.first_load = True

        # Load the logger
        self.logger = logging.getLogger()
        self.logger.debug("Entering MainWindow")

        # Load the settings
        self.settings = QtCore.QSettings()

        # shared.exp stores the data
        shared.exp = Experiment()

        # Load files in shared.exp
        Load(self, files_info)

        self.current_tab = None

        Progressbar()
        wl.widget_progressbar.set_label("Loading GUI ...")
        progressbar_range = misc_tools.get_hist_refresh_range(shared.exp, True)
        wl.widget_progressbar.set_range(0, progressbar_range + 31)

        self.default_font = QtWidgets.QApplication.font()
        self.default_font_name = self.default_font.defaultFamily()
        self.smallfont = QtWidgets.QApplication.font()  # Default font
        self.smallfont.setPointSize(10)
        self.canvas_resolution = 72
        self.update_results = False

        # Tabs
        self.tabs = QtWidgets.QTabWidget()

        # Create the 5 main widgets (will automatically be registered in wl)
        DataWidget(self)
        ComputeWidget(self)
        ResultsWidget(self)
        ResultsSingleWidget(self)
        ResultsGroupsWidget(self)
        ResultsExperimentWidget(self)

        self.tabs.addTab(wl.widget_data, "Data")
        self.tabs.addTab(wl.widget_compute, "Compute")
        self.tabs.addTab(wl.widget_results, "Results")
        self.tabs.addTab(wl.widget_results_single, "Plots")
        self.tabs.addTab(wl.widget_results_groups, "Grouped Plots")
        self.tabs.addTab(wl.widget_results_experiment, "Experiment Plots")

        # Create the main menu of pyAF
        self.list_actions = menu.create_menu(self)

        # Register the tabs_changed action after the creation of the menu
        self.tabs.currentChanged.connect(self.tabs_changed)
        self.tabs_changed()  # Call it to update the figure menu

        # Some automatic things to do at load if the options are set to
        # True in consts.py

        if consts.DEBUG:
            self.tabs.setCurrentIndex(consts.AUTOTAB)

        if consts.AUTOCALC_STIFFNESS:
            wl.widget_compute.button_clicked("get_stiffness")

        if consts.AUTOCALC_WORK:
            wl.widget_compute.button_clicked("get_work_and_rupture_force")

        if consts.AUTOCALC_EVENTS:
            shared.exp.list[0].nbr_events_found = 9
            wl.widget_compute.button_clicked("get_events")

        if consts.AUTO3D:
            wl.widget_results.button_clicked("button_display_in_3D")

            if consts.AUTO_LOAD_TIFF_IN_3D:
                wl.widget_vtk.menu.tab_Global.button_clicked("load")

            if consts.AUTO_LOAD_STACK_IN_3D:
                wl.widget_vtk.menu.tab_Global.button_clicked("load_stack")

        if consts.AUTOCALC_LOADING_RATE:
            wl.widget_compute.button_clicked("get_lr")

        if consts.OPENLOGSENDER:
            logwidget = LogSenderWidget(None)
            logwidget.show()

        self.setCentralWidget(self.tabs)
        id_selected = shared.exp.id_selected
        wl.widget_results.list_exp.setCurrentIndex(id_selected)
        wl.widget_compute.list_exp2.setCurrentIndex(id_selected)
        if wl.widget_progressbar is not None:
            wl.widget_progressbar.update()

        self.update_GUI("tabs")

        # Close progressbar and display message + beep
        if wl.widget_progressbar is not None:
            wl.widget_progressbar.close()
        if consts.AUTOLOAD is False and consts.UNIT_TESTING is False:
            QtWidgets.QApplication.beep()
            QtWidgets.QMessageBox.warning(self, "Information", "Loading done.")

        self.logger.debug("Opening done")

        self.setFocus()

        self.first_load = False

        # Do not create and link timers while running the unit tests,
        # this is often triggering segfaults.
        if consts.UNIT_TESTING is False:
            # Create a timer, calling every 5 minutes the method to display a
            # message to the user
            self.tray_enabled = True
            self.timer = QtCore.QTimer()
            self.timer.timeout.connect(self.display_message)
            self.timer.start(1000 * 60 * 5)

            # Register a shortuct to discard or enable the current curve
            seq = QtGui.QKeySequence("Ctrl+D")
            QtWidgets.QShortcut(seq, self, discard_shortcut)

            self.current_tab = self.tabs.currentIndex()

            # Refresh blitting
            if self.tabs.currentIndex() == 0:
                mesh = wl.widget_data.meshgrid
                self.timer2 = QtCore.QTimer.singleShot(
                    0, mesh.refresh_blitting)
            elif self.tabs.currentIndex() == 1:
                mesh = wl.widget_compute.MPL_meshgrid
                self.timer3 = QtCore.QTimer.singleShot(
                    0, mesh.refresh_blitting)
            elif self.tabs.currentIndex() == 2:
                mesh = wl.widget_results.MPL_canvas1
                self.timer4 = QtCore.QTimer.singleShot(
                    0, mesh.refresh_blitting)

    def display_message(self):
        """Displays a message telling the user he needs to save.

        Display the message every 5 minutes. Will work only on an OS
        supporting this kind of messages. Can be disabled in the
        preferences.
        """
        enabled = self.settings.value("RememberMeToSave", True)
        isavailable = QtWidgets.QSystemTrayIcon.isSystemTrayAvailable()

        if isavailable and self.tray_enabled and enabled:
            tray = QtWidgets.QSystemTrayIcon(self)
            icon = QtGui.QIcon("images/icon.png")
            tray.setIcon(icon)

            # On Windows, allow the creation of a new tray only if the user
            # has clicked on the old one and removed it.
            # Will not be triggered on OS X
            tray.messageClicked.connect(self.enable_tray)

            tray.show()
            tray.showMessage(
                "pyAF", "Do not forget to save your file !",
                QtWidgets.QSystemTrayIcon.Information, msecs=1000 * 10)

            # Disable the creation of a new tray on Windows and Linux
            if sys.platform != "darwin":
                self.tray_enabled = False

    def enable_tray(self):
        """Re-enable the tray message when the user clicks on it."""
        self.tray_enabled = True

        # Restart the timer with the same time interval
        self.timer.start()

    def closeEvent(self, event):
        """Override the PyQT closeEvent.

        Overriding the PyQT closeEvent will allow to reset the focus on
        the MainWindow, which will prevent some signals to be triggered
        during the closing (for example some editingFinished signals).

        Besides, it will hide the closing file message from pytables,
        and delete the temporary hdf5 file.
        """
        self.logger.debug("Closing pyAF")

        # Setting focus on main window, so no editingFinished signals will be
        # send from the inputs
        self.setFocus()

        # Define values for Qsettings (Reopen them, because they have been
        # closed BEFORE the closevent, and no value will be saved ...)
        QtCore.QCoreApplication.setOrganizationName("CMIP")
        QtCore.QCoreApplication.setOrganizationDomain("cmip.cnrs.fr")
        QtCore.QCoreApplication.setApplicationName("pyAF")
        settings = QtCore.QSettings()

        # Save current path for the next opening of files
        path = misc_tools.get_user_path(settings)
        misc_tools.set_user_path(settings, path)

        # Hide temp file closing message from terminal
        temp_file.display_close_file_message(False)

        # Delete temp file
        shared.exp.temp_file.delete_temp_file()

        # Reset the main widget to None
        wl.widget_main = None

        # Close the app
        event.accept()

    def keyPressEvent(self, event):
        """Catches key events with the arrow keys.

        Specific to the events related to the curve changing. Allows
        also to discard/re-enable curves with the d key.
        """
        if event.type() == QtCore.QEvent.KeyPress:
            key = event.key()
            if key == QtCore.Qt.Key_Left:
                self.change_curve_by_key("left")
            elif key == QtCore.Qt.Key_Right:
                self.change_curve_by_key("right")
            elif key == QtCore.Qt.Key_Up:
                self.change_curve_by_key("up")
            elif key == QtCore.Qt.Key_Down:
                self.change_curve_by_key("down")

        event.accept()

    def open_pref_widget(self):
        """Open the preferences Window."""
        if wl.widget_preferences is None:
            # Create a new widget
            PreferencesWidget(self)
            wl.widget_preferences.setWindowTitle("Preferences")
            wl.widget_preferences.show()
        else:
            # Bring to front
            wl.widget_preferences.activateWindow()
            wl.widget_preferences.raise_()

    def load_more(self, path=None):
        """Opem a new entry widget, which allows to load more files.

        The files_info gotten by the entry widget is passed to the
        loading class. Can be called by the menu or the + button in the
        data tab.
        """
        self.logger.debug("Load more called")

        # Save the currently defined parameters if apply to all is checked
        apply_to_all.apply_to_all("all")

        # Get the path from the settings, to be in the same folder as before.
        settings_path = misc_tools.get_user_path(self.settings)
        if consts.AUTOLOAD is False:
            load_dir_path = settings_path
        else:
            load_dir_path = consts.TEST_PATH

        # Open a file dialog to open files
        if not consts.UNIT_TESTING:
            path = load_dir_path
            # Ask a first time for files
            files, _ = QtWidgets.QFileDialog.getOpenFileNames(directory=path)

            if files:
                value = True
                # Loop as long as we are not ready.
                while value:
                    value, files_info, _ = \
                        open_entry_widget("load_more", files)
                    files = add_more(value, path, files)
            else:
                # Cancel was clicked or no file was selected, do nothing
                return False

        else:
            # Path is passed as argument
            files = []
            files.append(path)
            files_info = CheckFiles(files).infos

        if files_info is not None:
            # Load the files
            Load(self, files_info, load_more=True)

            # Close progressbar
            wl.widget_progressbar.close()

    def open_plugin_window(self):
        """Opens the plugins widget.

        Is called from the menu. See the plugins widget for more
        details (widgets/plugins.py).
        """
        self.logger.debug("Opening plugin window")

        if wl.widget_plugins is None:
            PluginsWidget(self)
            wl.widget_plugins.show()
        else:
            wl.widget_plugins.activateWindow()
            wl.widget_plugins.raise_()

    def change_curve_by_key(self, key):
        """Changes the curve upon a call by the arrow keys."""
        self.logger.debug("Changing curve by key press (%s)", key)

        current_xpos = shared.exp.meshgrid_click_xpos
        current_ypos = shared.exp.meshgrid_click_ypos

        data = shared.exp.current_data

        change = False
        if key == "left":
            if current_xpos - 1 != -1:
                shared.exp.meshgrid_click_xpos -= 1
                change = True
        elif key == "right":
            if current_xpos + 1 != data.nbr_pixels_x:
                shared.exp.meshgrid_click_xpos += 1
                change = True
        elif key == "up":
            if current_ypos + 1 != data.nbr_pixels_y:
                shared.exp.meshgrid_click_ypos += 1
                change = True
        elif key == "down":
            if current_ypos - 1 != -1:
                shared.exp.meshgrid_click_ypos -= 1
                change = True

        if change:
            self.change_curve(None, None, None)

    def change_curve(self, xpos, ypos, rand=None,
                     save_rand=False, save_in_last_ten=True):
        """Displays the asked curve (or a random curve).

        xpos, ypos are the new positions on the map, in user coordinates
        (from 1 to nbr_pixels_x)

        Can be called from the input fields to go to a specific curve,
        or from the random button. GUI updates and plot updates are done
        after the new xpos and ypos are saved as
        shared.exp.meshgrid_click_xpos and
        shared.exp.meshgrid_click_ypos.

        Can also be called from the change_curve_by_key method, then
        xpos and ypos are set to None because the positions were
        previously set by the change_curve_by_key function.

        If enabled in the preferences, the random part will save already
        seen curves and will no more display them, until all the curves
        have been seen.
        """
        string = "Changing for curve %s, %s, rand = %s"
        self.logger.debug(string, str(xpos), str(ypos), str(rand))

        # Reset zoomed in position in compute tab
        shared.zoomed_in_element = 0

        data = shared.exp.current_data

        # Allow or disable random control
        random_control = self.settings.value("RandomControl", False)

        # If rand is given, chose a new xpos and ypos
        if rand == "rand":
            nbr_pixels_x = data.nbr_pixels_x
            nbr_pixels_y = data.nbr_pixels_y

            xpos = random.randint(1, nbr_pixels_x)
            ypos = random.randint(1, nbr_pixels_y)

            string = "Go to rand curve %s, %s"
            self.logger.debug(string, str(xpos), str(ypos))

            if random_control:
                # Check if the position has already been parsed by the
                # random
                list_curves = shared.exp.parsed_random_curves

                if mt.in_list(list_curves, [xpos, ypos]):
                    # If there are only a limited number of curves left
                    # go trough them one by one and do not try a new
                    # random, this could take very long until the script
                    # finds the last curves ...
                    nbr_of_curves_total = nbr_pixels_x * nbr_pixels_y
                    limit = int((nbr_of_curves_total / 2.0))

                    if len(list_curves) == nbr_of_curves_total:
                        shared.exp.parsed_random_curves = []
                        self.change_curve(1, 1)

                    elif len(list_curves) > nbr_of_curves_total - limit:
                        # Find a curve not already seen
                        for i in range(nbr_pixels_x):
                            for j in range(nbr_pixels_y):
                                if mt.in_list(list_curves, [i, j]) is False:
                                    self.change_curve(i, j, save_rand=True)
                                    return

                    else:
                        # If already seen, search for a new valid position
                        notfound = True
                        while notfound:
                            xpos = random.randint(1, nbr_pixels_x)
                            ypos = random.randint(1, nbr_pixels_y)

                            if mt.in_list(list_curves, [xpos, ypos]) is False:
                                notfound = False

                        self.change_curve(xpos, ypos, rand="rand")

        if rand == "rand" or save_rand and random_control:
            shared.exp.parsed_random_curves.append([xpos, ypos])

        # Save the curve in the last ten curves list only if asked
        if save_in_last_ten:
            if len(shared.exp.last_ten_curves) < 10:
                last = shared.exp.last_ten_curves
                shared.exp.pos_in_last_ten_curves = len(last)
                shared.exp.last_ten_curves.append([xpos, ypos])
            else:
                # Remove first element
                shared.exp.last_ten_curves.pop(0)
                # Add new position
                shared.exp.last_ten_curves.append([xpos, ypos])
                shared.exp.pos_in_last_ten_curves = 9

        # xpos and ypos are given in user coordinates. Substract 1 to go to
        # array coordinates
        if xpos is not None and ypos is not None:
            shared.exp.meshgrid_click_xpos = xpos - 1
            shared.exp.meshgrid_click_ypos = ypos - 1

        # Update the curve selector
        if self.tabs.currentIndex() == 0:
            wl.widget_curve_selector_data.update_widget()
        if self.tabs.currentIndex() == 1:
            wl.widget_curve_selector_compute.update_widget()
        if self.tabs.currentIndex() == 2:
            wl.widget_curve_selector_results.update_widget()

        # Update the values in details widget if it exists
        if wl.widget_details is not None:
            wl.widget_details.update_widget()

        # Update the values in curve tilt widget if it exists
        if wl.widget_curve_mod is not None:
            wl.widget_curve_mod.update_MPL("MPL_canvas")
            wl.widget_curve_mod.IN_curve_x.changeValue(xpos)
            wl.widget_curve_mod.IN_curve_y.changeValue(ypos)

        # Update the position of the red square on the data meshgrid
        if self.tabs.currentIndex() == 0:
            wl.widget_data.meshgrid.canvas.update_blit("red_square")

        # Update the position of the red square on the compute meshgrid
        if self.tabs.currentIndex() == 1:
            wl.widget_compute.MPL_meshgrid.canvas.update_blit("red_square")

        # Update the position of the red square on the meshgrid
        if self.tabs.currentIndex() == 2:
            wl.widget_results.MPL_canvas1.canvas.update_blit("red_square")

        # Tilt check meshgrid
        if wl.widget_tilt_check is not None:
            wl.widget_tilt_check.MPL_canvas.canvas.update_blit("red_square")

        # Update the curves
        if self.tabs.currentIndex() == 0:
            wl.widget_data.curve.update_plot()
            if data.file_type in ("JPK (Single File)", "JPK (Force Map)", "JPK (QI)"):
                wl.widget_data.curve_1.update_plot()
                wl.widget_data.curve_2.update_plot()
        if self.tabs.currentIndex() == 1:
            wl.widget_compute.update_MPL("MPL_canvas")
        if self.tabs.currentIndex() == 2:
            wl.widget_results.update_MPL("MPL_canvas2")

    def file_changed(self, option=None, option2=None):
        """Method triggered after the user changed from file.

        Option2 is here for the case where you click on the
        multimeshgrid plots. We don't want to update the multimeshgrid
        widget in this case.
        """
        string = "Changing file, option = %s"
        self.logger.debug(string, str(option))

        # This list lets you chose the dataset
        # List1 is the list in the results tab
        # List2 is in the compute tab
        # Box is in the data tab

        # Check if we need to apply the display/compute options to all
        apply_to_all.apply_to_all("all", option)

        if option == "box":
            # Get the id of the selected file, has already been updated by
            # the BoxData class

            # Update the other list
            wl.widget_compute.list_exp2.setCurrentIndex(shared.exp.id_selected)
            wl.widget_results.list_exp.setCurrentIndex(shared.exp.id_selected)

        elif option == "list1":
            # Get the id of the selected file
            shared.exp.id_selected = wl.widget_results.list_exp.currentIndex()

            # Update the other list
            wl.widget_compute.list_exp2.setCurrentIndex(shared.exp.id_selected)
            wl.widget_data.update_GUI("list_widget")

        elif option == "list2":
            # Get the id of the selected file
            shared.exp.id_selected = wl.widget_compute.list_exp2.currentIndex()

            # Update the other list
            wl.widget_results.list_exp.setCurrentIndex(shared.exp.id_selected)
            wl.widget_data.update_GUI("list_widget")

        else:
            # Needed after update by load (when loading more data)
            wl.widget_results.list_exp.setCurrentIndex(shared.exp.id_selected)
            wl.widget_compute.list_exp2.setCurrentIndex(shared.exp.id_selected)

        string = "Changing file to id: %s"
        self.logger.debug(string, str(shared.exp.id_selected))

        # Get the data and update GUI elements
        data = shared.exp.list[shared.exp.id_selected]

        # Reset the list counting the already seen curves with the random
        # button
        shared.exp.parsed_random_curves = []

        # Check for the selected xpos, ypos, and set it to 0, 0
        # (=1, 1 for the user) if its out of bonds
        # The check has to be done for each dimension because there can also
        # be non square maps loaded
        new_pos = None
        if shared.exp.meshgrid_click_xpos >= data.nbr_pixels_x:
            new_pos = [1, shared.exp.meshgrid_click_ypos + 1]

        if shared.exp.meshgrid_click_ypos >= data.nbr_pixels_y:
            if new_pos is None:
                new_pos = [shared.exp.meshgrid_click_xpos + 1, 1]
            else:
                new_pos[1] = 1

        # Change the position if needed
        if new_pos is not None:
            self.change_curve(new_pos[0], new_pos[1])

        # Reset the last ten curves list, it could happen that the user wants
        # to go back to a curve out of bounds if scan sizes are different
        shared.exp.last_ten_curves = [[0, 0]]
        shared.exp.pos_in_last_ten_curves = 0

        self.update_GUI("single_or_forcevolume")
        self.update_GUI("tabs")  # For load_more

        # Change to defl-ext if there is no force computed
        data = shared.exp.current_data
        tr_ret = data.display_trace_retrace

        # If we are on "defl_ext" or "force_ext" we can stay on it; else we
        # have to check if the corresponding data has been computed
        if data.display_curve_type == "force":
            if tr_ret == 2 or tr_ret == 0:
                if data.stiffness_calculated is False:
                    data.display_curve_type = "defl_ext"

            if tr_ret == 2 or tr_ret == 1:
                if data.work_and_rupture_force1_calculated is False:
                    data.display_curve_type = "defl_ext"

        elif data.display_curve_type == "indentation":
            if data.stiffness_calculated is False:
                data.display_curve_type = "defl_ext"

        elif data.display_curve_type == "results_force_events":
            if data.events_calculated is False:
                data.display_curve_type = "defl_ext"

        # Uncheck some options if they are not aviable
        if data.stiffness_calculated is False:
            if data.display_poc:
                data.display_poc = False
            if data.display_fit_poc:
                data.display_fit_poc = False
            if data.display_segments:
                data.display_segments = False
            if data.display_fits_stiffness:
                data.display_fits_stiffness = False

        if data.work_and_rupture_force1_calculated is False:
            if data.display_joc:
                data.display_joc = False
            if data.display_fit_joc:
                data.display_fit_joc = False
            if data.display_surface:
                data.display_surface = False
            if data.display_force:
                data.display_force = False

        if data.events_calculated is False:
            if data.events_results_display_joc:
                data.events_results_display_joc = False
            if data.events_display_results_filter_dist:
                data.events_display_results_filter_dist = False
            if data.events_results_display_annotations != 0:
                data.events_results_display_annotations = 0

        # Apply to all the new options (no load_more option here,
        # we will use the current values.)
        apply_to_all.apply_to_all("all")

        # Change the meshgrid type if needed (going to a slice meshgrid
        # when there is only one slice ...)
        if data.stiffness_calculated:
            cond1 = data.meshgrid_type == "stiffness_slice"
            cond2 = data.indentation_step == 0
            if cond1 and cond2:
                data.meshgrid_type = "stiffness"
        if data.stiffness_corrected:
            cond1 = data.meshgrid_type == "stiffness_slice_corr"
            cond2 = data.indentation_step == 0
            if cond1 and cond2:
                data.meshgrid_type = "stiffness_corr"

        # Data widget
        if self.tabs.currentIndex() == 0:
            wl.widget_data.update_widget()

        # Compute Widget
        if self.tabs.currentIndex() == 1:
            wl.widget_compute.update_widget()

        # Results widget
        if self.tabs.currentIndex() == 2:
            wl.widget_results.update_widget()

        # Update widgets if they are opened :
        if wl.widget_info is not None:
            wl.widget_info.update_widget()
        if wl.widget_flatten is not None:
            wl.widget_flatten.update_widget()
        if wl.widget_slices is not None:
            wl.widget_slices.update_widget()
        if wl.widget_curve_mod is not None:
            wl.widget_curve_mod.update_widget()
        if wl.widget_details is not None:
            wl.widget_details.update_widget()
        if wl.widget_clean_up is not None:
            wl.widget_clean_up.update_widget()
        if wl.widget_meshgrid_options is not None:
            wl.widget_meshgrid_options.update_widget()
        if wl.widget_multimeshgrids is not None and option2 != "multimeshgrid":
            wl.widget_multimeshgrids.update_widget()
        if wl.widget_roi_manager is not None:
            wl.widget_roi_manager.update_widget()
        if wl.widget_indentation is not None:
            wl.widget_indentation.update_widget()
        if wl.widget_used_parameters is not None:
            wl.widget_used_parameters.update_widget()
        if wl.widget_tilt_check is not None:
            wl.widget_tilt_check.update_widget()

    def update_names(self, force_id=None):
        """Update all or one filename(s) in the GUI.

        If a force_id value is passed, only one name will be updated.
        Else, all the names will be updated.
        """
        files = []
        if force_id is not None:
            files.append(shared.exp.list[force_id])
        else:
            files = shared.exp.list

        for data in files:
            wl.widget_data.list_boxes[
                data.unique_id].label_name.setText(data.filename)
            wl.widget_results.list_exp.setItemText(
                data.unique_id, data.filename)
            wl.widget_compute.list_exp2.setItemText(
                data.unique_id, data.filename)
            if wl.widget_multimeshgrids is not None:
                wl.widget_multimeshgrids.update_widget()

    def update_GUI(self, element):
        """Method to update different GUI elements."""
        if element == "tabs":
            found = False
            for item in shared.exp.list:
                cond1 = item.stiffness_calculated
                cond2 = item.work_and_rupture_force1_calculated
                cond3 = item.events_calculated
                cond4 = item.events_calculated

                if cond1 or cond2 or cond3 or cond4:
                    found = True

            if found:
                self.tabs.setTabEnabled(2, True)
                self.tabs.setTabEnabled(3, True)
                self.tabs.setTabEnabled(4, True)
            else:
                self.tabs.setTabEnabled(2, False)
                self.tabs.setTabEnabled(3, False)
                self.tabs.setTabEnabled(4, False)

            # Update the back and forward buttons
            wl.widget_results.update_GUI("bts_back_forward")
            wl.widget_compute.update_GUI("bts_back_forward")

    def tabs_changed(self):
        """Method called whenever a tab is updated.

        Disable/enable some menu actions for the exporting of the
        figures. The menu actions are stored in self.list_actions
        """
        self.logger.debug("Going to tab: %s", str(self.tabs.currentIndex()))
        self.current_tab = self.tabs.currentIndex()

        if self.tabs.currentIndex() == 0:
            # Data tab
            self.list_actions["open_meshgrid"].setEnabled(True)
            self.list_actions["open_curve"].setEnabled(True)
            self.list_actions["open_results"].setEnabled(False)
            self.list_actions["row_color"].setEnabled(False)
            self.list_actions["export_selected"].setEnabled(False)
            self.list_actions["export_selected_r"].setEnabled(False)
            self.list_actions["rename_row"].setEnabled(False)
            self.list_actions["duplicate_row"].setEnabled(False)
            self.list_actions["remove_row"].setEnabled(False)
            self.list_actions["row_color"].setEnabled(False)

            wl.widget_data.update_widget()

        if self.tabs.currentIndex() == 1:
            # Compute tab
            self.list_actions["open_meshgrid"].setEnabled(False)
            self.list_actions["open_curve"].setEnabled(True)
            self.list_actions["open_results"].setEnabled(False)
            self.list_actions["export_selected"].setEnabled(False)
            self.list_actions["export_selected_r"].setEnabled(False)
            self.list_actions["rename_row"].setEnabled(False)
            self.list_actions["duplicate_row"].setEnabled(False)
            self.list_actions["remove_row"].setEnabled(False)
            self.list_actions["row_color"].setEnabled(False)

            wl.widget_compute.update_widget()

        elif self.tabs.currentIndex() == 2:
            # Results tab
            self.list_actions["open_meshgrid"].setEnabled(True)
            self.list_actions["open_curve"].setEnabled(True)
            self.list_actions["open_results"].setEnabled(False)
            self.list_actions["export_selected"].setEnabled(False)
            self.list_actions["export_selected_r"].setEnabled(False)
            self.list_actions["rename_row"].setEnabled(False)
            self.list_actions["duplicate_row"].setEnabled(False)
            self.list_actions["remove_row"].setEnabled(False)
            self.list_actions["row_color"].setEnabled(False)

            wl.widget_results.update_widget()

        elif self.tabs.currentIndex() == 3:
            # Plots tab
            self.list_actions["open_meshgrid"].setEnabled(False)
            self.list_actions["open_curve"].setEnabled(False)
            self.list_actions["open_results"].setEnabled(True)
            self.list_actions["export_selected"].setEnabled(True)
            self.list_actions["export_selected_r"].setEnabled(True)
            self.list_actions["rename_row"].setEnabled(True)
            self.list_actions["duplicate_row"].setEnabled(True)
            self.list_actions["remove_row"].setEnabled(True)
            self.list_actions["row_color"].setEnabled(True)
            text = "Change color of row"
            self.list_actions["row_color"].setText(text)

        elif self.tabs.currentIndex() == 4:
            # Plots (Grouped) tab
            self.list_actions["open_meshgrid"].setEnabled(False)
            self.list_actions["open_curve"].setEnabled(False)
            self.list_actions["open_results"].setEnabled(True)
            self.list_actions["export_selected"].setEnabled(True)
            self.list_actions["export_selected_r"].setEnabled(False)
            self.list_actions["rename_row"].setEnabled(False)
            self.list_actions["duplicate_row"].setEnabled(False)
            self.list_actions["remove_row"].setEnabled(False)
            self.list_actions["row_color"].setEnabled(True)
            text = "Change color of row (Groups)"
            self.list_actions["row_color"].setText(text)

        # Results update. Called from ROI manager
        if self.tabs.currentIndex() == 3 or self.tabs.currentIndex() == 4:
            if self.update_results:
                wl.widget_results_single.update_widget()
                wl.widget_results_groups.update_widget()
                self.update_results = False

        self.setFocus()

    def open_about(self):
        """Open about pyAF widget."""
        if wl.widget_about is None:
            # Create new widget
            AboutWidget(self)
            wl.widget_about.resize(500, 250)
            wl.widget_about.show()
        else:
            # Bring to front
            wl.widget_about.activateWindow()
            wl.widget_about.raise_()


def discard_shortcut():
    """ Discards or re-enables the current curve.

    Called by the CTRL-D shortcut.
    """
    discarding.change_single_curve_status()


def my_excepthook(etype, value, tback):
    """Override the python exception hook.

    This is needed because the Qt4Agg backend sends exception errors due to
    the matplotlib.close() function. This happens when you save figures without
    using the matplotlib.show() function. (Sometimes also on the Qt4MplCanvas).
    Perhaps this monkeypatching will not be necessary the day this bug is
    corrected by the PyQT/matplotlib community.

    The error messages can be displayed in DEBUG mode.

    If an error occurs, the LogSenderWidget will be opened. This will allow the
    users to fill in an issue in the bugtracker.

    The shared.error_triggered flag is set to true. When the user will try to
    save his file, he will get an error message telling him his data may have
    been corruped.
    """
    logger = logging.getLogger(__name__)

    if str(value) ==\
            "wrapped C/C++ object of type FigureCanvasQTAgg has been deleted":
        if consts.DEBUG:
            print("Info : wrapped C/C++ object" + \
                "of type FigureCanvasQTAgg has been deleted")
    elif str(value) ==\
            "wrapped C/C++ object of type Qt4MplCanvas has been deleted":
        if consts.DEBUG:
            print("Info : wrapped C/C++ object" + \
                "of type Qt4MplCanvas has been deleted")
    else:
        shared.error_triggered = True

        # Close the progressbar, because if there is an error, closing the
        # progressbar by yourself will make pyAF crash.
        if wl.widget_progressbar is not None:
            wl.widget_progressbar.close()

        # Log the error in the log file
        logger.error("Logging exception", exc_info=(etype, value, tback))

        if consts.DISABLE_ERROR_SENDER is False:
            # Display a widget to send the error directly to the bugtracker
            # While working with plugins or during debugging it can be useful
            # to set consts.DISABLE_ERROR_SENDER to True if you don't want to
            # be annoyed by this widget.

            logwidget = LogSenderWidget(None, (etype, value, tback))
            logwidget.show()

        sys.__excepthook__(etype, value, tback)


def open_entry_widget(mode, files, new_tmp_path=None):
    """Function used to launch the entry widget."""

    window_title = misc_tools.get_base_window_title()

    entry_window = EntryWidget(
        files=files, mode=mode, new_tmp_path=new_tmp_path)
    entry_window.center()
    entry_window.resize(600, 300)
    entry_window.setWindowTitle(window_title)
    entry_window.activateWindow()
    entry_window.show()
    # PyQt specific problem in OS X
    # Application window is not automatically brought to the front on
    # launch
    if sys.platform == "darwin":
        entry_window.raise_()

    # Set focus from outside (wont work inside)
    # Is important for load button callable by Enter Key press
    entry_window.setFocus()

    if not consts.UNIT_TESTING:
        entry_window.exec_()
    else:
        entry_window.close()

    # Get the results back
    files_info = entry_window.files_info
    value = entry_window.add_more
    ask_for_path = entry_window.ask_for_new_tmp_path

    return value, files_info, ask_for_path


def add_more(value, path, files):
    """Checks if the user wants to add more files to the entry widget.

    If the user adds files, the old files list and the new one are
    merged in a new one and returned. Else, nothing is done and the
    untouched files list is returned.
    """
    if value:
        # This means the user wants to add more data to the list
        # and clicked on the add more button.
        add_more_files, _ = QtWidgets.QFileDialog.getOpenFileNames(
            directory=path)
        if add_more_files:
            # Make a new list a files by merging the two lists
            new_files = []
            for old_file in files:
                new_files.append(old_file)
            for old_file in add_more_files:
                new_files.append(old_file)
            return new_files
        else:
            return files
    else:
        return files


def open_ask_for_path():
    """Ask the user to choose a folder for the tmp file.

    Open a dialog.
    """
    folder_name = QtWidgets.QFileDialog.getExistingDirectory()

    newpath = None
    if folder_name:
        newpath = os.path.normpath(str(folder_name) + "/")

    # Make sure we return an unicode path for python 2
    return str(newpath)


def main(app=None, files=None):
    """Main function, starting point of pyAF.

    To launch pyAF, you can either run this script in a terminal window
    with the following command : python main.py, or launch it with the
    bash script named pyAF, at the root folder of the application.

    In this function the QApplication is created and the entry window is
    displayed. Once the user has chosen the file(s) to load, the
    MainWindow is called.
    """
    logger = misc_tools.create_new_logger()
    logger.debug("Entering main function")

    # Redefine the exception hook, see "def my_excepthook()" above
    sys.excepthook = my_excepthook

    # Define values for Qsettings
    QtCore.QCoreApplication.setOrganizationName("CMPI")
    QtCore.QCoreApplication.setOrganizationDomain("cmpi.cnrs.fr")
    QtCore.QCoreApplication.setApplicationName("pyAF")

    # Init the GUI for the main window
    if app is None:
        app = QtWidgets.QApplication(["pyAF"])
        theme_handler("Dark")

    if not consts.AUTOLOAD:
        value = True
        new_tmp_path = None
        while value:
            value, files_info, ask_for_path = \
                open_entry_widget("first_load", files, new_tmp_path)
            if ask_for_path:
                new_tmp_path = open_ask_for_path()
                value = True

        settings = QtCore.QSettings()

        # Set the path for matplotlib at first load
        settings_path = misc_tools.get_user_path(settings, True)
        matplotlib.rcParams["savefig.directory"] = settings_path

        if not consts.AUTOLOAD:
            load_dir_path = settings_path
        else:
            load_dir_path = consts.TEST_PATH

        # Open a file dialog to select files
        if not consts.UNIT_TESTING:
            path = load_dir_path
            files, _ = QtWidgets.QFileDialog.getOpenFileNames(directory=path)

    if not consts.UNIT_TESTING:
        if not consts.AUTOLOAD and not files:
            # If the user clicks on "cancel" we quit the software
            sys.exit()

        elif not consts.AUTOLOAD and files:
            value = True
            # Loop as long as we are not ready. Allows to add more files
            # to the entry widget.
            new_tmp_path = None
            while value:
                value, files_info, ask_for_path = \
                    open_entry_widget("second_load", files, new_tmp_path)
                if ask_for_path:
                    new_tmp_path = open_ask_for_path()
                    value = True
                else:
                    files = add_more(value, path, files)

        else:
            # In autoload mode we load the file defined in consts.py
            files_info = CheckFiles(
                [consts.TEST_PATH + consts.LOAD_FILE]).infos

    else:
        # For unit testing, get the files variable which was passed
        # through and get the files_info
        files_info = CheckFiles(files).infos

    # Create main window once the loading window is closed
    main_window = MainWindow(files_info)
    main_window.move(0, 0)
    main_window.show()
    if sys.platform == "darwin":
        main_window.raise_()

    # Update the title if needed
    window_title = misc_tools.get_base_window_title()
    label = shared.exp.main_title_label
    if label:
        main_window.setWindowTitle(window_title + " - " + label)
    else:
        main_window.setWindowTitle(window_title)

    if consts.AUTO3D:
        # When autoloading with 3D, the 3D widget stays behind the mainwindow
        # whatever is done. To be sure to have direct access to the 3D widget
        # while debuggin, just move the main_window away.
        main_window.move(1200, 0)

    if consts.UNIT_TESTING is False:
        app.exec_()
