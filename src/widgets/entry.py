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

"""First widget displayed by the app."""

import os
import psutil
import sys
import logging
from PyQt5 import QtCore, QtWidgets
from .. import consts
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFComboBox
from ..tools import misc_tools
from ..tools import gui_tools
import tempfile
from PIL import Image
import numpy
from ..load_and_save.pre_check import CheckFiles
from ..tools.misc_tools import get_resource_path


class EntryWidget(QtWidgets.QDialog):
    """Widget displayed at launch.

    The entry widget is there to display a picture introducing the app (with
    name, authors, logos ...
    There is a checkbox to disable the histogram refreshing at load
    (loading is faster).
    The button will start the loading procedure of AFM files.

    There are two layouts, depending if the widgets is opened the first time
    or the second time.
    """

    def __init__(self, files=None, mode=None, new_tmp_path=None):
        super().__init__()

        if mode == "load_more":
            load_more = True
            mode = "second_load"
        else:
            load_more = False

        self.logger = logging.getLogger()
        self.logger.debug("Opening EntryWidget, load_more = %s", load_more)

        self.mode = mode

        self.smallfont = QtWidgets.QApplication.font()  # Default font
        self.smallfont.setPointSize(9)

        # Get the settings
        self.settings = QtCore.QSettings()

        # Define some variables
        self.label_path = None
        self.BT_change_path = None
        self.triangle_right = None
        self.label_mode = None
        self.label_triangle = None
        self.tableWidget = None
        self.BT_reset_path = None
        self.triangle_down = None
        self.GL_advanced = None
        self.checkbox = None
        self.BT_load = None

        # Value for closeEvent
        self.ask_for_new_tmp_path = False
        self.exit_from_button = False
        self.load_more = load_more
        self.add_more = False

        # The second time this widget is displayed, get the list of files
        self.files = files
        if self.mode == "second_load" or load_more:
            self.files_info = CheckFiles(self.files).infos
        else:
            self.files_info = None

        if consts.UNIT_TESTING and self.mode == "first_load":
            # In case of unit tests, fill the list with the passed files
            self.files_info = CheckFiles(self.files).infos

        # Is the advanced mode part displayed or not ?
        self.display_advanced = False

        temp_path = os.path.normpath(tempfile.gettempdir() + "/")
        self.tempPath = self.settings.value("tempPath", temp_path)

        # Check if path exists (could have been deleted or moved)
        # If it does not exist, fall back to tmp folder
        if os.path.isdir(str(self.tempPath).encode("utf8")) is False:
            self.tempPath = temp_path
            self.settings.setValue("tempPath", temp_path)

        # Get the path of the entry image
        path = os.path.normpath(misc_tools.get_app_path() + "/images")
        self.im_path = os.path.normpath(path)

        # List of images
        self.image_list = []
        self.load_images()

        self.VL_entry = QtWidgets.QVBoxLayout()

        self.HL_image = QtWidgets.QHBoxLayout()
        self.HL_image.setSpacing(0)
        self.labelfixed = QtWidgets.QLabel()
        self.labelchanging = QtWidgets.QLabel()
        self.HL_image.addWidget(self.labelfixed)
        self.HL_image.addWidget(self.labelchanging)
        self.labelchanging.setStyleSheet("background-color: white")

        # Set fixed image
        #path = os.path.normpath(self.im_path + "/fix.png")
        im = numpy.array(Image.open(get_resource_path("fix.png")).convert("RGBA"))
        self.labelfixed.setPixmap(misc_tools.get_pixmap(im))

        self.VL_entry.addLayout(self.HL_image)

        # Advanced mode
        self.label_triangle = QtWidgets.QLabel()
        #path = os.path.normpath(self.im_path + "/triangle_right.png")
        self.triangle_right = numpy.array(Image.open(get_resource_path("triangle_right.png")).convert("RGBA"))
        #path = os.path.normpath(self.im_path + "/triangle_down.png")
        self.triangle_down = numpy.array(Image.open(get_resource_path("triangle_down.png")).convert("RGBA"))
        pixmap = misc_tools.get_pixmap(self.triangle_right)
        self.label_triangle.setPixmap(pixmap)

        self.label_mode = QtWidgets.QLabel()
        self.label_mode.setText("Advanced")
        self.label_triangle.mousePressEvent = lambda event: self.change_mode()
        self.label_mode.mousePressEvent = lambda event: self.change_mode()

        # Define checkbox for advanced mode
        label = "Refresh histograms at load (slower)"
        self.checkbox = PYAFCheckBox(self, "refresh", label)
        if self.settings.contains("refreshHistsAtLoad"):
            # Get the value
            if self.settings.value("refreshHistsAtLoad") == "true":
                self.refresh = True
            elif self.settings.value("refreshHistsAtLoad") == "false":
                self.refresh = False
            else:
                self.refresh = self.settings.value("refreshHistsAtLoad")

        else:
            # Default value
            self.refresh = True
        self.checkbox.setChecked(self.refresh)

        # Path and buttons for advanced modes
        self.label_path = QtWidgets.QLabel()
        self.label_path.setFont(self.smallfont)
        self.label_path.setText(self.tempPath)
        self.BT_reset_path = PYAFButton(self, "reset_path", "Reset")
        self.BT_reset_path.setFont(self.smallfont)
        self.BT_change_path = PYAFButton(self, "change_path", "Change path")
        self.BT_change_path.setFont(self.smallfont)

        # Define advanced layout
        self.GL_advanced = QtWidgets.QGridLayout()

        if self.mode == "first_load":
            # The first time the widget is displayed we display the advanced
            # mode triangle and options
            self.init_first_UI()
        elif self.mode == "second_load" or load_more:
            # The second time we display the list with the files
            self.init_second_UI()

        self.setLayout(self.VL_entry)
        self.layout().setSizeConstraint(QtWidgets.QLayout.SetFixedSize)

        # Save new tmp path
        if new_tmp_path is not None:
            self.settings.setValue("tempPath", new_tmp_path)
            self.label_path.setText(new_tmp_path)
            self.logger.debug("Changed path to = %s", new_tmp_path)

        # Create a QTimer
        self.timer = QtCore.QTimer()
        self.alpha = 255
        self.pos = 0
        self.sign = -1
        # Connect it to fade
        self.timer.timeout.connect(self.fade)
        # Call every 0.1 seconds
        self.timer.start(50)

        if consts.UNIT_TESTING:
            if self.mode == "first_load":
                self.button_clicked("first_load")
            elif self.mode == "second_load":
                self.button_clicked("second_load")

    def init_first_UI(self):
        """First layout of the entry widget.

        Is called at launch.
        """
        self.logger.debug("EntryWidget init first UI")

        # Define loading button and mode
        HL_mode = QtWidgets.QHBoxLayout()

        self.BT_load = PYAFButton(self, "first_load", "Load")
        self.BT_load.setDefault(True)

        HL_mode.addWidget(self.label_triangle)
        HL_mode.addWidget(self.label_mode)
        HL_mode.addStretch(1)
        HL_mode.addWidget(self.BT_load)

        # Add to layout
        self.VL_entry.addLayout(HL_mode)
        self.VL_entry.addLayout(self.GL_advanced)

    def init_second_UI(self):
        """Second layer of the entry widget."""
        self.logger.debug("EntryWidget init second UI")

        self.tableWidget = QtWidgets.QTableWidget()
        behavior = QtWidgets.QAbstractItemView.SelectRows
        self.tableWidget.setSelectionBehavior(behavior)
        self.tableWidget.horizontalHeader().setHighlightSections(False)
        policy = QtCore.Qt.ScrollBarAlwaysOff
        self.tableWidget.setHorizontalScrollBarPolicy(policy)
        self.tableWidget.verticalHeader().hide()

        header_labels = ["Load", "Name", "Type", "Error"]
        self.tableWidget.setColumnCount(4)
        self.tableWidget.setHorizontalHeaderLabels(header_labels)

        for i in range(len(self.files)):
            self.add_line_in_tablewidget(i)

        # Set the width of the columns
        self.tableWidget.setColumnWidth(0, 40)
        self.tableWidget.setColumnWidth(1, 240)
        self.tableWidget.setColumnWidth(2, 100)
        self.tableWidget.setColumnWidth(3, 203)

        # Define check all checkbox
        checkbox = QtWidgets.QCheckBox("Select all")
        checkbox.setChecked(True)
        checkbox.clicked.connect(lambda: self.checkbox_clicked("check_all"))

        HL = QtWidgets.QHBoxLayout()
        # Add the button to load the selected data
        self.BT_load = PYAFButton(self, "second_load", "Load", 150)
        self.BT_load.setDefault(True)

        # Add the button to load the selected data
        button_add_more = PYAFButton(self, "add_more", "Add more", 150)

        HL_mode = QtWidgets.QHBoxLayout()
        HL_mode.addWidget(self.label_triangle)
        HL_mode.addWidget(self.label_mode)
        HL_mode.addStretch(1)

        if self.load_more is False:
            HL.addLayout(HL_mode)
        HL.addStretch(1)
        HL.addWidget(button_add_more)
        HL.addWidget(self.BT_load)

        self.VL_entry.addWidget(self.tableWidget)
        self.VL_entry.addWidget(checkbox)
        self.VL_entry.addLayout(HL)
        self.VL_entry.addLayout(self.GL_advanced)

    def change_mode(self):
        """Method to display or to hide the advanced mode part."""
        value = self.display_advanced
        self.logger.debug("Change mode, display advanced = %s", value)

        if self.display_advanced:
            self.display_advanced = False
            self.checkbox.setParent(None)
            self.BT_reset_path.setParent(None)
            self.BT_change_path.setParent(None)
            self.label_path.setParent(None)
            self.label_mode.setText("Advanced")
            pixmap = misc_tools.get_pixmap(self.triangle_right)
            self.label_triangle.setPixmap(pixmap)
        else:
            self.display_advanced = True
            self.GL_advanced.addWidget(self.checkbox, 0, 0, 1, 0)
            self.GL_advanced.addWidget(self.BT_reset_path, 1, 0)
            self.GL_advanced.addWidget(self.BT_change_path, 1, 1)
            self.GL_advanced.addWidget(self.label_path, 1, 2)
            self.label_mode.setText("Hide")
            pixmap = misc_tools.get_pixmap(self.triangle_down)
            self.label_triangle.setPixmap(pixmap)

        self.layout().setSizeConstraint(QtWidgets.QLayout.SetFixedSize)

    def add_line_in_tablewidget(self, i):
        """Add a line on the tablewidget."""
        info = self.files_info[i]
        error = info["error"]

        if error == "Not a valid AFM file" or error == "Corrupted file" \
                or error == "Permission denied" or error == "IOError":
            color = ("color : red")

        elif error == "Nanoscope version error" or error == "Incomplete file"\
                or error == "Decode error":
            color = ("color : orange")

        else:
            color = ("color : green")

        self.tableWidget.insertRow(i)

        # Default row height is 30, but it is too big
        self.tableWidget.setRowHeight(i, 30)

        # Checkbox (hide or display all)
        cb = gui_tools.CenteredCellCheckbox(self, "display", str(i))
        self.tableWidget.setCellWidget(i, 0, cb)
        self.tableWidget.cellWidget(i, 0).checkbox.setChecked(info["checked"])
        self.tableWidget.cellWidget(i, 0).checkbox.setEnabled(info["enabled"])

        # Name
        label = QtWidgets.QLabel()
        label.setFont(self.smallfont)
        label.setText(info["filename"])
        label.setStyleSheet(color)
        lb = gui_tools.CenteredCellwidget(label)
        self.tableWidget.setCellWidget(i, 1, lb)

        # Type
        label = QtWidgets.QLabel()
        label.setFont(self.smallfont)
        label.setText(info["file_type"])
        label.setStyleSheet(color)
        lb = gui_tools.CenteredCellwidget(label)
        self.tableWidget.setCellWidget(i, 2, lb)

        # Error message
        if info["error"] == "Nanoscope version error":
            content = PYAFComboBox(self, i)
            content.addItem("Chose a version")
            content.addItem("5.x")
            content.addItem("6.x")
            content.addItem("7.x")
            content.addItem("8.x")
            content.setFont(self.smallfont)
            content.setStyleSheet(color)
            lb = gui_tools.CenteredCellComboBox(content)
            self.tableWidget.setCellWidget(i, 3, lb)
        else:
            content = QtWidgets.QLabel()
            content.setFont(self.smallfont)
            content.setText(info["error"])
            content.setStyleSheet(color)
            lb = gui_tools.CenteredCellwidget(content)
            self.tableWidget.setCellWidget(i, 3, lb)

        # Define the size of the table to hide the white space at the bottom
        if self.tableWidget.verticalHeader().length() >= 200:
            table_height = 200
        else:
            table_height = self.tableWidget.verticalHeader().length()
        self.tableWidget.setFixedSize(600, table_height + 30)

    def fade(self):
        """Method which changes the alpha value of an image.

        The image is stored in a numpy array, the slicing allows a fast change
        of the alpha value. Depending on self.sign, the alpha value will
        increase or decrease, resulting in a fade in or a fade out.
        """
        outpos = self.pos

        if self.pos + 1 == len(self.image_list):
            outpos = self.pos

        fade_out = self.image_list[outpos].astype('float')

        self.alpha += 3 * self.sign

        fade_out[:, :, 3] = self.alpha
        self.labelchanging.setPixmap(misc_tools.get_pixmap(fade_out))

        QtWidgets.QApplication.processEvents()

        if self.alpha == 0 or self.alpha == 255:
            if self.pos + 1 == len(self.image_list):
                self.pos = 0
            else:
                self.pos += 1
            # Reset values
            if self.sign == 1:
                self.sign = -1
            else:
                self.sign = 1

    def load_images(self):
        """Load the images."""
        path = self.im_path

        # Display image in a label
        im1 = Image.open(get_resource_path("01.png")).convert("RGBA")
        im2 = Image.open(get_resource_path("02.png")).convert("RGBA")
        im3 = Image.open(get_resource_path("03.png")).convert("RGBA")
        im4 = Image.open(get_resource_path("04.png")).convert("RGBA")
        im5 = Image.open(get_resource_path("05.png")).convert("RGBA")
        im6 = Image.open(get_resource_path("06.png")).convert("RGBA")
        im7 = Image.open(get_resource_path("07.png")).convert("RGBA")
        self.image_list.append(numpy.array(im1))
        self.image_list.append(numpy.array(im2))
        self.image_list.append(numpy.array(im2))
        self.image_list.append(numpy.array(im3))
        self.image_list.append(numpy.array(im3))
        self.image_list.append(numpy.array(im4))
        self.image_list.append(numpy.array(im4))
        self.image_list.append(numpy.array(im5))
        self.image_list.append(numpy.array(im5))
        self.image_list.append(numpy.array(im6))
        self.image_list.append(numpy.array(im6))
        self.image_list.append(numpy.array(im7))
        self.image_list.append(numpy.array(im7))
        self.image_list.append(numpy.array(im1))

    def center(self):
        """Called by main.py to center the Widget in the middle of the screen."""
        geometry = self.frameGeometry()
        wcenter = QtWidgets.QDesktopWidget().availableGeometry().center()
        geometry.moveCenter(wcenter)
        self.move(geometry.topLeft())

    def button_clicked(self, name):
        """Method called when the a button is clicked."""
        self.logger.debug("Button clicked (%s)", name)

        if name == "first_load":
            self.exit_from_button = True
            self.close()

        elif name == "change_path":
            # For this action, we need to close the widget and call the
            # file opening window from pyaf.py
            self.ask_for_new_tmp_path = True
            self.exit_from_button = True
            self.close()

        elif name == "reset_path":
            temp_path = os.path.normpath(tempfile.gettempdir() + "/")
            self.settings.setValue("tempPath", temp_path)
            self.label_path.setText(temp_path)

        elif name == "second_load":
            # Check the number of files to load
            found = False
            for i in range(len(self.files_info)):
                if self.files_info[i]["checked"]:
                    found = True

            if not found:
                # Create a message box if no file was found
                text = "You should load at least one file"
                msg = QtWidgets.QMessageBox()
                msg.setText("Error (No file selected)")
                msg.setInformativeText(text)
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
                msg.exec_()

                # Set focus back
                self.setFocus()

            else:
                # Check for the folder's size
                temp_path = os.path.normpath(tempfile.gettempdir() + "/")
                path = self.settings.value("tempPath", temp_path)

                # get_free_space_mb uses os.statvfs, which is not compatible
                # with unicode. In this case try to encode the path to make
                # it work. Solution: move to python3
                path = str(path).encode("utf8")

                # Deprecated use of get_free_soace_mb in favour of psutil.disk_usage
                # available = misc_tools.get_free_space_mb(path)

                # Get free disk space in path in mb
                available = psutil.disk_usage(path).free / 1024.0 / 1024.0

                total_size = 0
                for file_info in self.files_info:
                    val = os.path.getsize(str(file_info["path"]))
                    # Will hold around 4*val on disk + 4 for the results.
                    # total_size += 8 * val / 1024.0 / 1024.0
                    total_size += 2 * val / 1024.0 / 1024.0

                if available < total_size:
                    end = "You can read more about this problem in the " + \
                        "user guide."

                    # There is not enough space
                    if self.load_more is False:
                        text = "There is not enough space in the " + \
                            "temporary folder. Please select another one. " + \
                            end
                    else:
                        text = "There is not enough space in the " + \
                            "temporary folder. Save your data, close the " + \
                            "software and reopen your .pyaf file. " + \
                            "Make sure you select a new temporary " + \
                            "folder. " + end

                    msg = QtWidgets.QMessageBox()
                    msg.setText("Error (Not enough space on disk)")
                    msg.setInformativeText(text)
                    msg.setIcon(QtWidgets.QMessageBox.Critical)
                    msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
                    msg.exec_()

                    # Set focus back
                    self.setFocus()

                else:
                    # Save current path (first file loaded ...)
                    settings = QtCore.QSettings()
                    path = os.path.dirname(str(self.files_info[0]["path"]))
                    misc_tools.set_user_path(settings, path)

                    self.exit_from_button = True

                    # Close widget, will load the file
                    self.close()

        elif name == "add_more":
            self.add_more = True
            self.exit_from_button = True

            # Close widget, the code in pyAF.py will ask the user to add
            # new files, and reopen this widget.
            self.close()

    def checkbox_clicked(self, name, row=None):
        """Method called upon click on a checkbox."""
        if name == "display":
            row = int(row)
            val = self.tableWidget.cellWidget(row, 0).checkbox.isChecked()
            self.files_info[row]["checked"] = val

        elif name == "check_all":
            val = self.checkbox.isChecked()

            for i in range(len(self.files_info)):
                if self.files_info[i]["enabled"]:
                    self.tableWidget.cellWidget(i, 0).checkbox.setChecked(val)
                    self.files_info[i]["checked"] = val

        elif name == "refresh":
            val = self.checkbox.isChecked()
            self.settings.setValue("refreshHistsAtLoad", val)

    def list_updated(self, i):
        """Method called when a list is clicked."""
        value = self.tableWidget.cellWidget(i, 3).list.currentText()
        self.logger.debug("List updated (%s)", value)

        # Count how many files there are with a version error
        count = 0
        for j in range(len(self.files_info)):
            if self.files_info[j]["error"] == "Nanoscope version error":
                count += 1

        if value != "Chose a version" and count <= 1:
            self.files_info[i]["version"] = value
            self.files_info[i]["enabled"] = True
            self.files_info[i]["checked"] = True
            self.tableWidget.cellWidget(i, 0).checkbox.setChecked(True)
            self.tableWidget.cellWidget(i, 0).checkbox.setEnabled(True)

        elif value != "Chose a version" and count > 1:
            text = "Do you want to apply this version number to all the \
                    files with missing version numbers ?"

            # Create a message box
            msg_corr = QtWidgets.QMessageBox()
            msg_corr.setText("Apply on all")
            msg_corr.setInformativeText(text)
            msg_corr.setIcon(QtWidgets.QMessageBox.Critical)
            msg_corr.addButton(QtWidgets.QMessageBox.Yes)
            msg_corr.addButton(QtWidgets.QMessageBox.No)
            msg_corr.addButton(QtWidgets.QMessageBox.Abort)
            ret = msg_corr.exec_()

            if ret == QtWidgets.QMessageBox.Yes:
                # Apply to all

                if value == "5.x":
                    index = 1
                if value == "6.x":
                    index = 2
                if value == "7.x":
                    index = 3
                if value == "8.x":
                    index = 4

                err = "Nanoscope version error"
                tw = self.tableWidget

                for j in range(len(self.files_info)):
                    if self.files_info[j]["error"] == err:
                        self.files_info[j]["version"] = value
                        self.files_info[j]["enabled"] = True
                        self.files_info[j]["checked"] = True
                        tw.cellWidget(j, 3).list.setCurrentIndex(index)
                        tw.cellWidget(j, 0).checkbox.setChecked(True)
                        tw.cellWidget(j, 0).checkbox.setEnabled(True)

            elif ret == QtWidgets.QMessageBox.No:
                # Apply only to the selected file
                self.files_info[i]["version"] = value
                self.files_info[i]["enabled"] = True
                self.files_info[i]["checked"] = True
                self.tableWidget.cellWidget(i, 0).checkbox.setChecked(True)
                self.tableWidget.cellWidget(i, 0).checkbox.setEnabled(True)

            elif ret == QtWidgets.QMessageBox.Abort:
                # Abort and reset combobox
                self.tableWidget.cellWidget(i, 3).list.setCurrentIndex(0)

    def closeEvent(self, event):
        """Reimplement the closeEvent.

        When closing the entry widget we don't want the app to be closed,
        we want to continue to go through the loading procedure.
        """
        # Stop the timer
        self.timer.stop()

        if self.exit_from_button:
            self.logger.debug("Closing and continue loading")
            # If a load button was called, close normally
            event.accept()

        else:
            # Else quit completely the app (command Q for example)
            # Do this only if we are not in load_more mode, else set files_info
            # to None, so nothing will be loaded
            if not self.load_more:
                sys.exit()
            else:
                self.logger.debug("Closing in load_more mode, doing nothing")
                self.files_info = None
                event.accept()
