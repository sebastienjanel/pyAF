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

"""Allows the user to launch home made plugins."""

import os
import logging
from PyQt5 import QtWidgets
from ..tools import misc_tools
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFButton
from ..tools.PYAFWidget import PYAFWidget
from ..tools.plugin_loader import Plugins


class PluginsWidget(PYAFWidget):
    """Widget with the list of plugins in the plugin folder and load button.

    Lets the user load and execute a plugin at runtime. The user can change the
    location of the plugin folder. If PYAF was installed from source some
    demo plugins are provided in the default plugin folder.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_plugins")

        self.parent = parent

        # Load logger
        self.logger = logging.getLogger()
        self.logger.debug("Opening Plugins Widget")

        self.resize(400, 100)

        # Get the path
        default_folder = misc_tools.get_app_path() + "/../../plugins"
        if os.path.isdir(default_folder):
            plugin_folder_path = default_folder
        else:
            plugin_folder_path = None

        self.plugin_loader = Plugins(plugin_folder_path)
        if self.plugin_loader.plugin_folder is not None:
            self.plugin_loader.get_plugins()

        self.VL = QtWidgets.QVBoxLayout()

        self.list = PYAFComboBox(self, None)
        self.BT_run = PYAFButton(self, "run", "Run", 120)
        self.BT_refresh_list = PYAFButton(self, "refresh", "Refresh list", 120)

        self.BT_change_folder = PYAFButton(
            self, "change", "Choose folder", 160)
        self.label_folder = QtWidgets.QLabel()

        self.VL.addWidget(self.list)
        self.VL.addWidget(self.BT_refresh_list)
        self.VL.addWidget(self.BT_run)
        self.VL.addWidget(self.BT_change_folder)
        self.VL.addWidget(self.label_folder)
        self.setLayout(self.VL)

        self.update_widget()

    def update_widget(self):
        """Updates the content of the widget."""
        self.logger.debug("Updating the widget")

        # Get new list of plugins
        if self.plugin_loader.plugin_folder is not None:
            self.plugin_loader.get_plugins()

        # The the label for the folder path
        self.label_folder.setText(str(self.plugin_loader.plugin_folder))

        # Remove all items from the plugins list
        self.list.clear()

        # Add new plugins
        for item in self.plugin_loader.plugins:
            self.list.addItem(item["name"])

        if self.plugin_loader.plugin_folder is None:
            self.list.addItem("No plugin folder")
            self.list.setEnabled(False)
            self.BT_run.setEnabled(False)
            self.BT_refresh_list.setEnabled(False)
        else:
            self.list.setEnabled(True)
            self.BT_run.setEnabled(True)
            self.BT_refresh_list.setEnabled(True)

        # If no plugin was found
        if len(self.plugin_loader.plugins) == 0 and \
                self.plugin_loader.plugin_folder is not None:
            self.list.addItem("No plugins found")
            self.BT_run.setEnabled(False)

    def button_clicked(self, name):
        """Called whenever a button is clicked."""
        if name == "run":
            # Check if no new plugin was added
            self.plugin_loader.get_plugins()
            if self.list.count() != len(self.plugin_loader.plugins):
                text = "Plugins were added or deleted. Please reselect the " \
                    + "plugin you want to run."

                msg_box = QtWidgets.QMessageBox()
                msg_box.setText("Error")
                msg_box.setInformativeText(text)
                msg_box.setIcon(QtWidgets.QMessageBox.Critical)
                msg_box.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
                msg_box.exec_()

                self.update_widget()

                return False

            # Get the id from the list
            listid = self.list.currentIndex()

            # Load the plugin
            plugin = self.plugin_loader.load_plugin(listid)
            plugin.run(self.parent)

        elif name == "change":
            foldername = QtWidgets.QFileDialog.getExistingDirectory()

            if not foldername:
                pass
            else:
                newpath = str(foldername) + "/"
                self.logger.debug("Changing plugins folder to %s", newpath)
                self.plugin_loader.plugin_folder = newpath
                self.plugin_loader.get_plugins()
                self.update_widget()

        elif name == "refresh":
            self.update_widget()
