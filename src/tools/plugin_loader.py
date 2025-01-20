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

"""Module used to make a list of plugins and load the wanted one."""

import os
import importlib.util

class Plugins:
    """Gets the list of plugins and allows to load a plugin at runtime.

    Inspired from :
    http://lkubuntu.wordpress.com/2012/10/02/writing-a-python-plugin-api/
    """

    def __init__(self, path):
        self.plugin_folder = path
        self.plugins = []

    def get_plugins(self):
        """Gets the list of plugins inside the plugins folder.

        A plugin is found if there is a plugin.py and an __init__.py file
        inside a folder, which should have as name the plugin's name.
        """
        self.plugins = []
        possible_plugins = os.listdir(self.plugin_folder)

        for i in possible_plugins:
            location = os.path.join(self.plugin_folder, i)
            is_dir = os.path.isdir(location)
            if not is_dir or ("plugin.py" not in os.listdir(location)):
                continue
            plugin_path = os.path.join(location, "plugin.py")
            spec = importlib.util.spec_from_file_location("plugin", plugin_path)

            if spec:
                self.plugins.append({"name": i, "spec": spec})

    def load_plugin(self, plugin_id):
        """Loads (runs) a plugin at runtime.

        The wanted plugin (given by the id) is loaded at runtime with the help
        of the importlib library. The loaded module is returned so that it can be
        used wherever you want.
        """
        # Recreate a fresh list of plugins (will "open" new plugins)
        self.get_plugins()

        # Get the current plugin from the list
        plugin_spec = self.plugins[plugin_id]["spec"]

        # Load the plugin module
        module = importlib.util.module_from_spec(plugin_spec)
        plugin_spec.loader.exec_module(module)

        return module