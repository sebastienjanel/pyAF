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


"""
Python application for fitting Hertz model to different force curves
"""
from qtpy.QtWidgets import QApplication
import os
import sys
from qtpy.QtCore import QDir, QEvent, QFile, QObject, QSettings, \
    QTextStream
from qtpy.QtWidgets import QGroupBox
from ..tools.misc_tools import get_resource_path
import matplotlib.pyplot as plt
from . import breeze_resources  # Do not delete, important for making the style themes work.


# Theme applied by theme_handler ("Dark" or "Light")
current_theme = "Dark"

# Colours of the curves: the light theme needs darker colours on the white
# plot background
CURVE_COLORS = {
    "Dark": {"approach": "#00FFFF", "retraction": "#FFA500",
             "pause": "#FF0000", "modulation": "#00FF80",
             "breakthrough": "#FF6B6B"},
    "Light": {"approach": "#007C99", "retraction": "#C25E00",
              "pause": "#B00000", "modulation": "#00884A",
              "breakthrough": "#D00000"}}


def curve_color(name):
    """Colour of a curve (approach, retraction, pause or modulation) or of
    the breakthroughs for the current theme."""
    return CURVE_COLORS.get(current_theme, CURVE_COLORS["Dark"])[name]


def saved_theme():
    """Theme chosen in the preferences (kept between the launches of pyAF),
    dark at the first launch."""
    theme = QSettings().value("theme", "Dark")
    return theme if theme in ("Light", "Dark") else "Dark"


def set_theme(theme):
    """Apply a theme and keep it for the next launches."""
    theme_handler(theme)
    QSettings().setValue("theme", theme)


def theme_handler(theme):
    global current_theme
    current_theme = theme
    # Get paths of the stylesheets and custom plot styles.

    # Apply theme selected by the user.
    if theme == "Light":
        switch_theme(get_resource_path("light.qss"))
        plt.style.use("default")
    elif theme == "Dark":
        switch_theme(get_resource_path("dark.qss"))
        plt.style.use(get_resource_path("dark_plot.mplstyle"))

class GroupBoxTitles(QObject):
    """Displays the titles of the group boxes in uppercase (section
    headers, distinct from the labels; the style sheets can not change the
    case of a text)."""

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Polish and isinstance(obj, QGroupBox):
            title = obj.title()
            if title and title != title.upper():
                obj.setTitle(title.upper())
        return False


def install_group_box_titles(app):
    """Install the GroupBoxTitles filter on the application (once)."""
    if getattr(app, "_pyaf_group_box_titles", None) is None:
        app._pyaf_group_box_titles = GroupBoxTitles(app)
        app.installEventFilter(app._pyaf_group_box_titles)


def get_icons_path():
    """Folder of the icons used in the stylesheets (url(pyaf:name.svg))."""
    if hasattr(sys, "_MEIPASS"):
        # Bundled application (see main_macOS.spec)
        return os.path.join(sys._MEIPASS, "gui_styles", "icons")
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "icons")


def switch_theme(path):
    # get the QApplication instance
    app = QApplication.instance()
    if app is None:
        raise RuntimeError("No Qt Application running.")

    QDir.setSearchPaths("pyaf", [get_icons_path()])
    install_group_box_titles(app)

    file = QFile(path)
    file.open(QFile.ReadOnly | QFile.Text)
    stream = QTextStream(file)
    app.setStyleSheet(stream.readAll())
