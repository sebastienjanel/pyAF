"""
Python application for fitting Hertz model to different force curves
"""
from qtpy.QtWidgets import QApplication
import os
import sys
from qtpy.QtCore import QDir, QEvent, QFile, QObject, QTextStream
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
             "pause": "#FF0000", "modulation": "#00FF80"},
    "Light": {"approach": "#007C99", "retraction": "#C25E00",
              "pause": "#B00000", "modulation": "#00884A"}}


def curve_color(name):
    """Colour of a curve (approach, retraction, pause or modulation) for
    the current theme."""
    return CURVE_COLORS.get(current_theme, CURVE_COLORS["Dark"])[name]


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
