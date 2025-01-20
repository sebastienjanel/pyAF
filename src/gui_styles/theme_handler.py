"""
Python application for fitting Hertz model to different force curves
"""
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QFile, QTextStream
from ..tools.misc_tools import get_resource_path
import matplotlib.pyplot as plt
from . import breeze_resources  # Do not delete, important for making the style themes work.


def theme_handler(theme):
    # Get paths of the stylesheets and custom plot styles.

    # Apply theme selected by the user.
    if theme == "Light":
        switch_theme(get_resource_path("light.qss"))
        plt.style.use("default")
    elif theme == "Dark":
        switch_theme(get_resource_path("dark.qss"))
        plt.style.use(get_resource_path("dark_plot.mplstyle"))

def switch_theme(path):
    # get the QApplication instance
    app = QApplication.instance()
    if app is None:
        raise RuntimeError("No Qt Application running.")

    file = QFile(path)
    file.open(QFile.ReadOnly | QFile.Text)
    stream = QTextStream(file)
    app.setStyleSheet(stream.readAll())
