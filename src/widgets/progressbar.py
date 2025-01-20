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

"""A progressbar."""

from PyQt5 import QtCore, QtWidgets
from ..tools.PYAFWidget import PYAFWidget


class Progressbar(PYAFWidget):
    """Creates a progressbar which stays on top."""

    def __init__(self, title="", label=""):
        super().__init__(None, "widget_progressbar")

        self.label = QtWidgets.QLabel(label)
        # Block the usage of the main app and let the progressbar stay on top.
        self.setWindowModality(QtCore.Qt.ApplicationModal)

        self.setWindowTitle(title)
        self.widget_VL = QtWidgets.QVBoxLayout()
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.progressbar = QtWidgets.QProgressBar()
        self.widget_VL.addWidget(self.label)
        self.widget_VL.addWidget(self.progressbar)
        self.widget_VL.addStretch(1)
        self.setLayout(self.widget_VL)
        self.counter = 0
        self.progressbar.setRange(0, 0)
        self.show()

    def update(self):
        """Update by one step."""
        self.progressbar.setValue(self.counter)
        self.counter += 1
        QtCore.QCoreApplication.processEvents()

    def set_label(self, text):
        """Update the label's text."""
        self.label.setText(text)
        QtCore.QCoreApplication.processEvents()

    def set_range(self, range1, range2):
        """Set the minimum and the maximum."""
        self.progressbar.setRange(range1, range2)

    def reset(self):
        """Reset the progressbar to 0."""
        self.counter = 0
        self.progressbar.reset()
        QtCore.QCoreApplication.processEvents()
