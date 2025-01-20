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

"""Window which opens when an error occurs."""

import logging
import webbrowser
import traceback
from PyQt5 import QtWidgets
from ..tools import misc_tools
from ..tools.gui_tools import PYAFButton
from ..tools.PYAFWidget import PYAFWidget


class LogSenderWidget(PYAFWidget):
    """Widget which displays the error."""

    def __init__(self, parent, tback):
        super().__init__(parent, "widget_log_sender")

        self.logger = logging.getLogger()

        self.VL = QtWidgets.QVBoxLayout()

        self.label = QtWidgets.QLabel()
        text = ("An error has occured. To get this bug solved and to further "
                "improve PYAF, you can create a new issue in our bugtracker."
                "<br><br>Feel free to modify the text if there are critical "
                "things for your research that should not be sent to us (for "
                "example the names of the files you opened).<br><br>You can "
                "also add a description to the text.<br><br>The text will be "
                "copied to your clipbord, so just use CTRL-V to paste it on "
                "the website.")
        self.label.setText(text)

        self.error = QtWidgets.QTextEdit()

        path_to_file = misc_tools.find_logger_filename(self.logger)
        if path_to_file is not None:
            thefile = open(path_to_file, "r")
            text = thefile.read()
            thefile.close()
        else:
            # The logging was disabled because the log folder is not writable
            # Just show the error message
            text = "".join(traceback.format_exception(*tback))

        self.error.setText(text)

        self.HL_buttons = QtWidgets.QHBoxLayout()
        self.BT_create = PYAFButton(self, "create_issue", "Create issue")
        self.BT_cancel = PYAFButton(self, "cancel", "Cancel")

        self.HL_buttons.addWidget(self.BT_cancel)
        self.HL_buttons.addStretch(1)
        self.HL_buttons.addWidget(self.BT_create)

        self.VL.addWidget(self.label)
        self.VL.addWidget(self.error)
        self.VL.addLayout(self.HL_buttons)

        self.setLayout(self.VL)

    def button_clicked(self, button):
        """Called when a button is clicked."""

        if button == "cancel":
            self.close()

        elif button == "create_issue":
            clipboard = QtWidgets.QApplication.clipboard()

            # Format the text for the bitbucket markdown
            text = self.error.toPlainText()
            text = text.replace("\n", "\r\n\r\n")
            text = "```\r\n\r\n#!python\r\n\r\n" + text + "```"

            clipboard.setText(str(text))

            webbrowser.open("https://bitbucket.org/adujardin/pyaf/issues/new")

            self.close()
