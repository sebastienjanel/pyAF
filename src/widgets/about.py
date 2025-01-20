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

"""About widget."""

import logging
from PyQt5 import QtWidgets
from .. import consts
from ..tools.PYAFWidget import PYAFWidget


class AboutWidget(PYAFWidget):
    """Widget displaying the version number and some informations."""

    def __init__(self, parent):
        super().__init__(parent, "widget_about")

        self.logger = logging.getLogger()
        self.logger.debug("Opening About widget")

        self.VL = QtWidgets.QVBoxLayout()

        label_version = QtWidgets.QLabel("pyAF version " + consts.VERSION)
        label_authors1 = QtWidgets.QLabel("<b>Michka Popoff<b/>")
        label_authors2 = QtWidgets.QLabel(
            "Sebastien Janel, Simone Bovio, Yann Ciczora, Vincent Dupres, Antoine Dujardin, Javier Lopez-Alonso, Nuno Duarte")
        label_authors3 = QtWidgets.QLabel("Frank Lafont")
        label_website = QtWidgets.QLabel("https://bitbucket.org/cmip/pyaf-advanced/")

        self.VL.addWidget(label_version)
        self.VL.addWidget(label_authors1)
        self.VL.addWidget(label_authors2)
        self.VL.addWidget(label_authors3)
        self.VL.addWidget(label_website)
        self.VL.addStretch(1)

        self.setLayout(self.VL)
