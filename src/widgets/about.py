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

"""About widget."""

import logging
from qtpy import QtWidgets
from .. import consts
from ..tools.PYAFWidget import PYAFWidget


class AboutWidget(PYAFWidget):
    """Widget displaying the version number and some informations."""

    def __init__(self, parent):
        super().__init__(parent, "widget_about")

        self.logger = logging.getLogger()
        self.logger.debug("Opening About widget")

        self.VL = QtWidgets.QVBoxLayout()

        link = consts.REPOSITORY.replace("https://", "")
        label = QtWidgets.QLabel(
            "<h3>pyAF " + consts.VERSION + "</h3>"
            "<p>AFM force curve analysis software dedicated to biologists.<br>"
            "CMPI lab (<a style='color: #3b8eea' href='https://cmpi.cnrs.fr'>cmpi.cnrs.fr</a>)</p>"
            "<p><b>Michka Popoff</b> (main initial coding), Antoine Dujardin, "
            "Javier Lopez-Alonso,<br>Simone Bovio, Nuno Duarte, "
            "<b>S\u00e9bastien Janel</b> (maintainer)<br>"
            "Contributors: Yann Ciczora, Vincent Dupres, Frank Lafont</p>"
            "<p>When using pyAF, cite: Stiffness tomography of eukaryotic "
            "intracellular compartments<br>by atomic force microscopy, "
            "Nanoscale, 2019, 11, 10320-10328, "
            "<a style='color: #3b8eea' href='https://doi.org/10.1039/C8NR08955H'>"
            "doi:10.1039/C8NR08955H</a></p>"
            "<p>Lab-developed software, not a commercial one: check its "
            "results against other software.<br>CeCILL license "
            "(<a style='color: #3b8eea' href='http://www.cecill.info'>cecill.info</a>). Source code "
            "and bug reports:<br><a style='color: #3b8eea' href='" + consts.REPOSITORY + "'>" + link +
            "</a></p>")
        label.setOpenExternalLinks(True)
        self.VL.addWidget(label)
        self.VL.addStretch(1)

        self.setLayout(self.VL)
