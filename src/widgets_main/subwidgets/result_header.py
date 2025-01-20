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

"""
Custom table header class allowing to embed a checkbox on the header first column.
This checkbox is used to allow the user to select all the rows at once.
"""

from PyQt5.QtCore import Qt, QRect, pyqtSignal
from PyQt5.QtWidgets import QTableWidget, QApplication, QHeaderView, QStyleOptionButton, QStyle


class QCheckableHeader(QHeaderView):

    # Used to handle row selection
    select_all = False
    all_selected = pyqtSignal()

    def __init__(self, orientation, parent=None):
        QHeaderView.__init__(self, orientation, parent)

    def paintSection(self, painter, rect, logicalIndex):
        painter.save()
        QHeaderView.paintSection(self, painter, rect, logicalIndex)
        painter.restore()

        if logicalIndex == 0:
            option = QStyleOptionButton()
            checkbox_rect = QRect(3, 10, 16, 16)
            checkbox_rect.moveCenter(rect.center())
            option.rect = checkbox_rect

            if self.select_all:
                option.state = QStyle.State_On

            else:
                option.state = QStyle.State_Off
            self.style().drawControl(QStyle.CE_CheckBox, option, painter)

    def mousePressEvent(self, event):
        if self.logicalIndexAt(event.pos()) == 0:
            self.select_all = not self.select_all
            self.all_selected.emit()
            self.updateSection(0)
            QHeaderView.mousePressEvent(self, event)

    def changeState(self, state):
        self.select_all = state
        self.updateSection(0)
