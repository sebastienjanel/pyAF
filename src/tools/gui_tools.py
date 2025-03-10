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

"""Some tools for the graphical interface of PYAF."""

import logging
import qtwidgets
from .. import shared
from .. import widgets_list
from PyQt5 import QtCore, QtWidgets
from ..plots.PYAFPlot import PYAFPlot


class BoxData(QtWidgets.QWidget):
    """A box with the piezo meshgrid and minimal information for the data tab."""

    def __init__(self, i):
        super().__init__()

        self.setFixedSize(290, 72)

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        self.setSizePolicy(sizePolicy)

        self.data_id = i
        data = shared.exp.list[self.data_id]
        self.smallfont = QtWidgets.QApplication.font()
        self.smallfont.setPointSize(10)
        self.smallfontbold = QtWidgets.QApplication.font()
        self.smallfontbold.setPointSize(10)
        self.smallfontbold.setBold(True)

        self.VL = QtWidgets.QVBoxLayout()
        self.VL.setContentsMargins(0, 0, 0, 0)

        self.widget = QtWidgets.QWidget()
        self.widget.setObjectName("data_" + str(i))

        self.widget.mousePressEvent = lambda event: self.change_id()

        self.HL_widget = QtWidgets.QHBoxLayout()
        self.HL_widget.setContentsMargins(0, 0, 0, 0)

        self.VL_labels = QtWidgets.QVBoxLayout()

        self.label_name = QtWidgets.QLabel()
        self.label_type_of_file = QtWidgets.QLabel()
        self.label_date = QtWidgets.QLabel()

        self.VL_labels.addStretch(1)
        self.VL_labels.addWidget(self.label_name)
        self.VL_labels.addWidget(self.label_type_of_file)
        self.VL_labels.addWidget(self.label_date)
        self.VL_labels.addStretch(1)

        # Canvas Meshgrid
        self.VL_canvas = QtWidgets.QVBoxLayout()
        self.canvas = QtWidgets.QWidget()
        self.canvas.setFixedSize(72, 72)

        if data.is_single is False:
            self.meshgrid = PYAFPlot(
                self, "small_meshgrid", self.canvas,
                [72.0 / 72.0, 72.0 / 72.0, 72.0],
                data_id=self.data_id)
        else:
            self.meshgrid = None

        self.VL_canvas.addWidget(self.canvas)
        self.VL_canvas.addStretch(1)

        self.HL_widget.addLayout(self.VL_canvas)
        self.HL_widget.addLayout(self.VL_labels)

        self.widget.setLayout(self.HL_widget)

        self.VL.addWidget(self.widget)

        self.setLayout(self.VL)

        self.update_widget()

    def setStyleSheet(self, styleSheet):
        """Set the styleSheet. (Used to change the background color)"""
        self.widget.setStyleSheet(styleSheet)

    def update_widget(self):
        """Called to update the widget."""
        data = shared.exp.list[self.data_id]

        if self.meshgrid is not None:
            self.meshgrid.update_plot()

        self.label_name.setText(str(data.filename))
        self.label_type_of_file.setText(data.file_type)
        self.label_date.setText(data.date)

        # Reupdate the fonts each time or it will not work here, I don't
        # really know why I can't directly do self.setFont() ...
        self.label_name.setFont(self.smallfontbold)
        self.label_type_of_file.setFont(self.smallfont)
        self.label_date.setFont(self.smallfont)

    def change_id(self):
        """When there is a click on the box, go to the selected file."""
        # Set the id of the selected file
        shared.exp.id_selected = self.data_id

        # Update the two other lists
        # The main widget is not directly registered in widgets_list,
        # but we can access it trough the parent of one of the tabs.
        widgets_list.widget_data.parent.file_changed("box")


class PYAFTableWidget(QtWidgets.QTableWidget):
    """Custom tablewidget which allows to redefine the widget's size."""

    def __init__(self, width, height, length, nbr):
        super().__init__(length, nbr)

        self.width = width
        self.height = height

    def sizeHint(self):
        """Redefine the sizehint."""
        return QtCore.QSize(self.width, self.height)


class PYAFButtonGroup(QtWidgets.QButtonGroup):
    """Custom button group."""

    def __init__(self, parent, name):
        super().__init__()

        self.name = name
        self.parent = parent

        self.logger = logging.getLogger()

        self.buttonClicked.connect(self.button_clicked)

    def button_clicked(self):
        """Called when the button is clicked.

        Calls the parent button_clicked method and logs the button which was
        clicked.
        """
        self.logger.debug("GroupButton clicked : %s", self.name)

        self.parent.button_clicked(self.name)


class PYAFInput(QtWidgets.QWidget):
    """An input field with a label on the left."""

    def __init__(self, parent, name, label_text, width=50, height=30):
        super().__init__()

        self.layout = QtWidgets.QHBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)

        self.name = name

        self.logger = logging.getLogger()

        self.input = QtWidgets.QLineEdit()
        self.input.setFixedSize(width, height)
        if "_input_updated" in dir(parent):
            # If the parent class is a PYAFWidget (that's what the dir() check
            # is for), then it has an _input_updated method. This method
            # catches the editingFinished signal and sends it to input_updated
            # method defined in the "real" parent class which is you widget.
            # See in PYAFWidget how this works. In your code, just use
            # input_updated without the first underscore.
            self.input.editingFinished.connect(
                lambda: parent._input_updated(name))
        else:
            # In case the parent is not a PYAFWidget but for example a QWidget,
            # it has not _input_updated method (or it has one if you implented)
            # one, but normally you just have an input_updated method, so
            # that's what is called here.
            self.input.editingFinished.connect(
                lambda: parent.input_updated(name))
        self.label = QtWidgets.QLabel()
        # Add label text plus some space between label and Input
        self.label.setText("  " + label_text)

        self.layout.addWidget(self.input)
        self.layout.addWidget(self.label)
        self.layout.addStretch(1)
        self.setLayout(self.layout)

    def changeValue(self, value):
        """Change the value in the input field."""
        self.logger.debug("Input updated: %s, value: %s", self.name, value)

        self.input.setText(str(value))

    def setEnabled(self, value):
        """Enable or disable the input."""
        self.input.setEnabled(value)

    def setFont(self, value):
        """Set the font of the input and of the label."""
        self.input.setFont(value)
        self.label.setFont(value)

    def get_float_value(self):
        """Returns the value in the input as a python float."""
        return float(str(self.input.text()).replace(",", "."))

    def get_str_value(self, use_replace=False):
        """Returns the value in the input as a python string.

        If the use_replace argument is set to True, the commas are replaced
        by dots (useful for transforming it in float afterwards).
        """
        if use_replace:
            val = str(self.input.text()).replace(",", ".")
        else:
            val = str(self.input.text())

        return val

    def get_int_value(self):
        """Returns the value in the input as a python int."""
        return int(str(self.input.text()))


class QLineEditForTableWidet(QtWidgets.QLineEdit):
    """This is a subclass of a QtWidgets.QLineEdit.

    It adds directly the connected signal, and passes the row/column.
    """

    def __init__(self, parent, col, row):
        super().__init__()

        self.editingFinished.connect(lambda: parent.input_updated(col, row))


class CenteredCellwidget(QtWidgets.QWidget):
    """A centered widget, can be used for tablewidgets."""

    def __init__(self, item):
        super().__init__()

        self.item = item
        layout = QtWidgets.QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignHCenter)
        layout.addWidget(item)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)


class CenteredCellwidget2(QtWidgets.QWidget):
    """A centered widget, can be used for tablewidgets."""

    def __init__(self, items):
        super().__init__()

        self.items = items
        layout = QtWidgets.QHBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignHCenter)
        for item in items:
            layout.addWidget(item)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)


class CenteredCellCheckbox(QtWidgets.QWidget):
    """A centered checkbox for the tablewidgets."""

    def __init__(self, parent, name, cell_id):
        super().__init__()

        self.cell_id = cell_id
        layout = QtWidgets.QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignHCenter)
        self.checkbox = QtWidgets.QCheckBox()
        layout.addWidget(self.checkbox)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)

        self.checkbox.clicked.connect(
            lambda: parent.checkbox_clicked(name, cell_id))



class CenteredCellInput(QtWidgets.QWidget):
    """A centered input field for the tablewidgets."""

    def __init__(self, parent, name, cell_id):
        super().__init__()

        layout = QtWidgets.QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignHCenter)
        self.input = QtWidgets.QLineEdit()
        layout.addWidget(self.input)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)

        self.input.editingFinished.connect(
            lambda: parent.input_updated(name, cell_id))


class CenteredCellRadioButton(QtWidgets.QWidget):
    """A centered radiobutton for the tablewidgets."""

    def __init__(self):
        super().__init__()

        layout = QtWidgets.QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignHCenter)
        self.radio_button = QtWidgets.QRadioButton()
        layout.addWidget(self.radio_button)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)


class CenteredCellComboBox(QtWidgets.QWidget):
    """A centered drop down list for the tablewidgets."""

    def __init__(self, combobox):
        super().__init__()

        self.list = combobox
        layout = QtWidgets.QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignHCenter)
        layout.addWidget(self.list)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)


class PYAFButton(QtWidgets.QPushButton):
    """Reimplementation of the PyQT QPushButton.

    Will create a QPushButton, connected to the parent's button_clicked method.
    """

    def __init__(self, parent, name, text, size=None):
        super().__init__()

        self.parent = parent
        self.name = name

        # Load the logger
        self.logger = logging.getLogger()

        if size is None:
            size = 150
        self.setFixedWidth(size)
        self.setFixedHeight(40)
        self.setText(text)
        self.clicked.connect(self.button_clicked)

    def button_clicked(self):
        """Called when the button is clicked.

        Calls the parent button_clicked method and logs the button which was
        clicked.
        """
        self.logger.debug("Button clicked : %s", self.name)

        self.parent.button_clicked(self.name)


class PYAFCheckBox(QtWidgets.QCheckBox):
    """Reimplementation of the PyQT QCheckBox.

    Will create a QCheckBox, connected to the parent's checkbox_clicked method.
    """

    def __init__(self, parent, name, text):
        super().__init__()

        self.parent = parent
        self.name = name

        # Load the logger
        self.logger = logging.getLogger()

        self.setText(text)
        self.clicked.connect(self.checkbox_clicked)

    def checkbox_clicked(self):
        """Called when the checkbox is clicked.

        Calls the parent checkbox_clicked method and logs the checkbox which
        was clicked.
        """
        self.logger.debug("Checkbox clicked : %s", self.name)

        self.parent.checkbox_clicked(self.name)


class PYAFComboBox(QtWidgets.QComboBox):
    """Reimplementation of the PyQT QComboBox.

    The difference here is that scrolling with the mouse through the list is
    disabled, because this is a little bugged under OS X. Sometimes you would
    scroll through the list without wanting to do it, which can be confusing.

    The activated signal is connected to the parent's list_updated method.
    """

    def __init__(self, parent, name):
        super().__init__()

        self.name = name
        self.parent = parent

        # Load the logger
        self.logger = logging.getLogger()

        if name is not None:
            self.activated.connect(self.list_updated)

        self.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToContents)

    def list_updated(self):
        """Catch the list updated signal, and trigger the parent call.

        Adds a log input for debugging.
        """
        self.logger.debug("List updated : %s", self.name)

        self.parent.list_updated(self.name)

    def wheelEvent(self, event):
        """Redefine wheelEvent. It will be disabled."""
        pass


class PYAFToggle(QtWidgets.QWidget):
    """Reimplementation of the PyQT Toggle.

    Will create a Toggle, connected to the parent's toggle_clicked method.
    """

    def __init__(self, parent, name, text, size=None):
        super().__init__()

        self.parent = parent
        self.name = name

        self.layout = QtWidgets.QHBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)

        # Load the logger
        self.logger = logging.getLogger()

        if size is None:
            size = 100

        self.toggle = qtwidgets.Toggle(checked_color="#00FF44")
        self.toggle.setFixedWidth(70)
        self.toggle.clicked.connect(self.toggle_clicked)

        self.label = QtWidgets.QLabel()
        # Add label text plus some space between label and Toggle
        self.label.setText(text + "  ")

        self.layout.addWidget(self.label)
        self.layout.addWidget(self.toggle)
        self.layout.addStretch(1)
        self.setLayout(self.layout)

    def toggle_clicked(self):
        """Called when the toggle is clicked.

        Calls the parent toggle_clicked method and logs the toggle which was
        clicked.
        """
        self.logger.debug("Toggle clicked : %s", self.name)

        self.parent.toggle_clicked(self.name)


class ClearWidgetsFromLayout:
    """This class lets you remove all widgets from a layout."""

    def __init__(self, layout):
        while layout.count() > 0:
            item = layout.takeAt(0)
            if not item:
                continue
            widget = item.widget()
            if widget:
                widget.deleteLater()


def disable_combobox_item(combobox, position):
    """Function used to disable a row in a qcombobox."""
    # http://theworldwideinternet.blogspot.fr/
    # 2011/01/disabling-qcombobox-items.html

    index = combobox.model().index(position, 0)
    combobox.model().setData(index, QtCore.QVariant(0), QtCore.Qt.UserRole - 1)


def enable_combobox_item(combobox, position):
    """Function used to enable a row in a qcombobox."""
    index = combobox.model().index(position, 0)
    combobox.model().setData(index,
                             QtCore.QVariant(1 | 32), QtCore.Qt.UserRole - 1)
