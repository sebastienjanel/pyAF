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

"""Subclasses QWidget."""

from PyQt5 import QtWidgets
from .. import widgets_list


class PYAFWidget(QtWidgets.QWidget):
    """This class is used as parent class for all the widgets in PYAF.

    The widget is automatically registered in widgets_list. If the widget is
    closed his value in *widgets_list* is set to *None*, else it is set to
    *self*, so that the widget can be accessed through *widgets_list* from
    everywhere in PYAF.

    The parent argument is saved as self.parent, so that it can be used
    in the widgets to access the widget's parent.
    """

    def __init__(self, parent, full_name):
        super().__init__()
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Minimum,
            QtWidgets.QSizePolicy.Minimum)

        self.open = True
        self.parent = parent
        self.name = extract_widget_name(full_name)

        # Register the widget itself in the widgets_list
        self.modify_state(self)

    def modify_state(self, value):
        """Allows to register or to unregister the widget in the widgets_list."""
        # Get the list of widgets in widgets_list
        widgets = dir(widgets_list)

        # Parse the list of widgets to find the widgets
        for widget in widgets:
            # Get the name of the wigdet
            widget_name = extract_widget_name(widget)

            # If the name corresponds, set the value
            if widget_name == self.name:
                name = "widget_" + self.name
                setattr(widgets_list, name, value)

    def closeEvent(self, event):
        """Method called when widget is closed.

        Set *self.open* to *False*, so that the input_updated method can no
        more be called. (see *input_updated* method for explanations).

        Parse the *widgets_list* and reset the widget's attribute to *None*, so
        that it can no more be called from outside.
        """
        self.open = False

        # Set the widget to None in widgets_list
        self.modify_state(None)

        # Close the widget
        event.accept()

    def _input_updated(self, val, i=None):
        """Method which redirects the input signals to the user defined method.

        Sometimes, editingFinished is called at close by PyQt (it's called
        before the closing of the widget). This would call input_updated, but
        as the tmp file is closed (and shared.exp destroyed), this will result
        in an error. To prevent this, input_updated can only be called when
        *self.open* is set to *True*. (*self.open* is set to *False* in
        *closeEvent*, before the destruction.
        (Note : this is a hack, perhaps QT/PyQT will solve the problem on their
        side one day).
        """
        # Can be called only if the widget still exists
        if self.open:
            if i is None:
                # In the normal case, there is only the value to pass
                self.input_updated(val)
            else:
                # When the inputs are in a tablewidget I pass the row
                self.input_updated(val, i)


def extract_widget_name(full_name):
    """Extract the name of the widget from the full name."""
    # Get the name of the widget
    splitted = full_name.split("_")
    widget_name = splitted[1]
    for i in range(len(splitted)):
        if i > 1:
            widget_name = widget_name + "_" + splitted[i]
    return widget_name
