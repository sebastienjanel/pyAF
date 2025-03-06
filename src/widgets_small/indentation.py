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

"""Small widget with a list used to change the indentation depth."""

from PyQt5 import QtCore, QtWidgets
from ..tools.gui_tools import PYAFComboBox
from .. import widgets_list
from ..tools import apply_to_all
from ..tools.PYAFWidget import PYAFWidget
from .. import shared


class WidgetIndentation(PYAFWidget):
    """Small widget with a list used to change the indentation depth."""

    def __init__(self, parent):
        super().__init__(parent, "widget_indentation")

        # Load the settings
        self.settings = QtCore.QSettings()

        self.data_id = None

        self.smallfont = QtWidgets.QApplication.font()
        self.smallfont.setPointSize(10)

        self.VL = QtWidgets.QVBoxLayout()

        self.list_ind = PYAFComboBox(self, "indentation")
        self.list_ind.setFont(self.smallfont)

        self.VL.addWidget(self.list_ind)
        self.setLayout(self.VL)

        self.update_widget()

    def list_updated(self, name):
        """Called when the list is updated by the user."""
        settings_name = "refreshBackgroundInOpenGL"
        refresh = self.settings.value(settings_name, True)

        w_mesh = widgets_list.widget_multimeshgrids
        w_opengl = widgets_list.widget_vtk

        if name == "indentation":
            if self.data_id is None:
                data = shared.exp.current_data
            else:
                data = shared.exp.list[self.data_id]

            data.stiffness_depth_view = self.list_ind.currentIndex()
            apply_to_all.apply_to_all("display_options", oldid=self.data_id)

            # Refresh meshgrid
            if w_opengl is None or (w_opengl is not None and refresh):
                widgets_list.widget_results.update_MPL("MPL_canvas1")

            # Refresh multimeshgrids
            cond1 = (w_mesh is not None and refresh and w_opengl is not None)
            if cond1 or (w_mesh is not None and w_opengl is None):
                widgets_list.widget_multimeshgrids.update_widget()

            # Refresh vtk
            if widgets_list.widget_vtk is not None:
                widgets_list.widget_vtk.update_indentation()

    def update_widget(self, data_id=None):
        """Updates the indentation widget."""
        # Once the value is set change it only if updated (for vtk widget)
        if data_id is not None:
            self.data_id = data_id

        if self.data_id is None:
            data = shared.exp.current_data
        else:
            data = shared.exp.list[self.data_id]

        self.list_ind.clear()
        if shared.exp.apply_to_all_data:
            length = shared.exp.global_max_indentation_index
        else:
            length = data.max_indentation_index

        for i in range(length):
            if shared.exp.apply_to_all_data:
                # Display the current slice for the current dataset, and the
                # slice number (for the whole experiment, in case there are
                # different sizes of segements which were used)

                text = get_segment_text(data, i)

                text = "Seg. nbr " + str(i + 1) + text

            else:
                text = get_segment_text(data, i)

            self.list_ind.addItem(text)

        self.list_ind.setCurrentIndex(data.stiffness_depth_view)


def get_segment_text(data, i):
    """Returns the value to be displayed in the indentation list.

    If the stiffness is not calculated for the current file, just return
    an empty string.
    """
    if data.stiffness_calculated:
        if not data.used_tomography:
            start = str(0)
        else:
            start = str(i * data.used_indentation_step)

        step = data.used_indentation_step
        end = str(i * step + step)

        return " . Current file : " + start + " - " + end + " nm"

    else:
        return ""
