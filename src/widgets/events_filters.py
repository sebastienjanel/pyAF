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

"""Small widget allowing to change the distance filter."""

from PyQt5 import QtWidgets
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFButton
from .. import widgets_list
from .. import shared
from ..tools import apply_to_all
from ..tools import events_refresher


class ChangeFilterWidget(QtWidgets.QDialog):
    """Small widget for the changing of the distance filter.

    The distance filter can be changed on the results tab. It is used to
    exclude events between the jump of contact and a specific distance.
    """

    def __init__(self, parent, filt):
        super().__init__()

        self.open = True
        self.parent = parent
        self.filt = filt

        dt = shared.exp.current_data

        LAY = QtWidgets.QHBoxLayout()

        if self.filt == "slope" or self.filt == "slope_batch":
            self.IN_slope_min = PYAFInput(self, "input_min", "Min", width=100)
            self.IN_slope_max = PYAFInput(self, "input_max", "Max", width=100)

            self.IN_slope_min.changeValue(dt.events_results_filter_dist_left)
            self.IN_slope_max.changeValue(dt.events_results_filter_dist_right)

            LAY.addWidget(self.IN_slope_min)
            LAY.addWidget(self.IN_slope_max)

        elif self.filt == "dist" or self.filt == "dist_batch":
            LAY = QtWidgets.QHBoxLayout()

            self.IN_left = PYAFInput(self, "input_left", "Left", width=100)
            self.IN_right = PYAFInput(self, "input_right", "Right", width=100)
            self.CB_middle = PYAFCheckBox(self, "cb", "Keep middle")

            self.IN_left.changeValue(dt.events_results_filter_dist_left)
            self.IN_right.changeValue(dt.events_results_filter_dist_right)
            self.CB_middle.setChecked(
                dt.events_results_filter_dist_keep_middle)

            LAY.addWidget(self.IN_left)
            LAY.addWidget(self.IN_right)
            LAY.addWidget(self.CB_middle)

        if self.filt == "dist_batch" or self.filt == "slope_batch":
            BT_ok = PYAFButton(self, "ok", "OK")
            BT_ok.setDefault(True)
            LAY.addWidget(BT_ok)

        self.setLayout(LAY)

    def closeEvent(self, _):
        """Called when closing the widget."""
        self.open = False

    def _input_updated(self, val):
        """Called when the input is updated."""
        # Can be called only if the widget still exists
        if self.open:
            self.input_updated(val)

    def input_updated(self, name):
        """Called when the input is updated."""
        data = shared.exp.current_data
        filt = self.filt
        update = False

        if filt == "dist_batch":
            # Do nothing, the values are validated when the user clicks on OK.
            return False

        if name == "input_min":
            val = self.IN_slope_min.get_float_value()
            if val != data.events_results_filter_slope_min:
                update = True
                data.events_results_filter_slope_min = val

        elif name == "input_max":
            val = self.IN_slope_max.get_float_value()
            if val != data.events_results_filter_slope_max:
                update = True
                data.events_results_filter_slope_max = val

        elif name == "input_left":
            val_left = self.IN_left.get_float_value()
            if val_left != data.events_results_filter_dist_left:
                update = True
                data.events_results_filter_dist_left = val_left

        elif name == "input_right":
            val_right = self.IN_right.get_float_value()
            if val_right != data.events_results_filter_dist_right:
                update = True
                data.events_results_filter_dist_right = val_right

        # Update parent plots if we have been called from the results tab
        # and not from the plots tab. Check also if a value has changed
        # to prevent too many unwanted auto-updates.
        if filt != "dist_batch" and filt != "slope_batch" and update:
            self.update_parents()

    def checkbox_clicked(self, name):
        """Called when the checkbox is clicked."""
        data = shared.exp.current_data

        if self.filt == "dist_batch":
            # Do nothing, the values are validated when the user clicks on OK.
            return False

        if name == "cb":
            val_cb = self.CB_middle.isChecked()
            data.events_results_filter_dist_keep_middle = val_cb

        # Update parent plots if we have been called from the results tab
        # and not from the plots tab.
        if self.filt != "dist_batch" and self.filt != "slope_batch":
            self.update_parents()

    def button_clicked(self, name):
        """Called when a button is clicked."""
        if name == "ok":
            if self.filt == "dist_batch":
                val_left = self.IN_left.get_float_value()
                val_right = self.IN_right.get_float_value()
                val_cb = self.CB_middle.isChecked()
            elif self.filt == "slope_batch":
                val_slope_min = self.IN_slope_min.get_float_value()
                val_slope_max = self.IN_slope_max.get_float_value()

            wl = widgets_list.widget_results_single
            for item in wl.tableWidget.selectionModel().selectedRows():
                index = item.row()
                if self.filt == "dist_batch":
                    col = wl.columns.get_column("Dist").position
                    if col is not None:
                        cell = wl.tableWidget.cellWidget(index, col)
                        cell.items[0].setText(str(val_left))
                        cell.items[1].setText(str(val_right))
                        cell.items[2].setChecked(val_cb)
                elif self.filt == "slope_batch":
                    col = wl.columns.get_column("Slope").position
                    if col is not None:
                        cell = wl.tableWidget.cellWidget(index, col)
                        cell.items[0].setText(str(val_slope_min))
                        cell.items[1].setText(str(val_slope_max))

            self.close()

    def update_parents(self):
        """Update parents plots and values.

        The events per curve and rupture force 2 arrays are updated here, and
        not in tools.apply_to_all, because this is a very slow process which
        needs a specific progress bar.
        """
        # Apply the updated values to all the datasets
        apply_to_all.apply_to_all("display_options")

        # Update the events array (with a progressbar)
        events_refresher.update_events()

        # Update the plots
        widgets_list.widget_results.update_MPL("MPL_canvas1")
        widgets_list.widget_results.update_MPL("MPL_canvas2")

        if widgets_list.widget_multimeshgrids is not None:
            widgets_list.widget_multimeshgrids.update_widget()
