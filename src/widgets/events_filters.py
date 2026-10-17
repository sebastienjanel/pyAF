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

"""Small window allowing to change the distance or the slope filter of the
events (see tools.event_filter_tools)."""

from qtpy import QtCore, QtWidgets
from ..tools.gui_tools import PYAFButton, PYAFComboBox
from .. import widgets_list
from .. import shared
from ..tools import apply_to_all
from ..tools import events_refresher
from ..tools import event_filter_tools

DIST_TOOLTIP = (
    "Distances from the jump of contact, in positive nm.\n"
    "Remove: the events between the two distances are not counted.\n"
    "Keep only: only the events between the two distances are counted.\n"
    "0 to 0: no filter.")

SLOPE_TOOLTIP = (
    "Slope of the line fitted on the curve next to each event, in nN/nm.\n"
    "Only the events with a slope between the two values are counted.\n"
    "0 to 0: no filter.")


class ChangeFilterWidget(QtWidgets.QDialog):
    """Small window to change the distance or the slope filter.

    filt is "dist" or "slope" for the filter of the file displayed in the
    Results tab, "dist_batch" or "slope_batch" for the filters of the rows
    selected in the table of the Plots tab.
    """

    def __init__(self, parent, filt):
        super().__init__()

        self.parent = parent
        self.filt = filt
        self.is_dist = filt in ("dist", "dist_batch")
        self.batch = filt in ("dist_batch", "slope_batch")

        if self.is_dist:
            self.setWindowTitle("Change distance filter")
        else:
            self.setWindowTitle("Change slope filter")

        dist_min, dist_max, dist_keep, slope_min, slope_max = \
            self.get_current_values()

        self.IN_min = QtWidgets.QLineEdit(
            event_filter_tools.format_number(
                dist_min if self.is_dist else slope_min))
        self.IN_max = QtWidgets.QLineEdit(
            event_filter_tools.format_number(
                dist_max if self.is_dist else slope_max))
        for item in (self.IN_min, self.IN_max):
            item.setFixedWidth(70)
            item.setAlignment(QtCore.Qt.AlignRight)

        LAY_values = QtWidgets.QHBoxLayout()
        if self.is_dist:
            self.CB_mode = PYAFComboBox(self, None)
            self.CB_mode.addItem("Remove")
            self.CB_mode.addItem("Keep only")
            self.CB_mode.fit_to_contents()
            self.CB_mode.setCurrentIndex(1 if dist_keep else 0)
            LAY_values.addWidget(self.CB_mode)
            LAY_values.addWidget(QtWidgets.QLabel("events from"))
        else:
            LAY_values.addWidget(
                QtWidgets.QLabel("Keep only events with a slope from"))
        LAY_values.addWidget(self.IN_min)
        LAY_values.addWidget(QtWidgets.QLabel("to"))
        LAY_values.addWidget(self.IN_max)
        if self.is_dist:
            LAY_values.addWidget(
                QtWidgets.QLabel("nm from the jump of contact"))
        else:
            LAY_values.addWidget(QtWidgets.QLabel("nN/nm"))
        LAY_values.addStretch(1)

        self.setToolTip(DIST_TOOLTIP if self.is_dist else SLOPE_TOOLTIP)

        if self.batch:
            nbr = len(self.get_selected_rows())
            where = "the " + str(nbr) + " selected row" + \
                ("s" if nbr > 1 else "") + " of the table"
        elif shared.exp.apply_to_all_data:
            # Results tab with "apply to all" (see update_parents)
            where = "all the files"
        else:
            where = shared.exp.current_data.filename
        self.LB_info = QtWidgets.QLabel("Applies to " + where + ".")

        BT_cancel = PYAFButton(self, "cancel", "Cancel")
        BT_ok = PYAFButton(self, "ok", "OK")
        BT_ok.setDefault(True)
        LAY_buttons = QtWidgets.QHBoxLayout()
        LAY_buttons.addWidget(self.LB_info)
        LAY_buttons.addStretch(1)
        LAY_buttons.addWidget(BT_cancel)
        LAY_buttons.addWidget(BT_ok)

        LAY = QtWidgets.QVBoxLayout()
        LAY.addLayout(LAY_values)
        LAY.addLayout(LAY_buttons)
        self.setLayout(LAY)

    def get_selected_rows(self):
        """Rows selected in the table of the Plots tab."""
        table = widgets_list.widget_results_single.tableWidget
        return sorted(item.row() for item in
                      table.selectionModel().selectedRows())

    def get_current_values(self):
        """Values of the filters displayed when opening the window: the
        ones of the first selected row of the table, or of the file of the
        Results tab."""
        if self.batch:
            rows = self.get_selected_rows()
            if rows:
                result = shared.exp.results_list[rows[0]]
                return (result.filter_dist_min, result.filter_dist_max,
                        result.filter_dist_keep, result.filter_slope_min,
                        result.filter_slope_max)

        data = shared.exp.current_data
        return (data.events_results_filter_dist_min,
                data.events_results_filter_dist_max,
                data.events_results_filter_dist_keep,
                data.events_results_filter_slope_min,
                data.events_results_filter_slope_max)

    def button_clicked(self, name):
        """Called when a button is clicked."""
        if name == "cancel":
            self.close()
            return

        try:
            val_min = float(self.IN_min.text().replace(",", "."))
            val_max = float(self.IN_max.text().replace(",", "."))
        except ValueError:
            self.LB_info.setText("The limits must be numbers.")
            return
        keep = self.is_dist and self.CB_mode.currentIndex() == 1

        if self.batch:
            wl = widgets_list.widget_results_single
            for row in self.get_selected_rows():
                result = shared.exp.results_list[row]
                if self.is_dist:
                    result.filter_dist_min = val_min
                    result.filter_dist_max = val_max
                    result.filter_dist_keep = keep
                else:
                    result.filter_slope_min = val_min
                    result.filter_slope_max = val_max
                wl.update_filter_cells(row)

            # Recompute the values with the new filters
            wl.schedule_refresh()

        else:
            data = shared.exp.current_data
            if self.is_dist:
                data.events_results_filter_dist_min = val_min
                data.events_results_filter_dist_max = val_max
                data.events_results_filter_dist_keep = keep
            else:
                data.events_results_filter_slope_min = val_min
                data.events_results_filter_slope_max = val_max
            self.update_parents()
            self.update_plots_rows()

        self.close()

    def update_plots_rows(self):
        """Give the filter of the file (and of the other files if it was
        applied to all) to its rows in the table of the Plots tab."""
        if shared.exp.apply_to_all_data:
            data_ids = range(len(shared.exp.list))
        else:
            data_ids = [shared.exp.id_selected]

        wl = widgets_list.widget_results_single
        for row, result in enumerate(shared.exp.results_list):
            if result.data_id not in data_ids:
                continue
            data = shared.exp.list[result.data_id]
            if self.is_dist:
                result.filter_dist_min = data.events_results_filter_dist_min
                result.filter_dist_max = data.events_results_filter_dist_max
                result.filter_dist_keep = \
                    data.events_results_filter_dist_keep
            else:
                result.filter_slope_min = \
                    data.events_results_filter_slope_min
                result.filter_slope_max = \
                    data.events_results_filter_slope_max
            if wl is not None:
                wl.update_filter_cells(row)

        # The Plots tab is refreshed when it is displayed
        if wl is not None:
            wl.schedule_refresh()

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
