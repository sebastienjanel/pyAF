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

"""Widget used to go through the curves."""

from qtpy import QtCore, QtWidgets

from ... import widgets_list
from ...tools.gui_tools import PYAFButton
from ...tools.gui_tools import PYAFComboBox
from ...tools.gui_tools import PYAFToggle
from ...tools.PYAFWidget import PYAFWidget


class PlotSelectorWidget(PYAFWidget):
    """The plot selector allows to go through different plot types."""

    def __init__(self, parent, name):
        super().__init__(parent, name)

        if self.name == "plot_selector_data":
            self.parent_widget = widgets_list.widget_data

        # Same order as the plots of the data tab. The time plots are only
        # available for the files with a time (JPK, Nanoscope)
        self.plot_types = ("Deflection vs Z position", "Force vs distance",
                           "Deflection vs Time", "Z position vs Time")
        self.time_plots = (2, 3)

        self.BT_back = PYAFButton(self, "back", "<", size=60)
        self.BT_forward = PYAFButton(self, "forward", ">", size=60)

        self.list_plot_types = PYAFComboBox(self, None)

        for plot_type in self.plot_types:
            self.list_plot_types.addItem(plot_type)
        self.list_plot_types.setItemData(
            1, "Force vs tip-sample distance (Z position - deflection). On a "
            "hard substrate\n(glass), the contact is vertical when the "
            "deflection sensitivity is right.", QtCore.Qt.ToolTipRole)
        # Wide enough for the longest plot type
        self.list_plot_types.fit_to_contents()

        self.list_plot_types.activated.connect(lambda: self.list_updated("plot_types"))

        # Zoom on the contact of the force vs distance plot
        self.TG_zoom = PYAFToggle(self, "zoom", "Zoom on contact")
        self.TG_zoom.setToolTip(
            "Force vs distance: zoom on the contact of the approach curve, to "
            "see if it is vertical\non a hard substrate (deflection "
            "sensitivity).")
        zoom = QtCore.QSettings().value("ForceDistanceZoom", False)
        self.TG_zoom.setChecked(zoom not in ("false", False, None))

        self.HL_select_curve = QtWidgets.QHBoxLayout()
        self.HL_select_curve.addStretch(1)
        self.HL_select_curve.addWidget(self.BT_back)
        self.HL_select_curve.addWidget(self.list_plot_types)
        self.HL_select_curve.addWidget(self.BT_forward)
        self.HL_select_curve.addStretch(1)

        HL_zoom = QtWidgets.QHBoxLayout()
        HL_zoom.addStretch(1)
        HL_zoom.addWidget(self.TG_zoom)
        HL_zoom.addStretch(1)
        VL = QtWidgets.QVBoxLayout()
        VL.setContentsMargins(0, 0, 0, 0)
        VL.addLayout(self.HL_select_curve)
        VL.addLayout(HL_zoom)
        self.setLayout(VL)

        self.update_widget()

    def toggle_clicked(self, _name):
        """Zoom on the contact of the force vs distance plot (or not)."""
        QtCore.QSettings().setValue("ForceDistanceZoom",
                                    self.TG_zoom.isChecked())
        self.parent_widget.curve_3.update_plot()

    def available(self):
        """Indexes of the available plot types."""
        model = self.list_plot_types.model()
        return [i for i in range(len(self.plot_types))
                if model.item(i).isEnabled()]

    def set_time_plots(self, enabled):
        """Enable or disable the time plots (files with a time)."""
        model = self.list_plot_types.model()
        for i in self.time_plots:
            model.item(i).setEnabled(enabled)
        if self.list_plot_types.currentIndex() not in self.available():
            self.list_plot_types.setCurrentIndex(0)
        self.update_widget()

    def update_widget(self):
        """Update the widget."""
        available = self.available()
        index = self.list_plot_types.currentIndex()
        # The zoom is only for the force vs distance plot
        self.TG_zoom.setEnabled(index == 1)
        self.BT_back.setEnabled(any(i < index for i in available))
        self.BT_forward.setEnabled(any(i > index for i in available))

    def button_clicked(self, button):
        """Called when a button is clicked."""
        index = self.list_plot_types.currentIndex()

        available = self.available()
        before = [i for i in available if i < index]
        after = [i for i in available if i > index]
        if (button == "back" and not before) or \
                (button == "forward" and not after):
            return
        if button == "back":
            new_index = max(before)
            self.list_plot_types.setCurrentIndex(new_index)
            self.parent_widget.canvas_stack.setCurrentIndex(new_index)
            self.update_widget()

        if button == "forward":
            new_index = min(after)
            self.list_plot_types.setCurrentIndex(new_index)
            self.parent_widget.canvas_stack.setCurrentIndex(new_index)
            self.update_widget()

    def list_updated(self, name):
        """Called when a list is updated."""

        if name == "plot_types":
            self.parent_widget.canvas_stack.setCurrentIndex(self.list_plot_types.currentIndex())
            self.update_widget()
