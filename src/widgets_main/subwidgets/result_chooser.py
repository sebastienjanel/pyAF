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

"""Widget which allows to change the type of result to display."""


import logging
from ... import shared
from ... import widgets_list
from qtpy import QtCore, QtWidgets
from ...tools.gui_tools import PYAFComboBox
from ...tools.PYAFWidget import PYAFWidget
from ...tools import misc_tools
from ...tools.results_sorting import BREAKTHROUGH_TYPES, BRUSH_TYPES, \
    MICRORHEO_TYPES, VISCO_TYPES


def groups():
    """Result types by nature, in the order of the maps of the results
    tab: (title, ((type, name), ...))."""
    return (
        ("Elasticity", (("stiffness", "Elasticity"),) + tuple(
            (name, values["name"]) for name, values in BRUSH_TYPES.items())),
        ("Viscoelasticity", tuple(
            (name, values["name"]) for name, values in VISCO_TYPES.items())),
        ("Microrheology", tuple(
            (name, values["name"])
            for name, values in MICRORHEO_TYPES.items())),
        ("Adhesion", (("work", "Detachment work"),
                      ("rupture_force", "Detachment force"),
                      ("events_forces", "Events force"),
                      ("events_per_curve", "Events per curve"),
                      ("events_rupture_force", "Event max. force"),
                      ("loading_rates", "Loading rates"),
                      ("events_distance", "Distance to Joc"))),
        ("Breakthrough", tuple(
            (name, values["name"])
            for name, values in BREAKTHROUGH_TYPES.items())))


def available(result_type):
    """True if the result type is computed for at least one file."""
    datas = shared.exp.list
    if result_type == "stiffness":
        return any(item.stiffness_array is not None for item in datas)
    if result_type in ("work", "rupture_force"):
        return any(item.jocs1_indices is not None for item in datas)
    if result_type == "loading_rates":
        return any(item.loading_rates_calculated for item in datas)
    if result_type.startswith("events_"):
        return any(item.events_calculated for item in datas)
    if result_type in VISCO_TYPES:
        return any(item.visco_calculated for item in datas)
    if result_type in BREAKTHROUGH_TYPES:
        return any(item.breakthrough_calculated for item in datas)
    if result_type in MICRORHEO_TYPES:
        return any(item.microrheo_calculated for item in datas)
    if result_type in BRUSH_TYPES:
        return any(item.brush_calculated for item in datas)
    return False


class ResultsChooserWidget(PYAFWidget):
    """List of the types of result to display, by nature (the types not
    computed are greyed out)."""

    def __init__(self, parent, name):
        super().__init__(parent, name)

        # Load the logger
        self.logger = logging.getLogger()

        self.VL = QtWidgets.QVBoxLayout()

        self.list_choser = PYAFComboBox(self, "list_choser")
        self.list_choser.setMaxVisibleItems(30)
        model = self.list_choser.model()
        # Index of each result type in the list
        self.items = {}
        for title, types in groups():
            # Title of the group: not selectable, with the text colour of the
            # theme (the disabled items are the types not computed)
            self.list_choser.addItem(title)
            item = model.item(self.list_choser.count() - 1)
            item.setFlags(QtCore.Qt.ItemIsEnabled)
            font = item.font()
            font.setBold(True)
            item.setFont(font)
            for result_type, label in types:
                self.items[result_type] = self.list_choser.count()
                self.list_choser.addItem("    " + label, result_type)

        # Wide enough for the longest result type
        self.list_choser.fit_to_contents()

        self.VL.addWidget(self.list_choser)

        self.setLayout(self.VL)

        self.update_widget()

    def update_widget(self):
        """Update the widget."""
        self.update_GUI("all")

    def list_updated(self, _):
        """Called when a list is updated."""
        result_type = self.list_choser.currentData()
        if result_type is None:
            # Title of a group (keyboard): back to the displayed type
            self.update_GUI("all")
            return
        shared.exp.results_type = result_type
        # Reset the columns which are displayed or hidden
        self.parent.tw_labels_is_checked = None
        # Refresh histograms and tables completely (will also refresh groups)
        widgets_list.widget_results_single.reset_table()
        widgets_list.widget_results_single.button_clicked("refresh_hist")

    def update_GUI(self, what):
        """Update the GUI."""
        if what == "button_type_hist" or what == "all":
            # Stiffness when computed by Tomography
            self.list_choser.setItemText(
                self.items["stiffness"],
                "    " + misc_tools.elasticity_name())
            self.list_choser.fit_to_contents()

            model = self.list_choser.model()
            for result_type, index in self.items.items():
                model.item(index).setEnabled(available(result_type))

            if shared.exp.results_type in self.items:
                self.list_choser.setCurrentIndex(
                    self.items[shared.exp.results_type])
