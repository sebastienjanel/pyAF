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

"""Preferences widget, called from PYAF's menu."""

import logging
from .. import shared
import matplotlib
from packaging.version import Version
from PyQt5 import QtCore, QtWidgets
from ..tools import misc_tools
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFComboBox
from ..gui_styles.theme_handler import theme_handler
from .. import widgets_list


class PreferencesWidget(PYAFWidget):
    """Widget with options for PYAF."""

    def __init__(self, parent):
        super().__init__(parent, "widget_preferences")

        self.parent = parent

        # Load logger
        self.logger = logging.getLogger()
        self.logger.debug("Opening Preferences Widget")

        self.resize(400, 400)

        # Load the settings
        self.settings = QtCore.QSettings()

        self.VL = QtWidgets.QVBoxLayout()

        text = "Random Curve Control"
        tip = "The random button will only display curves once. Can be slow."
        self.CB_rand_ctrl = PYAFCheckBox(self, "rand_ctrl", text)
        self.CB_rand_ctrl.setToolTip(tip)

        tip = ("The defl-ext and force-ext curves will have their first values"
               " set to 0 deflection/force.")
        self.CB_zero_defl = PYAFCheckBox(self, "zero_defl", "Zero deflection")
        self.CB_zero_defl.setToolTip(tip)

        text = "Display raw data on curves"
        tip = "Will only be applied to the curve on the data tab."
        self.CB_raw_data = PYAFCheckBox(self, "raw_data", text)
        self.CB_raw_data.setToolTip(tip)

        text = "Display fit noise on curves"
        tip = "The fit noises are computed on the fitting segments."
        self.CB_fit_noise = PYAFCheckBox(self, "fit_noise", text)
        self.CB_fit_noise.setToolTip(tip)

        text = "Remember to save popup"
        tip = "Display or hide the popup with the remember to save message."
        self.CB_remember_to_save = PYAFCheckBox(self, "remember_to_save", text)
        self.CB_remember_to_save.setToolTip(tip)

        text = "Use classical meshgrid style"
        tip = "Changes the style of the meshgrid"
        self.CB_mesh_style = PYAFCheckBox(self, "mesh_style", text)
        self.CB_mesh_style.setToolTip(tip)
        if Version(matplotlib.__version__) < Version("1.4.0"):
            self.CB_mesh_style.setEnabled(False)

        self.IN_main_label = \
            PYAFInput(self, "main_label", "Title label ", width=200)

        text = "Plots font size "
        tip = "Is set to 12 by default"
        self.IN_labels_font_size = PYAFInput(self, "labels", text, width=40)

        text = "Tick labels font size "
        tip = "Set this value to 12 points for printing. Default is 8"
        self.IN_tick_labels_font_size = \
            PYAFInput(self, "tick_labels", text, width=40)

        self.list_theme = PYAFComboBox(self, "list_theme")
        self.list_theme.setFixedSize(200, 26)
        for theme in shared.exp.theme_list:
            self.list_theme.addItem(theme)

        self.VL.addWidget(self.CB_rand_ctrl)
        self.VL.addWidget(self.CB_zero_defl)
        self.VL.addWidget(self.CB_raw_data)
        self.VL.addWidget(self.CB_fit_noise)
        self.VL.addWidget(self.CB_remember_to_save)
        self.VL.addWidget(self.CB_mesh_style)
        self.VL.addWidget(self.IN_main_label)
        self.VL.addWidget(self.IN_labels_font_size)
        self.VL.addWidget(self.IN_tick_labels_font_size)
        self.VL.addWidget(self.list_theme)

        self.VL.addStretch(1)
        self.setLayout(self.VL)

        self.update_widget()

    def update_widget(self):
        """Updates the content of the widget."""
        random_control = self.settings.value("RandomControl", False)
        self.CB_rand_ctrl.setChecked(random_control)

        zero_defl = self.settings.value("ZeroDefl", True)
        self.CB_zero_defl.setChecked(zero_defl)

        display_raw_data = self.settings.value(
            "DisplayRawData",
            False)
        self.CB_raw_data.setChecked(display_raw_data)

        remember = self.settings.value("RememberMeToSave", True)
        self.CB_remember_to_save.setChecked(remember)

        st = self.settings.value("DisplayFitNoise", False)
        self.CB_fit_noise.setChecked(st)

        st = self.settings.value("ClassicMeshgridStyle", True)
        self.CB_mesh_style.setChecked(st)

        self.IN_main_label.changeValue(shared.exp.main_title_label)

        font_size = shared.exp.mpl_labels_font_size
        self.IN_labels_font_size.changeValue(font_size)

        font_size = shared.exp.mpl_tick_labels_font_size
        self.IN_tick_labels_font_size.changeValue(font_size)

        self.list_theme.setCurrentIndex(shared.exp.theme_list.index(shared.exp.current_theme))

    def checkbox_clicked(self, name):
        """Called whenever a checkbox is clicked."""
        if name == "rand_ctrl":
            # Get the value
            value = self.CB_rand_ctrl.isChecked()

            # Reset the list
            shared.exp.parsed_random_curves = []

            # Save the value
            self.settings.setValue("RandomControl", value)

        elif name == "zero_defl":
            # Get the value
            value = self.CB_zero_defl.isChecked()

            # Save the value
            self.settings.setValue("ZeroDefl", value)

            # Update the plots
            widgets_list.widget_data.curve.update_plot()
            widgets_list.widget_results.update_MPL("MPL_canvas2")
            widgets_list.widget_compute.update_MPL("MPL_canvas")
            if widgets_list.widget_curve_mod is not None:
                widgets_list.widget_curve_mod.update_MPL("MPL_canvas")

        elif name == "raw_data":
            # Get the value
            value = self.CB_raw_data.isChecked()

            # Save the value
            self.settings.setValue("DisplayRawData", value)

            widgets_list.widget_data.curve.update_plot()
            widgets_list.widget_results.update_MPL("MPL_canvas2")

        elif name == "fit_noise":
            # Get the value
            value = self.CB_fit_noise.isChecked()

            # Save the value
            self.settings.setValue("DisplayFitNoise", value)

            widgets_list.widget_compute.update_MPL("MPL_canvas")

        elif name == "remember_to_save":
            # Get the value
            value = self.CB_remember_to_save.isChecked()

            # Save the value
            self.settings.setValue("RememberMeToSave", value)

        elif name == "mesh_style":
            # Get the value
            value = self.CB_mesh_style.isChecked()

            # Save the value
            self.settings.setValue("ClassicMeshgridStyle", value)

            # Update the plots
            widgets_list.widget_data.meshgrid.update_plot()
            widgets_list.widget_results.update_MPL("MPL_canvas1")

    def input_updated(self, name):
        """Called when an input is updated."""
        update_plots = False

        if name == "main_label":
            val = self.IN_main_label.get_str_value(False)
            shared.exp.main_title_label = val
            label = shared.exp.main_title_label

            window_title = misc_tools.get_base_window_title()
            if label != "":
                self.parent.setWindowTitle(window_title + " - " + label)
            else:
                self.parent.setWindowTitle(window_title)

        elif name == "tick_labels":
            value = self.IN_tick_labels_font_size.get_str_value(False)
            shared.exp.mpl_tick_labels_font_size = value

            update_plots = True

        elif name == "labels":
            value = self.IN_labels_font_size.get_str_value(False)
            shared.exp.mpl_labels_font_size = value

            update_plots = True

        if update_plots:
            self.update_plots()

    def list_updated(self, name):
        if name == "list_theme":
            # Get current theme
            selected_theme = shared.exp.theme_list[self.list_theme.currentIndex()]
            # Set selected theme
            shared.exp.current_theme = selected_theme
            # Apply selected theme
            theme_handler(selected_theme)
            # Update plots to apply plot theme
            self.update_plots()


    def update_plots(self):
        widgets_list.widget_data.meshgrid.update_plot()
        widgets_list.widget_data.curve.update_plot()
        widgets_list.widget_compute.update_MPL("MPL_meshgrid")
        widgets_list.widget_compute.update_MPL("MPL_canvas")
        widgets_list.widget_results.update_MPL("MPL_canvas1")
        widgets_list.widget_results.update_MPL("MPL_canvas2")
        widgets_list.widget_results_single.update_MPL()
        widgets_list.widget_results_groups.update_MPL()

        if widgets_list.widget_curve_mod is not None:
            widgets_list.widget_curve_mod.update_widget()
        if widgets_list.widget_slices is not None:
            widgets_list.widget_slices.update_widget()
        if widgets_list.widget_profiles is not None:
            widgets_list.widget_profiles.update_widget()
