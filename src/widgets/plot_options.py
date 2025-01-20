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

"""Options for the meshgrids."""

from PyQt5 import QtCore, QtWidgets
from ..tools import apply_to_all
from ..tools import misc_tools
from ..tools.colortables import ColorTables
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFCheckBox
from .. import widgets_list
from .. import shared
from ..tools.PYAFWidget import PYAFWidget


class PlotOptionsWidget(PYAFWidget):
    """Widget used to change the colorscales of the meshgrid."""

    def __init__(self, parent, mode=None):
        super().__init__(parent,
                                                "widget_meshgrid_options")

        # Can be None or "vtk"
        self.mode = mode

        # Load the settings
        self.settings = QtCore.QSettings()

        self.data_id = None

        self.VL = QtWidgets.QVBoxLayout()

        # Color scales --------------------------------------------------
        self.box_colors = QtWidgets.QGroupBox("Color scale")

        self.gridlayout_box_colors = QtWidgets.QGridLayout()

        self.label_max = QtWidgets.QLabel()

        self.LS_colortable = PYAFComboBox(self, "list_colortable_chooser")
        self.LS_colortable.setSizeAdjustPolicy(
            QtWidgets.QComboBox.AdjustToContents)
        self.colortables = ColorTables()
        for name in self.colortables.colortables_list:
            self.LS_colortable.addItem(name)

        self.input_max = PYAFInput(self, "input_max", "Max", 150)
        self.input_max.input.setValidator(misc_tools.validator("UF"))
        self.input_middle = PYAFInput(self, "input_middle", "Middle", 150)
        self.input_middle.input.setValidator(misc_tools.validator("UF"))
        self.input_min = PYAFInput(self, "input_min", "Min", 150)
        self.input_min.input.setValidator(misc_tools.validator("UF"))

        BT_autoscale = PYAFButton(self, "autoscale", "Auto Scale")
        BT_sat_color = PYAFButton(self, "sat_color", "Saturation color")
        BT_neg_color = PYAFButton(self, "neg_color", "Negative color")
        BT_nan_color = PYAFButton(self, "nan_color", "Nan color")

        # Will be hidden for normal meshgrids; and displayed only in vtk
        # See update_widget method
        label = "Transparent Nan values"
        self.CB_nan_transp = PYAFCheckBox(self, "nan_transp", label)

        BT_apply = PYAFButton(self, "button_apply", "Apply")

        self.gridlayout_box_colors.addWidget(self.LS_colortable, 0, 0)
        self.gridlayout_box_colors.addWidget(self.label_max, 1, 0)
        self.gridlayout_box_colors.addWidget(self.input_max, 2, 0)
        self.gridlayout_box_colors.addWidget(self.input_middle, 3, 0)
        self.gridlayout_box_colors.addWidget(self.input_min, 4, 0)
        self.gridlayout_box_colors.addWidget(BT_autoscale, 5, 0)
        self.gridlayout_box_colors.addWidget(BT_sat_color, 6, 0)
        self.gridlayout_box_colors.addWidget(BT_neg_color, 7, 0)
        self.gridlayout_box_colors.addWidget(BT_nan_color, 8, 0)
        if self.mode == "vtk":
            self.gridlayout_box_colors.addWidget(self.CB_nan_transp, 9, 0)
        self.gridlayout_box_colors.addWidget(BT_apply, 10, 0)

        self.box_colors.setLayout(self.gridlayout_box_colors)

        # Checkbox apply options on all ---------------------------------------

        name = "checkbox_apply_to_all"
        self.CB_apply_to_all = PYAFCheckBox(self, name, "Apply to all")

        # Set up layout -------------------------------------------------------

        self.VL.addWidget(self.box_colors)
        self.VL.addWidget(self.CB_apply_to_all)
        self.VL.addStretch(1)
        self.setLayout(self.VL)

        self.update_widget()

    def checkbox_clicked(self, what):
        """Method called whenever a checkbox is clicked."""
        if self.data_id is None:
            data = shared.exp.current_data
        else:
            data = shared.exp.list[self.data_id]

        if what == "checkbox_apply_to_all":
            shared.exp.apply_to_all_data = self.CB_apply_to_all.isChecked()

            # Apply on all
            apply_to_all.apply_to_all("display_options", oldid=self.data_id)

            # Update the plots only if asked
            self.update_parent_elements()

            self.update_widget()

        elif what == "nan_transp":
            data.nan_color_transparent = self.CB_nan_transp.isChecked()
            widgets_list.widget_vtk.update_AFM_colors()

    def input_updated(self, what):
        """Method called when an input is updated."""
        if self.data_id is None:
            data = shared.exp.current_data
        else:
            data = shared.exp.list[self.data_id]

        mt = data.meshgrid_type

        # Set focus on widget to validate the inputs
        self.setFocus()

        # Get and store the values
        if what == "input_max" or what == "all":
            value = self.input_max.get_str_value()
            if value != "":
                if float(value) <= data.colortable_min_value:
                    display_bad_value_message()
                    self.update_widget()
                    return False

                if mt == "stiffness" or mt == "stiffness_slice" or \
                        mt == "stiffness_corr" or mt == "stiffness_corr_slice":
                    if data.used_stiffness_model_selected == 3:
                        data.colortable_max_value = float(value)
                    else:
                        data.colortable_max_value = float(value) * 1e3
                elif mt == "work":
                    data.colortable_max_value = float(value) * 1e-9
                elif mt == "rupture_force":
                    data.colortable_max_value = float(value) * 1e-12
                elif mt == "events_rupture_force":
                    data.colortable_max_value = float(value) * 1e-12
                elif mt == "events_per_curve":
                    if float(value) < 2:
                        text = "The maximum color value can not be set to a \
                                value smaller than 2 !"

                        # Create a message box
                        msg = QtWidgets.QMessageBox()
                        msg.setText("Error")
                        msg.setInformativeText(text)
                        msg.setIcon(QtWidgets.QMessageBox.Critical)
                        msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
                        msg.exec_()

                        val = data.colortable_max_value - 1
                        self.input_max.changeValue(val)
                        return False
                    else:
                        data.colortable_max_value = float(value) + 1
                else:
                    data.colortable_max_value = float(value)

        elif what == "input_middle" or what == "all":
            value = self.input_middle.get_str_value()
            if value != "":
                if float(value) <= data.colortable_min_value or \
                        float(value) >= data.colortable_max_value:
                    display_bad_value_message(mode="middle")
                    self.update_widget()
                    return False

                if mt == "stiffness" or mt == "stiffness_slice" or \
                        mt == "stiffness_corr" or mt == "stiffness_corr_slice":
                    if data.used_stiffness_model_selected == 3:
                        data.colortable_middle_value = float(value)
                    else:
                        data.colortable_middle_value = float(value) * 1e3
                elif mt == "work":
                    data.colortable_middle_value = float(value) * 1e-9
                elif mt == "rupture_force":
                    data.colortable_middle_value = float(value) * 1e-12
                elif mt == "events_rupture_force":
                    data.colortable_middle_value = float(value) * 1e-12
                elif mt == "events_per_curve":
                    data.colortable_middle_value = float(value) + 1
                else:
                    data.colortable_middle_value = float(value)

        elif what == "input_min" or what == "all":
            value = self.input_min.get_str_value()
            if value != "":
                if float(value) >= data.colortable_max_value:
                    display_bad_value_message()
                    self.update_widget()
                    return False

                if mt == "stiffness" or mt == "stiffness_slice" or \
                        mt == "stiffness_corr" or mt == "stiffness_corr_slice":
                    if data.used_stiffness_model_selected == 3:
                        data.colortable_min_value = float(value)
                    else:
                        data.colortable_min_value = float(value) * 1e3
                elif mt == "work":
                    data.colortable_min_value = float(value) * 1e-9
                elif mt == "rupture_force":
                    data.colortable_min_value = float(value) * 1e-12
                elif mt == "events_rupture_force":
                    data.colortable_min_value = float(value) * 1e-12
                elif mt == "events_per_curve":
                    if float(value) < 1:
                        text = "The minimum color value can not be set to a \
                                value smaller than 1 !"

                        # Create a message box
                        msg = QtWidgets.QMessageBox()
                        msg.setText("Error")
                        msg.setInformativeText(text)
                        msg.setIcon(QtWidgets.QMessageBox.Critical)
                        msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
                        msg.exec_()

                        val = data.colortable_min_value + 1
                        self.input_min.changeValue(val)
                        return False
                    else:
                        data.colortable_min_value = float(value) - 1
                else:
                    data.colortable_min_value = float(value)

        return True

    def list_updated(self, what):
        """Called when a list is updated."""
        if self.data_id is None:
            data = shared.exp.current_data
        else:
            data = shared.exp.list[self.data_id]

        if what == "list_colortable_chooser":
            val_max = self.input_max.get_float_value()
            val_middle = self.input_middle.get_float_value()
            val_min = self.input_min.get_float_value()

            colortable_id = self.LS_colortable.currentIndex()
            if colortable_id == 1:
                # Check for the middle value, if its misconfigured
                # replace it by the half of the max value
                if val_middle <= val_min or val_middle >= val_max:
                    data.colortable_middle_value = val_max / 2.0

            data.colortableid = colortable_id

            # Apply to all if needed
            apply_to_all.apply_to_all("display_options", oldid=self.data_id)

            # Update
            self.update_widget()
            self.update_parent_elements()

    def button_clicked(self, button):
        """Method called when a button is clicked."""
        if self.data_id is None:
            data = shared.exp.current_data
        else:
            data = shared.exp.list[self.data_id]

        mt = data.meshgrid_type
        exp = shared.exp
        dt = data

        if button == "autoscale":
            if shared.exp.apply_to_all_data:
                if mt == "stiffness" or mt == "stiffness_slice" or \
                        mt == "stiffness_corr" or mt == "stiffness_corr_slice":
                    dt.colortable_max_value = exp.global_max_stiffness
                    dt.colortable_middle_value = exp.global_max_stiffness / 2.0
                elif mt == "work":
                    dt.colortable_max_value = exp.global_max_work
                    dt.colortable_middle_value = exp.global_max_work / 2.0
                elif mt == "rupture_force":
                    dt.colortable_max_value = exp.global_max_rupture_force1
                    dt.colortable_middle_value = \
                        exp.global_max_rupture_force1 / 2.0
                elif mt == "events_rupture_force":
                    dt.colortable_max_value = exp.global_max_rupture_force2
                    dt.colortable_middle_value = \
                        exp.global_max_rupture_force2 / 2.0
                elif mt == "events_per_curve":
                    dt.colortable_max_value = exp.global_max_nbr_events
                    dt.colortable_middle_value = exp.global_max_nbr_events / \
                        2.0
                elif mt == "piezo":
                    dt.colortable_max_value = exp.global_max_piezo
                    dt.colortable_middle_value = exp.global_max_piezo / 2.0
                elif mt == "topo":
                    dt.colortable_max_value = exp.global_max_topo
                    dt.colortable_middle_value = exp.global_max_topo / 2.0

            else:
                if mt == "stiffness" or mt == "stiffness_slice" or \
                        mt == "stiffness_corr" or mt == "stiffness_corr_slice":
                    data.colortable_max_value = data.max_stiffness
                    data.colortable_middle_value = data.max_stiffness / 2.0
                elif mt == "work":
                    data.colortable_max_value = data.max_work
                    data.colortable_middle_value = data.max_work / 2.0
                elif mt == "rupture_force":
                    data.colortable_max_value = data.max_rupture_force1
                    data.colortable_middle_value = data.max_rupture_force1 / \
                        2.0
                elif mt == "events_rupture_force":
                    data.colortable_max_value = data.max_rupture_force2
                    data.colortable_middle_value = data.max_rupture_force2 / \
                        2.0
                elif mt == "events_per_curve":
                    data.colortable_max_value = data.max_nbr_events
                    data.colortable_middle_value = data.max_nbr_events / 2.0
                elif mt == "piezo":
                    data.colortable_max_value = data.max_piezo
                    data.colortable_middle_value = data.max_piezo / 2.0
                elif mt == "topo":
                    data.colortable_max_value = data.max_topo
                    data.colortable_middle_value = data.max_topo / 2.0

            data.colortable_min_value = 0

            # Apply to all if needed
            apply_to_all.apply_to_all("display_options", oldid=self.data_id)

            self.update_widget()

            self.button_clicked("button_apply")

        if button == "sat_color":
            color = misc_tools.ask_user_for_color(data.color_saturation)
            if color is not False:
                data.color_saturation = color

        elif button == "neg_color":
            color = misc_tools.ask_user_for_color(data.color_negative)
            if color is not False:
                data.color_negative = color

        elif button == "nan_color":
            color = misc_tools.ask_user_for_color(data.color_nan)
            if color is not False:
                data.color_nan = color

        elif button == "button_apply":
            # Get the values
            result = self.input_updated("all")

            if result:
                # Apply to all if needed
                apply_to_all.apply_to_all("display_options",
                                          oldid=self.data_id)

                # Update the plots
                self.update_parent_elements()

            # Update parent checkbox
            widgets_list.widget_results.update_GUI("apply_to_all")

    def update_widget(self, data_id=None):
        """This method lets you update the content of the widget, for example when
        you change the dataset or when some modifications are done.
        """
        # Once the value is set change it only if updated (for vtk widget)
        if data_id is not None:
            self.data_id = data_id

        if self.data_id is None:
            data = shared.exp.current_data
        else:
            data = shared.exp.list[self.data_id]

        mt = data.meshgrid_type

        # Define max label
        if shared.exp.apply_to_all_data:
            if mt == "piezo":
                text = "Max. piezo height (All) : " + \
                    str(shared.exp.global_max_piezo) + " nm"
            if mt == "topo":
                text = "Max. topography height (All) : " + \
                    str(shared.exp.global_max_topo) + " nm"
            if mt == "stiffness" or mt == "stiffness_slice" or \
                    mt == "stiffness_corr" or mt == "stiffness_corr_slice":
                if data.used_stiffness_model_selected == 3:
                    text = "Max. slope (All) : " + \
                        str(shared.exp.global_max_stiffness) + " N/m"
                else:
                    val = round(shared.exp.global_max_stiffness * 1e-3, 2)
                    text = "Max. stiffness (All) : " + str(val) + " kPa"
            if mt == "work":
                val = round(shared.exp.global_max_work * 1e15, 2)
                text = "Max. work (All) : " + str(val) + " fJ"
            if mt == "rupture_force":
                val = round(shared.exp.global_max_rupture_force1 * 1e12, 2)
                text = "Max. rupture force (All) : " + str(val) + " pN"
            if mt == "events_rupture_force":
                val = round(shared.exp.global_max_rupture_force2 * 1e12, 2)
                text = "Max. rupture force (All) : " + str(val) + " pN"
            if mt == "events_per_curve":
                val = int(shared.exp.global_max_nbr_events - 1)
                text = "Max. number of events (All) : " + str(val) + " events"
        else:
            if mt == "piezo":
                text = "Max. piezo height : " + str(data.max_piezo) + " nm"
            if mt == "topo":
                text = "Max. topography height : " + str(data.max_topo) + " nm"
            if mt == "stiffness" or mt == "stiffness_slice" or \
                    mt == "stiffness_corr" or mt == "stiffness_corr_slice":
                if data.used_stiffness_model_selected == 3:
                    text = "Max. slope : " + \
                        str(data.max_stiffness) + " N/m"
                else:
                    text = "Max. stiffness : " + \
                        str(round(data.max_stiffness * 1e-3, 2)) + \
                        " kPa"
            if mt == "work":
                text = "Max. work : " + \
                    str(round(data.max_work * 1e15, 2)) + " fJ"
            if mt == "rupture_force":
                text = "Max. rupture force : " + \
                    str(round(data.max_rupture_force1 * 1e12)) + " pN"
            if mt == "events_rupture_force":
                text = "Max. rupture force : " + \
                    str(round(data.max_rupture_force2 * 1e12)) + " pN"
            if mt == "events_per_curve":
                text = "Max. number of events : " + \
                    str(int(data.max_nbr_events + 1)) + " events"

        self.label_max.setText(text)

        # Define chosen colortable
        self.LS_colortable.setCurrentIndex(data.colortableid)

        # Define input values
        if mt == "stiffness" or mt == "stiffness_slice" or \
                mt == "stiffness_corr" or mt == "stiffness_corr_slice":
            if data.used_stiffness_model_selected == 3:
                val = data.colortable_max_value
            else:
                val = round(data.colortable_max_value * 1e-3, 2)
            self.input_max.changeValue(val)
        elif mt == "work":
            val = round(data.colortable_max_value * 1e9, 2)
            self.input_max.changeValue(val)
        elif mt == "rupture_force":
            val = round(data.colortable_max_value * 1e12, 2)
            self.input_max.changeValue(val)
        elif mt == "events_rupture_force":
            val = round(data.colortable_max_value * 1e12, 2)
            self.input_max.changeValue(val)
        elif mt == "events_per_curve":
            self.input_max.changeValue(data.colortable_max_value - 1)
        else:
            self.input_max.changeValue(data.colortable_max_value)

        if mt == "stiffness" or mt == "stiffness_slice" or \
                mt == "stiffness_corr" or mt == "stiffness_corr_slice":
            if data.used_stiffness_model_selected == 3:
                val = data.colortable_middle_value
            else:
                val = round(data.colortable_middle_value * 1e-3, 2)
            self.input_middle.changeValue(val)
        elif mt == "work":
            val = round(data.colortable_middle_value * 1e9, 2)
            self.input_middle.changeValue(val)
        elif mt == "rupture_force":
            val = round(data.colortable_middle_value * 1e12, 2)
            self.input_middle.changeValue(val)
        elif mt == "events_rupture_force":
            val = round(data.colortable_middle_value * 1e12, 2)
            self.input_middle.changeValue(val)
        elif mt == "events_per_curve":
            self.input_middle.changeValue(data.colortable_middle_value - 1)
        else:
            self.input_middle.changeValue(data.colortable_middle_value)

        if mt == "stiffness" or mt == "stiffness_slice" or \
                mt == "stiffness_corr" or mt == "stiffness_corr_slice":
            if data.used_stiffness_model_selected == 3:
                val = data.colortable_min_value
            else:
                val = round(data.colortable_min_value * 1e-3, 2)
            self.input_min.changeValue(val)
        elif mt == "work":
            val = round(data.colortable_min_value * 1e9, 2)
            self.input_min.changeValue(val)
        elif mt == "rupture_force":
            val = round(data.colortable_min_value * 1e12, 2)
            self.input_min.changeValue(val)
        elif mt == "events_rupture_force":
            val = round(data.colortable_min_value * 1e12, 2)
            self.input_min.changeValue(val)
        elif mt == "events_per_curve":
            self.input_min.changeValue(data.colortable_min_value + 1)
        else:
            self.input_min.changeValue(data.colortable_min_value)

        # Check checkbox if needed
        self.CB_apply_to_all.setChecked(shared.exp.apply_to_all_data)

        # Remove or add middle input field (needed only for Rainbow 1 scale)
        if data.colortableid == 1:
            self.gridlayout_box_colors.addWidget(self.input_middle, 3, 0)
        else:
            self.input_middle.setParent(None)

        self.CB_nan_transp.setChecked(data.nan_color_transparent)

    def update_parent_elements(self):
        """This method is called to update some elements from the parents widget.

        Recalculate only the color for the 3D (which is lighter than
        reinitializing the whole 3D scenery)

        Do not refresh the matplotlib canvas if the user has uncheked the
        option in the 3D VTK menu.
        """
        ref = self.settings.value("refreshBackgroundInOpenGL", True)

        widg_slices = widgets_list.widget_slices
        widg_mesh = widgets_list.widget_multimeshgrids
        widg_opengl = widgets_list.widget_vtk

        # Update the plots in the main GUI
        if widg_opengl is None or (widg_opengl is not None and ref):
            widgets_list.widget_data.meshgrid.update_plot()
            widgets_list.widget_results.update_MPL("MPL_canvas1")
        # Update the slices widget if needed
        if (widg_slices is not None and ref and widg_opengl is not None) or \
                (widg_slices is not None and widg_opengl is None):
            widgets_list.widget_slices.update_MPL()
        # Update the 3D colors if needed
        if widgets_list.widget_vtk is not None:
            widgets_list.widget_vtk.update_AFM_colors()
        # Update the multiple meshgrids widget if needed
        if (widg_mesh is not None and ref and widg_opengl is not None) or \
                (widg_mesh is not None and widg_opengl is None):
            widgets_list.widget_multimeshgrids.update_widget()


def display_bad_value_message(mode=None):
    """Tells the user the inputed values are not well defined."""
    if mode == "middle":
        text = "The middle value has to be between the min and max value."
    else:
        text = "The max value can not be lower than the min value."
    msg = QtWidgets.QMessageBox()
    msg.setText("Error : bad values")
    msg.setInformativeText(text)
    msg.setIcon(QtWidgets.QMessageBox.Critical)
    msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
    msg.exec_()
