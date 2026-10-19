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

import numpy
from qtpy import QtCore, QtWidgets
from ..single import robust_max, robust_range
from ..plots.meshgrids.plot_meshgrid import MODULUS_MAPS
from ..tools import apply_to_all
from ..tools import misc_tools
from ..tools import planefit
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

        BT_autoscale = PYAFButton(self, "autoscale", "Autoscale")
        BT_sat_color = PYAFButton(self, "sat_color", "Saturation color")
        BT_neg_color = PYAFButton(self, "neg_color", "Negative color")
        BT_nan_color = PYAFButton(self, "nan_color", "NaN color")

        # Will be hidden for normal meshgrids; and displayed only in vtk
        # See update_widget method
        label = "Transparent Nan values"
        self.CB_nan_transp = PYAFCheckBox(self, "nan_transp", label)

        # Log color scale (elasticity maps)
        self.CB_log = PYAFCheckBox(self, "log_scale", "Log scale")
        self.CB_log.setToolTip(
            "Logarithmic color scale for the elasticity maps (from the "
            "lowest positive value\nwhen the minimum is 0): soft and stiff "
            "regions are both visible (e.g. a cell on glass).")

        BT_apply = PYAFButton(self, "button_apply", "Apply")

        self.gridlayout_box_colors.addWidget(self.LS_colortable, 0, 0)
        self.gridlayout_box_colors.addWidget(self.label_max, 1, 0)
        self.gridlayout_box_colors.addWidget(self.input_max, 2, 0)
        self.gridlayout_box_colors.addWidget(self.input_middle, 3, 0)
        self.gridlayout_box_colors.addWidget(self.input_min, 4, 0)
        self.gridlayout_box_colors.addWidget(self.CB_log, 5, 0)
        self.gridlayout_box_colors.addWidget(BT_autoscale, 6, 0)
        self.gridlayout_box_colors.addWidget(BT_sat_color, 7, 0)
        self.gridlayout_box_colors.addWidget(BT_neg_color, 8, 0)
        self.gridlayout_box_colors.addWidget(BT_nan_color, 9, 0)
        if self.mode == "vtk":
            self.gridlayout_box_colors.addWidget(self.CB_nan_transp, 10, 0)
        self.gridlayout_box_colors.addWidget(BT_apply, 11, 0)

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
            apply_to_all.apply_to_all("display_options", oldid=self.data_id,
                                      from_color_options=True)

            # Update the plots only if asked
            self.update_parent_elements()

            self.update_widget()

        elif what == "log_scale":
            data.stiffness_color_log = self.CB_log.isChecked()
            apply_to_all.apply_to_all("display_options", oldid=self.data_id,
                                      from_color_options=True)
            self.update_parent_elements()

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

        # Get the fields to update: a single input, or all of them when the
        # apply button is clicked. The middle value is only used by the
        # Rainbow 1 colortable.
        inputs = {"max": self.input_max,
                  "middle": self.input_middle,
                  "min": self.input_min}
        if what == "all":
            fields = ["max", "min"]
            if data.colortableid == 1:
                fields.append("middle")
        elif what in ("input_max", "input_middle", "input_min"):
            fields = [what[len("input_"):]]
        else:
            return True

        # Start from the stored values, converted to the displayed units, so
        # that the new and old values are compared in the same units
        values = {}
        for field in inputs:
            stored = getattr(data, "colortable_" + field + "_value")
            factor, offset = color_value_conversion(data, field)
            values[field] = (stored - offset) / factor

        # Only the values modified by the user are stored. The others keep
        # their exact value (the inputs display rounded values: re-reading
        # them would for example lower the max set by autoscale)
        new_fields = []
        for field in fields:
            value = inputs[field].get_str_value(use_replace=True)
            if value != "" and \
                    float(value) != displayed_color_value(data, field):
                values[field] = float(value)
                new_fields.append(field)

        if mt == "events_per_curve":
            if "max" in new_fields and values["max"] < 2:
                display_events_limit_message(
                    "The maximum color value can not be set to a value "
                    "smaller than 2 !")
                self.update_widget()
                return False
            if "min" in new_fields and values["min"] < 1:
                display_events_limit_message(
                    "The minimum color value can not be set to a value "
                    "smaller than 1 !")
                self.update_widget()
                return False

        if ("max" in new_fields or "min" in new_fields) and \
                values["max"] <= values["min"]:
            display_bad_value_message()
            self.update_widget()
            return False

        if "middle" in new_fields and \
                not values["min"] < values["middle"] < values["max"]:
            display_bad_value_message(mode="middle")
            self.update_widget()
            return False

        # Store the values, converted back to the internal units
        for field in new_fields:
            factor, offset = color_value_conversion(data, field)
            setattr(data, "colortable_" + field + "_value",
                    values[field] * factor + offset)

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
                    factor, offset = color_value_conversion(data, "middle")
                    data.colortable_middle_value = \
                        val_max / 2.0 * factor + offset

            data.colortableid = colortable_id

            # Apply to all if needed
            apply_to_all.apply_to_all("display_options", oldid=self.data_id,
                                      from_color_options=True)

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

        if button == "autoscale":
            # Autoscale on the selected file. With apply to all, its scale
            # is then copied to all the files (same scale in the multi
            # meshgrids), except for the elasticity: autoscale of each file
            # (cells, gels and substrates differ by orders of magnitude)
            # Bottom of the scale: 0, unless set below
            data.colortable_min_value = 0
            elasticity = mt in ("stiffness", "stiffness_slice",
                                "stiffness_corr", "stiffness_corr_slice")
            modulus = mt in MODULUS_MAPS
            if elasticity or modulus:
                calculated = "stiffness_calculated" if elasticity \
                    else "visco_calculated"
                files = [dt for dt in shared.exp.list
                         if getattr(dt, calculated)] \
                    if shared.exp.apply_to_all_data else [data]
                for dt in files:
                    autoscale_elasticity(dt, mt)
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
            elif mt == "visco_tan_delta":
                # Diverging scale centered on tan delta = 1
                data.colortable_min_value = 0.0
                data.colortable_max_value = 2.0
                data.colortable_middle_value = 1.0
            elif mt.startswith("visco_") and not modulus:
                high = robust_max(getattr(data, mt))
                if high is not None:
                    data.colortable_max_value = high
                    data.colortable_middle_value = high / 2.0
            elif mt == "fit_r2":
                # From the 5th percentile of the R^2 to 1
                values = numpy.asarray(data.fit_r2, dtype=float)
                values = values[numpy.isfinite(values)]
                if values.size:
                    low = min(float(numpy.percentile(values, 5)), 0.99)
                    data.colortable_min_value = low
                    data.colortable_max_value = 1.0
                    data.colortable_middle_value = (low + 1.0) / 2.0
            elif mt == "fit_poc_fixed":
                data.colortable_min_value = 0.0
                data.colortable_max_value = 1.0
                data.colortable_middle_value = 0.5
            elif mt == "piezo":
                data.colortable_max_value = data.max_piezo
                data.colortable_middle_value = data.max_piezo / 2.0
            elif mt == "topo":
                data.colortable_max_value = data.max_topo
                data.colortable_middle_value = data.max_topo / 2.0

            if mt in ("piezo", "topo") and \
                    planefit.get_planefit_roi(data) is not None:
                # Planefitted heights: their whole range (negative values
                # around the substrate)
                planefit.reset_color_scales(data)

            # Apply to all if needed (each file keeps its own elasticity
            # scale)
            apply_to_all.apply_to_all("display_options", oldid=self.data_id,
                                      from_color_options=not (elasticity or
                                                              modulus))

            self.update_widget()

            if elasticity or modulus:
                # As button_apply, without copying the scale to the files
                if self.input_updated("all"):
                    self.update_parent_elements()
                widgets_list.widget_results.update_GUI("apply_to_all")
            else:
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
                                          oldid=self.data_id,
                                          from_color_options=True)

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

        # Define max label: maximum of the selected file (used by autoscale)
        # and, when the scale is applied to all the files, of all the files
        exp = shared.exp
        if mt in ("stiffness", "stiffness_slice",
                  "stiffness_corr", "stiffness_corr_slice"):
            if data.used_stiffness_model_selected == 3:
                name, unit, factor = "Max. slope", "N/m", 1
            else:
                name, unit, factor = "Max. stiffness", "kPa", 1e-3
            val, val_all = data.max_stiffness, exp.global_max_stiffness
        elif mt == "piezo":
            name, unit, factor = "Max. piezo height", "nm", 1
            val, val_all = data.max_piezo, exp.global_max_piezo
        elif mt == "topo":
            name, unit, factor = "Max. topography height", "nm", 1
            val, val_all = data.max_topo, exp.global_max_topo
        elif mt == "work":
            name, unit, factor = "Max. work", "fJ", 1e15
            val, val_all = data.max_work, exp.global_max_work
        elif mt == "rupture_force":
            name, unit, factor = "Max. rupture force", "pN", 1e12
            val, val_all = data.max_rupture_force1, \
                exp.global_max_rupture_force1
        elif mt == "events_rupture_force":
            name, unit, factor = "Max. rupture force", "pN", 1e12
            val, val_all = data.max_rupture_force2, \
                exp.global_max_rupture_force2
        elif mt == "events_per_curve":
            name, unit, factor = "Max. number of events", "events", 1
            val, val_all = data.max_nbr_events, exp.global_max_nbr_events
        elif mt.startswith("visco_"):
            name, unit, factor = {
                "visco_e0": ("Max. E\u2080", "kPa", 1e-3),
                "visco_storage": ("Max. E\u2032", "kPa", 1e-3),
                "visco_loss": ("Max. E\u2033", "kPa", 1e-3),
                "visco_fluidity": ("Max. fluidity", "", 1),
                "visco_tan_delta": ("Max. tan \u03b4", "", 1)}[mt]
            values = numpy.asarray(getattr(data, mt), dtype=float)
            values = values[numpy.isfinite(values)]
            val = float(values.max()) if values.size else 0.0
            val_all = None
        elif mt == "fit_r2":
            name, unit, factor = "Median R\u00b2", "", 1
            values = numpy.asarray(data.fit_r2, dtype=float)
            values = values[numpy.isfinite(values)]
            val = float(numpy.median(values)) if values.size else 0.0
            val_all = None
        elif mt == "fit_poc_fixed":
            name, unit, factor = "Fixed point of contact", "% of the fits", \
                100
            values = numpy.asarray(data.fit_poc_fixed, dtype=float)
            values = values[numpy.isfinite(values)]
            val = float(values.mean()) if values.size else 0.0
            val_all = None

        def fmt(value):
            if mt == "events_per_curve":
                return str(int(value))
            if mt in ("visco_fluidity", "visco_tan_delta", "fit_r2"):
                return str(round(float(value) * factor, 3))
            if mt == "fit_poc_fixed":
                return str(round(float(value) * factor, 1))
            return str(round(float(value) * factor, 2))

        text = (name + " : " + fmt(val) + " " + unit).rstrip()
        if exp.apply_to_all_data and val_all is not None:
            text += " (all files : " + fmt(val_all) + " " + unit + ")"

        self.label_max.setText(text)

        # Define chosen colortable
        self.LS_colortable.setCurrentIndex(data.colortableid)

        # Define input values
        self.input_max.changeValue(displayed_color_value(data, "max"))
        self.input_middle.changeValue(displayed_color_value(data, "middle"))
        self.input_min.changeValue(displayed_color_value(data, "min"))

        # Check checkbox if needed
        self.CB_apply_to_all.setChecked(shared.exp.apply_to_all_data)

        # Remove or add middle input field (needed only for Rainbow 1 scale)
        if data.colortableid == 1:
            self.gridlayout_box_colors.addWidget(self.input_middle, 3, 0)
        else:
            self.input_middle.setParent(None)

        self.CB_nan_transp.setChecked(data.nan_color_transparent)

        # Log scale only for the elasticity and E0 maps
        self.CB_log.setChecked(getattr(data, "stiffness_color_log", False))
        self.CB_log.setEnabled(mt in ("stiffness", "stiffness_slice",
                                      "stiffness_corr",
                                      "stiffness_corr_slice", "visco_e0",
                                      "visco_storage", "visco_loss"))

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


def color_value_conversion(data, field):
    """Returns the (factor, offset) used to convert a color value displayed in
    the inputs to the value stored in data:
    stored = displayed * factor + offset.

    field is "max", "middle" or "min".
    """
    mt = data.meshgrid_type
    if mt in ("stiffness", "stiffness_slice",
              "stiffness_corr", "stiffness_corr_slice"):
        if data.used_stiffness_model_selected == 3:
            # Slope in N/m
            return 1.0, 0.0
        # kPa to Pa
        return 1e3, 0.0
    elif mt in ("visco_e0", "visco_storage", "visco_loss"):
        # kPa to Pa
        return 1e3, 0.0
    elif mt == "work":
        # fJ to J
        return 1e-15, 0.0
    elif mt in ("rupture_force", "events_rupture_force"):
        # pN to N
        return 1e-12, 0.0
    elif mt == "events_per_curve":
        # The events colortable is shifted by one
        if field == "min":
            return 1.0, -1.0
        return 1.0, 1.0
    return 1.0, 0.0


def displayed_color_value(data, field):
    """Returns the color value of field ("max", "middle" or "min") as
    displayed in the inputs (displayed units, rounded)."""
    stored = getattr(data, "colortable_" + field + "_value")
    factor, offset = color_value_conversion(data, field)
    value = float((stored - offset) / factor)
    mt = data.meshgrid_type
    if mt in ("stiffness", "stiffness_slice",
              "stiffness_corr", "stiffness_corr_slice") and \
            data.used_stiffness_model_selected != 3:
        return round(value, 2)
    elif mt == "work":
        # Work values can be small, keep significant digits
        return float(f"{value:.4g}")
    elif mt in ("rupture_force", "events_rupture_force", "visco_e0",
                "visco_storage", "visco_loss"):
        return round(value, 2)
    elif mt in ("visco_fluidity", "visco_tan_delta", "fit_r2",
                "fit_poc_fixed"):
        return round(value, 3)
    return float(f"{value:.7g}")


def display_events_limit_message(text):
    """Tells the user the events colortable limits are not respected."""
    msg = QtWidgets.QMessageBox()
    msg.setText("Error")
    msg.setInformativeText(text)
    msg.setIcon(QtWidgets.QMessageBox.Critical)
    msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
    msg.exec_()


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


def autoscale_elasticity(data, meshgrid_type):
    """Automatic color scale of an elasticity or viscoelastic modulus map
    of a file: robust limits of the values (of the sample for the
    elasticity: substrate excluded; values outside drawn with the first and
    saturation colors), see robust_range."""
    if meshgrid_type in MODULUS_MAPS:
        values = getattr(data, meshgrid_type)
    else:
        values = stiffness_values(data, meshgrid_type)
    limits = robust_range(values, data.stiffness_color_log)
    if limits is not None:
        opts = data.get_color_opts(meshgrid_type)
        opts[1], opts[2], opts[3] = limits[0], limits[1], sum(limits) / 2.0


def stiffness_values(data, meshgrid_type):
    """Values of the sample in the displayed elasticity map of the file
    (corrected for the B.E.C. maps), for the automatic color scale."""
    corrected = meshgrid_type.startswith("stiffness_corr") and \
        data.stiffness_corrected
    return data.sample_elasticity_values(corrected=corrected)
