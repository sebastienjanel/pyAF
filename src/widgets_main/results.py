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

"""Tab with the meshgrids and the results of the fits."""

import os
import math
from .. import widgets_list
from .. import shared
import logging
from .. import consts
from PyQt5 import QtGui, QtCore, QtWidgets
from ..tools import math_tools
from ..tools import apply_to_all
from ..tools import misc_tools
from ..utils import meshgrid_menu
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFButtonGroup
from ..widgets_main.subwidgets.curve_selector import ResultsSelectorWidget
from ..widgets.multi_meshgrids import MultiMeshgridsWidget
from ..widgets.slices import SlicesWidget
from ..widgets.profiles import ProfileWidget
from ..widgets.details import DetailsWidget
from ..widgets.roi_manager import RoiManagerWidget
from ..widgets.events_filters import ChangeFilterWidget
from ..widgets.used_parameters import UsedParametersWidget
from ..plots.plot_main import MainPlot
from ..plots.PYAFPlot import PYAFPlot
from ..tools.PYAFWidget import PYAFWidget
from ..tools import events_refresher

if consts.ALLOW_VTK:
    from ..widgets_vtk.entry import VTKEntryWidget


class ResultsWidget(PYAFWidget):
    """Results widget, displaying the maps and the force curves.

    It's the third tab of PYAF. It is called from the QMainWindow in main.py.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_results")

        self.logger = logging.getLogger()

        # Fixed height for the 2 top boxes, to keep them aligned
        self.options_boxes_height = 160

        self.cursor_mode = "select"
        self.parent = parent

        # Load the settings
        self.settings = QtCore.QSettings()

        self.current_profile_erase = None

        self.canvas_resolution = self.parent.canvas_resolution
        self.meshgrid_canvas_size = [300, 340]
        self.curves_canvas_size = [400, 340]
        # Canvas for the meshgrid
        self.canvas_meshgrid = QtWidgets.QWidget()
        self.canvas_meshgrid.setFixedSize(self.meshgrid_canvas_size[0],
                                          self.meshgrid_canvas_size[1])
        # Canvas for the curves
        self.canvas_curves = QtWidgets.QWidget()
        self.canvas_curves.setFixedSize(self.curves_canvas_size[0],
                                        self.curves_canvas_size[1])

        sizes = [self.meshgrid_canvas_size[0] / self.canvas_resolution,
                 self.meshgrid_canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]
        self.MPL_canvas1 = PYAFPlot(
            self, "meshgrid", self.canvas_meshgrid, sizes)
        self.MPL_canvas1.mpl_connect("button_press_event",
                                     self.meshgrid_canvas_press)
        self.MPL_canvas1.mpl_connect("button_release_event",
                                     self.canvas_release)
        self.MPL_canvas1.mpl_connect("motion_notify_event", self.canvas_motion)
        self.MPL_canvas1.tab_id = 2
        sizes = [self.curves_canvas_size[0] / self.canvas_resolution,
                 self.curves_canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]
        self.MPL_canvas2 = PYAFPlot(
            self, "curve_results", self.canvas_curves, sizes)
        self.MPL_canvas2.mpl_connect(
            "button_press_event", self.curves_canvas_press)

        widgets_list.widget_progressbar.update()

        # --- Main layout -----------------------------------------------------

        self.VL = QtWidgets.QVBoxLayout()

        # --- Experiment ------------------------------------------------------

        # List (has to be in a layout to adjust to content)
        self.HL_choser = QtWidgets.QHBoxLayout()
        self.list_exp = PYAFComboBox(self, "list1")
        self.list_exp.setFixedSize(350, 26)
        for fv_data in shared.exp.list:
            self.list_exp.addItem(os.path.basename(str(fv_data.filename)))
        self.HL_choser.addWidget(self.list_exp)

        widgets_list.widget_progressbar.update()

        self.CB_app_to_all = PYAFCheckBox(self, "apply_to_all", "Apply to all")

        name = "widget_curve_selector_results"
        self.curve_selector = ResultsSelectorWidget(self, name)
        # self.curve_selector is added to the layout in update_widget()

        widgets_list.widget_progressbar.update()

        # --- Meshgrid --------------------------------------------------------

        self.box_meshgrid = QtWidgets.QGroupBox()
        self.GL_meshgrid = QtWidgets.QGridLayout()
        # Options (meshgrid)
        self.box_meshgrid_options = QtWidgets.QGroupBox()
        self.box_meshgrid_options.setFixedHeight(self.options_boxes_height)

        self.GRP_mesh_type = PYAFButtonGroup(self, "BTG_meshgrid_type")
        self.RBT_mesh_piezo = QtWidgets.QRadioButton("Piezo height")
        self.RBT_mesh_topography = QtWidgets.QRadioButton("Topography")
        self.RBT_mesh_stiffness = QtWidgets.QRadioButton()
        self.RBT_mesh_stiffness_slice = QtWidgets.QRadioButton()
        self.RBT_mesh_stiffness_corr = QtWidgets.QRadioButton(
            "Elasticity (BEC)")
        self.RBT_mesh_stiffness_slice_corr = QtWidgets.QRadioButton(
            "Elasticity (Slice) (BEC)")
        self.RBT_mesh_work = QtWidgets.QRadioButton("Detachment work")
        self.RBT_mesh_rupture_force1 = QtWidgets.QRadioButton("Detachment force")
        self.RBT_mesh_events_per_curve = QtWidgets.QRadioButton("Events per curve")
        self.RBT_mesh_rupture_force2 = QtWidgets.QRadioButton(
            "Event max. force")
        self.GRP_mesh_type.addButton(self.RBT_mesh_piezo, 0)
        self.GRP_mesh_type.addButton(self.RBT_mesh_topography, 1)
        self.GRP_mesh_type.addButton(self.RBT_mesh_stiffness, 2)
        self.GRP_mesh_type.addButton(self.RBT_mesh_stiffness_slice, 3)
        self.GRP_mesh_type.addButton(self.RBT_mesh_stiffness_corr, 4)
        self.GRP_mesh_type.addButton(self.RBT_mesh_stiffness_slice_corr, 5)
        self.GRP_mesh_type.addButton(self.RBT_mesh_work, 6)
        self.GRP_mesh_type.addButton(self.RBT_mesh_rupture_force1, 7)
        self.GRP_mesh_type.addButton(self.RBT_mesh_events_per_curve, 8)
        self.GRP_mesh_type.addButton(self.RBT_mesh_rupture_force2, 9)

        self.GL_mesh_options = QtWidgets.QGridLayout()
        self.VL_meshgrid_options = QtWidgets.QVBoxLayout()
        self.VL_meshgrid_options.setContentsMargins(0, 0, 0, 0)

        self.GL_mesh_options.addWidget(self.RBT_mesh_piezo, 0, 0)
        self.GL_mesh_options.addWidget(self.RBT_mesh_topography, 1, 0)

        self.GL_mesh_options.addWidget(self.RBT_mesh_stiffness, 0, 1)
        self.GL_mesh_options.addWidget(self.RBT_mesh_stiffness_slice, 1, 1)
        self.GL_mesh_options.addWidget(self.RBT_mesh_stiffness_corr, 2, 1)
        self.GL_mesh_options.addWidget(
            self.RBT_mesh_stiffness_slice_corr, 3, 1)

        self.GL_mesh_options.addWidget(self.RBT_mesh_work, 0, 2)
        self.GL_mesh_options.addWidget(self.RBT_mesh_rupture_force1, 1, 2)

        self.GL_mesh_options.addWidget(self.RBT_mesh_events_per_curve, 2, 2)
        self.GL_mesh_options.addWidget(self.RBT_mesh_rupture_force2, 3, 2)

        self.VL_meshgrid_options.addLayout(self.GL_mesh_options)
        self.VL_meshgrid_options.addStretch(1)

        self.HL_meshgrid_options = QtWidgets.QHBoxLayout()
        self.HL_meshgrid_options.addLayout(self.VL_meshgrid_options)

        self.box_meshgrid_options.setLayout(self.HL_meshgrid_options)
        widgets_list.widget_progressbar.update()

        # Some tools
        self.BT_roi_manager = PYAFButton(
            self, "button_roi_manager", "ROI Manager")
        self.BT_profiles = PYAFButton(self, "button_profiles", "Profiles")
        self.BT_slices = PYAFButton(self, "button_slices", "Slices")
        self.BT_multi_meshgrids = PYAFButton(
            self, "button_multi_meshgrids", "Multi Meshgrids")
        self.BT_display_in_3D = PYAFButton(
            self, "button_display_in_3D", "Display in 3D")
        self.BT_used_parameters = PYAFButton(
            self, "used_parameters", "Used parameters")

        widgets_list.widget_progressbar.update()

        # --- Curves ----------------------------------------------------------

        # Options curves
        self.box_curves_options = QtWidgets.QGroupBox()
        self.box_curves_options.setFixedHeight(self.options_boxes_height)

        self.BT_group_type_of_curve = PYAFButtonGroup(
            self, "button_group_type_of_curve")
        self.RB_defl_ext = QtWidgets.QRadioButton("Deflection - Extension curve")
        self.RB_force_ext = QtWidgets.QRadioButton("Force - Extension curve")
        self.RB_force = QtWidgets.QRadioButton("Force - Distance curve")
        self.RB_force_events = QtWidgets.QRadioButton("Force curve (events)")
        self.RB_indentation_curve = QtWidgets.QRadioButton("Indentation curve")
        self.RB_residuals = QtWidgets.QRadioButton("Residuals plot")
        self.BT_group_type_of_curve.addButton(self.RB_defl_ext, 0)
        self.BT_group_type_of_curve.addButton(self.RB_force_ext, 1)
        self.BT_group_type_of_curve.addButton(self.RB_force, 2)
        self.BT_group_type_of_curve.addButton(self.RB_force_events, 3)
        self.BT_group_type_of_curve.addButton(self.RB_indentation_curve, 4)
        self.BT_group_type_of_curve.addButton(self.RB_residuals, 5)

        self.BT_group_trace_or_retrace = PYAFButtonGroup(
            self, "button_group_trace_or_retrace")
        self.radiobutton_trace = QtWidgets.QRadioButton("Approach curve")
        self.radiobutton_retrace = QtWidgets.QRadioButton("Retraction curve")
        self.radiobutton_both = QtWidgets.QRadioButton("Both")
        self.BT_group_trace_or_retrace.addButton(self.radiobutton_trace, 0)
        self.BT_group_trace_or_retrace.addButton(self.radiobutton_retrace, 1)
        self.BT_group_trace_or_retrace.addButton(self.radiobutton_both, 2)
        self.CB_poc = PYAFCheckBox(self, "checkbox_poc", "Point of contact")
        self.CB_fit_poc = PYAFCheckBox(
            self, "checkbox_fit_poc", "POC fit")
        self.CB_segments = PYAFCheckBox(
            self, "checkbox_segments", "Segments")
        self.CB_display_fits = PYAFCheckBox(
            self, "display_fits_stiffness", "Elasticity fit")

        self.CB_joc = PYAFCheckBox(self, "checkbox_joc", "Jump of contact")
        self.CB_fit_joc = PYAFCheckBox(
            self, "checkbox_fit_joc", "JOC fit")
        self.CB_surface = PYAFCheckBox(
            self, "checkbox_surface", "Detach. work")
        self.CB_force = PYAFCheckBox(self, "force", "Detach. force")

        # Events
        self.BT_group_events = PYAFButtonGroup(self, "button_group_events")
        self.RBT_events_none = QtWidgets.QRadioButton("None")
        self.RBT_events_annotations = QtWidgets.QRadioButton("Annotations")
        self.BT_group_events.addButton(self.RBT_events_none, 0)
        self.BT_group_events.addButton(self.RBT_events_annotations, 1)
        self.CB_events_filter_dist = PYAFCheckBox(
            self, "checkbox_events_filter_dist", "Distance filter")
        self.CB_events_display_joc = PYAFCheckBox(
            self, "checkbox_events_display_joc", "Display Joc")
        self.CB_events_display_fit_joc = PYAFCheckBox(
            self, "checkbox_events_display_fit_joc", "Display Joc fit")

        self.GL_curves_options = QtWidgets.QGridLayout()
        self.VL_curves_options = QtWidgets.QVBoxLayout()
        widgets_list.widget_progressbar.update()

        # Col 1
        self.GL_curves_options.addWidget(self.RB_defl_ext, 0, 0)
        self.GL_curves_options.addWidget(self.RB_force_ext, 1, 0)
        self.GL_curves_options.addWidget(self.RB_force, 2, 0)
        self.GL_curves_options.addWidget(self.RB_force_events, 3, 0)
        self.GL_curves_options.addWidget(self.RB_indentation_curve, 4, 0)
        self.GL_curves_options.addWidget(self.RB_residuals, 5, 0)
        # Col 2
        self.GL_curves_options.addWidget(self.radiobutton_both, 0, 1)
        self.GL_curves_options.addWidget(self.radiobutton_trace, 1, 1)
        self.GL_curves_options.addWidget(self.radiobutton_retrace, 2, 1)
        # Col 3
        self.GL_curves_options.addWidget(self.CB_poc, 0, 2)
        self.GL_curves_options.addWidget(self.CB_fit_poc, 1, 2)
        self.GL_curves_options.addWidget(self.CB_segments, 2, 2)
        self.GL_curves_options.addWidget(self.CB_display_fits, 3, 2)
        # Col 4
        self.GL_curves_options.addWidget(self.CB_joc, 0, 3)
        self.GL_curves_options.addWidget(self.CB_fit_joc, 1, 3)
        self.GL_curves_options.addWidget(self.CB_surface, 2, 3)
        self.GL_curves_options.addWidget(self.CB_force, 3, 3)
        # Col 5
        self.GL_curves_options.addWidget(self.RBT_events_none, 0, 4)
        self.GL_curves_options.addWidget(self.RBT_events_annotations, 1, 4)
        self.GL_curves_options.addWidget(self.CB_events_filter_dist, 2, 4)
        self.GL_curves_options.addWidget(self.CB_events_display_joc, 3, 4)
        self.GL_curves_options.addWidget(self.CB_events_display_fit_joc, 4, 4)

        self.VL_curves_options.addLayout(self.GL_curves_options)
        self.VL_curves_options.addStretch(1)

        self.HL_curves_options = QtWidgets.QHBoxLayout()
        self.HL_curves_options.addLayout(self.VL_curves_options)
        self.HL_curves_options.addStretch(1)
        self.box_curves_options.setLayout(self.HL_curves_options)
        widgets_list.widget_progressbar.update()

        # --- Set up page -----------------------------------------------------
        self.GL = QtWidgets.QGridLayout()

        HL_empty = QtWidgets.QHBoxLayout()
        HL_empty.addStretch(1)

        self.GL_choser = QtWidgets.QGridLayout()
        self.GL_choser.addWidget(self.list_exp, 0, 0)
        self.GL_choser.addWidget(self.CB_app_to_all, 0, 2)
        self.GL_choser.addLayout(HL_empty, 0, 3)

        self.HL_buttons = QtWidgets.QHBoxLayout()
        self.HL_buttons.addWidget(self.BT_multi_meshgrids)
        self.HL_buttons.addWidget(self.BT_profiles)
        self.HL_buttons.addWidget(self.BT_slices)
        self.HL_buttons.addWidget(self.BT_roi_manager)
        self.HL_buttons.addWidget(self.BT_display_in_3D)
        self.HL_buttons.addWidget(self.BT_used_parameters)
        self.HL_buttons.addStretch(1)

        self.GL.addWidget(self.box_meshgrid_options, 0, 0)
        self.GL.addWidget(self.box_curves_options, 0, 1)
        self.GL.addLayout(self.GL_choser, 1, 0, 1, 0)
        self.GL.addWidget(self.MPL_canvas1, 2, 0)
        self.GL.addWidget(self.MPL_canvas2, 2, 1)
        self.GL.addLayout(self.HL_buttons, 3, 0, 1, 0)

        widgets_list.widget_progressbar.update()

        self.VL.addLayout(self.GL)
        self.setLayout(self.VL)

        self.update_widget()

        widgets_list.widget_progressbar.update()

        self.box_meshgrid_options.setFont(self.parent.smallfont)
        self.box_curves_options.setFont(self.parent.smallfont)

    def update_widget(self):
        """Update the widget."""
        data = shared.exp.list[shared.exp.id_selected]

        # For single files don't display the curve selector
        if data.is_single:
            self.curve_selector.setParent(None)
        else:
            self.GL_choser.addWidget(self.curve_selector, 0, 1)
            self.curve_selector.update_widget()

        self.update_GUI("all")
        self.update_MPL("MPL_canvas1")
        self.update_MPL("MPL_canvas2")

    def update_GUI(self, element):
        """Method to update different GUI elements."""
        data = shared.exp.list[shared.exp.id_selected]

        if element == "BTG_meshgrid_type" or element == "all":
            # Refresh the radiobuttons for the choice of the meshgrid
            # Enable or disable buttons
            val = data.stiffness_calculated
            self.RBT_mesh_topography.setEnabled(val)
            self.RBT_mesh_stiffness.setEnabled(val)
            if val:
                if data.indentation_step == 0:
                    self.RBT_mesh_stiffness_slice.setEnabled(False)
                else:
                    self.RBT_mesh_stiffness_slice.setEnabled(True)
            else:
                self.RBT_mesh_stiffness_slice.setEnabled(val)
            val = data.stiffness_corrected
            self.RBT_mesh_stiffness_corr.setEnabled(val)
            if val:
                if data.indentation_step == 0:
                    self.RBT_mesh_stiffness_slice_corr.setEnabled(False)
                else:
                    self.RBT_mesh_stiffness_slice_corr.setEnabled(True)
            else:
                self.RBT_mesh_stiffness_slice_corr.setEnabled(val)
            val = data.work_and_rupture_force1_calculated
            self.RBT_mesh_work.setEnabled(val)
            self.RBT_mesh_rupture_force1.setEnabled(val)
            val = data.events_calculated
            self.RBT_mesh_events_per_curve.setEnabled(val)
            self.RBT_mesh_rupture_force2.setEnabled(val)

            # Check the button
            index = misc_tools.get_meshgrid_type_by_id(data.meshgrid_type)
            self.GRP_mesh_type.button(index).setChecked(True)

        if element == "apply_to_all" or element == "all":
            self.CB_app_to_all.setChecked(shared.exp.apply_to_all_data)

        if element == "radio_buttons" or element == "all":
            if data.used_stiffness_model_selected == 3:
                self.RBT_mesh_stiffness.setText("Slope (Indentation)")
                self.RBT_mesh_stiffness_slice.setText("Slope (Slice)")
            else:
                self.RBT_mesh_stiffness.setText("Elasticity")
                self.RBT_mesh_stiffness_slice.setText("Elasticity (Slice)")

        if element == "button_display_in_3D" or element == "all":
            # 3D is only allowed if vtk is imported
            if consts.ALLOW_VTK:
                if data.stiffness_array is not None and data.scan_size_x != 0:
                    self.BT_display_in_3D.setEnabled(True)
                else:
                    self.BT_display_in_3D.setEnabled(False)
            else:
                self.BT_display_in_3D.setEnabled(False)

        if element == "button_profiles_and_slices" or element == "all":
            if not data.is_single:
                length = len(data.profile_list)
                if data.stiffness_array is not None and length > 0:
                    self.BT_profiles.setEnabled(True)
                else:
                    self.BT_profiles.setEnabled(False)
                if data.stiffness_array is not None:
                    self.BT_slices.setEnabled(True)
                else:
                    self.BT_slices.setEnabled(False)

            else:
                self.BT_profiles.setEnabled(False)
                self.BT_slices.setEnabled(False)

        if element == "buttons_options_curves" or element == "all":
            # Enable force curve option
            if data.stiffness_calculated or \
                    data.work_and_rupture_force1_calculated:
                self.RB_force.setEnabled(True)
            else:
                self.RB_force.setEnabled(False)

            # Enable force curve option (events)
            self.RB_force_events.setEnabled(data.events_calculated)

            # Enable indentation option
            self.RB_indentation_curve.setEnabled(data.stiffness_calculated)

            # Enable residuals option
            self.RB_residuals.setEnabled(data.perform_fit)

            # Check
            if data.display_curve_type == "defl_ext":
                self.RB_defl_ext.setChecked(True)
            elif data.display_curve_type == "force_ext":
                self.RB_force_ext.setChecked(True)
            elif data.display_curve_type == "force":
                self.RB_force.setChecked(True)
            elif data.display_curve_type == "results_force_events":
                self.RB_force_events.setChecked(True)
            elif data.display_curve_type == "indentation":
                self.RB_indentation_curve.setChecked(True)
            elif data.display_curve_type == "residuals":
                self.RB_residuals.setChecked(True)

            # Enable or disable radiobuttons for approach or retraction
            if data.display_curve_type == "force":
                self.radiobutton_trace.setEnabled(data.stiffness_calculated)

                val = data.work_and_rupture_force1_calculated
                self.radiobutton_retrace.setEnabled(val)

                val = data.stiffness_calculated and \
                    data.work_and_rupture_force1_calculated
                self.radiobutton_both.setEnabled(val)

            elif data.display_curve_type == "results_force_events":
                self.radiobutton_trace.setEnabled(False)
                self.radiobutton_retrace.setEnabled(True)
                self.radiobutton_both.setEnabled(False)

            elif data.display_curve_type == "defl_ext" or \
                    data.display_curve_type == "force_ext":
                self.radiobutton_trace.setEnabled(True)
                self.radiobutton_retrace.setEnabled(True)
                self.radiobutton_both.setEnabled(True)

            elif data.display_curve_type == "indentation":
                self.radiobutton_trace.setEnabled(True)
                self.radiobutton_retrace.setEnabled(False)
                self.radiobutton_both.setEnabled(False)

            elif data.display_curve_type == "residuals":
                self.radiobutton_trace.setEnabled(True)
                self.radiobutton_retrace.setEnabled(False)
                self.radiobutton_both.setEnabled(False)

            # Approach/Retraction/Both
            if data.display_trace_retrace == 0:
                self.radiobutton_trace.setChecked(True)
            elif data.display_trace_retrace == 1:
                self.radiobutton_retrace.setChecked(True)
            elif data.display_trace_retrace == 2:
                self.radiobutton_both.setChecked(True)

            # Check option to display POC and fit
            self.CB_poc.setChecked(data.display_poc)
            self.CB_fit_poc.setChecked(data.display_fit_poc)

            # Check option to display the segments for the stiffness
            self.CB_segments.setChecked(data.display_segments)

            # Check the option to display the stiffness fits
            self.CB_display_fits.setChecked(data.display_fits_stiffness)

            # Check option to display JOC, fit and work (surface)
            self.CB_joc.setChecked(data.display_joc)
            self.CB_fit_joc.setChecked(data.display_fit_joc)
            self.CB_surface.setChecked(data.display_surface)
            self.CB_force.setChecked(data.display_force)

            # Check option to display events filter
            self.CB_events_filter_dist.setChecked(
                data.events_display_results_filter_dist)

            # Check option to display events joc
            self.CB_events_display_joc.setChecked(
                data.events_results_display_joc)
            val = data.events_results_display_fit_joc
            self.CB_events_display_fit_joc.setChecked(
                data.events_results_display_fit_joc)

            # Display some results for the events
            if data.events_results_display_annotations == 0:
                self.RBT_events_none.setChecked(True)
            elif data.events_results_display_annotations == 1:
                self.RBT_events_annotations.setChecked(True)

            tp = data.display_curve_type
            dtr = data.display_trace_retrace

            # Approach curve box
            if data.stiffness_calculated:
                if (tp == "defl_ext" or tp == "force_ext") and \
                        (dtr == 0 or dtr == 2):
                    # Deflexion - extension curve or Deflexion - Force
                    # (trace or both)
                    self.CB_poc.setEnabled(True)
                    self.CB_fit_poc.setEnabled(True)
                    self.CB_segments.setEnabled(False)
                    self.CB_display_fits.setEnabled(False)
                elif tp == "force" and (dtr == 0 or dtr == 2):
                    # Force - distance curve (trace or both)
                    self.CB_poc.setEnabled(True)
                    self.CB_fit_poc.setEnabled(True)
                    self.CB_segments.setEnabled(True)
                    self.CB_display_fits.setEnabled(True)
                elif tp == "indentation" and (dtr == 0 or dtr == 2):
                    # Force - Indentation (trace or both)
                    self.CB_poc.setEnabled(False)
                    self.CB_fit_poc.setEnabled(False)
                    self.CB_segments.setEnabled(True)
                    self.CB_display_fits.setEnabled(True)
                elif tp == "residuals" and (dtr == 0 or dtr == 2):
                    self.CB_poc.setEnabled(False)
                    self.CB_fit_poc.setEnabled(False)
                    self.CB_segments.setEnabled(False)
                    self.CB_display_fits.setEnabled(False)
                else:
                    self.CB_poc.setEnabled(False)
                    self.CB_fit_poc.setEnabled(False)
                    self.CB_segments.setEnabled(False)
                    self.CB_display_fits.setEnabled(False)
            else:
                self.CB_poc.setEnabled(False)
                self.CB_fit_poc.setEnabled(False)
                self.CB_segments.setEnabled(False)
                self.CB_display_fits.setEnabled(False)

            # Retraction curve box
            if data.work_and_rupture_force1_calculated or \
                    data.events_calculated:
                if (tp == "defl_ext" or tp == "force_ext") and \
                        (dtr == 1 or dtr == 2):
                    if data.work_and_rupture_force1_calculated:
                        self.CB_joc.setEnabled(True)
                        self.CB_fit_joc.setEnabled(True)
                        self.CB_surface.setEnabled(True)
                        self.CB_force.setEnabled(True)
                    if data.events_calculated:
                        self.RBT_events_none.setEnabled(True)
                        self.RBT_events_annotations.setEnabled(True)
                        self.CB_events_filter_dist.setEnabled(True)
                        self.CB_events_display_joc.setEnabled(False)
                        self.CB_events_display_fit_joc.setEnabled(False)
                elif tp == "force" and (dtr == 1 or dtr == 2):
                    if data.work_and_rupture_force1_calculated:
                        self.CB_joc.setEnabled(True)
                        self.CB_fit_joc.setEnabled(True)
                        self.CB_surface.setEnabled(True)
                        self.CB_force.setEnabled(True)
                    if data.events_calculated:
                        self.RBT_events_none.setEnabled(True)
                        self.RBT_events_annotations.setEnabled(True)
                        self.CB_events_filter_dist.setEnabled(True)
                        self.CB_events_display_joc.setEnabled(False)
                        self.CB_events_display_fit_joc.setEnabled(False)
                elif data.display_curve_type == "results_force_events" and \
                        (dtr == 1 or dtr == 2):
                    if data.work_and_rupture_force1_calculated:
                        self.CB_joc.setEnabled(False)
                        self.CB_fit_joc.setEnabled(False)
                        self.CB_surface.setEnabled(False)
                        self.CB_force.setEnabled(False)
                    if data.events_calculated:
                        self.RBT_events_none.setEnabled(True)
                        self.RBT_events_annotations.setEnabled(True)
                        self.CB_events_filter_dist.setEnabled(True)
                        self.CB_events_display_joc.setEnabled(True)
                        self.CB_events_display_fit_joc.setEnabled(True)
                else:
                    if data.work_and_rupture_force1_calculated:
                        self.CB_joc.setEnabled(False)
                        self.CB_fit_joc.setEnabled(False)
                        self.CB_surface.setEnabled(False)
                        self.CB_force.setEnabled(False)
                    if data.events_calculated:
                        self.RBT_events_none.setEnabled(False)
                        self.RBT_events_annotations.setEnabled(False)
                        self.CB_events_filter_dist.setEnabled(False)
                        self.CB_events_display_joc.setEnabled(False)
                        self.CB_events_display_fit_joc.setEnabled(False)

            if not data.work_and_rupture_force1_calculated:
                self.CB_joc.setEnabled(False)
                self.CB_fit_joc.setEnabled(False)
                self.CB_surface.setEnabled(False)
                self.CB_force.setEnabled(False)
            if not data.events_calculated:
                self.RBT_events_none.setEnabled(False)
                self.RBT_events_annotations.setEnabled(False)
                self.CB_events_filter_dist.setEnabled(False)
                self.CB_events_display_joc.setEnabled(False)
                self.CB_events_display_fit_joc.setEnabled(False)

    def button_clicked(self, button):
        """Method called whenever a button is clicked."""
        data = shared.exp.list[shared.exp.id_selected]
        if button == "button_rand":
            widgets_list.widget_main.change_curve(0, 0, "rand")

        elif button == "back":
            pos = shared.exp.pos_in_last_ten_curves
            pos -= 1
            pos = shared.exp.pos_in_last_ten_curves
            posx = shared.exp.last_ten_curves[pos][0]
            posy = shared.exp.last_ten_curves[pos][1]

            widgets_list.widget_main.change_curve(
                posx, posy, save_in_last_ten=False)
            self.update_GUI("bts_back_forward")

        elif button == "forward":
            pos = shared.exp.pos_in_last_ten_curves
            pos += 1
            posx = shared.exp.last_ten_curves[pos][0]
            posy = shared.exp.last_ten_curves[pos][1]
            widgets_list.widget_main.change_curve(
                posx, posy, save_in_last_ten=False)
            self.update_GUI("bts_back_forward")

        elif button == "BTG_meshgrid_type":
            # Get the data :
            value = self.GRP_mesh_type.checkedId()
            data.meshgrid_type = misc_tools.get_meshgrid_type_as_string(value)
            self.logger.debug("Changing meshgrid type, %s", data.meshgrid_type)
            apply_to_all.apply_to_all("display_options")
            # Update the GUI elements
            self.update_MPL("MPL_canvas1")
            # Update the list of colortables and the min/max values
            if widgets_list.widget_meshgrid_options is not None:
                widgets_list.widget_meshgrid_options.update_widget()
            if widgets_list.widget_multimeshgrids is not None:
                widgets_list.widget_multimeshgrids.update_widget()

        elif button == "button_multi_meshgrids":
            if widgets_list.widget_multimeshgrids is None:
                # Create new widget
                MultiMeshgridsWidget(self)
                widgets_list.widget_multimeshgrids.resize(500, 250)
                widgets_list.widget_multimeshgrids.show()
            else:
                # Bring to front
                widgets_list.widget_multimeshgrids.activateWindow()
                widgets_list.widget_multimeshgrids.raise_()

        elif button == "button_profiles":
            if widgets_list.widget_profiles is None:
                # Create new widget
                ProfileWidget(self)
                widgets_list.widget_profiles.resize(800, 550)
                widgets_list.widget_profiles.setWindowTitle("Profiles")
                widgets_list.widget_profiles.show()
            else:
                # Bring to front
                widgets_list.widget_profiles.activateWindow()
                widgets_list.widget_profiles.raise_()

        elif button == "button_slices":
            if widgets_list.widget_slices is None:
                # Create new widget
                SlicesWidget(self)
                widgets_list.widget_slices.resize(550, 550)
                widgets_list.widget_slices.setWindowTitle("Slices")
                widgets_list.widget_slices.show()
                # Update the profiles
                self.MPL_canvas1.canvas.update_blit("all")
            else:
                # Bring to front
                widgets_list.widget_slices.activateWindow()
                widgets_list.widget_slices.raise_()

        elif button == "get_details":
            if widgets_list.widget_details is None:
                DetailsWidget(self)
                widgets_list.widget_details.resize(250, 500)
                widgets_list.widget_details.show()
            else:
                widgets_list.widget_details.activateWindow()
                widgets_list.widget_details.raise_()

        elif button == "button_roi_manager":
            if widgets_list.widget_roi_manager is None:
                RoiManagerWidget(self)
                widgets_list.widget_roi_manager.resize(350, 500)
                widgets_list.widget_roi_manager.show()
            else:
                widgets_list.widget_roi_manager.activateWindow()
                widgets_list.widget_roi_manager.raise_()

        elif button == "button_group_type_of_curve":
            # Save the curve type
            curve_type = self.BT_group_type_of_curve.checkedId()
            if curve_type == 0:
                data.display_curve_type = "defl_ext"
            elif curve_type == 1:
                data.display_curve_type = "force_ext"
            elif curve_type == 2:
                data.display_curve_type = "force"
            elif curve_type == 3:
                data.display_curve_type = "results_force_events"
            elif curve_type == 4:
                data.display_curve_type = "indentation"
            elif curve_type == 5:
                data.display_curve_type = "residuals"

            tp = data.display_curve_type
            self.logger.debug("Changing curve type, %s", tp)

            dtr = data.display_trace_retrace

            # Check if we can display trace or retrace
            if data.display_curve_type == "force":
                if not data.stiffness_calculated and (dtr == 0 or dtr == 2):
                    data.display_trace_retrace = 1
                elif not data.work_and_rupture_force1_calculated and \
                        (dtr == 1 or dtr == 2):
                    data.display_trace_retrace = 0
            elif data.display_curve_type == "indentation":
                if dtr == 1 or dtr == 2:
                    # Display approach only
                    data.display_trace_retrace = 0
            elif data.display_curve_type == "results_force_events":
                if dtr == 0 or dtr == 2:
                    # Display retraction only
                    data.display_trace_retrace = 1

            val = str(data.display_trace_retrace)
            self.logger.debug("Curve mode, %s", val)

            # Update plot and GUI
            self.update_GUI("buttons_options_curves")
            self.update_MPL("MPL_canvas2")

        elif button == "button_group_trace_or_retrace":
            val = self.BT_group_trace_or_retrace.checkedId()
            data.display_trace_retrace = val
            self.update_GUI("buttons_options_curves")
            self.update_MPL("MPL_canvas2")

        elif button == "button_group_events":
            val = self.BT_group_events.checkedId()
            data.events_results_display_annotations = val
            self.update_MPL("MPL_canvas2")

        elif button == "button_display_in_3D":
            # Go to stiffness meshgrid
            self.RBT_mesh_stiffness.setChecked(True)
            self.button_clicked("BTG_meshgrid_type")

            if widgets_list.widget_vtk is not None:
                widgets_list.widget_vtk.activateWindow()
                widgets_list.widget_vtk.raise_()

            else:
                # Make sure the indentation widget is closed
                if widgets_list.widget_indentation is not None:
                    widgets_list.widget_indentation.close()

                # Ask the user for which tomography has to be loaded
                self.open_vtk_entry_widget()

        elif button == "button_show_slices":
            if widgets_list.widget_slices is None:
                SlicesWidget(self.parent)
                widgets_list.widget_slices.resize(500, 500)
                widgets_list.widget_slices.show()
            else:
                widgets_list.widget_slices.activateWindow()
                widgets_list.widget_slices.raise_()

        elif button == "used_parameters":
            if widgets_list.widget_used_parameters is None:
                UsedParametersWidget(self)
                widgets_list.widget_used_parameters.show()
            else:
                widgets_list.widget_used_parameters.activateWindow()
                widgets_list.widget_used_parameters.raise_()

    def list_updated(self, name):
        """Called when a list is updated."""
        if name == "list1":
            self.parent.file_changed("list1")

    def checkbox_clicked(self, name):
        """Called when clicking on a checkbox."""
        data = shared.exp.current_data

        if name == "apply_to_all":
            shared.exp.apply_to_all_data = self.CB_app_to_all.isChecked()
            apply_to_all.apply_to_all("display_options")
            if shared.exp.apply_to_all_data:
                # Update the events array (with a progressbar)
                events_refresher.update_events()

        elif name == "checkbox_poc":
            data.display_poc = self.CB_poc.isChecked()
            self.update_MPL("MPL_canvas2")

        elif name == "checkbox_fit_poc":
            data.display_fit_poc = self.CB_fit_poc.isChecked()
            self.update_MPL("MPL_canvas2")

        elif name == "checkbox_segments":
            data.display_segments = self.CB_segments.isChecked()
            self.update_MPL("MPL_canvas2")

        elif name == "display_fits_stiffness":
            data.display_fits_stiffness = self.CB_display_fits.isChecked()
            self.update_MPL("MPL_canvas2")

        elif name == "checkbox_joc":
            data.display_joc = self.CB_joc.isChecked()
            self.update_MPL("MPL_canvas2")

        elif name == "checkbox_fit_joc":
            data.display_fit_joc = self.CB_fit_joc.isChecked()
            self.update_MPL("MPL_canvas2")

        elif name == "checkbox_surface":
            data.display_surface = self.CB_surface.isChecked()
            self.update_MPL("MPL_canvas2")

        elif name == "force":
            data.display_force = self.CB_force.isChecked()
            self.update_MPL("MPL_canvas2")

        elif name == "checkbox_events_filter_dist":
            val = self.CB_events_filter_dist.isChecked()
            data.events_display_results_filter_dist = val
            self.update_MPL("MPL_canvas2")

        elif name == "checkbox_events_display_joc":
            val = self.CB_events_display_joc.isChecked()
            data.events_results_display_joc = val
            self.update_MPL("MPL_canvas2")

        elif name == "checkbox_events_display_fit_joc":
            val = self.CB_events_display_fit_joc.isChecked()
            data.events_results_display_fit_joc = val
            self.update_MPL("MPL_canvas2")

    def open_vtk_entry_widget(self, mode=None):
        """Opens a widget asking the user which AFM tomograpy he wants to add."""
        if widgets_list.widget_vtk_entry is None:
            # The first time, ask the user which tomographies he wants
            # to load in VTK
            VTKEntryWidget(self, mode)
            if not consts.AUTO3D and consts.UNIT_TESTING is False:
                widgets_list.widget_vtk_entry.show()

        else:
            widgets_list.widget_vtk_entry.activateWindow()
            widgets_list.widget_vtk_entry.raise_()

    def update_MPL(self, what):
        """Update the plots."""
        if what == "MPL_canvas1":
            self.MPL_canvas1.update_plot()

        elif what == "MPL_canvas2":
            self.MPL_canvas2.update_plot()

    def meshgrid_canvas_press(self, ev):
        """Matplotlib method, called when the user clicks on the meshgrid."""
        data = shared.exp.list[shared.exp.id_selected]

        if not isinstance(widgets_list.widget_results.MPL_canvas1.canvas,
                          MainPlot):
            return  # Nothing to do if no real canvas

        # ev.button == 1 = left click
        axes = widgets_list.widget_results.MPL_canvas1.canvas.axes
        if ev.button == 1 and ev.inaxes == axes:
            xpos, ypos = misc_tools.get_position_on_meshgrid(ev)

            mode = self.cursor_mode

            if mode == "select":
                widgets_list.widget_main.change_curve(xpos, ypos)

            elif mode == "profile_add":
                data.profile_starts.append([xpos - 1, ypos - 1])
                data.profile_list.append([xpos - 1, ypos - 1])

            elif mode == "roi_add" or mode == "roi_erase":
                self.draw_roi(xpos, ypos)

            elif mode == "profile_erase":
                # Search in profile_list if there is a profile to delete
                if data.profile_list != []:
                    count = 0
                    for profile in data.profile_list:
                        for i in range(len(profile[0])):
                            pos = [profile[0][i], profile[1][i]]
                            if pos == [xpos - 1, ypos - 1]:
                                del data.profile_list[count]
                                del data.profile_starts[count]
                                data.profile_count -= 1
                                # We have found the profile, so we stop
                                break
                        count += 1

                # Close profile widget if needed
                w_profiles = widgets_list.widget_profiles
                if data.profile_count == 0 and w_profiles is not None:
                    widgets_list.widget_profiles.close()

                # Reset profiles in slices widget if there is no more profile
                # (and a custome profile was selected)
                direct = data.slice_view_from_direction
                if direct == 4 or direct == 5 and data.profile_list == []:
                    # Go back to first direction
                    data.slice_view_from_direction = 1

                    # If the slices widget is open, update it
                    if widgets_list.widget_slices is not None:
                        widgets_list.widget_slices.update_widget()
                        widgets_list.widget_slices.update_MPL()

                # Update button
                self.update_GUI("button_profiles_and_slices")
                # Finalize plotting
                self.MPL_canvas1.canvas.update_blit("all")

                # Reset mouse cursor
                self.cursor_mode = "select"
                QtWidgets.QApplication.restoreOverrideCursor()

        elif (ev.button == 3 and
              isinstance(self.MPL_canvas1.canvas, MainPlot) and
              ev.inaxes == self.MPL_canvas1.canvas.axes):
            # Get position in screen coordinates (top left position of canvas)
            globalpos = self.MPL_canvas1.mapToGlobal(QtCore.QPoint(0, 0))
            # Recalculate positions with canvas position
            x = globalpos.x() + ev.x
            y = globalpos.y() + (self.meshgrid_canvas_size[1] - ev.y)
            pos = QtCore.QPoint(x, y)
            meshgrid_menu.create_meshgrid_menu(self, "results", pos)

    def canvas_release(self, event):
        """The canvas_release is called whenever the mouse click is released."""
        data = shared.exp.list[shared.exp.id_selected]

        inaxes = event.inaxes
        button = event.button

        if not isinstance(widgets_list.widget_results.MPL_canvas1.canvas,
                          MainPlot):
            return  # Nothing to do if no real canvas

        # Release for profile, only inside canvas
        axes = widgets_list.widget_results.MPL_canvas1.canvas.axes
        mode = self.cursor_mode
        if inaxes == axes and button == 1:
            ls = data.profile_list
            if mode == "profile_add":
                length = len(ls[data.profile_count][0])

            if mode == "profile_add" and length >= 3:
                data.profile_count = data.profile_count + 1
                self.update_GUI("button_profiles_and_slices")

                # Finalize plotting
                self.MPL_canvas1.canvas.update_blit("all")

            elif mode == "profile_add" and length < 3 and ls != []:
                # Nothing has happened
                del data.profile_starts[data.profile_count]
                del data.profile_list[data.profile_count]

        # Releasing for the ROI can be done outside the canvas
        if self.cursor_mode == "roi_add" or self.cursor_mode == "roi_erase":
            # Finalize plotting
            if widgets_list.widget_roi_manager is not None:
                widgets_list.widget_roi_manager.update_widget()
            self.MPL_canvas1.canvas.update_blit("all")

            # Update the ROIs if displayed on multimeshgrids widget
            if widgets_list.widget_multimeshgrids is not None:
                if shared.exp.display_roi_in_multimeshgrid:
                    widgets_list.widget_multimeshgrids.update_MPL()

        # Reset cursor mode
        self.cursor_mode = "select"
        QtWidgets.QApplication.restoreOverrideCursor()

    def canvas_motion(self, event):
        """The canvas_motion is called when the mouse moves over the plots.

        Depending on the action, a profile or a ROI can be drawn with
        this method.
        """
        data = shared.exp.list[shared.exp.id_selected]

        inaxes = event.inaxes
        button = event.button

        if not isinstance(self.MPL_canvas1.canvas, MainPlot):
            return  # Nothing to do if no real canvas

        if inaxes == self.MPL_canvas1.canvas.axes and button == 1:
            # Mouse motion over meshgrid, with left click

            # Get the current position
            xpos, ypos = misc_tools.get_position_on_meshgrid(event)

            mode = self.cursor_mode

            if mode == "profile_add":
                # Delete current profile
                if data.profile_list != []:
                    del data.profile_list[data.profile_count]

                # Add profile to list
                data.profile_list.append(
                    misc_tools.make_path(
                        data.profile_starts[data.profile_count],
                        [xpos - 1, ypos - 1]))

                # Draw profile
                cv = self.MPL_canvas1.canvas
                prof = data.profile_list[data.profile_count]
                cv.update_blit("profile_preview", profile=prof)

            elif mode == "roi_add" or mode == "roi_erase":
                # Draw a ROI
                self.draw_roi(xpos, ypos)

        elif (inaxes == self.MPL_canvas1.canvas.axes and
              (self.cursor_mode == "roi_add" or
               self.cursor_mode == "roi_erase")):
            # Mouse motion over meshgrid, no click

            # Get the current position
            xpos, ypos = misc_tools.get_position_on_meshgrid(event)

            # Predraw the ROI
            self.pre_draw_roi(xpos, ypos)

        elif (inaxes == self.MPL_canvas1.canvas.axes
              and self.cursor_mode == "profile_erase"):
            # Get the current position
            xpos, ypos = misc_tools.get_position_on_meshgrid(event)

            # Search in profile_list if there is a profile to delete
            if data.profile_list != []:
                count = 0
                for profile in data.profile_list:
                    found = False
                    for i in range(len(profile[0])):
                        prof = [profile[0][i], profile[1][i]]
                        if prof == [xpos - 1, ypos - 1]:
                            # We have found the profile, so we stop
                            found = True
                            break
                    if not found:
                        count = count + 1

                if found and self.current_profile_erase is None:
                    self.current_profile_erase = count
                    cv = self.MPL_canvas1.canvas
                    prof = data.profile_list[count]
                    cv.update_blit("profile_preview_erase", profile=prof)
                elif found is False and self.current_profile_erase is not None:
                    self.MPL_canvas1.canvas.update_blit("all")
                    self.current_profile_erase = None

    def pre_draw_roi(self, xpos, ypos):
        """Pre draw the ROI on the meshgrid.

        Used when moving the mouse without clicking above the meshgrid.
        """
        data = shared.exp.list[shared.exp.id_selected]
        local_roi, local_roi_mask = prepare_roi(xpos, ypos)
        local_roi, local_roi_mask = \
            check_roi_boundaries(data, local_roi, local_roi_mask)
        cv = self.MPL_canvas1.canvas
        if self.cursor_mode == "roi_add":
            cv.update_blit("preview_add", local_roi, local_roi_mask)
        elif self.cursor_mode == "roi_erase":
            cv.update_blit("preview_erase", local_roi, local_roi_mask)

    def draw_roi(self, xpos, ypos):
        """Draw the ROI on the meshgrid."""
        data = shared.exp.list[shared.exp.id_selected]
        id_roi = data.roi_selected_row
        local_roi, local_roi_mask = prepare_roi(xpos, ypos)
        local_roi, local_roi_mask = \
            self.check_roi_values(data, local_roi, local_roi_mask, id_roi)

        cv = self.MPL_canvas1.canvas

        if self.cursor_mode == "roi_add":
            for i in range(len(local_roi)):
                if local_roi_mask[i] and local_roi[i] is not None:
                    data.roi_list[id_roi].values.extend([local_roi[i]])
            cv.update_blit("local_add", local_roi, local_roi_mask)

        elif self.cursor_mode == "roi_erase":
            for i in range(len(local_roi) - 1, -1, -1):
                vals = data.roi_list[id_roi].values
                inlist = math_tools.in_list(vals, local_roi[i])
                if local_roi_mask[i] and inlist:
                    ind = local_roi[i]
                    ls = data.roi_list[id_roi]
                    del ls.values[data.roi_list[id_roi].values.index(ind)]
            cv.update_blit("local_erase", local_roi, local_roi_mask)

        if widgets_list.widget_roi_manager is not None:
            widgets_list.widget_roi_manager.update_widget()

    def check_roi_values(self, data, values, mask, id_roi):
        """Return well formatted array."""
        new_values = []
        new_mask = []
        for i in range(len(values)):
            position = values[i]

            cond1 = position[0] >= 0 and position[1] >= 0
            cond2 = position[0] < data.nbr_pixels_x
            cond3 = position[1] < data.nbr_pixels_y

            if cond1 and cond2 and cond3:
                if self.cursor_mode == "roi_add":
                    vals = data.roi_list[id_roi].values
                    inlist = math_tools.in_list(vals, position)
                    if not inlist:
                        new_values.append(position)
                        new_mask.append(mask[i])
                    else:
                        new_values.append(None)
                        new_mask.append(False)
                elif self.cursor_mode == "roi_erase":
                    new_values.append(position)
                    new_mask.append(mask[i])

        return new_values, new_mask

    def curves_canvas_press(self, event):
        """Called when there is a click on the curve."""
        if event.button == 3 and event.inaxes == self.MPL_canvas2.canvas.axes:
            # Get position in screen coordinates (top left position of canvas)
            globalpos = self.MPL_canvas2.mapToGlobal(QtCore.QPoint(0, 0))
            # Recalculate positions with canvas position
            x = globalpos.x() + event.x
            y = globalpos.y() + (self.curves_canvas_size[1] - event.y)
            pos = QtCore.QPoint(x, y)
            self.popUpMenu(pos, "curve_results")

    def change_cursor_mode(self, mode):
        """Method to change the cursor mode.

        The value of self.cursor_mode is set, and the cursor shape is changed
        with the PyQT setOverrideCursor method.
        """
        self.cursor_mode = mode
        data = shared.exp.current_data
        if mode == "profile_add":
            curs = QtGui.QCursor(QtCore.Qt.CrossCursor)
            QtWidgets.QApplication.setOverrideCursor(curs)

        elif mode == "roi_add":
            curs = QtGui.QCursor(QtCore.Qt.CrossCursor)
            QtWidgets.QApplication.setOverrideCursor(curs)
            self.MPL_canvas1.canvas.first_blit = True

        elif mode == "roi_erase":
            curs = QtGui.QCursor(QtCore.Qt.CrossCursor)
            QtWidgets.QApplication.setOverrideCursor(curs)
            self.MPL_canvas1.canvas.first_blit = True
            # Hide all the other roi's
            for index in range(len(data.roi_list)):
                if index != data.roi_selected_row:
                    data.roi_list[index].display = False
            self.MPL_canvas1.canvas.update_blit("all")

        elif mode == "profile_erase":
            curs = QtGui.QCursor(QtCore.Qt.ForbiddenCursor)
            QtWidgets.QApplication.setOverrideCursor(curs)

    def popUpMenu(self, pos, plot_type):
        """Pop Up menu opened when right-clicking on the curve or the meshgrid."""
        data = shared.exp.current_data

        # Create menu
        menu = QtWidgets.QMenu()

        if plot_type == "curve_results":
            # Show details widget
            action = QtWidgets.QAction("Data (results)", self)
            action.triggered.connect(
                lambda: self.button_clicked("get_details"))
            menu.addAction(action)

            # Hide/display approach/retraction legend
            menu.addSeparator()

            if shared.exp.display_curve_legend:
                txt = "Hide legend"
            else:
                txt = "Display legend"

            action = QtWidgets.QAction(txt, self)
            action.triggered.connect(self.display_legend)
            menu.addAction(action)

            if data.events_calculated:
                menu.addSeparator()
                action = QtWidgets.QAction("Change dist filter", self)
                action.triggered.connect(lambda: self.change_filter("dist"))
                menu.addAction(action)

            # Units for the distance
            menu.addSeparator()
            if shared.exp.current_data.curve_distance_units == "nm":
                txt = "um"
            elif shared.exp.current_data.curve_distance_units == "um":
                txt = "nm"

            action = QtWidgets.QAction("Set distance units to " + txt, self)
            action.triggered.connect(self.change_distance_units)
            menu.addAction(action)

            tp = data.display_curve_type

            if tp == "force" or tp == "force_ext" or \
                    tp == "results_force_events" or tp == "indentation" or tp == "residuals":
                # Units for the force
                if shared.exp.current_data.curve_force_units == "nN":
                    txt = "pN"
                elif shared.exp.current_data.curve_force_units == "pN":
                    txt = "nN"

                action = QtWidgets.QAction("Set force units to " + txt, self)
                action.triggered.connect(self.change_force_units)
                menu.addAction(action)

        # Display menu
        menu.popup(pos, menu.menuAction())
        menu.exec_()

    def change_force_units(self):
        """Method called to change the units of the force in the curve plot.

        The units can be either nN or pN.
        """
        if shared.exp.current_data.curve_force_units == "nN":
            shared.exp.current_data.curve_force_units = "pN"
        elif shared.exp.current_data.curve_force_units == "pN":
            shared.exp.current_data.curve_force_units = "nN"

        self.update_MPL("MPL_canvas2")

    def change_distance_units(self):
        """Method called to change the units of the distances in the curve plot.

        The units can be either nm or um.
        """
        if shared.exp.current_data.curve_distance_units == "nm":
            shared.exp.current_data.curve_distance_units = "um"
        elif shared.exp.current_data.curve_distance_units == "um":
            shared.exp.current_data.curve_distance_units = "nm"

        self.update_MPL("MPL_canvas2")

    def display_legend(self):
        """Method called to display or to hide the legend on the curve plot.

        The legend displays the color and the name of the approach and the
        retraction curve.
        """
        if shared.exp.display_curve_legend:
            shared.exp.display_curve_legend = False
        else:
            shared.exp.display_curve_legend = True

        self.update_MPL("MPL_canvas2")

    def change_filter(self, filter_value):
        """Update the value of the filter for the events."""
        widget = ChangeFilterWidget(self, filter_value)
        if filter_value == "slope":
            widget.setWindowTitle("Change slope filter")
        elif filter_value == "dist":
            widget.setWindowTitle("Change dist filter")
        widget.resize(200, 100)
        widget.activateWindow()
        widget.show()


def prepare_roi(xpos, ypos):
    """Prepare the ROI.

    Returns an array with ids on the meshgrid, depending on the chosen
    pencil.
    """
    shape = shared.exp.roi_cursor_shape
    dot_size = shared.exp.roi_cursor_size
    local_roi = []
    local_roi_mask = []
    if shape == "square":
        for i in range(-dot_size - 1, dot_size):
            for j in range(-dot_size - 1, dot_size):
                local_roi.append([xpos + i, ypos + j])
                local_roi_mask.append(True)
    elif shape == "dot":
        local_roi.append([xpos - 1, ypos - 1])
        local_roi_mask.append(True)
    elif shape == "circle":
        circle_limit = math.sqrt((dot_size - 0.5) ** 2 + (dot_size - 0.5) ** 2)
        for i in range(-dot_size, dot_size + 1):
            for j in range(-dot_size, dot_size + 1):
                dist_to_center = math.sqrt(i ** 2 + j ** 2)
                local_roi.append([xpos + i - 1, ypos + j - 1])
                if dist_to_center < circle_limit:
                    local_roi_mask.append(True)
                else:
                    local_roi_mask.append(False)

    return local_roi, local_roi_mask


def check_roi_boundaries(data, values, mask):
    """Check if a ROI value is inside the meshgrid."""
    new_values = []
    new_mask = []
    for i in range(len(values)):
        position = values[i]
        if position is not None:
            cond1 = position[0] >= 0 and position[1] >= 0
            cond2 = position[0] < data.nbr_pixels_x
            cond3 = position[1] < data.nbr_pixels_y

            if cond1 and cond2 and cond3:
                new_values.append(position)
                new_mask.append(mask[i])
            else:
                new_values.append(None)
                new_mask.append(False)

    return new_values, new_mask
