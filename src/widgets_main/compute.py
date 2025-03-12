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

"""Main compute widget, used to compute elasticity, events, forces, ..."""

from .. import shared
from .. import widgets_list
from .. import experiment
import logging
from .. import consts
from PyQt5 import QtWidgets
from ..tools import misc_tools
from ..tools import apply_to_all
from ..tools import stat_tools
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFButtonGroup
from ..widgets.clean_up import CleanUpWidget
from ..widgets.flatten import FlattenWidget
from ..widgets_main.subwidgets.events_parameters import EventsParamsWidget
from ..widgets_main.subwidgets.stiffness_parameters import StiffnessParamsWidget
from ..widgets_main.subwidgets.curve_selector import ResultsSelectorWidget
from ..widgets_main.subwidgets.fit_parameters import FitParamsWidget
from ..widgets.curve_modify import CurveModWidget
from ..plots.PYAFPlot import PYAFPlot
from ..compute_tools.stiffness.main import Stiffness
from ..compute_tools.work_and_rupture_force.main import WorkAndRuptureForce
from ..compute_tools.stiffness_correction import StiffnessCorrection
from ..compute_tools.loading_rates import GetLoadingRate
from ..compute_tools.events import GetEvents
from ..tools.PYAFWidget import PYAFWidget


class ComputeWidget(PYAFWidget):
    """Main compute widget.

    The main compute widget is in the second tab of PYAF. It is called from the
    QMainWindow in main.py.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_compute")

        self.parent = parent

        self.logger = logging.getLogger()

        # Objects for the computations
        self.compute_stiffness = Stiffness()
        self.compute_work_and_rupture_force = WorkAndRuptureForce()
        self.stiffness_correction = StiffnessCorrection()
        self.compute_events = GetEvents()
        self.get_lr = GetLoadingRate()

        self.canvas_resolution = self.parent.canvas_resolution
        self.curve_canvas_size = [800, 350]
        self.meshgrid_canvas_size = [300, 300]

        # --- Main layout -----------------------------------------------------

        self.VL = QtWidgets.QVBoxLayout()

        # --- Experiment ------------------------------------------------------
        self.GL_choosers = QtWidgets.QGridLayout()
        self.GL_choosers.setContentsMargins(0, 0, 0, 0)

        self.HL_choser = QtWidgets.QHBoxLayout()
        self.list_exp2 = PYAFComboBox(self, "list2")
        self.list_exp2.setFixedSize(350, 26)
        for data in shared.exp.list:
            self.list_exp2.addItem(data.filename)
        self.HL_choser.addWidget(self.list_exp2)

        name = "widget_curve_selector_compute"
        self.curve_selector = ResultsSelectorWidget(self, name)
        # self.curve_selector is added to the layout in update_widget().

        HL_empty = QtWidgets.QHBoxLayout()
        HL_empty.addStretch(1)

        self.GL_choosers.addLayout(self.HL_choser, 0, 0)
        self.GL_choosers.addLayout(HL_empty, 0, 2)

        widgets_list.widget_progressbar.update()

        self.CB_zoom_in = PYAFCheckBox(self, "zoom_in", "Zoom in")
        self.IN_zoom_in_factor = PYAFInput(
            self, "zoom_in_factor", "Zoom factor")
        self.BT_next_zoom_element = PYAFButton(
            self, "next_zoom_element", ">", size=40)

        widgets_list.widget_progressbar.update()

        # Canvas Meshgrid
        sizes = [self.meshgrid_canvas_size[0] // self.canvas_resolution,
                 self.meshgrid_canvas_size[1] // self.canvas_resolution,
                 self.canvas_resolution]

        self.W_mesh = QtWidgets.QWidget()
        label = "small_meshgrid_compute"
        self.MPL_meshgrid = PYAFPlot(self, label, self.W_mesh, sizes)
        self.MPL_meshgrid.mpl_connect("button_press_event",
                                      self.mesh_canvas_press)
        self.MPL_meshgrid.tab_id = 1

        # Canvas curves
        sizes = [self.curve_canvas_size[0] // self.canvas_resolution,
                 self.curve_canvas_size[1] // self.canvas_resolution,
                 self.canvas_resolution]
        self.W_curve = QtWidgets.QWidget()
        self.MPL_curve = PYAFPlot(self, "fit_preview", self.W_curve, sizes)

        widgets_list.widget_progressbar.update()

        # Calc All
        self.CB_compute_all = PYAFCheckBox(self, "compute_all", "Compute all")
        label = "Apply to all"
        self.CB_apply_to_all = PYAFCheckBox(self, "apply_to_all", label)

        self.HL_options = QtWidgets.QHBoxLayout()
        self.BT_flatten = PYAFButton(self, "flatten", "Flatten")
        self.BT_curve_mod = PYAFButton(self, "curve_mod", "Curve modifications")
        self.BT_clean_up = PYAFButton(self, "clean_up", "Clean Up")
        self.HL_options.addWidget(self.BT_flatten)
        self.HL_options.addWidget(self.BT_curve_mod)
        self.HL_options.addWidget(self.BT_clean_up)

        self.HL_options.addWidget(self.CB_compute_all)
        self.HL_options.addWidget(self.CB_apply_to_all)
        self.HL_options.addWidget(self.CB_zoom_in)
        self.HL_options.addWidget(self.IN_zoom_in_factor)
        self.HL_options.addWidget(self.BT_next_zoom_element)
        self.HL_options.addStretch(1)

        widgets_list.widget_progressbar.update()

        # Tabs
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.currentChanged.connect(self.tabs_changed)

        self.widget_app = QtWidgets.QWidget()
        self.widget_ret = QtWidgets.QWidget()
        self.widget_events = QtWidgets.QWidget()
        self.widget_loading_rates = QtWidgets.QWidget()
        self.widget_stiff_corr = QtWidgets.QWidget()
        self.HL_widget_app = QtWidgets.QHBoxLayout()
        self.HL_widget_ret = QtWidgets.QHBoxLayout()
        self.VL_widget_events = QtWidgets.QVBoxLayout()

        # --- Stiffness widget ------------------------------------------------

        W_stiff_exp_params = QtWidgets.QWidget()
        W_stiff_segments = QtWidgets.QWidget()
        self.Box_stiff_exp_params = QtWidgets.QGroupBox("Model parameters")
        self.Box_stiff_segments = QtWidgets.QGroupBox("Computation method")

        # Options
        label = "Tip radius [nm]"
        self.IN_tip_radius = PYAFInput(self, "tip_radius", label, 60)
        label = "Tip half-opening angle (" + "\u00b0" + ")"
        self.IN_tip_angle = PYAFInput(self, "tip_angle", label)
        self.IN_poisson = PYAFInput(self, "poisson", "Poisson ratio")
        self.IN_tip_radius.input.setValidator(misc_tools.validator("UF"))
        self.IN_tip_angle.input.setValidator(misc_tools.validator("UF"))
        self.IN_poisson.input.setValidator(misc_tools.validator("UF"))

        # Model
        self.HL_list_stiff = QtWidgets.QHBoxLayout()
        self.label_models_list = QtWidgets.QLabel("Model : ")
        self.list_models = PYAFComboBox(self, "list_models")
        self.models = [
            "Hertz (Paraboloid)", "Sneddon (Cone)", "Bilodeau (Pyramid)", "Stiffness tomography", "Flat punch"]
        for model in self.models:
            self.list_models.addItem(str(model))

        self.HL_list_stiff.addWidget(self.label_models_list)
        self.HL_list_stiff.addWidget(self.list_models)
        self.HL_list_stiff.addStretch(1)

        VL_empty1 = QtWidgets.QVBoxLayout()
        VL_empty1.addStretch(1)

        self.GL_stiffness_params = QtWidgets.QGridLayout()
        self.GL_stiffness_params.addLayout(self.HL_list_stiff, 0, 0)
        self.GL_stiffness_params.addWidget(self.IN_poisson, 1, 0)
        # Position 2, defined in update_widget
        self.GL_stiffness_params.addLayout(VL_empty1, 3, 0)

        self.Box_stiff_exp_params.setLayout(self.GL_stiffness_params)
        VL_stiff_exp_params = QtWidgets.QVBoxLayout()
        VL_stiff_exp_params.addWidget(self.Box_stiff_exp_params)
        W_stiff_exp_params.setLayout(VL_stiff_exp_params)

        W_fit_params_pocs = FitParamsWidget(self, "widget_fit_param_pocs")

        # Fit parameters stiffness
        self.VL_fit_params_stiff = QtWidgets.QVBoxLayout()

        self.GRP_stiffness_params = PYAFButtonGroup(self, "stiffness_params")
        self.RBT_calculus_params = QtWidgets.QRadioButton("Calculus/Tomography")
        self.RBT_fitting_params = QtWidgets.QRadioButton("Fitting")
        self.GRP_stiffness_params.addButton(self.RBT_calculus_params, 0)
        self.GRP_stiffness_params.addButton(self.RBT_fitting_params, 1)

        BOX_stiffness_mode = QtWidgets.QGroupBox("Computation method")
        HL_stiffness_mode = QtWidgets.QHBoxLayout()
        HL_stiffness_mode.addStretch(1)
        HL_stiffness_mode.addWidget(self.RBT_calculus_params)
        HL_stiffness_mode.addStretch(1)
        HL_stiffness_mode.addWidget(self.RBT_fitting_params)
        HL_stiffness_mode.addStretch(1)
        BOX_stiffness_mode.setLayout(HL_stiffness_mode)

        # Stacked widget
        self.STK_stiffness_mode = QtWidgets.QStackedWidget(self)

        self.widget_calculus = StiffnessParamsWidget(self, "widget_calculus")
        self.widget_fitting = StiffnessParamsWidget(self, "widget_fitting")

        self.STK_stiffness_mode.addWidget(self.widget_calculus)
        self.STK_stiffness_mode.addWidget(self.widget_fitting)

        self.VL_fit_params_stiff.addWidget(BOX_stiffness_mode)
        self.VL_fit_params_stiff.addWidget(self.STK_stiffness_mode)

        # self.Box_stiff_segments.setLayout(self.GL_params_stiff)
        # VL_stiff_segments = QtWidgets.QVBoxLayout()
        # VL_stiff_segments.addWidget(self.Box_stiff_segments)
        W_stiff_segments.setLayout(self.VL_fit_params_stiff)

        VL_button_compute_stiffness = QtWidgets.QVBoxLayout()
        self.BT_compute_stiffness = PYAFButton(
            self, "get_stiffness", "Compute")
        self.BT_compute_stiffness.setFixedSize(100, 100)

        VL_button_compute_stiffness.addStretch(1)
        VL_button_compute_stiffness.addWidget(self.BT_compute_stiffness)
        VL_button_compute_stiffness.addStretch(1)

        self.HL_widget_app.addWidget(W_stiff_exp_params)
        self.HL_widget_app.addWidget(W_fit_params_pocs)
        self.HL_widget_app.addWidget(W_stiff_segments)
        self.HL_widget_app.addLayout(VL_button_compute_stiffness)

        self.widget_app.setLayout(self.HL_widget_app)

        widgets_list.widget_progressbar.update()

        # --- Work and rupture force widget -----------------------------------

        W_fit_params_jocs = FitParamsWidget(self, "widget_fit_param_jocs")

        self.BT_compute_work_and_rupture_force = PYAFButton(
            self, "get_work_and_rupture_force", "Compute")
        self.BT_compute_work_and_rupture_force.setFixedSize(100, 100)

        VL_button_compute_work = QtWidgets.QVBoxLayout()
        VL_button_compute_work.addStretch(1)
        VL_button_compute_work.addWidget(
            self.BT_compute_work_and_rupture_force)
        VL_button_compute_work.addStretch(1)

        self.HL_widget_ret.addWidget(W_fit_params_jocs)
        self.HL_widget_ret.addStretch(1)
        self.HL_widget_ret.addLayout(VL_button_compute_work)

        self.widget_ret.setLayout(self.HL_widget_ret)

        widgets_list.widget_progressbar.update()

        # Events --------------------------------------------------------------

        HL_events = QtWidgets.QHBoxLayout()

        self.widget_kernel = EventsParamsWidget(self, "widget_kernel")
        self.widget_msf = EventsParamsWidget(self, "widget_msf")

        # Algorithm type
        self.GRP_algorithms = PYAFButtonGroup(self, "events_algorithms")
        self.RBT_kernel = QtWidgets.QRadioButton("Kernel")
        self.RBT_msf = QtWidgets.QRadioButton("MSF")
        self.GRP_algorithms.addButton(self.RBT_kernel, 0)
        self.GRP_algorithms.addButton(self.RBT_msf, 1)

        BOX_algorithms = QtWidgets.QGroupBox("Methods")
        HL_algorithms = QtWidgets.QHBoxLayout()
        HL_algorithms.addStretch(1)
        HL_algorithms.addWidget(self.RBT_kernel)
        HL_algorithms.addStretch(1)
        HL_algorithms.addWidget(self.RBT_msf)
        HL_algorithms.addStretch(1)
        BOX_algorithms.setLayout(HL_algorithms)

        # Curve type
        self.GRP_events_curve_type = PYAFButtonGroup(self, "events_curve_type")
        self.RBT_curve = QtWidgets.QRadioButton("Curve")
        self.RBT_detection = QtWidgets.QRadioButton("Detection")
        self.GRP_events_curve_type.addButton(self.RBT_curve, 0)
        self.GRP_events_curve_type.addButton(self.RBT_detection, 1)

        BOX_events_curve_type = QtWidgets.QGroupBox("Display")
        HL_curve_type = QtWidgets.QHBoxLayout()
        HL_curve_type.addStretch(1)
        HL_curve_type.addWidget(self.RBT_curve)
        HL_curve_type.addStretch(1)
        HL_curve_type.addWidget(self.RBT_detection)
        HL_curve_type.addStretch(1)
        BOX_events_curve_type.setLayout(HL_curve_type)

        name = "widget_fit_param_jocs_events"
        W_fit_params_jocs_events = FitParamsWidget(self, name)

        self.BT_compute_events = PYAFButton(self, "get_events", "Compute")
        self.BT_compute_events.setFixedSize(100, 100)

        VL_button_compute_events = QtWidgets.QVBoxLayout()
        VL_button_compute_events.addStretch(1)
        VL_button_compute_events.addWidget(self.BT_compute_events)
        VL_button_compute_events.addStretch(1)

        HL_events_compute_type = QtWidgets.QHBoxLayout()
        HL_events_compute_type.addStretch(1)
        HL_events_compute_type.addWidget(BOX_events_curve_type)
        HL_events_compute_type.addWidget(BOX_algorithms)
        HL_events_compute_type.addStretch(1)

        self.GL_events = QtWidgets.QGridLayout()
        self.GL_events.addWidget(W_fit_params_jocs_events, 0, 1)

        VL_events_options = QtWidgets.QVBoxLayout()
        VL_events_options.addLayout(HL_events_compute_type)
        VL_events_options.addLayout(self.GL_events)

        HL_events.addLayout(VL_events_options)
        HL_events.addLayout(VL_button_compute_events)

        self.widget_events.setLayout(HL_events)

        widgets_list.widget_progressbar.update()

        # --- Loading rates ---------------------------------------------------

        self.BT_get_lr = PYAFButton(self, "get_lr", "Compute")
        self.BT_get_lr.setFixedSize(100, 100)

        label = "Loading Rate Coefficient"
        self.IN_lr_coef = PYAFInput(self, "input_lr_coef", label)
        self.IN_lr_coef.input.setValidator(misc_tools.validator("UF"))

        HL_loading_rates = QtWidgets.QHBoxLayout()

        VL_loading_rates = QtWidgets.QVBoxLayout()
        VL_loading_rates.addStretch(1)
        VL_loading_rates.addWidget(self.BT_get_lr)
        VL_loading_rates.addStretch(1)

        VL_coeff = QtWidgets.QVBoxLayout()
        VL_coeff.addWidget(self.IN_lr_coef)
        VL_coeff.addStretch(1)

        HL_loading_rates.addLayout(VL_coeff)
        HL_loading_rates.addStretch(1)
        HL_loading_rates.addLayout(VL_loading_rates)

        self.widget_loading_rates.setLayout(HL_loading_rates)

        # --- Stiffness height correction -------------------------------------

        self.BT_height_corr = PYAFButton(self, "height_correction", "Compute")
        self.BT_height_corr.setFixedSize(100, 100)

        # Instruction label
        self.label_info = QtWidgets.QLabel("Either select the substrate from the ROI manager (results tab) or input"
                                           " the known height of the sample")

        # ROI selection (Label + Dropdown)
        self.label_list_roi = QtWidgets.QLabel("ROI of substrate:")
        self.list_roi_glass = PYAFComboBox(self, "roi")
        self.list_roi_glass.setFixedSize(100, 26)  # Fixed width to prevent stretching

        # Noise floor input (Label + Input Box)
        self.label_noise_floor = QtWidgets.QLabel("Noise floor [nm]:")  # Label for noise floor
        self.IN_noise_floor = PYAFInput(self, "noise_floor", "", width=50)  # Input box for noise floor
        self.IN_noise_floor.input.setValidator(misc_tools.validator("UF"))  # Validation

        # Sample height input (Switched Label & Input Box)
        self.label_sample_height = QtWidgets.QLabel("Sample height [nm]:")  # Separate label
        self.IN_user_h = PYAFInput(self, "user_h", "", width=50)  # Empty label in PYAFInput
        self.IN_user_h.input.setValidator(misc_tools.validator("UF"))

        # Layout for ROI selection (aligns label & dropdown closely)
        self.stiff_corr_row = QtWidgets.QHBoxLayout()
        self.stiff_corr_row.addWidget(self.label_list_roi)
        self.stiff_corr_row.addWidget(self.list_roi_glass)
        self.stiff_corr_row.addStretch(1)  # Push everything to the left

        # Layout for Noise Floor (Label First)
        self.stiff_corr_noise_row = QtWidgets.QHBoxLayout()
        self.stiff_corr_noise_row.addWidget(self.label_noise_floor)  # Label
        self.stiff_corr_noise_row.addWidget(self.IN_noise_floor)  # Input box
        self.stiff_corr_noise_row.addStretch(1)  # Push to the left

        # Layout for Sample Height (Label First)
        self.stiff_corr_height_row = QtWidgets.QHBoxLayout()
        self.stiff_corr_height_row.addWidget(self.label_sample_height)  # Label first
        self.stiff_corr_height_row.addWidget(self.IN_user_h)  # Input box second
        self.stiff_corr_height_row.addStretch(1)  # Push everything to the left

        # Vertical layout to stack components
        self.stiff_corr_opts = QtWidgets.QVBoxLayout()
        self.stiff_corr_opts.addWidget(self.label_info)  # Instruction text
        self.stiff_corr_opts.addLayout(self.stiff_corr_row)  # ROI label + dropdown
        self.stiff_corr_opts.addLayout(self.stiff_corr_noise_row)  # Add noise floor layout
        self.stiff_corr_opts.addLayout(self.stiff_corr_height_row)  # Sample height label + input

        VL_stiff_corr_opts = QtWidgets.QVBoxLayout()
        VL_stiff_corr_opts.addLayout(self.stiff_corr_opts)
        VL_stiff_corr_opts.addStretch(1)

        VL_stiff_corr = QtWidgets.QVBoxLayout()
        VL_stiff_corr.addStretch(1)
        VL_stiff_corr.addWidget(self.BT_height_corr)
        VL_stiff_corr.addStretch(1)

        HL_stiff_corr = QtWidgets.QHBoxLayout()
        HL_stiff_corr.addLayout(VL_stiff_corr_opts)
        HL_stiff_corr.addStretch(1)
        HL_stiff_corr.addLayout(VL_stiff_corr)

        self.widget_stiff_corr.setLayout(HL_stiff_corr)

        # --- Set up page -----------------------------------------------------

        self.tabs.addTab(self.widget_app, "Elasticity")
        self.tabs.addTab(self.widget_stiff_corr, "B.E.C.")
        self.tabs.addTab(self.widget_ret, "Work, Rupture force")
        self.tabs.addTab(self.widget_events, "Unbinding events")
        self.tabs.addTab(self.widget_loading_rates, "Loading rates")

        self.GL_compute = QtWidgets.QGridLayout()
        self.GL_compute.addWidget(self.MPL_meshgrid, 0, 0)
        self.GL_compute.addWidget(self.MPL_curve, 0, 1)
        self.GL_compute.addLayout(self.GL_choosers, 1, 0, 1, 0)
        self.GL_compute.addLayout(self.HL_options, 2, 0, 1, 0)
        self.GL_compute.addWidget(self.tabs, 3, 0, 1, 0)

        self.VL.addLayout(self.GL_compute)

        self.setLayout(self.VL)

        self.update_widget()
        widgets_list.widget_progressbar.update()

    def update_widget(self):
        """Updates the content of the widget."""
        data = shared.exp.current_data

        # For single files don't display the curve selector
        if data.is_single:
            self.curve_selector.setParent(None)
        else:
            self.GL_choosers.addWidget(self.curve_selector, 0, 1)
            self.curve_selector.update_widget()

        if data.events_algorithm == "kernel":
            self.widget_msf.setParent(None)
            self.GL_events.addWidget(self.widget_kernel, 0, 0)

        elif data.events_algorithm == "msf":
            self.widget_kernel.setParent(None)
            self.GL_events.addWidget(self.widget_msf, 0, 0)

        self.STK_stiffness_mode.setCurrentIndex(data.stiffness_mode)

        widgets_list.widget_fit_param_pocs.update_widget()
        widgets_list.widget_fit_param_jocs.update_widget()
        widgets_list.widget_fit_param_jocs_events.update_widget()
        widgets_list.widget_kernel.update_widget()
        widgets_list.widget_msf.update_widget()
        self.update_GUI("all")
        self.update_MPL("MPL_meshgrid")
        self.update_MPL("MPL_canvas")

    def tabs_changed(self):
        """When changing a tab, change the computation type.

        This will update the fitting preview on the plot.
        """
        # Reset zoomed in position
        shared.zoomed_in_element = 0

        if self.tabs.currentIndex() == 0:
            shared.exp.compute_type = "stiffness"
            self.BT_next_zoom_element.setEnabled(False)
        elif self.tabs.currentIndex() == 1:
            shared.exp.compute_type = "stiffness_corr"
            self.BT_next_zoom_element.setEnabled(False)
        elif self.tabs.currentIndex() == 2:
            shared.exp.compute_type = "work_and_rupture_force"
            self.BT_next_zoom_element.setEnabled(True)
        elif self.tabs.currentIndex() == 3:
            shared.exp.compute_type = "events"
            self.BT_next_zoom_element.setEnabled(True)
        elif self.tabs.currentIndex() == 4:
            shared.exp.compute_type = "loading_rates"
            self.BT_next_zoom_element.setEnabled(False)

        self.update_MPL("MPL_canvas")

    def update_MPL(self, canvas):
        """Update the plots."""
        data = shared.exp.current_data

        if canvas == "MPL_canvas":
            self.MPL_curve.update_plot()

        elif canvas == "MPL_meshgrid":
            if not data.is_single:
                self.MPL_meshgrid.update_plot()

    def update_GUI(self, element):
        """Update GUI elements."""
        data = shared.exp.current_data

        if element == "single_or_forcevolume" or element == "all":
            if data.is_single:
                self.MPL_meshgrid.setParent(None)
                self.BT_clean_up.setEnabled(False)
            else:
                self.GL_compute.addWidget(self.MPL_meshgrid, 0, 0)
                self.BT_clean_up.setEnabled(True)

        if element == "apply_to_all" or element == "all":
            self.CB_apply_to_all.setChecked(shared.exp.apply_to_all_compute)

        if element == "stiffness_model_selected" or element == "all":
            self.list_models.setCurrentIndex(data.stiffness_model_selected)

            self.update_MPL("MPL_canvas")
            self.update_GUI("buttons_compute")

        if element == "BT_height_correction" or element == "all":
            # In the case of the model with the slope, do not allow the
            # height correction computation.

            model = data.used_stiffness_model_selected
            roi = data.roi_glass_id
            user_h = data.user_h
            if data.stiffness_calculated and (roi != 0 or (user_h is not None and user_h != 0)) and model != 3:
                self.BT_height_corr.setEnabled(True)
            else:
                self.BT_height_corr.setEnabled(False)

        if element == "roi_list_height_corr" or element == "all":
            # Remove old items
            self.list_roi_glass.clear()

            # Add rois
            self.list_roi_glass.addItem("None")
            count = 1
            for roi in shared.exp.current_data.roi_list:
                if roi.glass_coeffs is not None:
                    self.list_roi_glass.addItem(str(roi.roi_id + 1))
                    # Select
                    if data.roi_glass_id == roi.roi_id + 1:
                        self.list_roi_glass.setCurrentIndex(count)
                    count += 1

        if element == "IN_user_h" or element == "all":
            self.IN_user_h.input.setText("None")
            if data.user_h is not None:
                self.IN_user_h.input.setText(str(data.user_h))

        if element == "get_lr" or element == "all":
            self.BT_get_lr.setEnabled(data.events_calculated)

        if element == "grid_model_options" or element == "all":
            if data.stiffness_model_selected == 0:
                self.IN_tip_angle.setParent(None)
                self.IN_tip_radius.input.setText(str(data.tip_radius))
                self.GL_stiffness_params.addWidget(self.IN_tip_radius, 2, 0)
                self.GL_stiffness_params.addWidget(self.IN_poisson, 1, 0)
                # Enable fit options
                self.RBT_fitting_params.setChecked(True)
                self.RBT_fitting_params.setEnabled(True)
            elif data.stiffness_model_selected == 1:
                self.IN_tip_radius.setParent(None)
                self.IN_tip_angle.input.setText(str(data.tip_angle))
                self.GL_stiffness_params.addWidget(self.IN_tip_angle, 2, 0)
                self.GL_stiffness_params.addWidget(self.IN_poisson, 1, 0)
                # Enable fit options
                self.RBT_fitting_params.setChecked(True)
                self.RBT_fitting_params.setEnabled(True)
            elif data.stiffness_model_selected == 2:
                self.IN_tip_radius.setParent(None)
                self.IN_tip_angle.input.setText(str(data.tip_angle))
                self.GL_stiffness_params.addWidget(self.IN_tip_angle, 2, 0)
                self.GL_stiffness_params.addWidget(self.IN_poisson, 1, 0)
                # Enable fit options
                self.RBT_fitting_params.setChecked(True)
                self.RBT_fitting_params.setEnabled(True)
            elif data.stiffness_model_selected == 3:
                self.IN_tip_radius.setParent(None)
                self.IN_tip_angle.setParent(None)
                self.IN_poisson.setParent(None)
                # Disable fit options
                self.RBT_fitting_params.setChecked(False)
                self.RBT_fitting_params.setEnabled(False)
            if data.stiffness_model_selected == 4:
                self.IN_tip_angle.setParent(None)
                self.IN_tip_radius.input.setText(str(data.tip_radius))
                self.GL_stiffness_params.addWidget(self.IN_tip_radius, 2, 0)
                self.GL_stiffness_params.addWidget(self.IN_poisson, 1, 0)
                # Enable fit options
                self.RBT_fitting_params.setChecked(True)
                self.RBT_fitting_params.setEnabled(True)

            self.IN_poisson.input.setText(str(data.poisson_ratio))

        # if element == "input_segments" or element == "all":
            # self.IN_start.input.setText(str(data.indentation_start))
            # self.IN_stop.input.setText(str(data.indentation_stop))
            # self.IN_step.input.setText(str(data.indentation_step))

        if element == "zoom_in" or element == "all":
            self.CB_zoom_in.setChecked(data.auto_zoom_in)

        if element == "zoom_fit_preview_factor" or element == "all":
            self.IN_zoom_in_factor.changeValue(data.zoom_fit_preview_factor)

        if element == "compute_all" or element == "all":
            self.CB_compute_all.setChecked(shared.exp.calc_all)

        if element == "flatten" or element == "all":
            if data.scan_size_x != 0:
                self.BT_flatten.setEnabled(True)
            else:
                self.BT_flatten.setEnabled(False)

        # if element == "strict_stop" or element == "all":
            # self.CB_strict_stop.setChecked(data.strict_stop)
            # Enable or disable strict stop mode
            # if data.indentation_stop != 0:
                # self.CB_strict_stop.setEnabled(True)
            # else:
                # self.CB_strict_stop.setEnabled(False)

        # if element == "tomography" or element == "all":
            # self.CB_tomography.setChecked(data.tomography)

        if element == "input_lr_coef" or element == "all":
            self.IN_lr_coef.changeValue(data.lr_coef)

        if element == "events_algorithms" or element == "all":
            if data.events_algorithm == "kernel":
                self.RBT_kernel.setChecked(True)
            elif data.events_algorithm == "msf":
                self.RBT_msf.setChecked(True)

        if element == "stiffness_params" or element == "all":
            if data.stiffness_mode == 0:
                self.RBT_calculus_params.setChecked(True)
                data.perform_fit = False
            elif data.stiffness_mode == 1:
                self.RBT_fitting_params.setChecked(True)
                data.perform_fit = True

        if element == "events_curve_type" or element == "all":
            if data.events_curve_type == "curve":
                self.RBT_curve.setChecked(True)
            elif data.events_curve_type == "detection":
                self.RBT_detection.setChecked(True)

    def list_updated(self, name):
        """Called when a list is updated."""
        data = shared.exp.current_data

        if name == "list_models":
            # This list lets you chose the model for the stiffness computation.
            data.stiffness_model_selected = self.list_models.currentIndex()
            # Update the GUI
            self.update_GUI("grid_model_options")

        elif name == "roi":
            # Get the value from the text in the list, not the id, because
            # we can not rely on the id here.
            value = str(self.list_roi_glass.currentText())
            if value == "None":
                value = 0
                self.IN_noise_floor.input.setText("None")
                data.noise_floor = None
                self.IN_noise_floor.setEnabled(False)
                self.IN_user_h.setEnabled(True)
            else:
                value = int(value)
                self.IN_user_h.input.setText("None")
                data.user_h = None
                self.IN_user_h.setEnabled(False)
                self.IN_noise_floor.setEnabled(True)

            data.roi_glass_id = value

            self.update_GUI("BT_height_correction")

        elif name == "list2":
            self.parent.file_changed("list2")

    def input_updated(self, field):
        """Called when an input is updated."""
        data = shared.exp.list[shared.exp.id_selected]

        if field == "tip_radius":
            if self.IN_tip_radius.input.text() != "":
                data.tip_radius = self.IN_tip_radius.get_float_value()

        elif field == "tip_angle":
            if self.IN_tip_angle.input.text() != "":
                data.tip_angle = self.IN_tip_angle.get_float_value()

        elif field == "poisson":
            if self.IN_poisson.input.text() != "":
                data.poisson_ratio = self.IN_poisson.get_float_value()

        # elif field == "input_start":
            # if self.IN_start.input.text() != "":
                # data.indentation_start = self.IN_start.get_int_value()

        # elif field == "input_stop":
            # if self.IN_stop.input.text() != "":
                # data.indentation_stop = self.IN_stop.get_int_value()
                # self.update_GUI("strict_stop")

        # elif field == "input_step":
            # if self.IN_step.input.text() != "":
                # data.indentation_step = self.IN_step.get_int_value()

        elif field == "input_lr_coef":
            data.lr_coef = self.IN_lr_coef.get_float_value()

        elif field == "zoom_in_factor":
            val = self.IN_zoom_in_factor.get_float_value()
            data.zoom_fit_preview_factor = val

        elif field == "noise_floor":
            if self.IN_noise_floor.input.text() != "":
                data.noise_floor = self.IN_noise_floor.get_float_value()

        elif field == "user_h":
            if self.IN_user_h.input.text() != "":
                data.user_h = self.IN_user_h.get_float_value()

        self.update_MPL("MPL_canvas")
        self.update_GUI("BT_height_correction")

    def finish_computation(self, warnings=None):
        """At the end of the computation, display the errors and refresh the GUI."""
        if warnings is not None:
            stop_size_errors = warnings[0]
            no_segment_errors = warnings[1]
            smoothing_warnings = warnings[2]
            strict_stop_errors = warnings[3]

            found_error_list = []
            for i in range(len(stop_size_errors)):
                # True means no error, so some data was saved
                # Set to False at the begining, there should be at least one
                # True at the end
                found_error_list.append([False])

            list_stop_size_error = ""
            for i in range(len(stop_size_errors)):
                item = stop_size_errors[i]
                if item:
                    shared.exp.list[item[0]].stiffness_calculated = False
                    name = shared.exp.list[item[0]].filename
                    dist = item[1]
                    text = "<li>" + name + " : " + dist + " nm</li>"
                    list_stop_size_error += text
                else:
                    found_error_list[i] = True

            list_no_segment_errors = ""
            for i in range(len(no_segment_errors)):
                item = no_segment_errors[i]
                if item:
                    shared.exp.list[item[0]].stiffness_calculated = False
                    name = shared.exp.list[item[0]].filename
                    text = "<li>" + name + "</li>"
                    list_no_segment_errors += text
                else:
                    found_error_list[i] = True

            list_smoothing_warnings = ""
            text = ""
            for item in smoothing_warnings:
                if item:
                    name = shared.exp.list[item[0]].filename
                    for values in item[1]:
                        # values[0] is set to "smoothing_error"
                        # the curves are stored in values[1]
                        text += "<li>" + name + " : " + values[1] + "</li>"
                    list_smoothing_warnings += text

            list_strict_stop_error = ""
            for i in range(len(strict_stop_errors)):
                item = strict_stop_errors[i]
                if item:
                    shared.exp.list[item[0]].stiffness_calculated = False
                    name = shared.exp.list[item[0]].filename

                    text = "<li>" + name + "</li>"
                    list_strict_stop_error += text
                else:
                    found_error_list[i] = True

            if list_stop_size_error != "":
                text = "The 'stop' value is too big ! The maximum value you \
                        can use is : <ul>" + list_stop_size_error + "</ul>"

            if list_no_segment_errors != "":
                text = "Please chose a smaller step or check your \
                        start/stop values ! No segment was found. Files : \
                        <ul>" + list_no_segment_errors + "</ul>"

            if list_smoothing_warnings != "":
                text = "Some curves seem corrupted and the smoothing was not \
                        applied for these curves.<br><br>\
                        The following curves had a problem : <ul>" + \
                    list_smoothing_warnings + "</ul>"

            if list_strict_stop_error != "":
                text = "Strict stop was checked but you defined a stop value \
                        which goes beyond the last point of some curves.<br>\
                        <br>The following files had a problem : <ul>" + \
                    list_strict_stop_error + "</ul>"

            if (list_stop_size_error != "" or list_no_segment_errors != ""
                    or list_smoothing_warnings != ""
                    or list_strict_stop_error != "") \
                    and not consts.UNIT_TESTING:
                # Create a message box and display it
                msg = QtWidgets.QMessageBox()
                msg.setText("Error")
                msg.setInformativeText(text)
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
                msg.exec_()

        else:
            # No errors
            found_error_list = [True]

        # Check if we can update the data
        ok = False
        for item in found_error_list:
            if item:
                ok = True
                break

        # FINAL UPDATES (for the 3 types of calcs)
        if ok:
            # Update some values
            for data in shared.exp.list:
                data.update()

                # Check for special case where a stiffness was computed,
                # there was then an error during a second computation,
                # then we will display only the piezo meshgird.
                mtype = data.meshgrid_type == "topography" or \
                    data.meshgrid_type == "stiffness"
                if mtype and not data.stiffness_calculated:
                    data.meshgrid_type = "piezo"

            # Update some values
            shared.exp.update_global_values()
            # Propage color options to the other files if needed
            apply_to_all.apply_to_all("display_options")
            # Update Gui
            widgets_list.widget_fit_param_pocs.update_widget()
            widgets_list.widget_fit_param_jocs.update_widget()
            widgets_list.widget_fit_param_jocs_events.update_widget()
            self.update_GUI("BT_height_correction")
            self.update_GUI("get_lr")
            # Update the results tab (especially the plot, so that it is
            # initialized and not empty)
            widgets_list.widget_results.update_widget()
            # Update the table with the data
            stat_tools.get_values()
            widgets_list.widget_results_single.reset_table()
            stat_tools.fetch_group_data()
            stat_tools.fetch_conditions_data()
            widgets_list.widget_results_single.update_widget()
            widgets_list.widget_results_groups.update_widget()
            widgets_list.widget_results_experiment.update_widget()
            self.parent.update_GUI("tabs")
            # Update widgets
            if widgets_list.widget_details is not None:
                widgets_list.widget_details.update_widget()
            if widgets_list.widget_multimeshgrids is not None:
                widgets_list.widget_multimeshgrids.update_widget()
            if widgets_list.widget_meshgrid_options is not None:
                widgets_list.widget_meshgrid_options.update_widget()
            if widgets_list.widget_events_per_scan is not None:
                widgets_list.widget_events_per_scan.update_widget()
            if widgets_list.widget_used_parameters is not None:
                widgets_list.widget_used_parameters.update_widget()

            if not consts.UNIT_TESTING and not consts.AUTOCALC_STIFFNESS:
                QtWidgets.QApplication.beep()
                text = "Calculation done."
                QtWidgets.QMessageBox.warning(self, "Information", text)

    def button_clicked(self, button):
        """Method called whenever a button is clicked."""
        # Set the focus so that all the inputs which have been modified but
        # not validated will save their values.
        self.setFocus()

        data = shared.exp.list[shared.exp.id_selected]

        if button == "button_rand":
            widgets_list.widget_main.change_curve(0, 0, "rand")

        elif button == "back":
            lt = shared.exp.last_ten_curves
            lt -= 1
            posx = lt[shared.exp.pos_in_last_ten_curves][0]
            posy = lt[shared.exp.pos_in_last_ten_curves][1]

            widgets_list.widget_main.change_curve(
                posx, posy, save_in_last_ten=False)
            self.update_GUI("bts_back_forward")

        elif button == "forward":
            lt = shared.exp.last_ten_curves
            lt += 1
            posx = lt[shared.exp.pos_in_last_ten_curves][0]
            posy = lt[shared.exp.pos_in_last_ten_curves][1]
            widgets_list.widget_main.change_curve(
                posx, posy, save_in_last_ten=False)
            self.update_GUI("bts_back_forward")

        elif button == "get_stiffness":
            # Update all the parameters if needed
            apply_to_all.apply_to_all("calc_options")
            smoothing_warnings = []
            stop_size_errors = []
            no_segment_errors = []
            strict_stop_errors = []
            one_computed = False

            # Make a list of data to treat
            calc_list = []
            if shared.exp.calc_all:
                for i in range(len(shared.exp.list)):
                    calc_list.append(i)
            else:
                calc_list.append(shared.exp.id_selected)

            # Overwrite calc_list with a list defined externally
            if shared.force_calc_list is not None:
                calc_list = shared.force_calc_list

            for calc_id in calc_list:
                data = shared.exp.list[calc_id]
                self.logger.debug("Computing %s", data.filename)

                self.compute_stiffness.get_stiffness(calc_id)

                if self.compute_stiffness.error == "Error: Stop size":
                    val = str(int(self.compute_stiffness.error_value))
                    stop_size_errors.append([calc_id, val])
                    no_segment_errors.append([])
                    smoothing_warnings.append([])
                    strict_stop_errors.append([])

                elif self.compute_stiffness.error == "Error: no segments":
                    no_segment_errors.append([calc_id, None])
                    stop_size_errors.append([])
                    smoothing_warnings.append([])
                    strict_stop_errors.append([])

                elif self.compute_stiffness.error == "Error: strict stop":
                    strict_stop_errors.append([calc_id])
                    no_segment_errors.append([])
                    stop_size_errors.append([])
                    smoothing_warnings.append([])

                else:
                    # At least one file was computed without error
                    one_computed = True

                    stop_size_errors.append([])
                    no_segment_errors.append([])
                    strict_stop_errors.append([])

                    # Store smoothing warnings to display a message at the end
                    # of the computation
                    warnings = self.compute_stiffness.list_of_smoothing_errors
                    if warnings:
                        smoothing_warnings.append([calc_id, warnings])
                    else:
                        smoothing_warnings.append([])

                    data.max_indentation_index = \
                        self.compute_stiffness.max_indentation_index

                    # Define what is now displayed as meshgrid
                    data.display_topography = 1
                    data.display_stiffness = 0

                    # Add a result in the results tab if needed
                    name = shared.exp.list[calc_id].filename
                    experiment.add_result(calc_id, name)

                    # Update the indendation options
                    misc_tools.update_indentation_list(calc_id)
                    misc_tools.update_slice_position(calc_id)

                # Reset the error
                self.compute_stiffness.error = False

            if one_computed:
                if shared.exp.results_type is None:
                    shared.exp.results_type = "stiffness"

            warnings = [stop_size_errors,
                        no_segment_errors,
                        smoothing_warnings,
                        strict_stop_errors]

            self.finish_computation(warnings)

        elif button == "get_work_and_rupture_force":
            # Update all the parameters if needed
            apply_to_all.apply_to_all("calc_options")

            calc_list = []
            if shared.exp.calc_all:
                for i in range(len(shared.exp.list)):
                    calc_list.append(i)
            else:
                calc_list.append(shared.exp.id_selected)

            # Overwrite calc_list with a list defined externally
            if shared.force_calc_list is not None:
                calc_list = shared.force_calc_list

            for calc_id in calc_list:
                data = shared.exp.list[calc_id]
                self.logger.debug("Computing %s", data.filename)
                comp = self.compute_work_and_rupture_force
                comp.get_work_and_rupture_force(calc_id)

                # Add a result in the result's tab if needed
                name = shared.exp.list[calc_id].filename
                experiment.add_result(calc_id, name)

                if shared.exp.results_type is None:
                    shared.exp.results_type = "work"

            self.finish_computation()

        elif button == "get_events":
            # Update all the parameters if needed
            apply_to_all.apply_to_all("calc_options")

            # Create a list of the datasets
            calc_list = []
            if shared.exp.calc_all:
                for i in range(len(shared.exp.list)):
                    calc_list.append(i)
            else:
                calc_list.append(shared.exp.id_selected)

            # Overwrite calc_list with a list defined externally
            if shared.force_calc_list is not None:
                calc_list = shared.force_calc_list

            for calc_id in calc_list:
                data = shared.exp.list[calc_id]
                self.logger.debug("Computing %s", data.filename)
                self.compute_events.compute(calc_id)

                # Add a result to the result's tab if needed
                name = shared.exp.list[calc_id].filename
                experiment.add_result(calc_id, name)

                if shared.exp.results_type is None:
                    shared.exp.results_type = "events_forces"

            self.finish_computation()

        elif button == "height_correction":
            # Update all the parameters if needed
            apply_to_all.apply_to_all("calc_options")
            calc_id = shared.exp.current_data.unique_id

            if shared.exp.current_data.user_h is None or shared.exp.current_data.user_h == 0:
                # Height correction works only for a single file
                self.stiffness_correction.correct_stiffness(calc_id)

            else:
                # Make a list of data to treat
                calc_list = []
                if shared.exp.calc_all:
                    for i in range(len(shared.exp.list)):
                        calc_list.append(i)
                else:
                    calc_list.append(shared.exp.id_selected)

                # Overwrite calc_list with a list defined externally
                if shared.force_calc_list is not None:
                    calc_list = shared.force_calc_list

                for calc_id in calc_list:
                    data = shared.exp.list[calc_id]
                    if data.stiffness_calculated:
                        self.logger.debug("Correcting Elasticity %s", data.filename)
                        self.stiffness_correction.correct_stiffness(calc_id)

            self.finish_computation()

        elif button == "get_lr":
            # Update all the parameters if needed
            apply_to_all.apply_to_all("calc_options")

            calc_list = []
            if shared.exp.calc_all:
                for i in range(len(shared.exp.list)):
                    if shared.exp.list[i].events_calculated:
                        calc_list.append(i)
            else:
                calc_list.append(shared.exp.id_selected)

            # Overwrite calc_list with a list defined externally
            if shared.force_calc_list is not None:
                calc_list = shared.force_calc_list

            for calc_id in calc_list:
                data = shared.exp.list[calc_id]
                self.logger.debug("Computing %s", data.filename)
                self.get_lr.compute(calc_id)

            self.finish_computation()

        elif button == "curve_mod":
            if widgets_list.widget_curve_mod is None:
                # Create new widget
                CurveModWidget(self)
                widgets_list.widget_curve_mod.setWindowTitle("Modify curves")
                widgets_list.widget_curve_mod.show()
            else:
                # Bring to front
                widgets_list.widget_curve_mod.activateWindow()
                widgets_list.widget_curve_mod.raise_()

        elif button == "flatten":
            if widgets_list.widget_flatten is None:
                # Create new widget
                FlattenWidget(self)
                widgets_list.widget_flatten.setWindowTitle("Flatten")
                widgets_list.widget_flatten.show()
            else:
                # Bring to front
                widgets_list.widget_flatten.activateWindow()
                widgets_list.widget_flatten.raise_()

        elif button == "clean_up":
            if widgets_list.widget_clean_up is None:
                # Create new widget
                CleanUpWidget(self)
                widgets_list.widget_clean_up.setWindowTitle("Clean Up")
                widgets_list.widget_clean_up.show()
            else:
                # Bring to front
                widgets_list.widget_clean_up.activateWindow()
                widgets_list.widget_clean_up.raise_()

        elif button == "events_algorithms":
            value = self.GRP_algorithms.checkedId()
            if value == 0:
                data.events_algorithm = "kernel"
            elif value == 1:
                data.events_algorithm = "msf"
            self.update_widget()

        elif button == "events_curve_type":
            value = self.GRP_events_curve_type.checkedId()
            if value == 0:
                data.events_curve_type = "curve"
            elif value == 1:
                data.events_curve_type = "detection"
            self.update_widget()

        elif button == "stiffness_params":
            value = self.GRP_stiffness_params.checkedId()
            data.stiffness_mode = value
            if value == 0:
                data.perform_fit = False
            elif value == 1:
                data.perform_fit = True
            self.update_widget()

        elif button == "next_zoom_element":
            if shared.exp.compute_type == "stiffness":
                # Only one element, the POC
                val = None
            elif shared.exp.compute_type == "stiffness_corr":
                # No element to display here
                val = None
            elif shared.exp.compute_type == "work_and_rupture_force":
                # Two elements, joc1 and joc2
                val = 2
            elif shared.exp.compute_type == "events":
                # Multiple elements, + 1 (the joc)
                if data.events_curve_type == "curve":
                    val = len(self.MPL_curve.canvas.events_preview) + 1
                else:
                    val = None
            elif shared.exp.compute_type == "loading_rates":
                # No element to display here
                val = None

            if val is not None:
                if shared.zoomed_in_element + 1 < val:
                    # Go to next element
                    shared.zoomed_in_element += 1
                else:
                    # Go back to first element
                    shared.zoomed_in_element = 0

            self.update_MPL("MPL_canvas")

    def checkbox_clicked(self, checkbox):
        """Called when a checkbox is clicked."""
        data = shared.exp.list[shared.exp.id_selected]

        if checkbox == "compute_all":
            shared.exp.calc_all = self.CB_compute_all.isChecked()

        elif checkbox == "apply_to_all":
            # Ask the user if he is sure he wants to overwrite his parameters
            do = True
            value = self.CB_apply_to_all.isChecked()
            if value:
                text = "Are you sure you want to overwrite all the parameters \
                        with the current ones ?"
                msg = QtWidgets.QMessageBox()
                msg.setText("Apply parameters to all ?")
                msg.setInformativeText(text)
                msg.setIcon(QtWidgets.QMessageBox.Information)
                msg.addButton(QtWidgets.QMessageBox.Yes)
                msg.addButton(QtWidgets.QMessageBox.No)
                ret = msg.exec_()

                if ret == QtWidgets.QMessageBox.No:
                    self.CB_apply_to_all.setChecked(False)
                    do = False

            if do:
                # If we uncheck the checkbox, apply to all before unchecking
                # it, so that the last configuration is saved. Will do nothing
                # if we are checking the box as here apply_to_all_compute is
                # still False.
                apply_to_all.apply_to_all("calc_options")

                shared.exp.apply_to_all_compute = value

                # If we check the box, apply to all the current configuration
                apply_to_all.apply_to_all("calc_options")

        elif checkbox == "strict_stop":
            state = self.CB_strict_stop.isChecked()
            shared.exp.list[shared.exp.id_selected].strict_stop = state

        # elif checkbox == "tomography":
            # state = self.CB_tomography.isChecked()
            # shared.exp.list[shared.exp.id_selected].tomography = state

        elif checkbox == "zoom_in":
            data.auto_zoom_in = self.CB_zoom_in.isChecked()
            tp = shared.exp.compute_type
            if tp == "events" or tp == "work_and_rupture_force":
                self.BT_next_zoom_element.setEnabled(data.auto_zoom_in)
            self.update_MPL("MPL_canvas")

        elif checkbox == "perform_fit":
            data.perform_fit = self.CB_perform_fit.isChecked()
            if self.CB_perform_fit.isChecked():
                pass
                # Hide elasticity segment calculation
                # self.IN_start.setEnabled(False)
                # self.IN_stop.setEnabled(False)
                # self.IN_step.setEnabled(False)
                # self.CB_tomography.setChecked(False)
                # self.CB_tomography.setEnabled(False)
                # self.CB_strict_stop.setChecked(False)
                # self.CB_strict_stop.setEnabled(False)

            else:
                pass
                # self.IN_start.setEnabled(True)
                # self.IN_stop.setEnabled(True)
                # self.IN_step.setEnabled(True)
                # self.CB_tomography.setEnabled(True)
                # self.CB_strict_stop.setEnabled(True)



    def mesh_canvas_press(self, event):
        """Matplotlib method, called when the user clicks on the meshgrid."""
        # event.button == 1 = left click

        if event.button == 1 and event.inaxes == self.MPL_meshgrid.canvas.axes:
            xpos, ypos = misc_tools.get_position_on_meshgrid(event)

            widgets_list.widget_main.change_curve(xpos, ypos)
