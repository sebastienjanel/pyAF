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

import os
import numpy
from .. import shared
from .. import widgets_list
from .. import experiment
import logging
from .. import consts
from qtpy import QtCore, QtWidgets
from ..tools import misc_tools
from ..tools import gui_tools
from ..tools import apply_to_all
from ..tools import stat_tools
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFToggle
from ..widgets.load_parameters import load_parameters
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
from ..compute_tools import breakthrough
from ..compute_tools.stiffness import brush
from ..compute_tools.stiffness.main import Stiffness
from ..compute_tools.work_and_rupture_force.main import WorkAndRuptureForce
from ..compute_tools.stiffness_correction import StiffnessCorrection
from ..compute_tools.viscoelasticity import main as visco_main
from ..compute_tools import microrheology
from ..compute_tools import viscous_drag
from ..widgets.viscous_drag import ViscousDragDialog
from ..widgets.piezo_characterization import PiezoCharacterizationDialog
from ..compute_tools.stiffness import contact_models
from ..compute_tools.stiffness import bec
from ..widgets import plot_options
from ..compute_tools.loading_rates import GetLoadingRate
from ..compute_tools.events import GetEvents
from ..tools.PYAFWidget import PYAFWidget


def file_label(data):
    """Short name of a file in the messages (without the end added by the
    AFM software, as in the plots)."""
    return experiment.auto_short_name(os.path.basename(str(data.filename)))


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
        self.compute_visco = visco_main.ViscoCompute()
        self.compute_events = GetEvents()
        self.get_lr = GetLoadingRate()

        self.canvas_resolution = self.parent.canvas_resolution
        # Height of the plots: the height of the main window is set by
        # this tab (the boxes of the tabs below take the rest)
        self.curve_canvas_size = [800, 299]
        self.meshgrid_canvas_size = [300, 299]

        # --- Main layout -----------------------------------------------------

        self.VL = QtWidgets.QVBoxLayout()

        # --- Experiment ------------------------------------------------------
        self.GL_choosers = QtWidgets.QGridLayout()
        self.GL_choosers.setContentsMargins(8, 0, 8, 2)
        self.GL_choosers.setHorizontalSpacing(4)

        self.HL_choser = QtWidgets.QHBoxLayout()
        self.list_exp2 = PYAFComboBox(self, "list2")
        self.list_exp2.setFixedSize(374, 26)
        # The min-width of the application style sheet overrides the fixed
        # width: set it for this list (width of the text area, without the
        # padding and the drop-down arrow)
        self.list_exp2.setStyleSheet("QComboBox { min-width: 335px; }")
        for data in shared.exp.list:
            self.list_exp2.addItem(data.filename)
        self.HL_choser.addWidget(self.list_exp2)

        name = "widget_curve_selector_compute"
        self.curve_selector = ResultsSelectorWidget(self, name)
        # self.curve_selector is added to the layout in update_widget().

        self.GL_choosers.addLayout(self.HL_choser, 0, 0)


        widgets_list.widget_progressbar.update()

        self.TG_zoom_in = PYAFToggle(self, "zoom_in", "Zoom in")
        self.IN_zoom_in_factor = PYAFInput(
            self, "zoom_in_factor", "factor", 40)
        self.BT_next_zoom_element = PYAFButton(
            self, "next_zoom_element", ">", size=40)
        self.BT_next_zoom_element.setToolTip(
            "Zoom on the next element of the curve (jump of contact, events, "
            "breakthroughs)")

        widgets_list.widget_progressbar.update()

        # Canvas Meshgrid
        sizes = [self.meshgrid_canvas_size[0] // self.canvas_resolution,
                 self.meshgrid_canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]

        self.W_mesh = QtWidgets.QWidget()
        label = "small_meshgrid_compute"
        self.MPL_meshgrid = PYAFPlot(self, label, self.W_mesh, sizes)
        self.MPL_meshgrid.mpl_connect("button_press_event",
                                      self.mesh_canvas_press)
        self.MPL_meshgrid.tab_id = 1

        # Canvas curves
        sizes = [self.curve_canvas_size[0] // self.canvas_resolution,
                 self.curve_canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]
        self.W_curve = QtWidgets.QWidget()
        self.MPL_curve = PYAFPlot(self, "fit_preview", self.W_curve, sizes)

        widgets_list.widget_progressbar.update()

        # Calc All
        self.TG_compute_all = PYAFToggle(self, "compute_all", "Compute all")
        label = "Apply to all"
        self.TG_apply_to_all = PYAFToggle(self, "apply_to_all", label)

        self.BT_flatten = PYAFButton(self, "flatten", "Flatten", size=80)
        self.BT_curve_mod = PYAFButton(
            self, "curve_mod", "Tilting/Smoothing", size=150)
        self.BT_clean_up = PYAFButton(self, "clean_up", "Clean Up", size=95)

        # Map navigation: file, flatten, curve selector (added in
        # update_widget)
        self.GL_choosers.addWidget(self.BT_flatten, 0, 1)
        self.box_navigation = QtWidgets.QGroupBox("Map navigation")
        self.box_navigation.setLayout(self.GL_choosers)

        # Curves: modifications, clean up and zoom
        self.HL_curves = QtWidgets.QHBoxLayout()
        self.HL_curves.setContentsMargins(8, 0, 8, 2)
        self.HL_curves.addWidget(self.BT_clean_up)
        self.HL_curves.addWidget(self.BT_curve_mod)
        self.HL_curves.addSpacing(4)
        self.HL_curves.addWidget(self.TG_zoom_in)
        self.HL_curves.addWidget(self.IN_zoom_in_factor)
        self.HL_curves.addWidget(self.BT_next_zoom_element)
        self.box_curves = QtWidgets.QGroupBox("Curves")
        self.box_curves.setLayout(self.HL_curves)

        # Both boxes at the top of the tab
        self.HL_top = QtWidgets.QHBoxLayout()
        # (they share the width of the tab, aligned with the panels below)
        self.HL_top.addWidget(self.box_navigation)
        self.HL_top.addWidget(self.box_curves)

        # Options of the computation
        # Right aligned, above the compute buttons of the tabs
        # Load parameters on the left, below the maps, as high as the toggles
        # (no margin and vertical padding of the theme)
        self.BT_load_params = PYAFButton(
            self, "load_parameters", "Load parameters", size=120)
        self.BT_load_params.setToolTip(
            "Load parameters from previously exported JSON parameter file")
        self.BT_load_params.setStyleSheet(
            "QPushButton { margin: 0px; padding: 0px 1.2ex; }")
        metrics = self.BT_load_params.fontMetrics()
        self.BT_load_params.setFixedHeight(max(
            self.TG_apply_to_all.sizeHint().height(), metrics.height() + 8))
        self.HL_options = QtWidgets.QHBoxLayout()
        self.HL_options.addWidget(self.BT_load_params)
        self.HL_options.addStretch(1)
        self.HL_options.addWidget(self.TG_apply_to_all)
        self.HL_options.addWidget(self.TG_compute_all)

        widgets_list.widget_progressbar.update()

        # Tabs
        self.tabs = QtWidgets.QTabWidget()
        # Styled as pills by the stylesheets
        self.tabs.setObjectName("subTabs")
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
        # Same width for the three inputs of the model parameters
        self.IN_tip_angle = PYAFInput(self, "tip_angle", label, 60)
        self.IN_poisson = PYAFInput(self, "poisson", "Poisson ratio", 60)
        self.IN_shell_thickness = PYAFInput(
            self, "shell_thickness", "Shell thickness [nm]", 60)
        self.IN_shell_thickness.setToolTip(
            "Thickness h of the shell of the vesicle (lipid bilayer: about "
            "4 nm).\nThin shell (Reissner): F = 4 E h\u00b2 / (R "
            "\u221a(3 (1 - \u03bd\u00b2))) \u03b4, valid for an indentation "
            "small compared to the radius R of the vesicle\n(fit range: "
            "indentation up to about 0.1 R, Vorselen et al. 2020), with a "
            "tip sharper than the vesicle, on its top.")
        self.IN_tip_angle.setToolTip(
            "Half-angle between the axis of the tip and its face (for a "
            "pyramid: the face, not the edge).")
        self.IN_tip_radius.setToolTip(
            "Radius of the apex of the tip (sphere, paraboloid, blunted "
            "tips), of the flat punch, or of the cylinder of the needle "
            "and the nanowire.\nThin shell: radius of the vesicle (not of "
            "the tip).")
        # 0 is refused (division by zero in the models)
        self.IN_tip_radius.input.setValidator(misc_tools.validator("PF"))
        self.IN_tip_angle.input.setValidator(misc_tools.validator("PF"))
        self.IN_shell_thickness.input.setValidator(
            misc_tools.validator("PF"))
        self.IN_poisson.input.setValidator(misc_tools.validator("UF"))

        # Model
        self.HL_list_stiff = QtWidgets.QHBoxLayout()
        self.label_models_list = QtWidgets.QLabel("Model : ")
        self.list_models = PYAFComboBox(self, "list_models")
        self.models = contact_models.MODEL_NAMES
        for model in self.models:
            self.list_models.addItem(str(model))
        # Wide enough to display the longest model name. Set with a style
        # sheet: the min-width of the application style sheet would override
        # setMinimumWidth (here min-width is the width of the text area)
        metrics = self.list_models.fontMetrics()
        longest = max(metrics.horizontalAdvance(m) for m in self.models)
        self.list_models.setStyleSheet(
            "QComboBox { min-width: " + str(longest + 20) + "px; }")
        # The list that opens is also wide enough for the names
        self.list_models.view().setMinimumWidth(longest + 80)

        self.HL_list_stiff.addWidget(self.label_models_list)
        self.HL_list_stiff.addWidget(self.list_models)
        self.HL_list_stiff.addStretch(1)

        VL_empty1 = QtWidgets.QVBoxLayout()
        VL_empty1.addStretch(1)

        self.GL_stiffness_params = QtWidgets.QGridLayout()
        self.GL_stiffness_params.addLayout(self.HL_list_stiff, 0, 0)
        self.GL_stiffness_params.addWidget(self.IN_poisson, 1, 0)
        # Position 2, defined in update_widget
        # Brush model (paraboloid and sphere), rows 4 and 5, defined in
        # update_GUI
        self.TG_brush = PYAFToggle(self, "brush_model", "Brush model")
        self.TG_brush.setToolTip(
            "Brush model of the cells (Sokolov et al. 2007, Makarova et al. "
            "2023), with a spherical probe:\nthe cell body is fitted near "
            "the maximum load (brush squeezed), then the force of the brush "
            "below it,\nF = 100 kT R N^1.5 L exp(-2\u03c0 h / L) for 0.1 < "
            "h / L < 0.8: length L and grafting density N of the brush\n"
            "(effective chains per \u00b5m\u00b2). Computed with the "
            "elasticity (Fitting), "
            "maps in the Results tab. Zero of the force: baseline of the "
            "point of contact,\nwhose refit is set to 0 (it would move the "
            "baseline into the brush).")
        self.IN_brush_threshold_pct = PYAFInput(
            self, "brush_threshold_pct", "Body fit above [% max force]", 60)
        self.IN_brush_threshold_pct.setToolTip(
            "The brush is squeezed near the maximum load: the cell body is "
            "fitted above this percentage of the maximum force of each "
            "curve,\nthe brush below it (dotted line on the curve)")
        self.IN_brush_threshold_pct.input.setValidator(
            misc_tools.validator("RF", bottom=1, top=99, decimal=1))
        self.GL_stiffness_params.addLayout(VL_empty1, 6, 0)

        self.Box_stiff_exp_params.setLayout(self.GL_stiffness_params)
        VL_stiff_exp_params = QtWidgets.QVBoxLayout()
        VL_stiff_exp_params.addWidget(self.Box_stiff_exp_params)
        W_stiff_exp_params.setLayout(VL_stiff_exp_params)

        W_fit_params_pocs = FitParamsWidget(self, "widget_fit_param_pocs")

        # Fit parameters stiffness
        self.VL_fit_params_stiff = QtWidgets.QVBoxLayout()

        self.GRP_stiffness_params = PYAFButtonGroup(self, "stiffness_params")
        self.RBT_calculus_params = QtWidgets.QRadioButton("Tomography")
        self.RBT_calculus_params.setToolTip(
            "Stiffness tomography: slope of each segment of the indentation "
            "(Slope model)")
        self.RBT_fitting_params = QtWidgets.QRadioButton("Fitting")
        self.RBT_fitting_params.setToolTip(
            "Fit of the elasticity on the indentation (Hertz type models)")
        self.GRP_stiffness_params.addButton(self.RBT_calculus_params, 0)
        self.GRP_stiffness_params.addButton(self.RBT_fitting_params, 1)

        BOX_stiffness_mode = QtWidgets.QGroupBox("Computation method")
        HL_stiffness_mode = QtWidgets.QHBoxLayout()
        HL_stiffness_mode.addStretch(1)
        HL_stiffness_mode.addWidget(self.RBT_fitting_params)
        HL_stiffness_mode.addStretch(1)
        HL_stiffness_mode.addWidget(self.RBT_calculus_params)
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

        self.HL_widget_app.addWidget(W_fit_params_pocs)
        self.HL_widget_app.addWidget(W_stiff_exp_params)
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

        # Box centered in the tab (the spacing balances the compute button)
        self.HL_widget_ret.addSpacing(
            self.BT_compute_work_and_rupture_force.width())
        self.HL_widget_ret.addStretch(1)
        self.HL_widget_ret.addWidget(W_fit_params_jocs)
        self.HL_widget_ret.addStretch(1)
        self.HL_widget_ret.addLayout(VL_button_compute_work)

        self.widget_ret.setLayout(self.HL_widget_ret)

        widgets_list.widget_progressbar.update()

        # Events --------------------------------------------------------------

        HL_events = QtWidgets.QHBoxLayout()

        self.widget_kernel = EventsParamsWidget(self, "widget_kernel")
        self.widget_msf = EventsParamsWidget(self, "widget_msf")

        # Method of detection and display of the preview: selectors
        self.label_events_algorithm = QtWidgets.QLabel("Method :")
        self.list_events_algorithm = PYAFComboBox(self, "events_algorithm")
        self.list_events_algorithm.addItem("Kernel", "kernel")
        self.list_events_algorithm.addItem("MSF", "msf")
        self.list_events_algorithm.setToolTip(
            "Kernel (default): convolution of the force curve with a short "
            "ramp (local slope), threshold in multiples of the noise of the "
            "curve,\nfast, the same threshold works for files with "
            "different noise. Use the adaptive threshold on curves with "
            "sloped segments (polymer stretching).\n"
            "MSF: tests for a step (two parallel lines with an offset) in a "
            "window sliding along the curve, insensitive to sloped "
            "segments,\nthreshold in absolute units (scaled to a 10 pN "
            "step): adapt it when the noise changes between files.")
        self.label_events_curve_type = QtWidgets.QLabel("Display :")
        self.list_events_curve_type = PYAFComboBox(self, "events_curve_type")
        self.list_events_curve_type.addItem("Curve", "curve")
        self.list_events_curve_type.addItem("Detection", "detection")
        self.list_events_curve_type.setToolTip(
            "Curve: force curve with the detected events and their fits.\n"
            "Detection: signal of the detection (convolution or MSF) with "
            "the threshold (green line), to adjust the parameters:\nset "
            "the threshold just above the noise, so that only the jumps "
            "cross it, and check a few curves (Random curve).")
        for combo in (self.list_events_algorithm,
                      self.list_events_curve_type):
            combo.fit_to_contents()

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
        HL_events_compute_type.addWidget(self.label_events_algorithm)
        HL_events_compute_type.addWidget(self.list_events_algorithm)
        HL_events_compute_type.addSpacing(20)
        HL_events_compute_type.addWidget(self.label_events_curve_type)
        HL_events_compute_type.addWidget(self.list_events_curve_type)
        HL_events_compute_type.addStretch(1)

        # Jump of contact first (left), then the parameters of the
        # algorithm (column 1, see update_widget)
        self.GL_events = QtWidgets.QGridLayout()
        self.GL_events.addWidget(W_fit_params_jocs_events, 0, 0)

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

        # Description, small text below the coefficient (as in the
        # viscoelasticity tab)
        label_lr = QtWidgets.QLabel(
            "Effective loading rate of each event: slope of the force vs "
            "time (retraction velocity \u00d7 slope of the force vs "
            "distance),\nfitted on the part of the curve before the "
            "rupture, over this fraction of the rupture force (0.33: last "
            "third).\nEvans & Ritchie, Biophys. J. 1997; Ray et al., "
            "Langmuir 2007.")
        label_lr.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        label_lr.setFont(gui_tools.font_of_size(
            max(label_lr.font().pointSize() - 2, 8)))

        VL_coeff = QtWidgets.QVBoxLayout()
        VL_coeff.addStretch(1)
        VL_coeff.addWidget(self.IN_lr_coef)
        VL_coeff.addSpacing(6)
        VL_coeff.addWidget(label_lr)
        VL_coeff.addStretch(1)

        # Coefficient centered in the tab (the spacing balances the compute
        # button)
        HL_loading_rates.addSpacing(self.BT_get_lr.width())
        HL_loading_rates.addStretch(1)
        HL_loading_rates.addLayout(VL_coeff)
        HL_loading_rates.addStretch(1)
        HL_loading_rates.addLayout(VL_loading_rates)

        self.widget_loading_rates.setLayout(HL_loading_rates)

        # --- Stiffness height correction -------------------------------------

        self.BT_height_corr = PYAFButton(self, "height_correction", "Compute")
        self.BT_height_corr.setFixedSize(100, 100)

        # Main action buttons, highlighted by the stylesheet
        for button in (self.BT_compute_stiffness,
                       self.BT_compute_work_and_rupture_force,
                       self.BT_compute_events, self.BT_get_lr,
                       self.BT_height_corr):
            button.setObjectName("primaryButton")

        # Instruction label (models with a correction of Garcia & Garcia)
        names = [contact_models.MODEL_NAMES[m] for m in contact_models.HAS_BEC]
        self.label_info = QtWidgets.QLabel(
            "Bottom effect correction of a sample bonded to the substrate "
            "(Garcia & Garcia, Biophys. J. 2018).\nApplies to " +
            ", ".join(names[:-1]) + " and " + names[-1] +
            " models.\nEither select the substrate from the ROI manager "
            "(results tab) or input the known height of the sample")
        self.label_info.setTextInteractionFlags(
            QtCore.Qt.TextSelectableByMouse)
        # Small text, below the parameters, as in the viscoelasticity tab
        self.label_info.setFont(gui_tools.font_of_size(
            max(self.label_info.font().pointSize() - 2, 8)))

        # ROI selection (Label + Dropdown)
        self.label_list_roi = QtWidgets.QLabel("ROI of substrate:")
        self.list_roi_glass = PYAFComboBox(self, "roi")
        self.list_roi_glass.setFixedSize(100, 26)  # Fixed width to prevent stretching

        # Noise floor input (Label + Input Box)
        self.label_noise_floor = QtWidgets.QLabel(
            "Noise floor: discard the first")  # Label for noise floor
        self.IN_noise_floor = PYAFInput(self, "noise_floor", "[nm]", width=50)  # Input box for noise floor
        self.IN_noise_floor.input.setValidator(
            misc_tools.validator("UF_or_empty"))  # Validation

        # Sample height input (Switched Label & Input Box)
        self.label_sample_height = QtWidgets.QLabel("Sample height [nm]:")  # Separate label
        self.IN_user_h = PYAFInput(self, "user_h", "", width=50)  # Empty label in PYAFInput
        self.IN_user_h.input.setValidator(misc_tools.validator("UF_or_empty"))

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
        self.stiff_corr_opts.addLayout(self.stiff_corr_row)  # ROI label + dropdown
        self.stiff_corr_opts.addLayout(self.stiff_corr_noise_row)  # Add noise floor layout
        self.stiff_corr_opts.addLayout(self.stiff_corr_height_row)  # Sample height label + input
        # Instruction text, below the parameters (as in the viscoelasticity
        # tab)
        self.stiff_corr_opts.addSpacing(6)
        self.stiff_corr_opts.addWidget(self.label_info)

        VL_stiff_corr_opts = QtWidgets.QVBoxLayout()
        VL_stiff_corr_opts.addLayout(self.stiff_corr_opts)
        VL_stiff_corr_opts.addStretch(1)

        VL_stiff_corr = QtWidgets.QVBoxLayout()
        VL_stiff_corr.addStretch(1)
        VL_stiff_corr.addWidget(self.BT_height_corr)
        VL_stiff_corr.addStretch(1)

        # Options centered in the tab (the spacing balances the compute
        # button), as in the other tabs
        HL_stiff_corr = QtWidgets.QHBoxLayout()
        HL_stiff_corr.addSpacing(self.BT_height_corr.width())
        HL_stiff_corr.addStretch(1)
        HL_stiff_corr.addLayout(VL_stiff_corr_opts)
        HL_stiff_corr.addStretch(1)
        HL_stiff_corr.addLayout(VL_stiff_corr)

        self.widget_stiff_corr.setLayout(HL_stiff_corr)

        # --- Viscoelasticity -------------------------------------------------

        self.widget_visco = QtWidgets.QWidget()

        label = QtWidgets.QLabel(
            "Ting fit of the approach and the retraction (power-law "
            "rheology, Efremov et al. 2017; B.E.C.: Tejedor & Garcia "
            "2025).\nUses the point of contact, the baseline, the model and "
            "the modulus of the elasticity fit (Fitting).\n"
            "Model and fit procedure of PyFMRheo: "
            "github.com/DyNaMo-INSERM/PyFMlab_DyNaMo")
        label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        # Small text, below the parameters
        label.setFont(gui_tools.font_of_size(
            max(label.font().pointSize() - 2, 8)))

        self.Box_visco = QtWidgets.QGroupBox("Ting fit")
        GL_visco = QtWidgets.QGridLayout()
        GL_visco.setVerticalSpacing(4)

        self.label_visco_solver = QtWidgets.QLabel("Solver :")
        self.list_visco_solver = PYAFComboBox(self, "visco_solver")
        self.list_visco_solver.addItem("Analytical (constant speeds)")
        self.list_visco_solver.addItem("Numerical (measured indentation)")
        HL_solver = QtWidgets.QHBoxLayout()
        HL_solver.addWidget(self.label_visco_solver)
        HL_solver.addWidget(self.list_visco_solver)
        HL_solver.addStretch(1)

        self.CB_visco_bec = PYAFCheckBox(
            self, "visco_bec", "Bottom effect correction (B.E.C. tab)")
        self.CB_visco_bec.setToolTip(
            "Height of the sample defined in the B.E.C. tab (substrate ROI "
            "or sample height). Uses the numerical solver.")
        self.CB_visco_drag = PYAFCheckBox(
            self, "visco_drag_correction", "Correct viscous drag")
        self.CB_visco_drag.setToolTip(
            "The retract baseline is shifted by the viscous drag of the "
            "cantilever: its offset from the approach baseline is removed.")
        self.IN_visco_offset = PYAFInput(
            self, "visco_contact_offset", "Baseline before contact [nm]",
            60)
        self.IN_visco_offset.setToolTip(
            "Part of the curve before the point of contact included in the "
            "fit (as the contact_offset of PyFMRheo): the time of contact is "
            "fitted, it can move from the point of contact of the elasticity "
            "fit.")
        self.IN_visco_offset.input.setValidator(misc_tools.validator("UF"))
        self.IN_visco_points = PYAFInput(
            self, "visco_nbr_points", "Points (downsampling)", 60)
        self.IN_visco_points.setToolTip(
            "The curve is averaged down to this number of points before the "
            "fit (as pts_downsample in PyFMRheo): the numerical solver "
            "grows with the square of the number of points.")
        self.IN_visco_points.input.setValidator(
            misc_tools.validator("MI", bottom=20))
        self.CB_visco_beta_auto = PYAFCheckBox(
            self, "visco_beta_auto",
            "Automatic initial fluidity (0.05 if E > 10 kPa, else 0.25)")
        self.IN_visco_beta = PYAFInput(
            self, "visco_beta_init", "Initial fluidity", 60)
        self.IN_visco_beta.input.setValidator(
            misc_tools.validator("RF", bottom=0.01, top=0.99, decimal=3))
        self.IN_visco_frequency = PYAFInput(
            self, "visco_frequency", "Frequency of E\u2032, E\u2033 [Hz]",
            60)
        self.IN_visco_frequency.setToolTip(
            "Frequency of the storage and loss moduli maps (Results tab), "
            "computed from E\u2080 and the fluidity\n(power law, no new fit). "
            "Limited to the time scales of the curves, measured by the "
            "computation (medians of the curves):\nfrom 1 / contact time "
            "to 1 / (2 x time between two fitted points); outside, the "
            "values would be an\nextrapolation of the model. Set to 1 / "
            "contact time (time scale of the measurement) by the "
            "computation, unless a frequency\nsupported by the curves was "
            "entered. Each file has its own frequency; with apply to "
            "all, a new value is also given to the files whose curves "
            "support it.")
        self.IN_visco_frequency.input.setValidator(misc_tools.validator("PF"))

        GL_visco.addLayout(HL_solver, 0, 0)
        GL_visco.addWidget(self.CB_visco_bec, 1, 0)
        GL_visco.addWidget(self.CB_visco_drag, 2, 0)
        GL_visco.addWidget(self.IN_visco_offset, 0, 1)
        GL_visco.addWidget(self.IN_visco_points, 1, 1)
        GL_visco.addWidget(self.CB_visco_beta_auto, 2, 1)
        GL_visco.addWidget(self.IN_visco_beta, 3, 1)
        GL_visco.addWidget(self.IN_visco_frequency, 4, 1)
        GL_visco.setColumnStretch(2, 1)
        self.Box_visco.setLayout(GL_visco)

        VL_visco_opts = QtWidgets.QVBoxLayout()
        VL_visco_opts.addWidget(self.Box_visco)
        VL_visco_opts.addWidget(label)
        VL_visco_opts.addStretch(1)

        self.BT_visco = PYAFButton(self, "compute_visco", "Compute")
        self.BT_visco.setFixedSize(100, 100)
        self.BT_visco.setObjectName("primaryButton")
        VL_visco = QtWidgets.QVBoxLayout()
        VL_visco.addStretch(1)
        VL_visco.addWidget(self.BT_visco)
        VL_visco.addStretch(1)

        # Options centered in the tab (the spacing balances the compute
        # button), as in the other tabs
        HL_visco = QtWidgets.QHBoxLayout()
        HL_visco.addSpacing(self.BT_visco.width())
        HL_visco.addStretch(1)
        HL_visco.addLayout(VL_visco_opts)
        HL_visco.addStretch(1)
        HL_visco.addLayout(VL_visco)
        self.widget_visco.setLayout(HL_visco)

        # --- Microrheology ---------------------------------------------------

        self.widget_microrheo = QtWidgets.QWidget()

        label = QtWidgets.QLabel(
            "Complex shear modulus G* = G\u2032 + i G\u2033 from the sine "
            "modulations of the piezo in contact,"
            "\none value per frequency: F/\u03b4 at the frequency, "
            "linearized contact model of the elasticity fit at the working "
            "indentation (Alcaraz et al. 2003).\n"
            "Uses the point of contact, the baseline, the model and the tip "
            "of the elasticity fit (Fitting). Method of PyFMRheo: "
            "github.com/DyNaMo-INSERM/PyFMlab_DyNaMo")
        label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        # Small text, below the parameters
        label.setFont(gui_tools.font_of_size(
            max(label.font().pointSize() - 2, 8)))

        self.Box_microrheo = QtWidgets.QGroupBox("Microrheology")
        GL_microrheo = QtWidgets.QGridLayout()
        GL_microrheo.setVerticalSpacing(4)

        self.label_microrheo_method = QtWidgets.QLabel("Method :")
        self.list_microrheo_method = PYAFComboBox(self, "microrheo_method")
        for name, key in zip(microrheology.METHOD_NAMES,
                             microrheology.METHODS):
            self.list_microrheo_method.addItem(name, key)
        self.list_microrheo_method.setToolTip(
            "FFT: ratio of the Fourier transforms of the force and of the "
            "indentation at the frequency (transfer function).\n"
            "Sine fit: ratio of the sines fitted at the frequency on the "
            "force and on the indentation.\nThe slow drift of each "
            "modulation is removed first (rolling average over one "
            "period). Both give the same values for whole periods.")
        self.list_microrheo_method.fit_to_contents()
        HL_microrheo_method = QtWidgets.QHBoxLayout()
        HL_microrheo_method.addWidget(self.label_microrheo_method)
        HL_microrheo_method.addWidget(self.list_microrheo_method)
        HL_microrheo_method.addStretch(1)

        self.IN_microrheo_drag = PYAFInput(
            self, "microrheo_drag",
            "Viscous drag b(0) [\u00b5N\u00b7s/m]", 60)
        self.IN_microrheo_drag.setToolTip(
            "Viscous drag coefficient of the liquid on the cantilever at "
            "the surface, b(0): i 2\u03c0 f b(0) is removed from F/\u03b4 "
            "(Alcaraz et al. 2003).\nMeasured by oscillating the "
            "cantilever at several heights above the surface and "
            "extrapolating to contact (0: no correction).")
        self.IN_microrheo_drag.input.setValidator(
            misc_tools.validator("UF"))
        self.BT_microrheo_drag = PYAFButton(
            self, "microrheo_drag_file", "From file\u2026")
        self.BT_microrheo_drag.setToolTip(
            "b(0) from a viscous drag measurement loaded in pyAF: "
            "modulations of the cantilever at several heights above the "
            "surface (no contact).\nb(h) = Im(F/\u03b4) / (2\u03c0 f) at "
            "each height, extrapolated to the surface with "
            "b(h) = A / (h + h_eff) (Alcaraz et al. 2002).")
        HL_microrheo_drag = QtWidgets.QHBoxLayout()
        HL_microrheo_drag.addWidget(self.IN_microrheo_drag)
        HL_microrheo_drag.addWidget(self.BT_microrheo_drag)
        HL_microrheo_drag.addStretch(1)

        # Maximal frequency (head of the AFM: resonances of the head above
        # about 500 Hz, the deflection no longer follows the piezo)
        self.IN_microrheo_max_freq = PYAFInput(
            self, "microrheo_max_freq", "Max. frequency [Hz]", 60)
        self.IN_microrheo_max_freq.setToolTip(
            "The modulations above this frequency are not analysed (0: all "
            "the frequencies), as max_freq in PyFMRheo.\nAbove the "
            "bandwidth of the AFM head (resonances), the deflection does not "
            "follow the piezo, even on the substrate:\ncheck the piezo "
            "characterization (amplitude quotient far below 1, phase above "
            "about 45\u00b0).")
        self.IN_microrheo_max_freq.input.setValidator(
            misc_tools.validator("UF"))
        HL_microrheo_method.insertSpacing(HL_microrheo_method.count() - 1, 20)
        HL_microrheo_method.insertWidget(HL_microrheo_method.count() - 1,
                                         self.IN_microrheo_max_freq)
        # Preview: the oscillations (raw data) or G' and G''
        self.TG_microrheo_raw = PYAFToggle(
            self, "microrheo_raw", "Show the oscillations")
        self.TG_microrheo_raw.setToolTip(
            "Preview of the displayed curve: the oscillations (deflection vs "
            "time of each modulation) instead of G\u2032 and G\u2033 vs "
            "the frequency.")
        HL_microrheo_method.insertSpacing(HL_microrheo_method.count() - 1, 20)
        HL_microrheo_method.insertWidget(HL_microrheo_method.count() - 1,
                                         self.TG_microrheo_raw)
        GL_microrheo.addLayout(HL_microrheo_method, 0, 0)
        # Row 2: after the piezo correction (used by the drag measurement)
        GL_microrheo.addLayout(HL_microrheo_drag, 2, 0)

        # Piezo correction: lag and gain of the piezo measured on a curve of
        # a rigid sample (substrate)
        self.label_microrheo_piezo = QtWidgets.QLabel()
        self.BT_microrheo_piezo = PYAFButton(
            self, "microrheo_piezo", "From this curve")
        self.BT_microrheo_piezo.setToolTip(
            "Lag and gain of the piezo measured on the displayed curve, which "
            "must be on a rigid sample (substrate: the deflection follows the "
            "piezo),\nwith the same modulations: the height is corrected at "
            "each frequency (PyFMRheo, piezo characterization).")
        self.BT_microrheo_piezo_clear = PYAFButton(
            self, "microrheo_piezo_clear", "Clear")
        HL_microrheo_piezo = QtWidgets.QHBoxLayout()
        HL_microrheo_piezo.addWidget(QtWidgets.QLabel("Piezo correction : "))
        HL_microrheo_piezo.addWidget(self.label_microrheo_piezo)
        HL_microrheo_piezo.addWidget(self.BT_microrheo_piezo)
        HL_microrheo_piezo.addWidget(self.BT_microrheo_piezo_clear)
        self.CB_microrheo_piezo_amplitude = PYAFCheckBox(
            self, "microrheo_piezo_amplitude", "Correct the amplitude")
        self.CB_microrheo_piezo_amplitude.setToolTip(
            "Also correct the amplitude of the height by the amplitude "
            "quotient of the reference (as corr_amp in PyFMRheo).\nOff: "
            "the phase lag only, as in Alcaraz et al. 2002 and the "
            "JPK/Bruker note TN301.")
        HL_microrheo_piezo.addWidget(self.CB_microrheo_piezo_amplitude)
        HL_microrheo_piezo.addStretch(1)
        GL_microrheo.addLayout(HL_microrheo_piezo, 1, 0)
        self.Box_microrheo.setLayout(GL_microrheo)

        VL_microrheo_opts = QtWidgets.QVBoxLayout()
        VL_microrheo_opts.addWidget(self.Box_microrheo)
        VL_microrheo_opts.addWidget(label)
        VL_microrheo_opts.addStretch(1)

        self.BT_microrheo = PYAFButton(self, "compute_microrheo", "Compute")
        self.BT_microrheo.setFixedSize(100, 100)
        self.BT_microrheo.setObjectName("primaryButton")
        VL_microrheo = QtWidgets.QVBoxLayout()
        VL_microrheo.addStretch(1)
        VL_microrheo.addWidget(self.BT_microrheo)
        VL_microrheo.addStretch(1)

        # Options centered in the tab (the spacing balances the compute
        # button), as in the other tabs
        HL_microrheo = QtWidgets.QHBoxLayout()
        HL_microrheo.addSpacing(self.BT_microrheo.width())
        HL_microrheo.addStretch(1)
        HL_microrheo.addLayout(VL_microrheo_opts)
        HL_microrheo.addStretch(1)
        HL_microrheo.addLayout(VL_microrheo)
        self.widget_microrheo.setLayout(HL_microrheo)

        # --- Breakthroughs of supported lipid bilayers -----------------------

        self.widget_breakthrough = QtWidgets.QWidget()
        self.Box_breakthrough = QtWidgets.QGroupBox(
            "Breakthroughs (approach curves)")
        self.IN_bt_n_sigma = PYAFInput(
            self, "breakthrough_n_sigma", "Threshold (\u00d7 noise)", 60)
        self.IN_bt_n_sigma.setToolTip(
            "A breakthrough is a sudden jump of the separation (tip-sample "
            "distance) on the approach, with a drop of the force: the jump "
            "must exceed this number of times the noise of the separation "
            "of the curve")
        self.IN_bt_n_sigma.input.setValidator(misc_tools.validator("PF"))
        self.IN_bt_min_thickness = PYAFInput(
            self, "breakthrough_min_thickness", "Minimum jump [nm]", 60)
        self.IN_bt_min_thickness.setToolTip(
            "Minimum jump of the separation: the thickness of the bilayer "
            "compressed at the rupture (about 4-5 nm for a bilayer)")
        self.IN_bt_min_thickness.input.setValidator(
            misc_tools.validator("PF"))
        self.IN_bt_min_force = PYAFInput(
            self, "breakthrough_min_force", "Minimum force [nN]", 60)
        self.IN_bt_min_force.setToolTip(
            "Minimum breakthrough force (Fb): higher to skip the jumps at "
            "the contact")
        self.IN_bt_min_force.input.setValidator(misc_tools.validator("UF"))
        self.IN_bt_max_events = PYAFInput(
            self, "breakthrough_max_events", "Maximum per curve", 60)
        self.IN_bt_max_events.setToolTip(
            "Maximum number of breakthroughs per curve (stacked bilayers); "
            "the maps and plots use the first one")
        self.IN_bt_max_events.input.setValidator(
            misc_tools.validator("RI", bottom=1, top=20))
        GL_bt = QtWidgets.QGridLayout()
        GL_bt.addWidget(self.IN_bt_n_sigma, 0, 0)
        GL_bt.addWidget(self.IN_bt_min_thickness, 1, 0)
        GL_bt.addWidget(self.IN_bt_min_force, 0, 1)
        GL_bt.addWidget(self.IN_bt_max_events, 1, 1)
        GL_bt.setColumnStretch(2, 1)
        self.Box_breakthrough.setLayout(GL_bt)
        label_bt = QtWidgets.QLabel(
            "Breakthrough force Fb (force at the rupture of the bilayer) and "
            "jump of the separation (thickness of the compressed bilayer).\n"
            "Preview: approach force vs distance to the support, with the "
            "breakthroughs found with these settings.\n"
            "Nominal loading rate: spring constant \u00d7 approach velocity "
            "(Gum\u00ed-Audenis & Giannotti 2019).\n"
            "Correct the baseline tilt first (Tilting/Smoothing): a tilted "
            "baseline biases Fb and can hide the low breakthroughs.")
        label_bt.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        # Small text, below the parameters, as in the viscoelasticity tab
        label_bt.setFont(gui_tools.font_of_size(
            max(label_bt.font().pointSize() - 2, 8)))
        VL_bt_opts = QtWidgets.QVBoxLayout()
        VL_bt_opts.addWidget(self.Box_breakthrough)
        VL_bt_opts.addWidget(label_bt)
        VL_bt_opts.addStretch(1)

        self.BT_breakthrough = PYAFButton(
            self, "compute_breakthrough", "Compute")
        self.BT_breakthrough.setFixedSize(100, 100)
        self.BT_breakthrough.setObjectName("primaryButton")
        VL_bt = QtWidgets.QVBoxLayout()
        VL_bt.addStretch(1)
        VL_bt.addWidget(self.BT_breakthrough)
        VL_bt.addStretch(1)

        # Options centered in the tab, as in the other tabs
        HL_bt = QtWidgets.QHBoxLayout()
        HL_bt.addSpacing(self.BT_breakthrough.width())
        HL_bt.addStretch(1)
        HL_bt.addLayout(VL_bt_opts)
        HL_bt.addStretch(1)
        HL_bt.addLayout(VL_bt)
        self.widget_breakthrough.setLayout(HL_bt)

        # --- Set up page -----------------------------------------------------

        self.tabs.addTab(self.widget_app, "Elasticity")
        self.tabs.addTab(self.widget_stiff_corr, "B.E.C.")
        self.tabs.addTab(self.widget_visco, "Viscoelasticity")
        self.tabs.addTab(self.widget_microrheo, "Microrheology")
        self.tabs.addTab(self.widget_ret, "Work, Rupture force")
        self.tabs.addTab(self.widget_events, "Unbinding events")
        self.tabs.addTab(self.widget_loading_rates, "Loading rates")
        self.tabs.addTab(self.widget_breakthrough, "Breakthrough")

        # Margins of the pages, with the padding of the pane (no padding in
        # the stylesheet, see QTabWidget#subTabs::pane)
        for i in range(self.tabs.count()):
            layout = QtWidgets.QWidget.layout(self.tabs.widget(i))
            m = layout.contentsMargins()
            layout.setContentsMargins(m.left() + 4, m.top() + 4,
                                      m.right() + 4, m.bottom() + 4)

        # The boxes keep their height (no empty space at their bottom) when
        # the tabs are taller
        self.HL_widget_app.setAlignment(QtCore.Qt.AlignTop)
        self.HL_widget_ret.setAlignment(QtCore.Qt.AlignTop)

        # The events tab is the tallest: no vertical margins around the
        # parameter boxes (the tab and the boxes keep theirs), to leave more
        # height to the plots
        for widget in (self.widget_kernel, self.widget_msf):
            # (some widgets have an attribute named layout)
            layout = QtWidgets.QWidget.layout(widget)
            margins = layout.contentsMargins()
            layout.setContentsMargins(margins.left(), 0, margins.right(), 0)

        self.GL_compute = QtWidgets.QGridLayout()
        # Compact rows (plots, navigation, options, computation tabs)
        self.GL_compute.setVerticalSpacing(4)
        self.GL_compute.addLayout(self.HL_top, 0, 0, 1, 0)
        self.GL_compute.addWidget(self.MPL_meshgrid, 1, 0)
        self.GL_compute.addWidget(self.MPL_curve, 1, 1)
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
            self.GL_choosers.addWidget(self.curve_selector, 0, 2)
            self.curve_selector.update_widget()

        if data.events_algorithm == "kernel":
            self.widget_msf.setParent(None)
            self.GL_events.addWidget(self.widget_kernel, 0, 1)

        elif data.events_algorithm == "msf":
            self.widget_kernel.setParent(None)
            self.GL_events.addWidget(self.widget_msf, 0, 1)

        self.STK_stiffness_mode.setCurrentIndex(data.stiffness_mode)

        # Tomography and Fitting parameters (start, segment size,
        # ranges) of the selected file
        self.widget_calculus.update_widget()
        self.widget_fitting.update_widget()
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

        widget = self.tabs.currentWidget()
        if widget is self.widget_app:
            shared.exp.compute_type = "stiffness"
            self.BT_next_zoom_element.setEnabled(False)
        elif widget is self.widget_stiff_corr:
            shared.exp.compute_type = "stiffness_corr"
            self.BT_next_zoom_element.setEnabled(False)
        elif widget is self.widget_visco:
            shared.exp.compute_type = "viscoelasticity"
            self.BT_next_zoom_element.setEnabled(False)
        elif widget is self.widget_microrheo:
            shared.exp.compute_type = "microrheology"
            self.BT_next_zoom_element.setEnabled(False)
        elif widget is self.widget_breakthrough:
            shared.exp.compute_type = "breakthrough"
            self.BT_next_zoom_element.setEnabled(True)
        elif widget is self.widget_ret:
            shared.exp.compute_type = "work_and_rupture_force"
            self.BT_next_zoom_element.setEnabled(True)
        elif widget is self.widget_events:
            shared.exp.compute_type = "events"
            self.BT_next_zoom_element.setEnabled(True)
        elif widget is self.widget_loading_rates:
            shared.exp.compute_type = "loading_rates"
            self.BT_next_zoom_element.setEnabled(True)

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
                self.GL_compute.addWidget(self.MPL_meshgrid, 1, 0)
                self.BT_clean_up.setEnabled(True)

        if element == "apply_to_all" or element == "all":
            self.TG_apply_to_all.setChecked(shared.exp.apply_to_all_compute)

        if element == "stiffness_model_selected" or element == "all":
            # Tomography: Slope model only. Fitting: the other models (files
            # saved with the previous Calculus mode)
            tomography = not data.perform_fit
            if tomography and data.stiffness_model_selected != 3:
                self.set_model(data, 3)
                data.tomography = 1
            elif not tomography and data.stiffness_model_selected == 3:
                self.set_model(data, 0)
            for index in range(len(self.models)):
                if (index == 3) == tomography:
                    gui_tools.enable_combobox_item(self.list_models, index)
                else:
                    gui_tools.disable_combobox_item(self.list_models, index)
            self.list_models.setCurrentIndex(data.stiffness_model_selected)

            self.update_MPL("MPL_canvas")
            self.update_GUI("buttons_compute")

        if element == "BT_height_correction" or element == "all":
            # Bottom effect correction: elasticity, model with a correction,
            # and a height (substrate ROI or sample height)
            reason = bec.check_file(data)
            self.BT_height_corr.setEnabled(reason is None)
            self.BT_height_corr.setToolTip(
                "" if reason is None else "Not available: " + reason)

        if element == "breakthrough" or element == "all":
            self.IN_bt_n_sigma.changeValue(
                misc_tools.format_number(data.breakthrough_n_sigma))
            self.IN_bt_min_thickness.changeValue(
                misc_tools.format_number(data.breakthrough_min_thickness))
            self.IN_bt_min_force.changeValue(
                misc_tools.format_number(data.breakthrough_min_force))
            self.IN_bt_max_events.changeValue(
                int(data.breakthrough_max_events))

        if element == "visco" or element == "all":
            self.list_visco_solver.setCurrentIndex(
                1 if data.visco_bec else data.visco_solver)
            self.list_visco_solver.setEnabled(not data.visco_bec)
            self.CB_visco_bec.setChecked(data.visco_bec)
            self.CB_visco_drag.setChecked(data.visco_drag_correction)
            self.IN_visco_offset.changeValue(data.visco_contact_offset)
            self.IN_visco_points.changeValue(data.visco_nbr_points)
            self.CB_visco_beta_auto.setChecked(data.visco_beta_auto)
            self.IN_visco_beta.changeValue(data.visco_beta_init)
            # Frequency of E' and E'': within the time scales of the curves
            # (window measured by the computation)
            window = getattr(data, "visco_frequency_window", None)
            label = "Frequency of E\u2032, E\u2033 [Hz]"
            if window is not None:
                label += " (" + "\u2013".join(
                    "%.3g" % f if f < 1000 else "%.0f" % f
                    for f in window) + ")"
                self.IN_visco_frequency.input.setValidator(
                    misc_tools.validator("RF", bottom=window[0],
                                         top=window[1], decimal=3))
            else:
                self.IN_visco_frequency.input.setValidator(
                    misc_tools.validator("PF"))
            self.IN_visco_frequency.label.setText("  " + label)
            self.IN_visco_frequency.changeValue(
                float("%.3g" % getattr(data, "visco_frequency", 1.0)))
            self.IN_visco_beta.setEnabled(not data.visco_beta_auto)
            # The elasticity fit is needed (and a height for the B.E.C.)
            reason = visco_main.check_file(data)
            self.BT_visco.setEnabled(reason is None)
            self.BT_visco.setToolTip(
                "" if reason is None else "Not available: " + reason)

        if element == "microrheo" or element == "all":
            self.list_microrheo_method.setCurrentIndex(
                self.list_microrheo_method.findData(data.microrheo_method))
            self.IN_microrheo_drag.changeValue(
                misc_tools.format_number(data.microrheo_drag))
            self.IN_microrheo_max_freq.changeValue(
                misc_tools.format_number(data.microrheo_max_freq or 0))
            piezo = data.microrheo_piezo
            self.label_microrheo_piezo.setText(
                "None" if not piezo else piezo["file"] + " [" + ", ".join(
                    str(k) for k in piezo["curve"]) + "]")
            self.BT_microrheo_piezo.setEnabled(
                microrheology.available(data))
            self.BT_microrheo_piezo_clear.setEnabled(bool(piezo))
            self.CB_microrheo_piezo_amplitude.setChecked(
                bool(data.microrheo_piezo_amplitude))
            self.CB_microrheo_piezo_amplitude.setEnabled(bool(piezo))
            self.TG_microrheo_raw.setChecked(
                bool(getattr(data, "microrheo_preview_raw", False)))
            self.TG_microrheo_raw.setEnabled(microrheology.available(data))
            # Modulation segments and the elasticity fit are needed
            reason = microrheology.check_file(data)
            self.BT_microrheo.setEnabled(reason is None)
            self.BT_microrheo.setToolTip(
                "" if reason is None else "Not available: " + reason)

        if element == "roi_list_height_corr" or element == "all":
            # ROIs set as substrate in the ROI manager, numbered as in it
            # (position in the list; roi_glass_id: position + 1, 0: none)
            substrates = [index + 1 for index, roi in
                          enumerate(data.roi_list)
                          if roi.glass_coeffs is not None]
            if data.roi_glass_id not in substrates:
                # The substrate ROI was removed or unset
                data.roi_glass_id = 0
            if data.roi_glass_id == 0 and len(substrates) == 1 and \
                    not data.user_h:
                # A single substrate: chosen
                data.roi_glass_id = substrates[0]
            self.list_roi_glass.blockSignals(True)
            self.list_roi_glass.clear()
            self.list_roi_glass.addItem("None")
            for number in substrates:
                self.list_roi_glass.addItem(str(number))
            self.list_roi_glass.setCurrentText(
                str(data.roi_glass_id) if data.roi_glass_id else "None")
            self.list_roi_glass.blockSignals(False)
            # Inputs and button with the substrate (also after "all": the
            # button is updated above, before the substrate)
            self.update_GUI("IN_user_h")
            self.update_GUI("BT_height_correction")

        if element == "IN_user_h" or element == "all":
            # Empty input when no value is defined
            self.IN_user_h.changeValue(
                "" if data.user_h is None else data.user_h)
            self.IN_noise_floor.changeValue(
                "" if data.noise_floor is None else data.noise_floor)
            # The noise floor is used with a substrate ROI, the sample height
            # without one
            self.IN_noise_floor.setEnabled(data.roi_glass_id != 0)
            self.IN_user_h.setEnabled(data.roi_glass_id == 0)

        if element == "get_lr" or element == "all":
            self.BT_get_lr.setEnabled(data.events_calculated)

        if element == "grid_model_options" or element == "all":
            # Inputs of the tip of the model: radius and/or angle (none for
            # the slope)
            model = data.stiffness_model_selected
            self.IN_tip_radius.setParent(None)
            self.IN_tip_angle.setParent(None)
            self.IN_shell_thickness.setParent(None)
            self.IN_poisson.setParent(None)
            if model != contact_models.SLOPE:
                self.GL_stiffness_params.addWidget(self.IN_poisson, 1, 0)
            row = 2
            if model in contact_models.USES_RADIUS:
                self.IN_tip_radius.label.setText(
                    "  " + contact_models.RADIUS_LABELS.get(
                        model, contact_models.RADIUS_LABEL))
                self.IN_tip_radius.input.setText(str(data.tip_radius))
                self.GL_stiffness_params.addWidget(self.IN_tip_radius, row, 0)
                row += 1
            if model in contact_models.USES_ANGLE:
                self.IN_tip_angle.input.setText(str(data.tip_angle))
                self.GL_stiffness_params.addWidget(self.IN_tip_angle, row, 0)
            if model in contact_models.USES_THICKNESS:
                self.IN_shell_thickness.input.setText(
                    misc_tools.format_number(data.shell_thickness))
                self.GL_stiffness_params.addWidget(
                    self.IN_shell_thickness, row, 0)

            # Brush model: models of a sphere
            self.TG_brush.setParent(None)
            self.IN_brush_threshold_pct.setParent(None)
            if brush.available(data):
                self.TG_brush.setChecked(data.brush_model)
                self.IN_brush_threshold_pct.input.setText(
                    misc_tools.format_number(data.brush_threshold_pct))
                self.IN_brush_threshold_pct.setEnabled(data.brush_model)
                self.GL_stiffness_params.addWidget(self.TG_brush, 4, 0)
                self.GL_stiffness_params.addWidget(
                    self.IN_brush_threshold_pct, 5, 0)

            self.IN_poisson.input.setText(str(data.poisson_ratio))

            # Display the mode of the file (the slope model may have changed it)
            self.update_GUI("stiffness_params")

        # if element == "input_segments" or element == "all":
            # self.IN_start.input.setText(str(data.indentation_start))
            # self.IN_stop.input.setText(str(data.indentation_stop))
            # self.IN_step.input.setText(str(data.indentation_step))

        if element == "zoom_in" or element == "all":
            self.TG_zoom_in.setChecked(data.auto_zoom_in)

        if element == "zoom_fit_preview_factor" or element == "all":
            self.IN_zoom_in_factor.changeValue(data.zoom_fit_preview_factor)

        if element == "compute_all" or element == "all":
            self.TG_compute_all.setChecked(shared.exp.calc_all)

        if element == "flatten" or element == "all":
            if data.scan_size_x != 0:
                self.BT_flatten.setEnabled(True)
            else:
                self.BT_flatten.setEnabled(False)

        # if element == "tomography" or element == "all":
            # self.CB_tomography.setChecked(data.tomography)

        if element == "input_lr_coef" or element == "all":
            self.IN_lr_coef.changeValue(data.lr_coef)

        if element == "events_algorithms" or element == "all":
            self.list_events_algorithm.blockSignals(True)
            self.list_events_algorithm.setCurrentIndex(
                self.list_events_algorithm.findData(data.events_algorithm))
            self.list_events_algorithm.blockSignals(False)

        if element == "stiffness_params" or element == "all":
            if data.stiffness_mode == 0:
                self.RBT_calculus_params.setChecked(True)
            elif data.stiffness_mode == 1:
                self.RBT_fitting_params.setChecked(True)

        if element == "events_curve_type" or element == "all":
            self.list_events_curve_type.blockSignals(True)
            self.list_events_curve_type.setCurrentIndex(
                self.list_events_curve_type.findData(data.events_curve_type))
            self.list_events_curve_type.blockSignals(False)

    @staticmethod
    def set_model(data, model):
        """Selects the model of the elasticity. Each model keeps its own tip
        radius and angle (default values of the model until the user changes
        them)."""
        old = data.stiffness_model_selected
        if model != old:
            if old in contact_models.USES_RADIUS:
                data.tip_radii[old] = data.tip_radius
            if old in contact_models.USES_ANGLE:
                data.tip_angles[old] = data.tip_angle
            if model in contact_models.USES_RADIUS:
                data.tip_radius = data.tip_radii.get(
                    model, contact_models.DEFAULT_RADIUS[model])
            if model in contact_models.USES_ANGLE:
                data.tip_angle = data.tip_angles.get(
                    model, contact_models.DEFAULT_ANGLE[model])
        data.stiffness_model_selected = model

    def list_updated(self, name):
        """Called when a list is updated."""
        data = shared.exp.current_data

        if name == "list_models":
            # This list lets you chose the model for the stiffness computation.
            self.set_model(data, self.list_models.currentIndex())
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

        elif name == "microrheo_method":
            data.microrheo_method = self.list_microrheo_method.currentData()
            self.update_MPL("MPL_canvas")

        elif name == "events_algorithm":
            data.events_algorithm = self.list_events_algorithm.currentData()
            self.update_widget()

        elif name == "events_curve_type":
            data.events_curve_type = self.list_events_curve_type.currentData()
            self.update_widget()

        elif name == "visco_solver":
            data.visco_solver = self.list_visco_solver.currentIndex()
            self.update_MPL("MPL_canvas")

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

        elif field == "shell_thickness":
            if self.IN_shell_thickness.input.text() != "":
                data.shell_thickness = \
                    self.IN_shell_thickness.get_float_value()

        elif field == "poisson":
            if self.IN_poisson.input.text() != "":
                data.poisson_ratio = self.IN_poisson.get_float_value()

        # elif field == "input_start":
            # if self.IN_start.input.text() != "":
                # data.indentation_start = self.IN_start.get_int_value()

        # elif field == "input_stop":
            # if self.IN_stop.input.text() != "":
                # data.indentation_stop = self.IN_stop.get_int_value()

        # elif field == "input_step":
            # if self.IN_step.input.text() != "":
                # data.indentation_step = self.IN_step.get_int_value()

        elif field == "input_lr_coef":
            data.lr_coef = self.IN_lr_coef.get_float_value()

        elif field == "zoom_in_factor":
            val = self.IN_zoom_in_factor.get_float_value()
            data.zoom_fit_preview_factor = val

        elif field == "noise_floor":
            # An empty input removes the value
            if self.IN_noise_floor.input.text() != "":
                data.noise_floor = self.IN_noise_floor.get_float_value()
            else:
                data.noise_floor = None

        elif field == "user_h":
            if self.IN_user_h.input.text() != "":
                data.user_h = self.IN_user_h.get_float_value()
            else:
                data.user_h = None

        elif field == "visco_contact_offset":
            data.visco_contact_offset = self.IN_visco_offset.get_float_value()

        elif field == "brush_threshold_pct":
            data.brush_threshold_pct = \
                self.IN_brush_threshold_pct.get_float_value()

        elif field == "breakthrough_n_sigma":
            data.breakthrough_n_sigma = self.IN_bt_n_sigma.get_float_value()

        elif field == "breakthrough_min_thickness":
            data.breakthrough_min_thickness = \
                self.IN_bt_min_thickness.get_float_value()

        elif field == "breakthrough_min_force":
            data.breakthrough_min_force = \
                self.IN_bt_min_force.get_float_value()

        elif field == "breakthrough_max_events":
            data.breakthrough_max_events = \
                self.IN_bt_max_events.get_int_value()

        elif field == "visco_nbr_points":
            data.visco_nbr_points = self.IN_visco_points.get_int_value()

        elif field == "visco_beta_init":
            data.visco_beta_init = self.IN_visco_beta.get_float_value()

        elif field == "microrheo_max_freq":
            if self.IN_microrheo_max_freq.input.text() != "":
                data.microrheo_max_freq = \
                    self.IN_microrheo_max_freq.get_float_value()

        elif field == "microrheo_drag":
            if self.IN_microrheo_drag.input.text() != "":
                data.microrheo_drag = self.IN_microrheo_drag.get_float_value()

        elif field == "visco_frequency":
            if self.IN_visco_frequency.input.text() != "":
                value = self.IN_visco_frequency.get_float_value()
                data.visco_frequency = visco_main.clamp_frequency(data, value)
                # With apply to all: also the files whose curves support this
                # frequency (the others keep theirs)
                files = [data]
                if shared.exp.apply_to_all_compute:
                    files += [dt for dt in shared.exp.list if dt is not data
                              and visco_main.clamp_frequency(dt, value) ==
                              value]
                for dt in files:
                    dt.visco_frequency = data.visco_frequency
                    # Kept by the next computations (if supported)
                    dt.visco_frequency_user = True
                    # New values: automatic scales of the moduli maps
                    if dt.visco_calculated:
                        for mt in ("visco_storage", "visco_loss"):
                            plot_options.autoscale_elasticity(dt, mt)
                # Maps of the moduli at this frequency (no new fit)
                if data.meshgrid_type in ("visco_storage", "visco_loss"):
                    widgets_list.widget_results.update_MPL("MPL_canvas1")

        self.update_MPL("MPL_canvas")
        self.update_GUI("BT_height_correction")
        self.update_GUI("visco")
        self.update_GUI("microrheo")

    def piezo_from_curve(self):
        """Piezo characterization on the displayed curve (on a rigid
        sample), shown before it is used."""
        data = shared.exp.current_data
        i = shared.exp.meshgrid_click_xpos
        j = shared.exp.meshgrid_click_ypos
        result = microrheology.piezo_characterization(data, i, j)
        if result is None:
            return
        dialog = PiezoCharacterizationDialog(self, result)
        if dialog.exec_():
            files = shared.exp.list if shared.exp.apply_to_all_compute \
                else [data]
            for dt in files:
                dt.microrheo_piezo = dict(result)
            self.update_GUI("microrheo")
            self.update_MPL("MPL_canvas")

    def drag_from_file(self):
        """b(0) from a viscous drag measurement loaded in pyAF (chosen if
        there are several), shown with its fit before it is used."""
        files = viscous_drag.drag_files()
        if not files:
            gui_tools.show_success(
                self, "No viscous drag measurement loaded",
                details="Load the file of the viscous drag measurement in "
                "pyAF (Data tab): modulations of the cantilever at several "
                "heights above the surface, without contact (JPK force map "
                "with modulation segments separated by retract segments).",
                warning=True)
            return
        # Always chosen by the user (also when there is one)
        names = [file_label(dt) for dt in files]
        name, ok = QtWidgets.QInputDialog.getItem(
            self, "Viscous drag", "Viscous drag measurement :", names, 0,
            False)
        if not ok:
            return
        data = files[names.index(name)]
        # With the piezo correction of the sample file (set it first)
        sample = shared.exp.current_data
        piezo = sample if sample.microrheo_piezo else None
        result = viscous_drag.compute(data, piezo)
        dialog = ViscousDragDialog(self, data, result, sample=sample)
        if dialog.exec_() and numpy.isfinite(result["b0"]):
            value = float("%.4g" % (result["b0"] * 1e6))
            files = shared.exp.list if shared.exp.apply_to_all_compute \
                else [shared.exp.current_data]
            for dt in files:
                dt.microrheo_drag = value
            self.update_GUI("microrheo")
            self.update_MPL("MPL_canvas")

    @staticmethod
    def radius_notes(calc_list):
        """Message for the files indented deeper than the radius of the tip
        with the paraboloid (Hertz is valid for an indentation much smaller
        than the radius), None if there is none."""
        lines = ""
        for calc_id in calc_list:
            data = shared.exp.list[calc_id]
            if not data.stiffness_calculated or \
                    data.used_stiffness_model_selected not in \
                    contact_models.RADIUS_LIMITED:
                continue
            depth = numpy.asarray(data.indentation[:], dtype=float)
            depth = depth[numpy.isfinite(depth)]
            nbr = int(numpy.sum(depth > data.used_tip_radius))
            if nbr:
                lines += "<li>" + file_label(data) + " : " + \
                    str(nbr) + " pixels out of " + str(depth.size) + \
                    " (median indentation " + \
                    str(int(round(float(numpy.median(depth))))) + " nm, tip " \
                    "radius " + str(data.used_tip_radius) + " nm)</li>"
        if not lines:
            return None
        return "The indentation is larger than the tip radius: the " \
            "paraboloid (Hertz) is valid for an indentation much smaller " \
            "than the radius. Use the Sneddon (Sphere) model for a bead, or " \
            "a blunted cone or pyramid for a sharp tip.<ul>" + lines + "</ul>"

    def finish_computation(self, warnings=None, notes=None, info=None):
        """At the end of the computation, display the errors and refresh the GUI.

        notes: text added to the final message (displayed as a warning).
        info: text of the final message without warning.
        """
        # Indentation range not attained on some pixels (final message)
        range_text = None

        if warnings is not None:
            stop_size_errors = warnings[0]
            no_segment_errors = warnings[1]
            smoothing_warnings = warnings[2]
            range_warnings = warnings[3] if len(warnings) > 3 else []

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
                    dt = shared.exp.list[item[0]]
                    dt.stiffness_calculated = False
                    text = "<li>" + file_label(dt)
                    # Baseline window of the point of contact longer than
                    # the curve: no point of contact, no fit
                    window = (dt.fitparam_poc_skip_start or 0) + \
                        (dt.fitparam_poc_fit_length or 0)
                    if dt.ramp_size and window >= 0.8 * dt.ramp_size:
                        text += " : the baseline window of the point of " \
                            "contact (skip + length = " + \
                            misc_tools.format_number(window) + \
                            " nm) is longer than the curve (ramp " + \
                            misc_tools.format_number(dt.ramp_size) + \
                            " nm): reduce it in the Initial point of " \
                            "contact box"
                    list_no_segment_errors += text + "</li>"
                else:
                    found_error_list[i] = True

            list_smoothing_warnings = ""
            list_fit_warnings = ""
            text = ""
            for item in smoothing_warnings:
                if item:
                    name = shared.exp.list[item[0]].filename
                    for values in item[1]:
                        # values[0] is "smoothing_error" or "fit_error"
                        # the curves are stored in values[1]
                        line = "<li>" + name + " : " + values[1] + "</li>"
                        if values[0] == "fit_error":
                            list_fit_warnings += line
                        else:
                            list_smoothing_warnings += line

            if list_stop_size_error != "":
                text = "The 'stop' value is too big ! The maximum value you \
                        can use is : <ul>" + list_stop_size_error + "</ul>"

            if list_no_segment_errors != "":
                text = "No elasticity value was found (no point of " \
                    "contact, or no segment: check the start and stop " \
                    "values). Files :<ul>" + list_no_segment_errors + "</ul>"

            if list_smoothing_warnings != "":
                text = "Some curves seem corrupted and the smoothing was not \
                        applied for these curves.<br><br>\
                        The following curves had a problem : <ul>" + \
                    list_smoothing_warnings + "</ul>"

            if list_fit_warnings != "":
                fit_text = "The fit failed for some curves (no contact or \
                        not enough points in the fit range). No stiffness \
                        value was stored for these curves (displayed as \
                        NaN).<br><br>The following curves had a problem : \
                        <ul>" + list_fit_warnings + "</ul>"
                if list_smoothing_warnings != "":
                    text += "<br>" + fit_text
                else:
                    text = fit_text

            if (list_stop_size_error != "" or list_no_segment_errors != ""
                    or list_smoothing_warnings != ""
                    or list_fit_warnings != "") \
                    and not consts.UNIT_TESTING:
                # Create a message box and display it
                msg = QtWidgets.QMessageBox()
                msg.setText("Error")
                msg.setInformativeText(text)
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
                msg.exec_()

            # The pixels are computed up to their maximum indentation when the
            # selected range is not attained: the user is informed in the
            # final message
            if range_warnings:
                lines = ""
                for calc_id, (stop, nbr, total) in range_warnings:
                    name = shared.exp.list[calc_id].filename
                    lines += "<li>" + name + " : selected indentation range " \
                        + str(stop) + " nm was not attained on " + str(nbr) \
                        + " pixels out of " + str(total) + "</li>"
                range_text = "The stiffness was computed up to the maximum " \
                    "indentation for these pixels.<ul>" + lines + "</ul>"

        else:
            # No errors
            found_error_list = [True]

        if notes:
            range_text = notes if range_text is None else range_text + notes

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
            self.update_GUI("visco")
            self.update_GUI("microrheo")
            self.update_GUI("get_lr")
            # Update the results tab (especially the plot, so that it is
            # initialized and not empty)
            widgets_list.widget_results.update_widget()
            # Update the table with the data
            stat_tools.get_values()
            widgets_list.widget_results_single.reset_table()
            stat_tools.fetch_group_data()
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
                details = range_text
                if info:
                    details = info if details is None else details + info
                gui_tools.show_success(self, "Calculation complete",
                                       details,
                                       warning=range_text is not None)

    def button_clicked(self, button):
        """Method called whenever a button is clicked."""
        # Set the focus so that all the inputs which have been modified but
        # not validated will save their values.
        self.setFocus()

        data = shared.exp.list[shared.exp.id_selected]

        if button == "load_parameters":
            # Settings of an exported parameters file (see
            # widgets.load_parameters), shown in the inputs
            if load_parameters(self):
                self.update_widget()

        elif button == "button_rand":
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
            range_warnings = []
            stop_size_errors = []
            no_segment_errors = []
            one_computed = False
            brush_info = []
            refit_disabled = []

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
                if self.compute_stiffness.refit_disabled:
                    refit_disabled.append(
                        os.path.basename(str(data.filename)))

                if self.compute_stiffness.error == "Error: Stop size":
                    val = str(int(self.compute_stiffness.error_value))
                    stop_size_errors.append([calc_id, val])
                    no_segment_errors.append([])
                    smoothing_warnings.append([])

                elif self.compute_stiffness.error == "Error: no segments":
                    no_segment_errors.append([calc_id, None])
                    stop_size_errors.append([])
                    smoothing_warnings.append([])

                else:
                    # At least one file was computed without error
                    one_computed = True

                    stop_size_errors.append([])
                    no_segment_errors.append([])

                    # Store smoothing warnings to display a message at the end
                    # of the computation
                    warnings = self.compute_stiffness.list_of_smoothing_errors
                    if warnings:
                        smoothing_warnings.append([calc_id, warnings])
                    else:
                        smoothing_warnings.append([])

                    # Pixels where the selected indentation range was not
                    # attained
                    not_attained = self.compute_stiffness.range_not_attained
                    if not_attained:
                        range_warnings.append([calc_id, not_attained])

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

                    # Brush model, with the points of contact of this fit
                    if data.brush_model and data.used_perform_fit and \
                            brush.available(data, used=True):
                        found = brush.compute(calc_id)
                        brush_info.append(
                            file_label(data) + " : " +
                            str(found) + " of " +
                            str(data.nbr_pixels_x * data.nbr_pixels_y) +
                            " curves with a brush")

                # Reset the error
                self.compute_stiffness.error = False

            if one_computed:
                if shared.exp.results_type is None:
                    shared.exp.results_type = "stiffness"

            warnings = [stop_size_errors,
                        no_segment_errors,
                        smoothing_warnings,
                        range_warnings]

            if refit_disabled:
                # Settings changed for the brush model
                widgets_list.widget_fit_param_pocs.update_widget()
                brush_info.append(
                    "Refit of the point of contact set to 0 for the brush "
                    "model")

            self.finish_computation(warnings,
                                    notes=self.radius_notes(calc_list),
                                    info="<ul>" + "".join(
                                        "<li>" + text + "</li>"
                                        for text in brush_info) + "</ul>"
                                    if brush_info else None)

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

        elif button == "compute_visco":
            # Update all the parameters if needed
            apply_to_all.apply_to_all("calc_options")

            calc_list = list(range(len(shared.exp.list))) \
                if shared.exp.calc_all else [shared.exp.id_selected]
            # Overwrite calc_list with a list defined externally
            if shared.force_calc_list is not None:
                calc_list = shared.force_calc_list

            skipped = []
            for calc_id in calc_list:
                data = shared.exp.list[calc_id]
                reason = visco_main.check_file(data)
                if reason is not None:
                    skipped.append(file_label(data) + " : " +
                                   reason)
                    continue
                self.logger.debug("Viscoelasticity %s", data.filename)
                self.compute_visco.compute(calc_id)
                # Add a result to the result's tab if needed
                experiment.add_result(calc_id, data.filename)
                if shared.exp.results_type is None:
                    shared.exp.results_type = "visco_e0"
                if self.compute_visco.nbr_failed:
                    skipped.append(
                        file_label(data) + " : the fit failed "
                        "for " + str(self.compute_visco.nbr_failed) +
                        " curves (NaN)")

            notes = None
            if skipped:
                notes = "<ul>" + "".join("<li>" + text + "</li>"
                                         for text in skipped) + "</ul>"
            self.finish_computation(notes=notes)

        elif button == "microrheo_drag_file":
            self.drag_from_file()

        elif button == "microrheo_piezo":
            self.piezo_from_curve()

        elif button == "microrheo_piezo_clear":
            files = shared.exp.list if shared.exp.apply_to_all_compute \
                else [data]
            for dt in files:
                dt.microrheo_piezo = None
            self.update_GUI("microrheo")
            self.update_MPL("MPL_canvas")

        elif button == "compute_microrheo":
            # Update all the parameters if needed
            apply_to_all.apply_to_all("calc_options")

            calc_list = list(range(len(shared.exp.list))) \
                if shared.exp.calc_all else [shared.exp.id_selected]
            # Overwrite calc_list with a list defined externally
            if shared.force_calc_list is not None:
                calc_list = shared.force_calc_list

            skipped = []
            lines = []
            for calc_id in calc_list:
                data = shared.exp.list[calc_id]
                reason = microrheology.check_file(data)
                if reason is not None:
                    skipped.append(file_label(data) + " : " + reason)
                    continue
                self.logger.debug("Microrheology %s", data.filename)
                nbr, not_analysed = microrheology.compute(calc_id)
                experiment.add_result(calc_id, data.filename)
                if shared.exp.results_type is None:
                    shared.exp.results_type = "microrheo_storage"
                shown = len(microrheology.shown_labels(data))
                total = len(microrheology.frequencies(data))
                lines.append(file_label(data) + " : " + str(nbr) +
                             " curves, " + str(shown) + " frequencies" +
                             ("" if shown == total else " (" +
                              str(total - shown) + " above " +
                              format(data.microrheo_max_freq, "g") +
                              " Hz not analysed)"))
                # Curves or frequencies not analysed (NaN in the maps)
                for reason, count in not_analysed.items():
                    if ": " in reason:
                        # A frequency of the curves: "10 Hz: reason"
                        label, why = reason.split(": ", 1)
                        skipped.append(file_label(data) + " : " +
                                       label + " not analysed for " +
                                       str(count) + " curves (" + why + ")")
                    else:
                        skipped.append(file_label(data) + " : " +
                                       str(count) + " curves not analysed (" +
                                       reason + ")")
                # Piezo correction outside its frequencies
                outside = microrheology.piezo_outside(data)
                if outside:
                    piezo = data.used_microrheo_piezo
                    skipped.append(
                        file_label(data) + " : piezo correction of " +
                        ", ".join(format(f, "g") for f in outside) +
                        " Hz taken at the nearest frequency of the "
                        "characterization (" + format(min(
                            piezo["frequencies"]), "g") + "\u2013" +
                        format(max(piezo["frequencies"]), "g") + " Hz)")

            notes = None
            if skipped:
                notes = "<ul>" + "".join("<li>" + text + "</li>"
                                         for text in skipped) + "</ul>"
            info = None
            if lines and not skipped:
                info = "<ul>" + "".join("<li>" + text + "</li>"
                                        for text in lines) + "</ul>"
            elif lines:
                notes = "<ul>" + "".join("<li>" + text + "</li>"
                                         for text in lines + skipped) + \
                    "</ul>"
            self.finish_computation(notes=notes, info=info)

        elif button == "compute_breakthrough":
            # Update all the parameters if needed
            apply_to_all.apply_to_all("calc_options")

            calc_list = list(range(len(shared.exp.list))) \
                if shared.exp.calc_all else [shared.exp.id_selected]
            # Overwrite calc_list with a list defined externally
            if shared.force_calc_list is not None:
                calc_list = shared.force_calc_list

            lines = []
            for calc_id in calc_list:
                data = shared.exp.list[calc_id]
                if data.is_single:
                    continue
                self.logger.debug("Breakthroughs %s", data.filename)
                found = breakthrough.compute(calc_id)
                experiment.add_result(calc_id, data.filename)
                if shared.exp.results_type is None:
                    shared.exp.results_type = "breakthrough_force"
                total = int(numpy.sum(numpy.asarray(
                    data.discarded_curves) == 0))
                lines.append(file_label(data) + " : " +
                             str(found) + " of " + str(total) +
                             " curves with a breakthrough")

            info = None
            if lines:
                info = "<ul>" + "".join("<li>" + text + "</li>"
                                        for text in lines) + "</ul>"
            self.finish_computation(info=info)

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

        elif button == "stiffness_params":
            value = self.GRP_stiffness_params.checkedId()
            data.stiffness_mode = value
            # Tomography: slope of the segments. Fitting: Hertz by default
            if value == 0:
                self.set_model(data, 3)
                data.tomography = 1
            elif data.stiffness_model_selected == 3:
                self.set_model(data, 0)
            self.update_GUI("grid_model_options")
            self.update_GUI("stiffness_model_selected")
            self.update_widget()

        elif button == "next_zoom_element":
            if shared.exp.compute_type == "stiffness":
                # Only one element, the POC
                val = None
            elif shared.exp.compute_type in ("stiffness_corr",
                                             "viscoelasticity"):
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
                # The events and the joc
                points = self.MPL_curve.canvas.loading_rates_zoom_points
                val = len(points) if points else None
            elif shared.exp.compute_type == "breakthrough":
                # The breakthroughs
                events = self.MPL_curve.canvas.breakthrough_events
                val = len(events) if events else None

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
            shared.exp.calc_all = self.TG_compute_all.isChecked()

        elif checkbox == "brush_model":
            data.brush_model = self.TG_brush.isChecked()
            self.IN_brush_threshold_pct.setEnabled(data.brush_model)
            # No refit of the baseline of the point of contact (it would
            # move the baseline, the zero of the force, into the brush)
            if data.brush_model:
                data.fitparam_poc_refit_option = 0
                data.fitparam_poc_refit_times = 0
                widgets_list.widget_fit_param_pocs.update_widget()
            # Threshold on the preview
            self.update_MPL("MPL_canvas")

        elif checkbox == "apply_to_all":
            # Ask the user if he is sure he wants to overwrite his parameters
            do = True
            value = self.TG_apply_to_all.isChecked()
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
                    self.TG_apply_to_all.setChecked(False)
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

        # elif checkbox == "tomography":
            # state = self.CB_tomography.isChecked()
            # shared.exp.list[shared.exp.id_selected].tomography = state

        elif checkbox == "visco_bec":
            data.visco_bec = self.CB_visco_bec.isChecked()
            self.update_GUI("visco")
            self.update_GUI("microrheo")
            self.update_MPL("MPL_canvas")

        elif checkbox == "microrheo_raw":
            data.microrheo_preview_raw = self.TG_microrheo_raw.isChecked()
            self.update_MPL("MPL_canvas")

        elif checkbox == "microrheo_piezo_amplitude":
            data.microrheo_piezo_amplitude = \
                self.CB_microrheo_piezo_amplitude.isChecked()
            self.update_MPL("MPL_canvas")

        elif checkbox == "visco_drag_correction":
            data.visco_drag_correction = self.CB_visco_drag.isChecked()
            # Baseline of the retraction in the preview
            self.update_MPL("MPL_canvas")

        elif checkbox == "visco_beta_auto":
            data.visco_beta_auto = self.CB_visco_beta_auto.isChecked()
            self.update_GUI("visco")
            self.update_GUI("microrheo")
            self.update_MPL("MPL_canvas")

        elif checkbox == "zoom_in":
            data.auto_zoom_in = self.TG_zoom_in.isChecked()
            tp = shared.exp.compute_type
            if tp in ("events", "work_and_rupture_force", "loading_rates"):
                self.BT_next_zoom_element.setEnabled(True)
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

            else:
                pass
                # self.IN_start.setEnabled(True)
                # self.IN_stop.setEnabled(True)
                # self.IN_step.setEnabled(True)
                # self.CB_tomography.setEnabled(True)



    def mesh_canvas_press(self, event):
        """Matplotlib method, called when the user clicks on the meshgrid."""
        # event.button == 1 = left click

        if event.button == 1 and event.inaxes == self.MPL_meshgrid.canvas.axes:
            xpos, ypos = misc_tools.get_position_on_meshgrid(event)

            widgets_list.widget_main.change_curve(xpos, ypos)
