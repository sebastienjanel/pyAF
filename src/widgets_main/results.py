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

"""Tab with the meshgrids and the results of the fits."""

import os
import math
from .. import widgets_list
from .. import shared
import logging
from .. import consts
from qtpy import QtGui, QtCore, QtWidgets
from ..tools import math_tools
from ..tools import apply_to_all
from ..tools import misc_tools
from ..utils import meshgrid_menu
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import PYAFToggle
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
from ..plots.PYAFPlot import PYAFPlot, curve_type_available
from ..tools.PYAFWidget import PYAFWidget
from ..tools import events_refresher
from ..widgets import plot_options
from ..compute_tools import microrheology

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
        # Room for the title and the rows of options
        self.options_boxes_height = 185

        self.cursor_mode = "select"
        self.parent = parent

        # Load the settings
        self.settings = QtCore.QSettings()

        self.current_profile_erase = None
        # A profile is being drawn (the last one of the profile list)
        self.drawing_profile = False

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
        # As wide as the map
        self.list_exp.setSizeAdjustPolicy(
            QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.list_exp.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                    QtWidgets.QSizePolicy.Fixed)
        for fv_data in shared.exp.list:
            self.list_exp.addItem(os.path.basename(str(fv_data.filename)))
        self.HL_choser.addWidget(self.list_exp)

        widgets_list.widget_progressbar.update()

        self.TG_app_to_all = PYAFToggle(self, "apply_to_all", "Apply to all")
        self.TG_app_to_all.setToolTip(
            "Same display options for all the files: map, curve and its "
            "options, filters of the events, units and colour scales")

        name = "widget_curve_selector_results"
        self.curve_selector = ResultsSelectorWidget(self, name)
        # self.curve_selector is added to the layout in update_widget()

        widgets_list.widget_progressbar.update()

        # --- Meshgrid --------------------------------------------------------

        # Maps: one list, by nature of the map; the maps not computed for
        # the file are greyed out
        # Nature of each map (elastic, viscoelastic, adhesion) and limits
        derived = (
            "\nDerived from E\u2080 and the fluidity (power-law model), not "
            "measured at this frequency: valid only if the sample follows a "
            "single power law,\nand within the time scales of the curves "
            "(about 1 / contact time to 1 / time between two points); "
            "outside, it is an extrapolation.")
        tips = {
            "piezo": "Height of the piezo at the end of the "
            "approach (includes the indentation).",
            "topo": "Height of the sample at the point of "
            "contact (zero force): topography without the indentation.",
            "stiffness": "Elastic: Young's modulus fitted on the "
            "approach (contact model of the Compute tab).\nOn a "
            "viscoelastic sample it depends on the indentation speed.",
            "stiffness_slice": "Elastic: elasticity of one "
            "indentation segment (depth).",
            "stiffness_corr": "Elastic: elasticity corrected for "
            "the stiff substrate under thin samples (bottom effect, Garcia "
            "& Garcia).",
            "stiffness_corr_slice": "Elastic: elasticity of one "
            "indentation segment, corrected for the substrate.",
            "visco_e0": "Viscoelastic: modulus at 1 s of the "
            "power-law relaxation E(t) = E\u2080 (t / 1 s)^-\u03b2\n(Ting fit "
            "of the approach and the retraction).",
            "visco_fluidity": "Viscoelastic: exponent \u03b2 of "
            "the power law, 0 elastic solid, 1 viscous liquid (cells: about "
            "0.1-0.3).",
            "visco_tan_delta": "Viscoelastic: loss tangent "
            "E\u2033 / E\u2032 = tan(\u03c0 \u03b2 / 2), the fluidity "
            "in rheology terms:\nbelow 1 solid-like, above 1 liquid-like; "
            "the same at all frequencies for a power law.",
            "visco_storage": "Viscoelastic: storage (elastic) "
            "modulus E\u2032 at the frequency of the Viscoelasticity tab."
            + derived,
            "visco_loss": "Viscoelastic: loss (dissipative) "
            "modulus E\u2033 at the frequency of the Viscoelasticity tab."
            + derived,
            "work": "Adhesion: work to detach the tip (area "
            "of the retraction curve below zero force).",
            "rupture_force": "Adhesion: maximum force of "
            "detachment.",
            "events_per_curve": "Adhesion: number of unbinding "
            "events of each curve.",
            "events_rupture_force": "Adhesion: largest force of the "
            "unbinding events of each curve.",
            "breakthrough_force": "Breakthrough: force Fb at the rupture of "
            "the supported lipid bilayer (first breakthrough of the "
            "approach).",
            "breakthrough_thickness": "Breakthrough: jump of the separation "
            "at the first breakthrough (thickness of the compressed "
            "bilayer).",
            "breakthrough_count": "Breakthrough: number of breakthroughs of "
            "each curve (stacked bilayers).",
            "brush_modulus": "Elastic: modulus of the cell body of the "
            "brush model, fitted above the threshold (brush compressed).",
            "brush_length": "Elastic: length L of the pericellular brush "
            "(brush model): its force decays as exp(-2\u03c0 h / L).",
            "brush_density": "Elastic: grafting density N of the brush "
            "(brush model): effective number of chains per \u00b5m\u00b2, "
            "from the amplitude of its force.",
            "microrheo_storage": "Microrheology: storage modulus G\u2032 "
            "at the frequency of the maps (Microrheology tab of Compute).",
            "microrheo_loss": "Microrheology: loss modulus G\u2033 at the "
            "frequency of the maps (Microrheology tab of Compute).",
            "microrheo_tan_delta": "Microrheology: loss tangent "
            "G\u2033 / G\u2032 at the frequency of the maps.",
            "microrheo_crossover": "Microrheology: frequency where G\u2033 "
            "becomes larger than G\u2032 (log interpolation between the "
            "analysed frequencies; none: no crossing).",
            "fit_r2": "R\u00b2 of the elasticity fit, on the fitted points "
            "after the point of contact.",
            "fit_poc_fixed": "Pixels where the point of contact could not be "
            "fitted (too stiff sample, failed fit, or point of contact "
            "moved too far):\nelasticity fitted with the detected point of "
            "contact (1)."}
        groups = (
            ("Heights", (("piezo", "Piezo height"),
                         ("topo", "Topography (zero force)"))),
            ("Elasticity", (("stiffness", "Elasticity"),
                            ("stiffness_slice", "Elasticity (Slice)"),
                            ("stiffness_corr", "Elasticity (BEC)"),
                            ("stiffness_corr_slice",
                             "Elasticity (Slice) (BEC)"),
                            ("fit_r2", "Fit R\u00b2"),
                            ("fit_poc_fixed", "Fit: fixed point of contact"),
                            ("brush_modulus", "Elasticity (brush model)"),
                            ("brush_length", "Brush length"),
                            ("brush_density", "Grafting density"))),
            ("Viscoelasticity", (("visco_e0", "E\u2080"),
                                 ("visco_fluidity", "Fluidity"),
                                 ("visco_tan_delta", "Loss tangent"),
                                 ("visco_storage", "Storage modulus E\u2032"),
                                 ("visco_loss", "Loss modulus E\u2033"))),
            ("Microrheology", (("microrheo_storage", "Storage modulus G\u2032"),
                               ("microrheo_loss", "Loss modulus G\u2033"),
                               ("microrheo_tan_delta", "Loss tangent"),
                               ("microrheo_crossover",
                                "Crossover frequency G\u2032 = G\u2033"))),
            ("Adhesion", (("work", "Detachment work"),
                          ("rupture_force", "Detachment force"),
                          ("events_per_curve", "Events per curve"),
                          ("events_rupture_force", "Event max. force"))),
            ("Breakthrough", (
                ("breakthrough_force", "Breakthrough force"),
                ("breakthrough_thickness", "Breakthrough jump"),
                ("breakthrough_count", "Breakthroughs per curve"))))
        self.list_maps = PYAFComboBox(self, "list_maps")
        self.list_maps.setMaxVisibleItems(40)
        # As wide as its longest map (and its list, open over the map)
        self.list_maps.setSizePolicy(QtWidgets.QSizePolicy.Fixed,
                                     QtWidgets.QSizePolicy.Fixed)
        model = self.list_maps.model()
        # Index of each map in the list
        self.map_items = {}
        for group, maps in groups:
            # Title of the group: not selectable, but enabled to have the
            # text colour of the theme (the disabled items are the maps not
            # computed)
            self.list_maps.addItem(group)
            item = model.item(self.list_maps.count() - 1)
            item.setFlags(QtCore.Qt.ItemIsEnabled)
            font = item.font()
            font.setBold(True)
            item.setFont(font)
            for name, label in maps:
                self.map_items[name] = self.list_maps.count()
                self.list_maps.addItem("    " + label, name)
                self.list_maps.setItemData(self.list_maps.count() - 1,
                                           tips[name],
                                           QtCore.Qt.ToolTipRole)
        self.list_maps.fit_to_contents()
        self.label_maps = QtWidgets.QLabel("Map")
        # Frequency of the microrheology maps (shown with these maps only)
        self.label_microrheo_frequency = QtWidgets.QLabel("Frequency:")
        self.list_microrheo_frequency = PYAFComboBox(
            self, "list_microrheo_frequency")
        self.list_microrheo_frequency.setToolTip(
            "Modulation frequency of the microrheology maps [Hz]")
        widgets_list.widget_progressbar.update()
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
        self.box_curves_options = QtWidgets.QGroupBox("Force curves")
        self.box_curves_options.setFixedHeight(self.options_boxes_height)

        # Type of curve: one list (the curves not available for the file
        # are greyed out)
        self.curve_types = (
            ("defl_ext", "Deflection - Z curve", None),
            ("force_ext", "Force - Z curve", None),
            ("force", "Force - Distance curve", None),
            ("indentation", "Indentation curve", None),
            ("residuals", "Residuals plot", None),
            ("visco", "Viscoelastic fit", None),
            ("brush", "Brush fit", "Brush model: force - distance curve of "
             "the approach with the fits of the cell body and of the brush"),
            ("microrheology", "Microrheology G\u2032, G\u2033",
             "Microrheology: storage and loss moduli vs the frequency of "
             "the modulations"),
            ("microrheology_map", "Microrheology G\u2032, G\u2033 (average)",
             "Microrheology, average of the file: median of the analysed "
             "curves with the interquartile range (band) at each frequency"),
            ("results_force_events", "Force curve (events)", None))
        self.list_curves = PYAFComboBox(self, "list_curves")
        for name, label, tip in self.curve_types:
            self.list_curves.addItem(label, name)
            if tip:
                self.list_curves.setItemData(self.list_curves.count() - 1,
                                             tip, QtCore.Qt.ToolTipRole)
        # As wide as the longest type (same width for all the types)
        self.list_curves.fit_to_contents()

        self.BT_group_trace_or_retrace = PYAFButtonGroup(
            self, "button_group_trace_or_retrace")
        self.radiobutton_trace = QtWidgets.QRadioButton("Approach curve")
        self.radiobutton_retrace = QtWidgets.QRadioButton("Retract curve")
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
        self.CB_breakthrough = PYAFCheckBox(
            self, "checkbox_breakthrough", "Breakthroughs")
        self.CB_breakthrough.setToolTip(
            "Breakthroughs of the bilayer on the approach: force peak and "
            "end of the jump, with Fb and the jump (Breakthrough tab of the "
            "Compute tab)")

        # Events
        self.BT_group_events = PYAFButtonGroup(self, "button_group_events")
        self.RBT_events_none = QtWidgets.QRadioButton("None")
        self.RBT_events_annotations = QtWidgets.QRadioButton("Annotations")
        self.BT_group_events.addButton(self.RBT_events_none, 0)
        self.BT_group_events.addButton(self.RBT_events_annotations, 1)
        self.CB_events_filter_dist = PYAFCheckBox(
            self, "checkbox_events_filter_dist", "Distance filter")
        self.CB_events_filter_dist.setToolTip(
            "Shades in grey the events removed by the distance filter of "
            "this file.\nTo set up the filter, right-click on the curve "
            "and choose \"Change distance filter\".")
        self.CB_events_display_joc = PYAFCheckBox(
            self, "checkbox_events_display_joc", "Display Joc")
        self.CB_events_display_fit_joc = PYAFCheckBox(
            self, "checkbox_events_display_fit_joc", "Display Joc fit")

        self.GL_curves_options = QtWidgets.QGridLayout()
        # Same spacing between all the rows (the style's default spacing
        # depends on the neighbouring widgets)
        self.GL_curves_options.setVerticalSpacing(4)
        self.GL_curves_options.setHorizontalSpacing(24)
        self.VL_curves_options = QtWidgets.QVBoxLayout()
        widgets_list.widget_progressbar.update()

        # Col 1: type of curve, and the parts of the curve below
        self.GL_curves_options.addWidget(self.list_curves, 0, 0)
        self.GL_curves_options.addWidget(self.radiobutton_both, 1, 0)
        self.GL_curves_options.addWidget(self.radiobutton_trace, 2, 0)
        self.GL_curves_options.addWidget(self.radiobutton_retrace, 3, 0)
        # Col 2
        self.GL_curves_options.addWidget(self.CB_poc, 0, 1)
        self.GL_curves_options.addWidget(self.CB_fit_poc, 1, 1)
        self.GL_curves_options.addWidget(self.CB_segments, 2, 1)
        self.GL_curves_options.addWidget(self.CB_display_fits, 3, 1)
        self.GL_curves_options.addWidget(self.CB_breakthrough, 4, 1)
        # Col 3
        self.GL_curves_options.addWidget(self.CB_joc, 0, 2)
        self.GL_curves_options.addWidget(self.CB_fit_joc, 1, 2)
        self.GL_curves_options.addWidget(self.CB_surface, 2, 2)
        self.GL_curves_options.addWidget(self.CB_force, 3, 2)
        # Col 4
        self.GL_curves_options.addWidget(self.RBT_events_none, 0, 3)
        self.GL_curves_options.addWidget(self.RBT_events_annotations, 1, 3)
        self.GL_curves_options.addWidget(self.CB_events_filter_dist, 2, 3)
        self.GL_curves_options.addWidget(self.CB_events_display_joc, 3, 3)
        self.GL_curves_options.addWidget(self.CB_events_display_fit_joc, 4, 3)
        # Apply to all (display options of the maps and the curves), in the
        # lower right corner, aligned with the checkboxes above
        self.GL_curves_options.addWidget(self.TG_app_to_all, 6, 3,
                                         QtCore.Qt.AlignLeft)

        self.VL_curves_options.addLayout(self.GL_curves_options)
        self.VL_curves_options.addStretch(1)

        self.HL_curves_options = QtWidgets.QHBoxLayout()
        self.HL_curves_options.addLayout(self.VL_curves_options)
        self.HL_curves_options.addStretch(1)
        self.box_curves_options.setLayout(self.HL_curves_options)
        widgets_list.widget_progressbar.update()

        # --- Set up page -----------------------------------------------------
        self.GL = QtWidgets.QGridLayout()

        # Above the map: file, map, and the navigation between the curves
        # (centered on the map)
        self.GL_choser = QtWidgets.QGridLayout()
        self.GL_choser.setVerticalSpacing(8)
        self.GL_choser.addWidget(QtWidgets.QLabel("File"), 0, 0)
        self.GL_choser.addWidget(self.list_exp, 0, 1)
        self.GL_choser.addWidget(self.label_maps, 1, 0)
        HL_map = QtWidgets.QHBoxLayout()
        HL_map.addWidget(self.list_maps)
        HL_map.addSpacing(12)
        HL_map.addWidget(self.label_microrheo_frequency)
        HL_map.addWidget(self.list_microrheo_frequency)
        HL_map.addStretch(1)
        self.GL_choser.addLayout(HL_map, 1, 1)
        self.GL_choser.setColumnStretch(1, 1)
        self.HL_navigation = QtWidgets.QHBoxLayout()
        self.HL_navigation.addStretch(1)
        # self.curve_selector is added here in update_widget()
        self.HL_navigation.addStretch(1)
        self.box_maps = QtWidgets.QGroupBox("Maps")
        self.box_maps.setLayout(self.GL_choser)
        VL_choser = QtWidgets.QVBoxLayout()
        VL_choser.addWidget(self.box_maps)
        VL_choser.addStretch(1)
        VL_choser.addLayout(self.HL_navigation)

        self.HL_buttons = QtWidgets.QHBoxLayout()
        self.HL_buttons.addWidget(self.BT_multi_meshgrids)
        self.HL_buttons.addWidget(self.BT_roi_manager)
        self.HL_buttons.addWidget(self.BT_profiles)
        self.HL_buttons.addWidget(self.BT_slices)
        self.HL_buttons.addWidget(self.BT_display_in_3D)
        self.HL_buttons.addWidget(self.BT_used_parameters)
        self.HL_buttons.addStretch(1)

        self.GL.addLayout(VL_choser, 0, 0)
        self.GL.addWidget(self.box_curves_options, 0, 1)
        self.GL.addWidget(self.MPL_canvas1, 1, 0)
        self.GL.addWidget(self.MPL_canvas2, 1, 1)
        self.GL.addLayout(self.HL_buttons, 2, 0, 1, 0)
        # The plots take the height of the window
        self.GL.setRowStretch(1, 1)

        widgets_list.widget_progressbar.update()

        self.VL.addLayout(self.GL)
        self.setLayout(self.VL)

        self.update_widget()

        widgets_list.widget_progressbar.update()

        self.box_maps.setFont(self.parent.smallfont)
        self.box_curves_options.setFont(self.parent.smallfont)

    def update_widget(self):
        """Update the widget."""
        data = shared.exp.list[shared.exp.id_selected]

        # For single files don't display the curve selector
        if data.is_single:
            self.curve_selector.setParent(None)
        else:
            if self.HL_navigation.indexOf(self.curve_selector) == -1:
                self.HL_navigation.insertWidget(1, self.curve_selector)
            self.curve_selector.update_widget()

        self.update_GUI("all")
        self.update_MPL("MPL_canvas1")
        self.update_MPL("MPL_canvas2")

    @staticmethod
    def match_curve_to_map(data):
        """Displays the curve corresponding to the map type, with its
        options: the approach force curve with the point of contact and the
        fits for the elasticity, the retraction force curve with the jump of
        contact, the work and the force for the detachment, the force curve
        of the events with their annotations, distance filter and jumps of
        contact for the events, and the default deflection curve for the
        heights. Returns True if the curve changed."""
        mt = data.meshgrid_type
        if mt in ("stiffness", "stiffness_slice", "stiffness_corr",
                  "stiffness_corr_slice", "fit_r2", "fit_poc_fixed") and \
                data.stiffness_calculated:
            data.display_curve_type = "force"
            data.display_trace_retrace = 0
            data.display_poc = True
            data.display_fit_poc = True
            data.display_segments = True
            data.display_fits_stiffness = True
        elif mt in ("work", "rupture_force") and \
                data.work_and_rupture_force1_calculated:
            data.display_curve_type = "force"
            data.display_trace_retrace = 1
            data.display_joc = True
            data.display_fit_joc = True
            data.display_surface = True
            data.display_force = True
        elif mt in ("events_per_curve", "events_rupture_force") and \
                data.events_calculated:
            data.display_curve_type = "results_force_events"
            data.display_trace_retrace = 1
            data.events_results_display_annotations = 1
            data.events_display_results_filter_dist = True
            data.events_results_display_joc = True
            data.events_results_display_fit_joc = True
        elif mt.startswith("visco_") and data.visco_calculated:
            data.display_curve_type = "visco"
        elif mt.startswith("brush_") and data.brush_calculated:
            data.display_curve_type = "brush"
        elif mt.startswith("microrheo_") and data.microrheo_calculated:
            data.display_curve_type = "microrheology"
        elif mt.startswith("breakthrough_") and data.breakthrough_calculated:
            # Approach in force vs distance, with the breakthroughs
            data.display_curve_type = "force"
            data.display_trace_retrace = 0
            data.display_breakthrough = True
        elif mt in ("piezo", "topo"):
            # Default display of the curves, without the results
            data.display_curve_type = "defl_ext"
            data.display_trace_retrace = 2
            data.display_poc = False
            data.display_fit_poc = False
            data.display_segments = False
            data.display_fits_stiffness = False
            data.display_joc = False
            data.display_fit_joc = False
            data.display_surface = False
            data.display_force = False
            data.events_results_display_annotations = 0
            data.events_display_results_filter_dist = False
            data.events_results_display_joc = False
            data.events_results_display_fit_joc = False
        else:
            return False
        return True

    def update_GUI(self, element):
        """Method to update different GUI elements."""
        data = shared.exp.list[shared.exp.id_selected]

        if element == "BTG_meshgrid_type" or element == "all":
            # Maps computed for the file (the others are greyed out)
            stiff = data.stiffness_calculated
            corr = data.stiffness_corrected
            slices = (stiff or corr) and data.get_used_indentation()[1] != 0
            adhesion = data.work_and_rupture_force1_calculated
            available = {
                "piezo": True, "topo": stiff, "stiffness": stiff,
                "stiffness_slice": stiff and slices,
                "stiffness_corr": corr,
                "stiffness_corr_slice": corr and slices,
                "fit_r2": data.has_fit_quality(),
                "fit_poc_fixed": data.has_fit_quality(),
                "work": adhesion, "rupture_force": adhesion,
                "events_per_curve": data.events_calculated,
                "events_rupture_force": data.events_calculated}
            model = self.list_maps.model()
            for name, index in self.map_items.items():
                if name.startswith("visco_"):
                    enabled = data.visco_calculated
                elif name.startswith("breakthrough_"):
                    enabled = data.breakthrough_calculated
                elif name.startswith("brush_"):
                    enabled = data.brush_calculated
                elif name.startswith("microrheo_"):
                    enabled = data.microrheo_calculated
                else:
                    enabled = available[name]
                model.item(index).setEnabled(bool(enabled))

            # Select the map
            self.list_maps.setCurrentIndex(
                self.map_items[data.meshgrid_type])

            self.update_microrheo_frequency()

        if element == "apply_to_all" or element == "all":
            self.TG_app_to_all.setChecked(shared.exp.apply_to_all_data)

        if element == "radio_buttons" or element == "all":
            if data.used_stiffness_model_selected == 3:
                labels = ("Apparent stiffness", "Apparent stiffness (Z-slice)")
            else:
                labels = ("Elasticity", "Elasticity (Slice)")
            for name, text in zip(("stiffness", "stiffness_slice"), labels):
                self.list_maps.setItemText(self.map_items[name], "    " + text)
            self.list_maps.fit_to_contents()

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
                # Slices of the elasticity in depth: only with segments
                # (Tomography, or files computed with segments before)
                _, step, tomography = data.get_used_indentation()
                depth = data.stiffness_calculated and (
                    not data.used_perform_fit or (tomography and step))
                self.BT_slices.setEnabled(bool(depth))
                self.BT_slices.setToolTip(
                    "" if depth else "Slices of the elasticity in depth: "
                    "available when the elasticity is computed with "
                    "Tomography\n(one value per indentation segment; "
                    "Fitting gives one value per curve)")

            else:
                self.BT_profiles.setEnabled(False)
                self.BT_slices.setEnabled(False)

        if element == "buttons_options_curves" or element == "all":
            # A curve type applied from another file can be missing here
            if not curve_type_available(data, data.display_curve_type):
                data.display_curve_type = "defl_ext"

            # Curves available for the file, and the displayed one
            model = self.list_curves.model()
            for index, (name, _, _) in enumerate(self.curve_types):
                model.item(index).setEnabled(
                    curve_type_available(data, name))
            self.list_curves.setCurrentIndex(
                self.list_curves.findData(data.display_curve_type))

            # Enable or disable radiobuttons for approach or retraction
            if data.display_curve_type == "force":
                self.radiobutton_trace.setEnabled(
                    data.stiffness_calculated or data.breakthrough_calculated)

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

            elif data.display_curve_type in ("visco", "brush",
                                             "microrheology",
                                             "microrheology_map"):
                # Approach and retraction in time, approach of the brush
                self.radiobutton_trace.setEnabled(False)
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

            # Breakthroughs, on the approach of the deflection and force
            # curves
            self.CB_breakthrough.setChecked(data.display_breakthrough)
            self.CB_breakthrough.setEnabled(
                data.breakthrough_calculated and dtr != 1 and
                tp in ("defl_ext", "force_ext", "force"))

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
            data.meshgrid_type = self.list_maps.currentData()
            self.logger.debug("Changing meshgrid type, %s", data.meshgrid_type)
            # The curve corresponding to the map
            curve_changed = self.match_curve_to_map(data)
            apply_to_all.apply_to_all("display_options")
            self.update_microrheo_frequency()
            # Update the GUI elements
            self.update_MPL("MPL_canvas1")
            if curve_changed:
                self.update_GUI("buttons_options_curves")
                self.update_MPL("MPL_canvas2")
            # Update the list of colortables and the min/max values
            if widgets_list.widget_meshgrid_options is not None:
                widgets_list.widget_meshgrid_options.update_widget()
            if widgets_list.widget_multimeshgrids is not None:
                widgets_list.widget_multimeshgrids.update_widget()
            # Profiles of the values of the new map
            if widgets_list.widget_profiles is not None:
                widgets_list.widget_profiles.update_widget()

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
            data.display_curve_type = self.list_curves.currentData()

            tp = data.display_curve_type
            self.logger.debug("Changing curve type, %s", tp)

            dtr = data.display_trace_retrace

            # Check if we can display trace or retrace
            if data.display_curve_type == "force":
                # Approach without the elasticity: breakthroughs
                if not (data.stiffness_calculated or
                        data.breakthrough_calculated) and \
                        (dtr == 0 or dtr == 2):
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
            self.select_map("stiffness")

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
        elif name == "list_curves":
            self.button_clicked("button_group_type_of_curve")
        elif name == "list_microrheo_frequency":
            self.microrheo_frequency_changed()
        elif name == "list_maps":
            if self.list_maps.currentData() is None:
                # Title of a group (keyboard): back to the displayed map
                self.update_GUI("BTG_meshgrid_type")
                return
            self.button_clicked("BTG_meshgrid_type")

    def update_microrheo_frequency(self):
        """Frequency selector of the microrheology maps: shown with these
        maps only, with the frequencies of the file."""
        data = shared.exp.current_data
        # Not for the crossover frequency (all the frequencies)
        microrheo = data.meshgrid_type.startswith("microrheo_") and \
            data.meshgrid_type != "microrheo_crossover" and \
            data.microrheo_calculated
        self.label_microrheo_frequency.setVisible(microrheo)
        self.list_microrheo_frequency.setVisible(microrheo)
        if not microrheo:
            return
        names = microrheology.shown_labels(data)
        label = microrheology.map_label(data)
        self.list_microrheo_frequency.blockSignals(True)
        self.list_microrheo_frequency.clear()
        for name in names:
            self.list_microrheo_frequency.addItem(name, name)
        if label is not None:
            self.list_microrheo_frequency.setCurrentIndex(names.index(label))
        self.list_microrheo_frequency.fit_to_contents()
        self.list_microrheo_frequency.blockSignals(False)

    def microrheo_frequency_changed(self):
        """Frequency of the microrheology maps chosen: maps at this
        frequency (no new computation), for all the files with this
        frequency with Apply to all."""
        data = shared.exp.current_data
        label = self.list_microrheo_frequency.currentData()
        files = [data]
        if shared.exp.apply_to_all_data:
            files += [dt for dt in shared.exp.list if dt is not data and
                      label in microrheology.shown_labels(dt)]
        for dt in files:
            dt.microrheo_frequency = label
            # New values: automatic scales of the moduli maps
            if dt.microrheo_calculated:
                for mt in ("microrheo_storage", "microrheo_loss"):
                    plot_options.autoscale_elasticity(dt, mt)
        self.update_MPL("MPL_canvas1")
        if widgets_list.widget_multimeshgrids is not None:
            widgets_list.widget_multimeshgrids.update_widget()

    def select_map(self, map_type):
        """Displays a map (as chosen in the list of the maps)."""
        self.list_maps.setCurrentIndex(self.map_items[map_type])
        self.button_clicked("BTG_meshgrid_type")

    def checkbox_clicked(self, name):
        """Called when clicking on a checkbox."""
        data = shared.exp.current_data

        if name == "apply_to_all":
            shared.exp.apply_to_all_data = self.TG_app_to_all.isChecked()
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

        elif name == "checkbox_breakthrough":
            data.display_breakthrough = self.CB_breakthrough.isChecked()
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
                # New profile, at the end of the list until the release
                data.profile_starts.append([xpos - 1, ypos - 1])
                data.profile_list.append([[xpos - 1], [ypos - 1]])
                data.profile_count = len(data.profile_list) - 1
                self.drawing_profile = True

            elif mode == "roi_add" or mode == "roi_erase":
                self.draw_roi(xpos, ypos)

            elif mode == "profile_erase":
                # Delete the first profile passing by this pixel
                for count, profile in enumerate(data.profile_list):
                    if [xpos - 1, ypos - 1] in \
                            [list(p) for p in zip(profile[0], profile[1])]:
                        del data.profile_list[count]
                        if count < len(data.profile_starts):
                            del data.profile_starts[count]
                        break
                data.profile_count = len(data.profile_list)
                self.current_profile_erase = None

                # Close profile widget if needed, else remove the profile
                w_profiles = widgets_list.widget_profiles
                if data.profile_count == 0 and w_profiles is not None:
                    widgets_list.widget_profiles.close()
                elif w_profiles is not None:
                    w_profiles.update_widget()

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
            # Position of the mouse on screen (the matplotlib event position
            # is in physical pixels, twice the Qt pixels on Retina screens)
            pos = QtGui.QCursor.pos()
            meshgrid_menu.create_meshgrid_menu(self, "results", pos)

    def canvas_release(self, event):
        """The canvas_release is called whenever the mouse click is released."""
        data = shared.exp.list[shared.exp.id_selected]

        button = event.button

        if not isinstance(widgets_list.widget_results.MPL_canvas1.canvas,
                          MainPlot):
            return  # Nothing to do if no real canvas

        # End of a profile (also when the mouse is released outside of the
        # map): kept if it has at least 3 pixels
        if self.drawing_profile and button == 1:
            self.drawing_profile = False
            ls = data.profile_list
            if ls and len(ls[-1][0]) < 3:
                del ls[-1]
                if data.profile_starts:
                    del data.profile_starts[-1]
            data.profile_count = len(ls)
            self.update_GUI("button_profiles_and_slices")
            # Finalize plotting
            self.MPL_canvas1.canvas.update_blit("all")
            # New profile in the profiles window
            if widgets_list.widget_profiles is not None and ls:
                widgets_list.widget_profiles.update_widget()

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

            if mode == "profile_add" and self.drawing_profile:
                # The profile being drawn is the last one: from its start to
                # the mouse
                data.profile_list[-1] = misc_tools.make_path(
                    list(data.profile_starts[-1]), [xpos - 1, ypos - 1])

                # Draw profile
                cv = self.MPL_canvas1.canvas
                cv.update_blit("profile_preview",
                               profile=data.profile_list[-1])

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
            if data.profile_list:
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
            # Position of the mouse on screen (the matplotlib event position
            # is in physical pixels, twice the Qt pixels on Retina screens)
            pos = QtGui.QCursor.pos()
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
                action = QtWidgets.QAction("Change distance filter", self)
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
        widget.exec()


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
