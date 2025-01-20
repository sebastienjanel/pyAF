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

"""Main plotting module, called from the GUI to display the data."""

from .. import shared
from .. import widgets_list
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib
from .meshgrids.piezo_meshgrid import PlotPiezoMeshGrid
from .meshgrids.topo_meshgrid import PlotTopoMeshGrid
from .meshgrids.stiffness_meshgrid import PlotStiffnessMeshGrid
from .meshgrids.stiffness_corr_meshgrid import PlotStiffnessCorrMeshGrid
from .meshgrids.stiffness_slice_meshgrid import PlotStiffnessSliceMeshGrid
from .meshgrids.stiffness_slice_corr_meshgrid import \
    PlotStiffnessSliceCorrMeshGrid
from .meshgrids.work_meshgrid import PlotWorkMeshGrid
from .meshgrids.rupture_force1_meshgrid import PlotRuptureForce1MeshGrid
from .meshgrids.rupture_force2_meshgrid import PlotRuptureForce2MeshGrid
from .meshgrids.events_per_curve_meshgrid import PlotEventsPerCurveMeshGrid
from .meshgrids.small_meshgrid import PlotSmallMeshGrid
from .meshgrids.small_meshgrid_compute import PlotSmallMeshGridCompute
from .meshgrids.flatten_meshgrid import PlotFlattenMeshGrid
from .meshgrids.empty_meshgrid import PlotEmptyMeshGrid
from .meshgrids.tilt_check import PlotTiltCheckMeshGrid
from .results.single import PlotSingleResults
from .results.groups import PlotGroupedResults
from .results.experiment import PlotExperimentResults
from .results.stats import PlotExperimentStats
from .slices.plot_slices import PlotSlices
from .events_per_scan.plot_events_per_scan import PlotEventsPerScan
from .profiles.plot_profiles import PlotProfiles
from .curves.tilt_correction import PlotTiltPreview

# Data tab ====================================================================

from .curves.data.data import PlotDataCurve
from .curves.data.data_a import PlotDataCurveA
from .curves.data.data_b import PlotDataCurveB

# Compute tab =================================================================

# Stiffness
from .curves.compute.stiffness import ComputeStiffness
from .curves.compute.work_and_rupture_force import ComputeWorkAndRuptureForce

# Events
from .curves.compute.events_force_curve_kernel import ComputeEventsKernel
from .curves.compute.events_force_curve_msf import ComputeEventsMSF
from .curves.compute.events_detection_kernel import \
    ComputeEventsDetectionKernel
from .curves.compute.events_detection_msf import ComputeEventsDetectionMSF

# Loading rates
from .curves.compute.loading_rates import ComputeLoadingRates

# Results tab =================================================================

# Deflexion - extension
from .curves.results.defl_ext import ResultsDeflExt
# Force - extension
from .curves.results.force_ext import ResultsForceExt
# Force - distance
from .curves.results.force_dist import ResultsForceDist
# Indentation
from .curves.results.indentation import ResultsIndentation
# Residuals
from .curves.results.residuals import ResultsResiduals

# Events
from .curves.results.force_events import ResultsForceEvents


params = {"legend.fontsize": 10}
matplotlib.rcParams.update(params)


class PYAFPlot(FigureCanvasQTAgg):
    """Canvas widget for the plots."""

    def __init__(self, parent, plot_type, canvas=None, sizes=None,
                 save_path=None, noborder=False, data_id=None, min_size=None):
        """Create the canvas widget.

        If the canvas argument is None, the plot has to be opened in a
        separate window. As there is no size defined, I set some
        arbitrary size to create the FigureCanvasQTAgg. A new figure
        (self.fig) will be created in the MainPlot class.
    """

        self.parent = parent
        self.plot_type = plot_type
        self.save_path = save_path
        self.noborder = noborder
        self.data_id = data_id
        self.canvas = canvas

        self.curve_type = None
        self.mesh_type = None
        self.data = None

        if plot_type == "small_meshgrid_for_correction":
            self.window_color_level = 2
        elif plot_type == "multimeshgrid":
            self.window_color_level = 1
        else:
            self.window_color_level = 0

        if self.canvas is None:
            self.single = True
            sizes = [6, 6, 72]
        else:
            self.single = False

        # frameon = False = no background color
        # (BUG while False; let it on True ?)
        # Happens only on OS X I think
        # Perhaps this is fixed now ? Should be reported to the matplotlib
        # guys.

        self.fig = Figure(figsize=(sizes[0], sizes[1]), dpi = sizes[2], frameon = False)

        super().__init__(self.fig)

        # Refresh delay for the resizeEvent
        self.timeout = 150  # ms
        self.resize_timer = QtCore.QTimer(self)
        self.resize_timer.timeout.connect(self.delayed_update)
        self.timer_is_running = False
        self.is_not_first = False
        self.tab_id = None
        self.timer2 = None
        self.timer2a = None

        if min_size is not None:
            # Force the usage of a minimun size (used for multimeshgrids)
            self.setMinimumSize(min_size[0], min_size[1])

        # the widget can't shrink, but should be allowed
        # to grow as much as possible.
        super().setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding,
                              QtWidgets.QSizePolicy.MinimumExpanding)
        FigureCanvasQTAgg.updateGeometry(self)
        self.setParent(canvas)

        # Directly plot the figure in a new window
        if self.canvas is None:
            self.update_plot()

    def update_plot(self, mode=None):
        """Refreshes the plot."""
        if self.data_id is None:
            # Get the current data
            dt = shared.exp.current_data
        else:
            # Get the data associated with data_id
            dt = shared.exp.list[self.data_id]

        c_type = dt.display_curve_type

        # Called by the plotting classes
        self.curve_type = c_type
        # Replace the curve type by defl-ext in case of the fit previews
        # or the data curve
        pt = self.plot_type
        if pt == "fit_preview":
            self.curve_type = "defl_ext"

        self.data = dt

        # For multimeshgrid, use the multimeshgrid's type
        if self.plot_type == "multimeshgrid":
            self.mesh_type = shared.exp.multi_meshgrid_type
        else:
            self.mesh_type = dt.meshgrid_type

        mesh = self.mesh_type

        # Clear the figure
        self.fig.clf()

        p_type = self.plot_type
        plot_selected = shared.exp.hist_single_plot_selected

        # Plot the asked figure
        if p_type == "meshgrid_tilt_check":
            self.canvas = PlotTiltCheckMeshGrid(self)

        if (p_type in ("meshgrid", "multimeshgrid", "meshgrid_data") and
                dt.is_single is False):
            if mesh == "piezo" or p_type == "meshgrid_data":
                self.canvas = PlotPiezoMeshGrid(self)
            elif mesh == "topo" and dt.stiffness_calculated:
                self.canvas = PlotTopoMeshGrid(self)
            elif mesh == "stiffness" and dt.stiffness_calculated:
                self.canvas = PlotStiffnessMeshGrid(self)
            elif mesh == "work" and dt.work_and_rupture_force1_calculated:
                self.canvas = PlotWorkMeshGrid(self)
            elif mesh == "rupture_force" and dt.jocs1_indices is not None:
                self.canvas = PlotRuptureForce1MeshGrid(self)
            elif mesh == "events_rupture_force" and dt.events_calculated:
                self.canvas = PlotRuptureForce2MeshGrid(self)
            elif mesh == "events_per_curve" and dt.events_calculated:
                self.canvas = PlotEventsPerCurveMeshGrid(self)
            elif mesh == "stiffness_slice" and dt.stiffness_calculated:
                self.canvas = PlotStiffnessSliceMeshGrid(self)
            elif mesh == "stiffness_corr" and dt.stiffness_corrected:
                self.canvas = PlotStiffnessCorrMeshGrid(self)
            elif mesh == "stiffness_corr_slice" and dt.stiffness_corrected:
                self.canvas = PlotStiffnessSliceCorrMeshGrid(self)
            else:
                # Plot an empty meshgrid
                self.canvas = PlotEmptyMeshGrid(self)
        elif p_type == "small_meshgrid" and dt.is_single is False:
            self.canvas = PlotSmallMeshGrid(self)
        elif p_type == "small_meshgrid_compute" and dt.is_single is False:
            self.canvas = PlotSmallMeshGridCompute(self)
        elif (p_type in ("small_meshgrid_for_correction",
                         "small_meshgrid_for_correction_original")
              and dt.is_single is False):
            self.canvas = PlotFlattenMeshGrid(self)
        elif (p_type in ("curve_results", "fit_preview",
                         "tilt_curve_preview", "curve_data", "curve_dataa", "curve_datab")):
            ct = shared.exp.compute_type
            if p_type == "curve_data":
                self.canvas = PlotDataCurve(self)
            elif p_type == "curve_dataa":
                self.canvas = PlotDataCurveA(self)
            elif p_type == "curve_datab":
                self.canvas = PlotDataCurveB(self)
            elif p_type == "tilt_curve_preview":
                self.canvas = PlotTiltPreview(self)
            elif p_type == "curve_results":
                if c_type == "defl_ext":
                    self.canvas = ResultsDeflExt(self)
                elif c_type == "force_ext":
                    self.canvas = ResultsForceExt(self)
                elif c_type == "force":
                    self.canvas = ResultsForceDist(self)
                elif c_type == "indentation":
                    self.canvas = ResultsIndentation(self)
                elif c_type == "results_force_events":
                    self.canvas = ResultsForceEvents(self)
                elif c_type == "residuals":
                    self.canvas = ResultsResiduals(self)
            elif p_type == "fit_preview":
                if ct == "stiffness":
                    self.canvas = ComputeStiffness(self)
                elif ct == "work_and_rupture_force":
                    self.canvas = ComputeWorkAndRuptureForce(self)
                elif ct == "events":
                    if dt.events_curve_type == "curve":
                        if dt.events_algorithm == "kernel":
                            self.canvas = ComputeEventsKernel(self)
                        elif dt.events_algorithm == "msf":
                            self.canvas = ComputeEventsMSF(self)
                    elif dt.events_curve_type == "detection":
                        if dt.events_algorithm == "kernel":
                            self.canvas = ComputeEventsDetectionKernel(self)
                        elif dt.events_algorithm == "msf":
                            self.canvas = ComputeEventsDetectionMSF(self)
                elif ct == "loading_rates":
                    self.canvas = ComputeLoadingRates(self)
        elif p_type == "results_single":
                self.canvas = PlotSingleResults(self)
        elif p_type == "results_groups":
            self.canvas = PlotGroupedResults(self)
        elif p_type == "results_experiment":
            self.canvas = PlotExperimentResults(self)
        elif p_type == "results_stats":
            self.canvas = PlotExperimentStats(self)
        elif p_type == "meshgrid_slices":
            self.canvas = PlotSlices(self)
        elif p_type == "events_per_scan":
            self.canvas = PlotEventsPerScan(self)
        elif p_type == "profiles":
            self.canvas = PlotProfiles(self)
        else:
            pass

        # Call update on the widget with a delay to force a refresh.
        # Needed for Windows to prevent graphical bugs
        self.timer2a = QtCore.QTimer.singleShot(self.timeout, self.update)

        # When loading the GUI the first time, no need to refresh the blit,
        # this is done manually at the end of the loading process (see end
        # of __init__ method in the pyaf.py file)
        cond = widgets_list.widget_main.first_load is False

        if (p_type in ("meshgrid", "small_meshgrid_compute",
                       "meshgrid_data", "meshgrid_tilt_check")
                and cond and mode is None):
            self.refresh_blitting()

    def delayed_blit(self):
        """Save clean (empty) meshgrid background and refresh blit."""
        bbox = self.canvas.axes.bbox
        self.canvas.empty_bbox = self.copy_from_bbox(bbox)
        self.canvas.current_bbox = self.copy_from_bbox(bbox)

        self.canvas.update_blit("all")

    def refresh_blitting(self):
        """Save a clean background once a meshgrid has been updated.

        Put the saving and blit update methods in a single shot timer
        so that there is a slight delay (100 ms), which lets time to the PyQt
        GUI to be refreshed. If you don't do this the wrong data is buffered
        when saving the clean backgrounds. (especially when changing tabs).
        """
        if not shared.exp.current_data.is_single:
            self.timer2 = QtCore.QTimer.singleShot(100, self.delayed_blit)
            if self.is_not_first is False:
                self.is_not_first = True

    def resizeEvent(self, ev):
        """Reimplement resizeEvent to manage plot redrawing."""
        p_type = self.plot_type

        # Check if the resize is done on the current tab. If not it means we
        # changed from one tab to another, and the first thing which is called
        # (before the tabs_changed() call) is the resizeEvent here. So we have
        # to mak sure we don't update the plot here.
        cond = self.tab_id == widgets_list.widget_main.current_tab

        if (p_type in ("meshgrid", "small_meshgrid_compute",
                       "meshgrid_data", "meshgrid_tilt_check")
                and self.timer_is_running is False and
                self.is_not_first and cond):
            self.update_plot("no_blit")
            self.resize_timer.start(self.timeout)
            self.timer_is_running = True

        super().resizeEvent(ev)

    def delayed_update(self):
        """Update the plots after a certain delay.

        Used to avoid blitting bugs.
        """
        self.timer_is_running = False
        self.resize_timer.stop()
        self.update_plot()
