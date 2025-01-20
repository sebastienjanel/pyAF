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

"""
This module contains references to all the data and options of pyAF.

The Experiment class contains top level information about the files and global
options for pyAF. The most important variables are **list** and
**id_selected**.

There is also a getter to access the current data. You can access it like
this :
.. code-block:: python

    data = shared.exp.current_data
    data = shared.exp.list[shared.exp.id_selected] # Equivalent

"""

import numpy
import logging
from . import shared
from .single import SingleFile
from .tools import math_tools


class Experiment:
    """The experiment's data is stored in this object."""

    def __init__(self):
        # List of things wich are note saved in the pickle dump
        self.dontsave = ["list", "temp_file", "logger"]

        # The logger logs the user's action in a log file. See the
        # logger setup in the src/tools/misc_tools.py module.
        self.logger = logging.getLogger()
        self.logger.debug("Init Experiment")

        # A label that can be displayed near the title of the app.
        # Allow to give manually another title to the app, can be useful when
        # you open 2 or 3 apps, you can discriminate between them with this
        # label.
        self.main_title_label = ""

        # List is where each file (Data object, see Data class below) is saved.
        self.list = []
        self.id_selected = 0  # Id of the currently selected file

        # List of files with issues of missmatch between number of segments described
        # in the header and found while loading.
        self.segment_handling =[]


        # Reference to the currently opened hdf5 file.
        # Will be defined in src/load_and_save/load.py module
        self.temp_file = None

        # Number of cores, changed at load in src/load_and_save/load.py
        self.cores = 1

        # Compute all the files
        self.calc_all = False

        # Defines if options in pyAF are applied to all the files or a single
        # one. Find out more in the src/tools/apply_to_all.py module.
        self.apply_to_all_data = True
        self.apply_to_all_compute = True
        self.discard_corrupted_curves_applied_all = False

        # List containg Results and Groups objects (see class definitions
        # below).
        self.results_list = []
        self.groups_list = []
        self.conditions_list = []
        self.sample_groups_list = []

        # Selected statistical test
        self.selected_statistical_test = "Two-sample T-test"
        # List containing the pairs for performing statistical tests.
        self.sample_groups = []
        # Dictionary containing the results of statistical tests. Where -> {"statistical_test": [result]}
        self.statistical_test_results = {}

        # Add five empty groups by default
        for i in range(6):
            self.groups_list.append(ResultGroup(i))
            # i+1 to skip group 0, it's the None group
            if i != 0:
                self.groups_list[i].color = shared.colors_list[i - 1]

        # Add two empty conditions by default (Control + Treatment)
        for i in range(3):
            self.conditions_list.append(ConditionGroup(i))
            # i+1 to skip group 0, it's the None group
            if i != 0:
                self.conditions_list[i].color = shared.colors_list[i - 1]

        # Size of the matplotlibs fonts
        self.mpl_labels_font_size = 10
        self.mpl_tick_labels_font_size = 8

        # List of curves seen by the random algorithm
        self.parsed_random_curves = []

        # Last ten selected curves?
        self.last_ten_curves = [[0, 0]]
        self.pos_in_last_ten_curves = 0

        # Display approach/retraction legend on curve plot
        self.display_curve_legend = True

        # scatter or points
        self.hist_lr_single_display_mode = "scatter"
        self.hist_lr_groups_display_mode = "scatter"

        # Stat plot type
        self.stat_plot_type = "Boxplot"

        self.hist_single_bins_mode = "manual"
        self.hist_groups_bins_mode = "manual"
        self.hist_experiment_bins_mode = "manual"
        self.hist_single_bins = 50
        self.hist_groups_bins = 50
        self.hist_experiment_bins = 50

        self.hist_single_x_mode = "max"
        self.hist_single_y_mode = "auto"
        self.hist_groups_x_mode = "max"
        self.hist_groups_y_mode = "auto"
        self.hist_experiment_x_mode = "max"
        self.hist_experiment_y_mode = "auto"
        self.hist_single_max_x = 100
        self.hist_single_min_x = 0
        self.hist_single_max_y = 100
        self.hist_single_min_y = 0
        self.hist_groups_max_x = 100
        self.hist_groups_min_x = 0
        self.hist_groups_max_y = 100
        self.hist_groups_min_y = 0
        self.hist_experiment_max_x = 100
        self.hist_experiment_min_x = 0
        self.hist_experiment_max_y = 100
        self.hist_experiment_min_y = 0
        self.hist_remove_zeros = True
        self.hist_log = False
        self.hist_single_display_pdf_fit = True
        self.hist_groups_display_pdf_fit = True
        self.hist_values_or_frequencies = "values"
        self.hist_single_bw_mode = "auto"
        self.hist_groups_bw_mode = "auto"
        self.hist_experiment_bw_mode = "auto"
        self.hist_single_bw = 1
        self.hist_groups_bw = 1
        self.hist_experiment_bw = 1
        # Legends for histograms
        self.display_legend_hist_single = True
        self.display_legend_hist_groups = True
        self.display_legend_hist_experiment = True
        # Norm histograms
        self.norm_hist_single = True
        self.norm_hist_groups = True
        self.norm_hist_experiment = True

        self.hist_single_plot_selected = "hist"
        self.hist_groups_plot_selected = "hist"

        # Selected flatten order
        self.selected_flatten_order = 1

        # Data types (meshgrid types follow the same naming convention)
        # piezo
        # topography
        # stiffness
        # stiffness_corr
        # work
        # rupture_force
        # events_forces
        # events_per_curve
        # events_rupture_force
        # loading_rates
        # events_distance
        # stiffness_slice
        # stiffness_slice_corr

        self.results_type = None

        # Multi meshgrid
        self.multi_meshgrid_type = "piezo"
        self.missing_z_positions = []

        # Info selected
        self.info_selected = "spring"

        self.tilt_correction_all = False
        self.flatten_correction_all = False
        # Profiles
        self.show_profiles_bars = True
        self.profiles_bar_1 = 100
        self.profiles_bar_2 = 200
        # Color options
        self.all_stiffness_depth_view = 0
        self.meshgrid_click_xpos = 0
        self.meshgrid_click_ypos = 0
        self.roi_cursor_shape = "square"
        self.roi_cursor_size = 2
        self.meshgrid_display_red_square = True
        self.meshgrid_display_discarded = True
        self.meshgrid_display_missing_z = True
        self.display_line_stiffness_slice_at_height = True

        # Opengl
        self.opengl_background_color = [0.0, 0.0, 0.0, 0.0]
        self.opengl_screenshot_magn = 1

        # Type of selected computations box in compute tab
        # stiffness, stiffness_corr, work_and_rupture_force, events,
        # loading_rate
        self.compute_type = "stiffness"
        self.compute_mode = "approach"  # Approach or retraction

        # Global display values
        self.global_max_indentation_index = None
        self.global_max_piezo = None
        self.global_max_topo = None
        self.global_max_stiffness = None
        self.global_max_work = None
        self.global_max_rupture_force1 = None
        self.global_max_nbr_events = None
        self.global_max_rupture_force2 = None

        # DFS
        self.dfs_mode = "single"
        self.list_dfs_single_fits = []
        self.list_dfs_single_results = []
        self.list_dfs_groups_fits = []
        self.list_dfs_groups_results = []

        # Units
        self.infowidget_threshold_units = "nN"

        # Tilt or stretch is selected in curve modifcation widget
        self.display_curve_modif_mode = "tilt"
        self.display_stretch_mode = "approach"

        # Display or hide ROIs in multimeshgrids widget
        self.display_roi_in_multimeshgrid = False

        # Events
        # List of reference events for the convolutions
        self.list_events_ref_selected = 0
        self.list_events_ref = [{
            "unique_id": 0,
            "pos_x": 0,
            "pos_y": 0,
            "ref_start": 0,
            "ref_stop": 100}]

        # Available Themes
        self.theme_list = ["Light", "Dark"]
        # Light theme is the default
        self.current_theme = "Dark"

    @property
    def current_data(self):
        """Getter whichs return the currently selected file."""
        return self.list[self.id_selected]

    def addData(self, filename, unique_id):
        """Will create a SingleFile object and add it to the list of files."""
        string = "Add data to exp, name = %s, id = %s"
        self.logger.debug(string, filename, str(unique_id))

        self.list.append(SingleFile(self, filename, self.temp_file, unique_id))

    def update_global_values(self):
        """Updates values like maximum stiffness, maximum indentation ..."""
        self.logger.debug("Update global values")

        array_indentation = []
        array_piezo = []
        array_topo = []
        array_stiffness = []
        array_work = []
        array_rupture_force1 = []
        array_nbr_events = []
        array_rupture_force2 = []
        for i in range(len(self.list)):
            array_piezo.append(self.list[i].max_piezo)
            if self.list[i].max_topo is not None:
                array_topo.append(self.list[i].max_topo)
                array_stiffness.append(self.list[i].max_stiffness)
                array_indentation.append(self.list[i].max_indentation_index)
            if self.list[i].max_work is not None:
                array_work.append(self.list[i].max_work)
            if self.list[i].max_rupture_force1 is not None:
                array_rupture_force1.append(self.list[i].max_rupture_force1)
            if self.list[i].max_nbr_events is not None:
                array_nbr_events.append(self.list[i].max_nbr_events)
                array_rupture_force2.append(self.list[i].max_rupture_force2)

        self.global_max_piezo = numpy.amax(array_piezo)
        if array_topo != []:
            self.global_max_topo = numpy.amax(array_topo)
            self.global_max_stiffness = numpy.amax(array_stiffness)
            self.global_max_indentation_index = numpy.amax(array_indentation)
        if array_work != []:
            self.global_max_work = numpy.amax(array_work)
        if array_rupture_force1 != []:
            self.global_max_rupture_force1 = numpy.amax(array_rupture_force1)
        if array_nbr_events != []:
            self.global_max_nbr_events = numpy.amax(array_nbr_events)
            self.global_max_rupture_force2 = numpy.amax(array_rupture_force2)


class ResultGroup:
    """Object used to store settings for a result's line in the grouped results."""

    def __init__(self, group_id):
        super().__init__()

        self.group_id = group_id
        self.display = True
        if self.group_id == 0:
            self.name = "None"
        else:
            self.name = "Group " + str(group_id)
        self.color = None  # Color for the histograms

        self.condition = 0  # 0 = No condition


class ConditionGroup:
    """Object used to store settings for a result's line in the experiment results."""

    def __init__(self, condition_id):
        super().__init__()

        self.condition_id = condition_id
        self.display = True
        if self.condition_id == 0:
            self.name = "None"

        else:
            self.name = "Condition " + str(condition_id)
        self.color = None  # Color for the histograms


class SampleGroup:
    """Object used to store settings for a sample group line in the stats experiment results."""

    def __init__(self, sample_group_id):
        super().__init__()

        self.sample_group_id = sample_group_id
        self.display = True
        self.color = None  # Color for the histograms


def add_result(unique_id, name):
    """Appends a new Result object to the result's list.

    If the result is already in the list, do not add it (Can happen if called
    from the computations tab, if data is recalcultated).

    Defines also the color of the histogram
    """
    # Get the list of the old id's
    ids = []
    for result in shared.exp.results_list:
        ids.append(result.data_id)

    # If the result is not in the list, add it
    if not math_tools.in_list(ids, unique_id):
        result = Result(len(shared.exp.results_list), unique_id)
        shared.exp.results_list.append(result)
        result.name = name

        # Set the color. Go through the list and determine a color
        # periodically.
        pos = len(shared.exp.results_list) - 1
        while pos >= len(shared.colors_list):
            pos = pos - len(shared.colors_list)
        result.color = shared.colors_list[pos]


class Result:
    """Object used to store settings for a result's line in the results tab."""

    def __init__(self, result_id, data_id):
        super().__init__()

        self.result_id = result_id
        self.data_id = data_id  # The result id in exp.list (where the data is)
        self.display = False
        self.group = 0  # 0 = no group
        self.slice = 0  # 0 = All

        # Filtering of the events by distance from JOC, in nm
        self.filter_dist_left = 0
        self.filter_dist_right = 0
        self.filter_dist_keep_middle = False

        # Event filtering by slope of the fit on the right
        self.filter_slope_min = 0
        self.filter_slope_max = 0

        self.roi = 0  # 0 = no ROI
        self.dist_to_joc = 0
        self.sub_type = 0  # 0 = normal, 1 = corrected
        self.name = ""
        self.color = None  # Color for the histograms (defined by add_result)

        # Experimental condition
        self.condition = 0  # 0 = no condition


class ROI:
    """ROI class.

    Used to save all the informations about a ROI (color, size, pixels, ...).
    The ROIs are saved inside the roi_list for each dataset.

    Roi id 0 = No Roi
    Roi id 1 = First Roi
    etc..
    """

    def __init__(self, roi_id):
        self.roi_id = roi_id
        self.display = True
        self.current = False
        self.color = 0
        self.glass_coeffs = None
        self.values = []
