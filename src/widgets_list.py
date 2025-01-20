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
List of all the widgets in PYAF.

By default the widgets are set to None. Once a widget is opened the
corresponding variable contains a PYAFWidget (see PYAFWidget.py in tools
folder).

"""

# Main widget
widget_main = None

# Main widgets
widget_data = None
widget_compute = None
widget_results = None
widget_results_single = None
widget_results_groups = None
widget_results_experiment = None

# About widget
widget_about = None

# Sub widgets
widget_results_chooser_single = None
widget_results_hist_options_single = None
widget_results_chooser_groups = None
widget_results_hist_options_groups = None
widget_curve_selector_data = None
widget_curve_selector_compute = None
widget_curve_selector_results = None
widget_kernel = None
widget_msf = None

# Progressbar
widget_progressbar = None

# Other widgets
widget_preferences = None
widget_info = None
widget_flatten = None
widget_multimeshgrids = None
widget_events_per_scan = None
widget_meshgrid_options = None
widget_slices = None
widget_clean_up = None
widget_profiles = None
widget_details = None
widget_roi_manager = None
widget_rename_groups = None
widget_rename_conditions = None
widget_dfs = None
widget_indentation = None
widget_vtk = None
widget_curve_mod = None
widget_tilt_check = None
widget_plugins = None
widget_slicing = None
widget_fiducials = None
widget_vtk_entry = None
widget_used_parameters = None
widget_vtk_glass_choser = None
widget_vtk_binary_threshold = None
widget_advanced_options = None
widget_sample_grouping_dialog = None

widget_fit_param_pocs = None
widget_fit_param_jocs = None
widget_fit_param_jocs_events = None
