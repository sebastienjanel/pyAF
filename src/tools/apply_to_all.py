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

"""Used to apply preferences to all files in the experiment."""

import copy
from .. import shared
from .. import widgets_list


def apply_to_all(mode, option=None, oldid=None):
    """Used to apply preferences to all files in the experiment.

    Lets says you want to change the values of a colorscale for a meshgrid in
    PYAF. You would want to be able to change it for one single file, or for
    all the files. That's why there are apply to all checkboxes in the GUI.

    Technically, there are two possibilities to manage this :
    - Store a global value which is used when the checkbox is checked, and
      single values for each file when the checkbox unchecked.
      That means that each time you change a preferences, or you fetch a
      preference, you have to have an if/else statement checking if we want
      to use the global or single value.

    - Do not store a global value, but overwrite all the single values with
      the current value when the checkbox is checked. This is the solution
      which was chosen for PYAF. The apply_to_all function does this, and is
      called whenever needed in the software.
    """
    exp = shared.exp

    if oldid is not None:
        # Force the usage of a particular id
        oldid = oldid

    if option == "load_more":
        # Will take values from the first file. If apply to all is checked
        # this is ok, all the values are the same. If not, defining oldid
        # as 0 is not important as nothing will be done below.
        oldid = 0
    elif option != "load_more" and oldid is None:
        oldid = exp.id_selected

    if exp.apply_to_all_data and (mode == "display_options" or mode == "all"):
        for i in range(len(exp.list)):
            dt = exp.list[i]
            cur = exp.list[oldid]

            mt = cur.meshgrid_type
            st_calc = dt.stiffness_calculated
            dct = cur.display_curve_type
            wf_calc = dt.work_and_rupture_force1_calculated

            if (mt == "topo" or mt == "stiffness" or mt == "stiffness_slice") \
                    and st_calc:
                dt.meshgrid_type = mt
            if (mt == "stiffness_corr" or mt == "stiffness_corr_slice") and \
                    dt.stiffness_corrected:
                dt.meshgrid_type = mt
            elif (mt == "work" or mt == "rupture_force") and \
                    dt.work_and_rupture_force1_calculated:
                dt.meshgrid_type = mt
            elif (mt == "events_per_curve" or mt == "events_rupture_force") \
                    and dt.events_calculated:
                dt.meshgrid_type = mt
            elif mt == "piezo":
                dt.meshgrid_type = mt

            if (st_calc or dt.work_and_rupture_force1_calculated) \
                    and dct == "force":
                dt.display_curve_type = dct
            elif st_calc and dct == "indentation":
                dt.display_curve_type = dct
            elif dt.events_calculated and dct == "results_force_events":
                dt.display_curve_type = dct
            else:
                dt.display_curve_type = dct

            if cur.display_trace_retrace == 0 and \
                    (dct == "force" or dct == "indentation") and st_calc:
                dt.display_trace_retrace = cur.display_trace_retrace
            elif cur.display_trace_retrace == 1 and dct == "force" and wf_calc:
                dt.display_trace_retrace = cur.display_trace_retrace
            elif cur.display_trace_retrace == 2 and dct == "force" and \
                    st_calc and wf_calc:
                dt.display_trace_retrace = cur.display_trace_retrace
            elif dct == "defl_ext" or dct == "force_ext":
                dt.display_trace_retrace = cur.display_trace_retrace

            if st_calc:
                dt.display_poc = cur.display_poc
                dt.display_fit_poc = cur.display_fit_poc
                dt.display_segments = cur.display_segments
                dt.display_fits_stiffness = cur.display_fits_stiffness
                dt.stiffness_depth_view = cur.stiffness_depth_view
                dt.stiffness_slice_depth = cur.stiffness_slice_depth
            if dt.work_and_rupture_force1_calculated:
                dt.display_joc = cur.display_joc
                dt.display_fit_joc = cur.display_fit_joc
                dt.display_surface = cur.display_surface
                dt.display_force = cur.display_force
            if dt.events_calculated:
                dt.events_results_display_annotations = \
                    cur.events_results_display_annotations
                dt.events_results_filter_dist_left = \
                    cur.events_results_filter_dist_left
                dt.events_results_filter_dist_right = \
                    cur.events_results_filter_dist_right
                dt.events_results_filter_dist_keep_middle = \
                    cur.events_results_filter_dist_keep_middle
                dt.events_display_results_filter_dist = \
                    cur.events_display_results_filter_dist
                dt.events_results_filter_slope_min = \
                    cur.events_results_filter_slope_min
                dt.events_results_filter_slope_max = \
                    cur.events_results_filter_slope_max
                dt.events_results_display_joc = cur.events_results_display_joc
                dt.events_results_display_fit_joc = \
                    cur.events_results_display_fit_joc

            dt.curve_force_units = cur.curve_force_units
            dt.curve_distance_units = cur.curve_distance_units
            dt.meshgrid_units = cur.meshgrid_units

            # Colortable colors
            # (important to have deepcopy here because these are lists)

            if cur.is_single is False:
                # This one is only updated in case of a meshgrid. Else the max
                # value is 0 and all the meshgrids will have a colorscale
                # with a maximum set to 0, resulting in an error
                dt.color_opts_piezo = copy.deepcopy(cur.color_opts_piezo)

            dt.color_opts_topo = copy.deepcopy(cur.color_opts_topo)
            dt.color_opts_stiffness = copy.deepcopy(cur.color_opts_stiffness)
            dt.color_opts_work = copy.deepcopy(cur.color_opts_work)
            dt.color_opts_rupture_force1 = \
                copy.deepcopy(cur.color_opts_rupture_force1)
            dt.color_opts_events_per_curve = \
                copy.deepcopy(cur.color_opts_events_per_curve)
            dt.color_saturation = cur.color_saturation
            dt.color_negative = cur.color_negative
            dt.color_nan = cur.color_nan

            dt.nan_color_transparent = cur.nan_color_transparent

    if exp.apply_to_all_compute and (mode == "calc_options" or mode == "all"):
        for i in range(len(exp.list)):
            dt = exp.list[i]
            cur = exp.list[oldid]

            dt.sg_smoothing_enabled = cur.sg_smoothing_enabled
            dt.sg_smoothing_order = cur.sg_smoothing_order
            dt.sg_smoothing_width = cur.sg_smoothing_width
            dt.sg_smoothing_uniform = cur.sg_smoothing_uniform

            dt.check_tilt_option = cur.check_tilt_option

        if exp.compute_type == "stiffness" or mode == "all":
            for i in range(len(exp.list)):
                dt = exp.list[i]
                cur = exp.list[oldid]
                st_calc = dt.stiffness_calculated

                dt.stiffness_model_selected = cur.stiffness_model_selected
                dt.tomography = cur.tomography
                dt.tip_radius = cur.tip_radius
                dt.tip_angle = cur.tip_angle
                dt.poisson_ratio = cur.poisson_ratio
                dt.perform_fit = cur.perform_fit
                dt.fit_range_type = cur.fit_range_type
                dt.indentation_start = cur.indentation_start
                if dt.fit_range_type == 0:
                    dt.force_start = cur.force_start
                    dt.force_stop = cur.force_stop
                elif dt.fit_range_type == 1:
                    dt.indentation_stop = cur.indentation_stop
                    dt.indentation_step = cur.indentation_step
                dt.strict_stop = cur.strict_stop
                dt.fitparam_poc_skip_start = cur.fitparam_poc_skip_start
                dt.fitparam_poc_fit_length = cur.fitparam_poc_fit_length
                dt.fitparam_poc_refit_option = cur.fitparam_poc_refit_option
                dt.fitparam_poc_noise_multiplicator = \
                    cur.fitparam_poc_noise_multiplicator
                dt.fitparam_poc_refit_times = cur.fitparam_poc_refit_times
                dt.stretch_app_lim1 = cur.stretch_app_lim1
                dt.stretch_app_lim2 = cur.stretch_app_lim2
                dt.stretch_ret_lim1 = cur.stretch_ret_lim1
                dt.stretch_ret_lim2 = cur.stretch_ret_lim2
                dt.stretch_len_app = cur.stretch_len_app
                dt.stretch_len_ret = cur.stretch_len_ret
                dt.stretch_applied_app = cur.stretch_applied_app
                dt.stretch_applied_ret = cur.stretch_applied_ret

        if exp.compute_type == "stiffness_corr" or mode == "all":
            for i in range(len(exp.list)):
                dt = exp.list[i]
                cur = exp.list[oldid]
                dt.user_h = cur.user_h

        if exp.compute_type == "work_and_rupture_force" or mode == "all":
            for i in range(len(exp.list)):
                dt = exp.list[i]
                cur = exp.list[oldid]

                dt.fitparam_joc_skip_start = cur.fitparam_joc_skip_start
                dt.fitparam_joc_fit_length = cur.fitparam_joc_fit_length
                dt.fitparam_joc_refit_option = cur.fitparam_joc_refit_option
                dt.fitparam_joc_noise_multiplicator = \
                    cur.fitparam_joc_noise_multiplicator
                dt.fitparam_joc_refit_times = cur.fitparam_joc_refit_times

        if exp.compute_type == "events" or mode == "all":
            for i in range(len(exp.list)):
                dt = exp.list[i]
                cur = exp.list[oldid]

                dt.events_algorithm = cur.events_algorithm
                dt.events_curve_type = cur.events_curve_type
                dt.fitparam_events_joc_skip_start = \
                    cur.fitparam_events_joc_skip_start
                dt.fitparam_events_joc_fit_length = \
                    cur.fitparam_events_joc_fit_length
                dt.fitparam_events_joc_refit_option = \
                    cur.fitparam_events_joc_refit_option
                dt.fitparam_events_joc_noise_multiplicator = \
                    cur.fitparam_events_joc_noise_multiplicator
                dt.fitparam_events_joc_refit_times = \
                    cur.fitparam_events_joc_refit_times
                dt.display_fit_events_joc2_preview = \
                    cur.display_fit_events_joc2_preview

                dt.events_kernel_detection_threshold = \
                    cur.events_kernel_detection_threshold
                dt.events_msf_detection_threshold = \
                    cur.events_msf_detection_threshold
                dt.events_fit_size = cur.events_fit_size
                dt.events_msf_window_size = cur.events_msf_window_size
                dt.kernel_size = cur.kernel_size
                dt.adaptive_threshold_option = \
                    cur.adaptive_threshold_option
                dt.events_kernel_adaptive_threshold = \
                    cur.events_kernel_adaptive_threshold
                dt.adaptive_smoothing_window = \
                    cur.adaptive_smoothing_window
                dt.adaptive_smoothing_order = cur.adaptive_smoothing_order

                dt.display_fit_event = cur.display_fit_event
                dt.display_fit_event_seg = cur.display_fit_event_seg

        if exp.compute_type == "loading_rate" or mode == "all":
            for i in range(len(exp.list)):
                dt = exp.list[i]
                cur = exp.list[oldid]

                if dt.is_single is False:
                    dt.lr_coef = cur.lr_coef

    if exp.tilt_correction_all and (mode == "tilt_all" or mode == "all"):
        for i in range(len(exp.list)):
            dt = exp.list[i]
            cur = exp.list[oldid]
            dt.tilt_limit_1 = cur.tilt_limit_1
            dt.tilt_limit_2 = cur.tilt_limit_2

    # Update some widgets
    if widgets_list.widget_indentation is not None:
        widgets_list.widget_indentation.update_widget()

    if widgets_list.widget_meshgrid_options is not None:
        widgets_list.widget_meshgrid_options.CB_apply_to_all.\
            setChecked(exp.apply_to_all_data)

    return True
