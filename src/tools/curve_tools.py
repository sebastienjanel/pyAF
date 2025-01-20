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

"""Some important tools for the force curves."""

# pylint: disable=E1103, W0702

import numpy
from ..tools import math_tools


def get_curve(data, pos, approach=None, retraction=None, pause=None, modulation=None,
              mode=None, zero_defl=False, dt=None, get_time=False):
    """Return a curve.

    The curve is corrected for the tilt and/or flattened if needed.

    The full curve is returned (it is not truncated with the approach and
    retraction positions). Those are only needed for the smoothing and tilt
    correction."""
    # Should be saved as used in the future ? :
    # data.discarded_curves ?
    i, j = pos[0], pos[1]

    valid_modes = ["compute", "fit_preview", "tilt_curve_preview",
                   "curve_results", "curve_data", "curve_dataa", "curve_datab", "next_event", "save_event",
                   "details_widget", "results_sorting", "compute_no_shared"]

    if mode not in valid_modes:
        raise Exception("The mode in get curve is not valid : " + str(mode))

    # The curve can be set to 0 deflection if asked. Then the value of this
    # correction is returned. This is only used for the plots, in the
    # computations there is no need to correct the deflection to 0.
    zero = 0

    stretch_applied_app = False
    stretch_applied_ret = False

    curve_x_app = None
    curve_x_time_app = None
    curve_y_app = None

    curve_x_ret = None
    curve_x_time_ret = None
    curve_y_ret = None

    curve_x_pause = None
    curve_x_time_pause = None
    curve_y_pause = None

    curve_x_modulation = None
    curve_x_time_modulation = None
    curve_y_modulation = None

    if mode == "compute":
        # Used during the computations
        if approach is not None:
            curve_x_app = approach[0]
            curve_y_app = approach[1]
            if get_time:
                curve_x_time_app = approach[2]
        if retraction is not None:
            curve_x_ret = retraction[0]
            curve_y_ret = retraction[1]
            if get_time:
                curve_x_time_ret = retraction[2]
        if pause is not None:
            curve_x_pause = pause[0]
            curve_y_pause = pause[1]
            if get_time:
                curve_x_time_pause = pause[2]
        if modulation is not None:
            curve_x_modulation = modulation[0]
            curve_y_modulation = modulation[1]
            if get_time:
                curve_x_time_modulation = modulation[2]
        deflection_sensitivity = dt["deflection_sensitivity"]

    else:

        get_time_segments_flag = data.file_type in ("JPK (Single File)", "JPK (Force Map)", "JPK (QI)")

        # print(time_segments)

        # Used for the plots, or the details widget, results_sorting
        curve_x_app = data.curves_approach[i][j][0]
        curve_y_app = data.curves_approach[i][j][1]
        curve_x_ret = data.curves_retraction[i][j][0]
        curve_y_ret = data.curves_retraction[i][j][1]

        if get_time_segments_flag:
            curve_x_time_app = data.curves_approach[i][j][2]
            curve_x_time_ret = data.curves_retraction[i][j][2]

        if data.nbr_points_per_pause_curve is not None:
            curve_x_pause = data.curves_pause[i][j][0]
            curve_y_pause = data.curves_pause[i][j][1]

            if get_time_segments_flag:
                curve_x_time_pause = data.curves_pause[i][j][2]

        if data.nbr_points_per_modulation_curve is not None:
            curve_x_modulation = data.curves_modulation[i][j][0]
            curve_y_modulation = data.curves_modulation[i][j][1]

            if get_time_segments_flag:
                curve_x_time_modulation = data.curves_modulation[i][j][2]

        deflection_sensitivity = data.deflection_sensitivity

    if dt is None:
        is_discarded = data.discarded_curves[i][j]
    else:
        is_discarded = dt["is_discarded"]

    if (mode == "fit_preview" or mode == "tilt_curve_preview" or
            mode == "save_event" or mode == "next_event"):
        tilt_limit_1 = data.tilt_limit_1
        tilt_limit_2 = data.tilt_limit_2
        tilt_applied = data.tilt_applied
        sg_smoothing_enabled = data.sg_smoothing_enabled
        sg_smoothing_order = data.sg_smoothing_order
        sg_smoothing_width = data.sg_smoothing_width
        sg_smoothing_uniform = data.sg_smoothing_uniform

        if not is_discarded:
            stretch_applied_app = data.stretch_applied_app
            stretch_applied_ret = data.stretch_applied_ret

            if stretch_applied_app:
                stretch_limit1 = data.stretch_app_lim1
                stretch_limit2 = data.stretch_app_lim2
                stretch_length = data.stretch_len_app
            elif stretch_applied_ret:
                stretch_limit1 = data.stretch_ret_lim1
                stretch_limit2 = data.stretch_ret_lim2
                stretch_length = data.stretch_len_ret

    elif mode == "curve_results" or mode == "details_widget" or \
            mode == "results_sorting" or mode == "compute_no_shared":
        deflection_sensitivity = data.used_deflection_sensitivity
        if deflection_sensitivity is None:
            # If nothing has been computed, there is no used_defl_sens
            # but we still want to display the curve : use the normal one
            deflection_sensitivity = data.deflection_sensitivity

        tilt_limit_1 = data.used_tilt_limit_1
        tilt_limit_2 = data.used_tilt_limit_2
        tilt_applied = data.used_tilt_applied
        sg_smoothing_enabled = data.used_sg_smoothing_enabled
        sg_smoothing_order = data.used_sg_smoothing_order
        sg_smoothing_width = data.used_sg_smoothing_width
        sg_smoothing_uniform = data.used_sg_smoothing_uniform

        if is_discarded == 0:
            stretch_applied_app = data.used_stretch_applied_app
            stretch_applied_ret = data.used_stretch_applied_ret

            if stretch_applied_app:
                stretch_limit1 = data.used_stretch_app_lim1
                stretch_limit2 = data.used_stretch_app_lim2
                stretch_length = data.used_stretch_len_app
            elif stretch_applied_ret:
                stretch_limit1 = data.used_stretch_ret_lim1
                stretch_limit2 = data.used_stretch_ret_lim2
                stretch_length = data.used_stretch_len_ret

    elif mode == "compute":
        tilt_limit_1 = dt["tilt_limit_1"]
        tilt_limit_2 = dt["tilt_limit_2"]
        tilt_applied = dt["tilt_applied"]
        stretch_applied_app = dt["stretch_applied_app"]
        stretch_applied_ret = dt["stretch_applied_ret"]

        if stretch_applied_app:
            stretch_limit1 = dt["stretch_app_lim1"]
            stretch_limit2 = dt["stretch_app_lim2"]
            stretch_length = dt["stretch_len_app"]
        elif stretch_applied_ret:
            stretch_limit1 = dt["stretch_ret_lim1"]
            stretch_limit2 = dt["stretch_ret_lim2"]
            stretch_length = dt["stretch_len_ret"]

        sg_smoothing_enabled = dt["sg_smoothing_enabled"]
        sg_smoothing_order = dt["sg_smoothing_order"]
        sg_smoothing_width = dt["sg_smoothing_width"]
        sg_smoothing_uniform = dt["sg_smoothing_uniform"]

    elif mode == "curve_data" or mode == "curve_dataa" or mode == "curve_datab":
        tilt_applied = False
        sg_smoothing_enabled = False

    if mode == "compute":
        if curve_x_app is not None:
            ind_offset_app = dt["approach_position"]
        if curve_x_ret is not None:
            ind_offset_ret = dt["retraction_position"]
    else:
        ind_offset_app = data.approach_positions[i][j]
        ind_offset_ret = data.retraction_positions[i][j]

    app = None
    ret = None
    pau = None
    mod = None

    if curve_x_app is not None:
        app = [curve_x_app, curve_y_app * deflection_sensitivity, curve_x_time_app]
    if curve_x_ret is not None:
        ret = [curve_x_ret, curve_y_ret * deflection_sensitivity, curve_x_time_ret]
    if curve_x_pause is not None:
        pau = [curve_x_pause, curve_y_pause * deflection_sensitivity, curve_x_time_pause]
    if curve_x_modulation is not None:
        mod = [curve_x_modulation, curve_y_modulation * deflection_sensitivity, curve_x_time_modulation]

    # Smooth data if needed
    smoothing_error = False

    if sg_smoothing_enabled and is_discarded == 0:
        # Do not take into account the corrupted data from the curves
        # (only the beggining of the curve)

        try:
            if app is not None:
                seg = slice(ind_offset_app[0], len(app[0]))
                seg0 = slice(0, ind_offset_app[0])

                curve = math_tools.sg_filter(app[0][seg], app[1][seg],
                                             sg_smoothing_order,
                                             sg_smoothing_width,
                                             sg_smoothing_uniform)
                app[1] = numpy.concatenate([numpy.array(app[1][seg0]), curve])

            if ret is not None:
                seg = slice(ind_offset_ret[0], len(ret[0]))
                seg0 = slice(0, ind_offset_ret[0])

                ret[1] = math_tools.sg_filter(
                    ret[0], ret[1], sg_smoothing_order, sg_smoothing_width,
                    sg_smoothing_uniform)

        except:
            # Sometimes you will get trouble with corrupted curves
            # return an error
            smoothing_error = True

    # Find a 0 for the deflection (first value of the curve). Only if asked.
    if zero_defl and (
            mode == "curve_results" or mode == "curve_data" or mode == "curve_dataa" or mode == "curve_datab"):
        zero = app[1][int(data.approach_positions[i][j][0])]
        app[1] = app[1] - zero
        ret[1] = ret[1] - zero

        if pau is not None:
            pau[1] = pau[1] - zero

        if mod is not None:
            mod[1] = mod[1] - zero

        if tilt_applied:
            # If a tilt is applied, do not use the zero correction.
            # In the computation, no zero correction is used, and the fit for
            # the pocs and jocs is determined on the tilt corrected curve
            # (with F around 0 due to the correction)
            # Then, there is no meaning to correct for zero.
            zero = 0.0

    # Correct tilt if needed
    if (is_discarded == 0 and
            (tilt_applied == "trace" or tilt_applied == "retrace")):
        approach, retraction, _ = math_tools.correct_tilt(
            app,
            ret,
            tilt_limit_1,
            tilt_limit_2,
            tilt_applied,
            ind_offset_app,
            ind_offset_ret)  # Apply for modulation and pause curves.
    else:
        approach, retraction, pause, modulation = app, ret, pau, mod

    # Add stretch
    if stretch_applied_app:
        approach = apply_stretch(approach, stretch_limit1, stretch_limit2,
                                 ind_offset_app[0], stretch_length)

    if stretch_applied_ret:
        retraction = apply_stretch(retraction, stretch_limit1, stretch_limit2,
                                   ind_offset_ret[0], stretch_length)

    return approach, retraction, pause, modulation, zero, smoothing_error


def get_force_curves(curve_x, curve_y, center, spring_constant):
    """Returns a force curve from a deflection-extension curve.

    Center is a real position (x, y) given on the curve, which needs to
    be substracted to set the 0, 0 point (piezo = 0, force = 0).
    """
    # Set the center position to 0, 0 and get a force curve
    center_force_x = center[0] - center[1]
    center_force_y = center[1] * spring_constant

    return [numpy.array(curve_x - curve_y - center_force_x),
            numpy.array(curve_y * spring_constant - center_force_y)]


def apply_stretch(curve, limit1, limit2, ind_offset, length):
    """Add some values at the end of the curve.

    Can be helpfull for fits where the flat part of the curve is not
    long enough.
    """
    # The limits are garanteed to be > ind_offset (approach/retraction
    # position). Because the apply stretch button checks for this before
    # allowing the stretch to be applied.

    # Find the indice
    limit1_ind = 0
    limit2_ind = 0
    found1 = False

    for i in range(ind_offset, len(curve[0]), 1):
        if curve[0][i] > limit1 and found1 is False:
            limit1_ind = int(i)
            found1 = True
        if curve[0][i] > limit2:
            limit2_ind = int(i)
            break

    seg = slice(limit1_ind, limit2_ind)

    coeffs, _ = math_tools.fit_linear(curve[0][seg], curve[1][seg])

    # Get size :
    xsize = abs(curve[0][seg][0] - curve[0][seg][-1]) / len(curve[0][seg])

    # Replace
    for i in range(limit1_ind, -1, -1):
        curve[0][i] = xsize * i
        curve[1][i] = numpy.polyval([coeffs[0], coeffs[1]], xsize * i)

    # Add in front
    indtoadd = length / xsize - limit1_ind

    new_curve_x = []
    new_curve_y = []
    for i in range(int(indtoadd)):
        new_curve_x.append(-xsize * i)
        new_curve_y.append(numpy.polyval([coeffs[0], coeffs[1]], -xsize * i))

    curve_x = numpy.concatenate((new_curve_x, curve[0]), axis=0)
    curve_y = numpy.concatenate((new_curve_y, curve[1]), axis=0)

    return [curve_x, curve_y]
