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

"""Computation of a part of a map, run in the worker processes.

This module is imported by each worker process of the pool (see
multiproc.py): it only imports what is needed for the computations, not
the GUI.
"""

from .stiffness.stiffness_compute import main_compute_stiff
from .stiffness.stiffness_fit import main_fit_stiffness
from .work_and_rupture_force.work_and_rupture_force_compute \
    import main_compute_work


# Parameters copied for each curve in the dictionary given to the
# computation functions
CURVE_PARAMS = (
    "tilt_limit_1", "tilt_limit_2", "deflection_sensitivity",
    "spring_constant", "tilt_applied", "sg_smoothing_enabled",
    "sg_smoothing_order", "sg_smoothing_width", "sg_smoothing_uniform",
    "stretch_applied_app", "stretch_applied_ret", "stretch_app_lim1",
    "stretch_app_lim2", "stretch_len_app", "stretch_ret_lim1",
    "stretch_ret_lim2", "stretch_len_ret")

STIFFNESS_FIT_PARAMS = {
    "poc_skip_start": "poc_skip_start",
    "poc_fit_length": "poc_fit_length",
    "poc_refit_option": "poc_refit_option",
    "poc_noise_multiplicator": "poc_noise_multiplicator",
    "poc_refit_times": "poc_refit_times",
    "trig_threshold": "trig_threshold",
    "poisson_ratio": "poisson_ratio",
    "tip_radius": "tip_radius",
    "tip_angle": "tip_angle",
    "fit_range_type": "fit_range_type"}

STIFFNESS_OPTIONS = {
    "indentation_start": "indentation_start",
    "indentation_stop": "indentation_stop",
    "indentation_step": "indentation_step",
    "force_start": "force_start",
    "force_stop": "force_stop",
    "model_selected": "model_selected",
    "tomography": "tomography",
    "coeff": "coeff"}

WORK_FIT_PARAMS = {
    "joc_skip_start": "joc_skip_start",
    "joc_fit_length": "joc_fit_length",
    "joc_refit_option": "joc_refit_option",
    "joc_noise_multiplicator": "joc_noise_multiplicator",
    "joc_refit_times": "joc_refit_times"}


def compute_columns(task):
    """Computes the curves of some columns of a map.

    task is a dictionary with:
    - mode: "stiffness" or "work"
    - params: the parameters of the computation (see the stiffness and
      work main.py files), without the arrays of the whole map
    - first_col, first_row: position of the part in the map
    - app, ret: approach and retraction curves of the part, indexed as
      [col][row] (None if not needed)
    - app_positions, ret_positions: valid parts of the curves
    - discarded: discarded curves, indexed as [col][row]
    - piezo_image: piezo image of the part (None if not needed)

    Returns the list of (status, position, result) for each curve, column
    after column. status is "continue", "fit_error" or "smoothing_error",
    position is the position of the curve for the users ("[col, row]",
    from 1), and result is the result of the computation function of the
    curve.
    """
    mode = task["mode"]
    params = task["params"]
    first_col = task["first_col"]
    first_row = task["first_row"]
    app = task["app"]
    ret = task["ret"]
    nbr_cols = len(task["discarded"])
    nbr_rows = len(task["discarded"][0])

    dt = {name: params[name] for name in CURVE_PARAMS}
    if mode == "stiffness":
        fit_params = {key: params[name]
                      for key, name in STIFFNESS_FIT_PARAMS.items()}
        options = {key: params[name]
                   for key, name in STIFFNESS_OPTIONS.items()}
    else:
        fit_params = {key: params[name]
                      for key, name in WORK_FIT_PARAMS.items()}

    results = []
    for i in range(nbr_cols):
        col = first_col + i
        for j in range(nbr_rows):
            row = first_row + j
            dt["is_discarded"] = task["discarded"][i][j]
            dt["approach_position"] = None
            dt["retraction_position"] = None
            dt["piezo_image"] = None

            curve_app = None
            if app is not None:
                curve_app = app[i][j]
                dt["approach_position"] = task["app_positions"][i][j]
                if task["piezo_image"] is not None:
                    dt["piezo_image"] = task["piezo_image"][i][j]

            curve_ret = None
            if ret is not None:
                curve_ret = ret[i][j]
                dt["retraction_position"] = task["ret_positions"][i][j]

            if mode == "stiffness":
                if params["perform_fit"]:
                    result = main_fit_stiffness(col, row, curve_app,
                                                curve_ret, dt, fit_params,
                                                options)
                else:
                    result = main_compute_stiff(col, row, curve_app,
                                                curve_ret, dt, fit_params,
                                                options)
            else:
                result = main_compute_work(col, row, curve_app, curve_ret,
                                           dt, fit_params)

            # Errors of the curve (the computation continues)
            position = "[" + str(col + 1) + ", " + str(row + 1) + "]"
            status = "continue"
            if result is not None and mode == "stiffness":
                if params["perform_fit"] and len(result[5]) == 0:
                    # The fit failed, no stiffness value for this curve
                    status = "fit_error"
                elif result[7]:
                    status = "smoothing_error"

            results.append((status, position, result))

    return results
