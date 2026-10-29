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

"""Parameters of the exports (JSON): the version of pyAF, the exported
result and its options, and for each file of the export its calibration,
its acquisition and the parameters used for its computations (the used_*
values of the dataset, not the current settings of the interface).

Every data export writes them next to its file (name.params.json, see
write_sidecar), and in the description of the TIFF images. The "settings"
of each file are the raw values of pyAF's settings (see SETTINGS), reloaded
by Load parameters (Compute tab, see tools.params_import). The units are
in the names of the keys; null: not computed or not known. The file names
are given without their folder (no user name in the paths).
"""

import datetime
import json
import math
import os

import numpy

from .. import consts
from .. import shared
from ..compute_tools.stiffness import contact_models

# Version of the format of the parameters
SCHEMA = 1


def jsonable(value):
    """Value for JSON: numpy numbers and arrays as Python ones, NaN and
    infinities as null."""
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, numpy.ndarray)):
        return [jsonable(v) for v in value]
    if isinstance(value, (numpy.bool_, bool)):
        return bool(value)
    if isinstance(value, (numpy.integer, int)):
        return int(value)
    if isinstance(value, (numpy.floating, float)):
        value = float(value)
        return value if math.isfinite(value) else None
    if value is None or isinstance(value, str):
        return value
    return str(value)


def round_value(value, digits=6):
    """Float with 6 significant digits (no noise of the conversions)."""
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return value
    if not math.isfinite(value):
        return None
    return float("%.*g" % (digits, value))


def calibration(data):
    """Calibration used for the computations (else the current one)."""
    spring = data.used_spring_constant
    if spring is None:
        spring = data.spring_constant
    sensitivity = data.used_deflection_sensitivity
    if sensitivity is None:
        sensitivity = data.deflection_sensitivity
    return {"spring_constant_N_m": round_value(spring * 1e9)
            if spring is not None else None,
            "deflection_sensitivity_nm_V": round_value(sensitivity),
            "temperature_C": round_value(data.temperature - 273.15)}


def acquisition(data):
    """Acquisition of the file."""
    def velocity(value):
        return round_value(value / 1000.0) if value else None

    return {"file_type": data.file_type,
            "microscope": data.microscope_name,
            "date": data.date,
            "scan_size_um": [round_value((data.scan_size_x or 0) / 1000.0),
                             round_value((data.scan_size_y or 0) / 1000.0)],
            "pixels": [data.nbr_pixels_x, data.nbr_pixels_y],
            "ramp_size_nm": round_value(data.ramp_size),
            "approach_velocity_um_s": velocity(data.approach_velocity),
            "retraction_velocity_um_s": velocity(data.retraction_velocity),
            "z_from_height_sensor": data.zsensor_used}


def curve_processing(data):
    """Smoothing and corrections of the curves used for the
    computations."""
    smoothing = None
    if data.used_sg_smoothing_enabled:
        smoothing = {"width": data.used_sg_smoothing_width,
                     "order": data.used_sg_smoothing_order,
                     "uniform": data.used_sg_smoothing_uniform}
    tilt = None
    if data.used_tilt_applied:
        tilt = {"limit_1_nm": data.used_tilt_limit_1,
                "limit_2_nm": data.used_tilt_limit_2}
    return {"savitzky_golay_smoothing": smoothing, "tilt_correction": tilt}


def fit_of_contact(prefix, data):
    """Parameters of the detection of a point of contact (poc) or of a jump
    of contact (joc, events_joc) by fitting."""
    def get(name):
        return getattr(data, "used_fitparam_" + prefix + "_" + name, None)

    return {"skip_start_nm": get("skip_start"),
            "fit_length_nm": get("fit_length"),
            "noise_multiplicator": get("noise_multiplicator"),
            "refit_nm": get("refit_option"),
            "refit_times": get("refit_times")}


def elasticity(data):
    """Parameters of the elasticity, None if not computed."""
    if not data.stiffness_calculated:
        return None
    model = data.used_stiffness_model_selected
    name = contact_models.MODEL_NAMES[model] \
        if model is not None and 0 <= model < len(
            contact_models.MODEL_NAMES) else None
    start, step, tomography = data.get_used_indentation()
    if data.used_fit_range_type == 1:
        fit_range = {"type": "indentation",
                     "start_nm": data.used_indentation_start,
                     "stop_nm": data.used_indentation_stop}
    else:
        fit_range = {"type": "force", "start_pct": data.used_force_start,
                     "stop_pct": data.used_force_stop}
    return {
        "model": name,
        "tip_radius_nm": data.used_tip_radius
        if model in contact_models.USES_RADIUS and
        model != contact_models.THIN_SHELL else None,
        # Thin shell: radius of the vesicle
        "vesicle_radius_nm": data.used_tip_radius
        if model == contact_models.THIN_SHELL else None,
        "tip_half_angle_deg": data.used_tip_angle
        if model in contact_models.USES_ANGLE else None,
        "shell_thickness_nm": data.get_used_shell_thickness()
        if model in contact_models.USES_THICKNESS else None,
        "poisson_ratio": data.used_poisson_ratio,
        "method": "Fitting" if data.used_perform_fit else "Tomography",
        "fit_range": fit_range,
        "segments_nm": {"start": start, "step": step,
                        "stop": data.used_indentation_stop}
        if tomography else None,
        "point_of_contact_window_pct": data.used_fit_poc_window,
        "min_indentation_pct": data.used_fit_min_indentation,
        "point_of_contact": fit_of_contact("poc", data),
        "bottom_effect_correction": bool(data.stiffness_corrected),
        "brush_model": {"body_fit_above_pct_of_max_force":
                        data.used_brush_threshold_pct}
        if data.brush_calculated else None}


def viscoelasticity(data):
    """Parameters of the viscoelasticity, None if not computed."""
    if not data.visco_calculated:
        return None
    return {"solver": "numerical" if data.used_visco_bec or
            data.used_visco_solver == 1 else "analytical",
            "bottom_effect_correction": bool(data.used_visco_bec),
            "viscous_drag_corrected": bool(data.used_visco_drag_correction),
            "baseline_before_contact_nm": data.used_visco_contact_offset,
            "points": data.used_visco_nbr_points,
            "initial_fluidity": "auto" if data.used_visco_beta_auto
            else data.used_visco_beta_init}


def microrheology(data):
    """Parameters of the microrheology, None if not computed."""
    if not data.microrheo_calculated:
        return None
    return {"method": "FFT" if data.used_microrheo_method == "fft"
            else "sine fit",
            "viscous_drag_b0_uNs_per_m": data.used_microrheo_drag,
            "piezo_correction": data.used_microrheo_piezo,
            "piezo_amplitude_corrected": bool(
                data.used_microrheo_piezo_amplitude),
            "max_frequency_Hz": data.used_microrheo_max_freq or None,
            "frequencies_Hz": sorted(data.modulation_frequencies or []),
            "modulation_amplitudes_nm": data.modulation_amplitudes}


def work_and_rupture_force(data):
    """Parameters of the detachment work and force, None if not
    computed."""
    if not data.work_and_rupture_force1_calculated:
        return None
    return {"jump_of_contact": fit_of_contact("joc", data)}


def events(data):
    """Parameters of the detection of the events, None if not computed."""
    if not data.events_calculated:
        return None
    algorithm = data.used_events_algorithm
    detection = {"algorithm": algorithm,
                 "event_fit_size": data.used_events_fit_size,
                 "adaptive_threshold": data.used_adaptive_threshold_option}
    if algorithm == "kernel":
        detection.update({
            "threshold": data.used_events_kernel_detection_threshold,
            "kernel": data.used_kernel,
            "kernel_size": data.used_kernel_size})
    elif algorithm == "msf":
        detection.update({
            "threshold": data.used_events_msf_detection_threshold,
            "window_size": data.used_events_msf_window_size})
    if data.used_adaptive_threshold_option:
        detection.update({
            "adaptive_smoothing_window": data.used_adaptive_smoothing_window,
            "adaptive_smoothing_order": data.used_adaptive_smoothing_order})
    return {"detection": detection,
            "jump_of_contact": fit_of_contact("events_joc", data),
            "loading_rate_coefficient": data.used_lr_coef
            if data.loading_rates_calculated else None}


def breakthrough(data):
    """Parameters of the breakthroughs, None if not computed. The loading
    rate is the nominal one (spring constant x approach velocity)."""
    if not data.breakthrough_calculated:
        return None
    # k [N/m] x v [nm/s] = nN/s
    velocity = data.approach_velocity
    rate = data.spring_constant * 1e9 * velocity if velocity else None
    return {"n_sigma": data.used_breakthrough_n_sigma,
            "min_thickness_nm": data.used_breakthrough_min_thickness,
            "min_force_nN": data.used_breakthrough_min_force,
            "max_events_per_curve": data.used_breakthrough_max_events,
            "nominal_loading_rate_nN_per_s": rate,
            "thickness": "jump of the separation at the rupture (thickness "
                         "of the compressed bilayer)"}


# Settings of the computations, by part (names of the attributes of the
# datasets, whose used_* values are the ones used for the computations):
# reloaded by tools.params_import (Load parameters, Compute tab)
SETTINGS = {
    "elasticity": (
        "stiffness_model_selected", "stiffness_mode", "fit_range_type",
        "brush_model", "brush_threshold_pct",
        "tip_radius", "tip_angle", "shell_thickness", "poisson_ratio",
        "indentation_start",
        "indentation_stop", "indentation_step", "force_start", "force_stop",
        "fit_min_indentation", "fit_poc_window", "tomography",
        "fitparam_poc_skip_start", "fitparam_poc_fit_length",
        "fitparam_poc_refit_option", "fitparam_poc_noise_multiplicator",
        "fitparam_poc_refit_times"),
    "viscoelasticity": (
        "visco_solver", "visco_bec", "visco_contact_offset",
        "visco_nbr_points", "visco_beta_auto", "visco_beta_init",
        "visco_drag_correction"),
    "microrheology": ("microrheo_method", "microrheo_drag",
                      "microrheo_piezo", "microrheo_max_freq",
                      "microrheo_piezo_amplitude"),
    "work_and_rupture_force": (
        "fitparam_joc_skip_start", "fitparam_joc_fit_length",
        "fitparam_joc_refit_option", "fitparam_joc_noise_multiplicator",
        "fitparam_joc_refit_times"),
    "events": (
        "fitparam_events_joc_skip_start", "fitparam_events_joc_fit_length",
        "fitparam_events_joc_refit_option",
        "fitparam_events_joc_noise_multiplicator",
        "fitparam_events_joc_refit_times", "events_algorithm",
        "events_kernel_detection_threshold",
        "events_kernel_adaptive_threshold", "events_msf_detection_threshold",
        "events_msf_window_size", "events_fit_size", "kernel_size",
        "adaptive_threshold_option", "adaptive_smoothing_window",
        "adaptive_smoothing_order"),
    "loading_rates": ("lr_coef",),
    "breakthrough": (
        "breakthrough_n_sigma", "breakthrough_min_thickness",
        "breakthrough_min_force", "breakthrough_max_events"),
    "curves": (
        "sg_smoothing_enabled", "sg_smoothing_order", "sg_smoothing_width",
        "sg_smoothing_uniform", "tilt_applied", "tilt_limit_1",
        "tilt_limit_2", "stretch_app_lim1", "stretch_app_lim2",
        "stretch_len_app", "stretch_ret_lim1", "stretch_ret_lim2",
        "stretch_len_ret"),
    "calibration": ("spring_constant", "deflection_sensitivity"),
}

# Computations of the parts (their settings are given if computed)
COMPUTED = {"elasticity": "stiffness_calculated",
            "viscoelasticity": "visco_calculated",
            "microrheology": "microrheo_calculated",
            "work_and_rupture_force": "work_and_rupture_force1_calculated",
            "events": "events_calculated",
            "loading_rates": "loading_rates_calculated",
            "breakthrough": "breakthrough_calculated"}


def settings(data):
    """Settings used for the computations of a file, by part (None for the
    parts not computed); spring constant in N/m."""
    result = {}
    for part, names in SETTINGS.items():
        flag = COMPUTED.get(part)
        if flag is not None and not getattr(data, flag, False):
            result[part] = None
            continue
        values = {}
        for name in names:
            value = getattr(data, "used_" + name, None)
            if part == "calibration" and name == "spring_constant":
                # used_spring_constant is in N/nm (property)
                value = value * 1e9 if value is not None else None
            if value is not None:
                values[name] = value
        result[part] = values or None
    return result


def file_params(data):
    """Parameters of a file."""
    return {"name": os.path.basename(str(data.filename)),
            "calibration": calibration(data),
            "acquisition": acquisition(data),
            "curves": curve_processing(data),
            "elasticity": elasticity(data),
            "viscoelasticity": viscoelasticity(data),
            "microrheology": microrheology(data),
            "work_and_rupture_force": work_and_rupture_force(data),
            "events": events(data),
            "breakthrough": breakthrough(data),
            "settings": settings(data)}


def row_params(row):
    """A row of the Plots table: its file, replicate, condition, ROI,
    slice and filters of the events."""
    from ..experiment import get_short_name

    result = shared.exp.results_list[row]
    group = shared.exp.groups_list[result.group]
    condition = shared.exp.conditions_list[group.condition]
    experiment = getattr(group, "experiment", None)

    def event_filter(low, high, extra=None):
        if high <= low:
            return None
        value = {"min": low, "max": high}
        if extra:
            value.update(extra)
        return value

    return {"name": get_short_name(result),
            "file": os.path.basename(
                str(shared.exp.list[result.data_id].filename)),
            "replicate": group.name if result.group else None,
            "experiment": experiment if result.group else None,
            "condition": condition.name if group.condition else None,
            "roi": result.roi or None,
            "slice": result.slice or None,
            "distance_filter_nm": event_filter(
                result.filter_dist_min, result.filter_dist_max,
                {"mode": "keep" if result.filter_dist_keep else "remove"}),
            "slope_filter_nN_nm": event_filter(
                result.filter_slope_min, result.filter_slope_max)}


def build(export, file_name, data_ids=None, rows=None, options=None):
    """Parameters of an export: its name and file, the options of the
    result, the rows of the Plots table and the files (the files of the
    rows if data_ids is None, all the files if both are None)."""
    if rows is not None and data_ids is None:
        data_ids = sorted({shared.exp.results_list[row].data_id
                           for row in rows})
    if data_ids is None:
        data_ids = range(len(shared.exp.list))
    params = {
        "schema": SCHEMA,
        "pyaf": {"version": consts.VERSION,
                 "exported": datetime.datetime.now().astimezone().isoformat(
                     timespec="seconds"),
                 "export": export,
                 "file": os.path.basename(file_name) if file_name else None},
        "result": {"type": shared.exp.results_type,
                   "remove_values_le_0": bool(shared.exp.hist_remove_zeros)},
        "files": [file_params(shared.exp.list[i]) for i in data_ids]}
    if options:
        params["result"].update(options)
    if rows is not None:
        params["rows"] = [row_params(row) for row in rows]
    return jsonable(params)


def sidecar_path(path):
    """Path of the parameters of an exported file: name.params.json."""
    return os.path.splitext(path)[0] + ".params.json"


def write(path, params):
    """Write the parameters (indented JSON, UTF-8)."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(params, f, indent=2, ensure_ascii=False)
        f.write("\n")


def write_sidecar(path, export, data_ids=None, rows=None, options=None):
    """Write the parameters of an exported file next to it. Returns the
    path of the parameters."""
    json_path = sidecar_path(path)
    write(json_path, build(export, path, data_ids, rows, options))
    return json_path


def tiff_description(export, path, data_ids, options=None):
    """Parameters as text for the description of a TIFF image (compact
    JSON, shown by ImageJ in Image > Show Info)."""
    return json.dumps(build(export, path, data_ids, options=options),
                      ensure_ascii=False)
