"""Bottom effect correction (B.E.C.) of the elasticity, and the points of
the elasticity fit of a curve.

The B.E.C. is applied inside the fit model: the elasticity and the point
of contact are fitted again on the points of the elasticity fit with
F = E g(delta - delta0) phi(delta - delta0, h) (g: contact model, phi:
factor of Garcia & Garcia 2018, h: height of the sample). The point of
contact has to be fitted again: the fit without the correction moves it
into the sample to follow the steeper curve, and the elasticity fitted
with that point of contact would be overestimated.
"""

import numpy

from ...tools import curve_tools
from ...tools import fit_tools
from . import contact_models


def has_height(data):
    """True if a height is defined for the bottom effect correction (sample
    height, or substrate ROI with a plane fit)."""
    if data.user_h is not None and data.user_h != 0:
        return True
    return data.roi_glass_id != 0 and \
        data.roi_list[data.roi_glass_id - 1].glass_coeffs is not None


def check_file(data):
    """Returns None if the B.E.C. can be computed for the file, else the
    reason."""
    if not data.stiffness_calculated:
        return "the elasticity has not been computed"
    if data.used_stiffness_model_selected not in contact_models.HAS_BEC:
        return "no bottom effect correction for the model " + \
            contact_models.MODEL_NAMES[data.used_stiffness_model_selected] + \
            " (Garcia & Garcia 2018: paraboloid, cone, pyramid, flat punch)"
    if data.user_h is not None and data.user_h != 0:
        return None
    if data.roi_glass_id == 0:
        return "select a substrate ROI in the list (ROI set as substrate " \
            "in the ROI manager), or enter the sample height"
    if data.roi_list[data.roi_glass_id - 1].glass_coeffs is None:
        return "the selected ROI is no longer set as substrate"
    return None


def get_height_map(data):
    """Height of the sample [nm] for the bottom effect correction (nan if
    not corrected: below the noise floor)."""
    nx, ny = data.nbr_pixels_x, data.nbr_pixels_y
    if data.user_h is not None and data.user_h != 0:
        return numpy.full((nx, ny), float(data.user_h))
    d, a, b, _ = data.roi_list[data.roi_glass_id - 1].glass_coeffs
    x_size = data.scan_size_x / nx
    y_size = data.scan_size_y / ny
    cols = numpy.arange(nx)[:, None] * x_size
    rows = numpy.arange(ny)[None, :] * y_size
    height = numpy.asarray(data.topography, float) - (a * rows + b * cols + d)
    noise_floor = data.noise_floor if data.noise_floor is not None else 0
    height[~(height >= noise_floor)] = numpy.nan
    return height


def approach_points(data, i, j):
    """Indentation [m] and force [N] of the valid part of the approach,
    relative to the point of contact of the elasticity fit, with the
    baseline of the search of the point of contact removed (as in the
    fit), and the piezo position [nm]. None without point of contact."""
    poc = data.pocs_real[i][j]
    z_poc, d_poc = float(poc[0]), float(poc[1])
    if not (numpy.isfinite(z_poc) and numpy.isfinite(d_poc)):
        return None
    app, _, _, _, _, _ = curve_tools.get_curve(data, [i, j],
                                               mode="curve_results")
    start, end = data.approach_positions[i][j]
    z = numpy.asarray(app[0], float)[int(start):int(end)]
    d = numpy.asarray(app[1], float)[int(start):int(end)]
    if len(z) < 2:
        return None
    slope = float(data.fits_poc[i][j][0])
    if not numpy.isfinite(slope):
        slope = 0.0
    d = d - slope * (z - z_poc)
    k = data.used_spring_constant
    if k is None:
        k = data.spring_constant
    defl = d - d_poc
    return (z - z_poc - defl) * 1e-9, defl * k, z


def fit_range_mask(data, delta, force):
    """Points of the fit range of the elasticity fit (force or indentation
    range), delta in m and force in N."""
    mask = delta > 0
    range_type = data.get_used_fit_range_type()
    if range_type == 0:
        start = data.used_force_start if data.used_force_start is not None \
            else data.force_start
        stop = data.used_force_stop if data.used_force_stop is not None \
            else data.force_stop
        setpoint = data.trig_threshold * data.deflection_sensitivity * \
            data.spring_constant
        mask &= (force >= start / 100 * setpoint) & \
            (force <= stop / 100 * setpoint)
    elif range_type == 1:
        mask &= delta >= data.used_indentation_start * 1e-9
        if data.used_indentation_stop != 0:
            mask &= delta <= data.used_indentation_stop * 1e-9
    if not numpy.any(mask):
        mask = delta > 0
    return mask


def unit_force(data, delta, height=None):
    """Force per unit modulus of the model used for the elasticity, with the
    B.E.C. factor for a height [m] (none if None or nan)."""
    model = data.used_stiffness_model_selected
    radius = data.used_tip_radius * 1e-9
    g = contact_models.unit_force(model, radius, data.used_tip_angle,
                                  data.used_poisson_ratio, delta,
                                  data.get_used_shell_thickness() * 1e-9)
    if height is not None and numpy.isfinite(height):
        g = g * contact_models.bec_factor(delta, model, radius,
                                          data.used_tip_angle, height)
    return g


def fit_points(data, i, j):
    """Points of the elasticity fit of a curve: baseline before the contact
    (segment of the search of the point of contact) and contact points of
    the fit range. Indentation [m] and force [N], or None."""
    points = approach_points(data, i, j)
    if points is None:
        return None
    delta, force, z = points
    skip = data.used_fitparam_poc_skip_start or 0
    length = data.used_fitparam_poc_fit_length or 0
    baseline = (delta <= 0) & (z <= z[0] + skip + length)
    contact = fit_range_mask(data, delta, force)
    keep = baseline | contact
    if numpy.count_nonzero(contact) < 2:
        return None
    return delta[keep], force[keep]


def fit_elasticity(data, delta, force, height):
    """Elasticity fitted with the point of contact free (delta0, within the
    window of the elasticity fit), else with the point of contact fixed,
    for a height [m] of the B.E.C. (None: no correction)."""

    def model(x, delta0, e0):
        return e0 * unit_force(data, numpy.asarray(x) - delta0, height)

    def model_lmfit(params, x, y):
        return y - model(x, params["delta0_fit"].value,
                         params["e0_fit"].value)

    delta0, e0 = fit_tools.do_hertz_fit(model, model_lmfit, delta, force)
    window = (data.used_fit_poc_window or 0) / 100.0 * numpy.max(delta)
    if numpy.isfinite(e0) and abs(delta0) <= window:
        return e0
    # Point of contact of the elasticity fit (E = sum(F g) / sum(g g))
    g = unit_force(data, delta, height)
    norm = numpy.sum(g * g)
    return numpy.sum(force * g) / norm if norm > 0 else numpy.nan


def correction_factor(data, i, j, height):
    """Factor of the elasticity of a pixel for the B.E.C. (<= 1): ratio of
    the elasticities fitted with and without the B.E.C. on the points of
    the elasticity fit (point of contact fitted again).

    height in nm; 1 if there is no height (nan: below the noise floor), nan
    if the curve can not be fitted.
    """
    if not numpy.isfinite(height):
        return 1.0
    points = fit_points(data, i, j)
    if points is None:
        return numpy.nan
    delta, force = points
    e = fit_elasticity(data, delta, force, None)
    e_bec = fit_elasticity(data, delta, force, height * 1e-9)
    if not (numpy.isfinite(e) and numpy.isfinite(e_bec)) or e <= 0 or \
            e_bec <= 0:
        return numpy.nan
    return float(e_bec / e)
