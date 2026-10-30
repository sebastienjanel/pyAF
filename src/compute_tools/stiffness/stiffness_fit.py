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

"""Computes the stiffness by fitting the selected model to experimental data."""

import numpy
import pandas as pd
from ...tools import curve_tools
from ...tools import fit_tools
from . import stiffness_models as sm
from . import contact_models


def main_fit_stiffness(i, j, app, ret, dt, fit_params, o_params):
    """Gets the points of contact and computes the stiffness fitting
    Hertz model using Levenberg-Marquard(least squares) algorithm."""

    model_selected = o_params["model_selected"]
    # Exponent of the indentation in the model (for the fit with a fixed
    # point of contact)
    exponent = {0: 1.5, 1: 2.0, 2: 2.0, 4: 1.0,
                contact_models.THIN_SHELL: 1.0}.get(model_selected)

    if model_selected == 0:
        # Hertz with sphere
        fit_model = sm.hertz_sphere
        div = numpy.sqrt(fit_params["tip_radius"] * 1e-9)
        sm.coeff = (4.0 / 3.0) * div / (1.0 - fit_params["poisson_ratio"] ** 2)

    elif model_selected == 1:
        # Sneddon (cone)
        fit_model = sm.sneddon_cone
        div = numpy.tan(fit_params["tip_angle"] * numpy.pi / 180)
        sm.coeff = (2.0 / numpy.pi) * div / (1.0 - fit_params["poisson_ratio"] ** 2)

    elif model_selected == 2:
        # Barber & Billings (4-sided pyramid)
        fit_model = sm.barber_billings_pyramid
        div = numpy.tan(fit_params["tip_angle"] * numpy.pi / 180)
        sm.coeff = (1.0 / numpy.sqrt(2)) * div / (1.0 - fit_params["poisson_ratio"] ** 2)

    elif model_selected == 4:
        # Flat punch
        fit_model = sm.flat_punch
        div = fit_params["tip_radius"] * 1e-9
        sm.coeff = 2 * div / (1.0 - fit_params["poisson_ratio"] ** 2)

    elif model_selected == contact_models.THIN_SHELL:
        # Thin shell (Reissner): linear, as the flat punch, the radius being
        # the one of the vesicle
        fit_model = sm.flat_punch
        sm.coeff = contact_models.power_law(
            model_selected, fit_params["tip_radius"] * 1e-9, 0.0,
            fit_params["poisson_ratio"],
            fit_params["shell_thickness"] * 1e-9)[0]

    else:
        # Sphere (exact), blunted cone and pyramid, needle and nanowire:
        # force per unit modulus of contact_models
        fit_model = sm.tabulated
        radius = fit_params["tip_radius"] * 1e-9
        angle = fit_params["tip_angle"]
        nu = fit_params["poisson_ratio"]
        sm.unit = lambda d: contact_models.unit_force(
            model_selected, radius, angle, nu, d)

    # Get approach curve.
    app, _, _, _, _, smoothing_error = curve_tools.get_curve(None, [i, j], app, ret,
                                                             mode="compute", dt=dt, get_time=False)

    # Get data for the fit.
    delta, data_to_fit, final_fit, first, stop, poc_index_i, poc_pos_i, \
        poc_found = fit_tools.prepare_data_hertzfit(app, dt, fit_params)

    # Indentation at the end of the curve (from the detected point of
    # contact), and piezo travel after the contact [m]
    max_indentation = numpy.max(delta) if len(delta) else 0.0
    travel = app[0][-1] - app[0][poc_index_i] if poc_found else 0.0
    # Sample too stiff for the cantilever (e.g. glass): the deflection
    # follows the piezo, the indentation is too small to be fitted
    too_stiff = travel > 0 and max_indentation < \
        o_params["fit_min_indentation"] / 100.0 * travel

    if not poc_found or len(delta) <= 1 or len(data_to_fit) <= 1:
        # No initial POC or no data to fit: no stiffness value, keep the
        # initial POC
        poc_index_fit = poc_index_i
        real_poc = poc_pos_i
        stiffness = []
        indentation_val = 0
        fit_quality = [numpy.nan, numpy.nan, numpy.nan]

    else:

        # A setpoint that later will be used to select the data for the fit is defined
        setpoint = fit_params["trig_threshold"] * dt["deflection_sensitivity"] * dt["spring_constant"]

        # Fitted points: the baseline of the search of the point of contact,
        # and a contiguous part of the curve (index range, not a point by
        # point selection on the noisy force: that kept only the points
        # with a positive noise near the contact, which biased the point of
        # contact and the modulus). From the start of the range (force or
        # indentation) to the first point beyond its stop; with a range
        # starting at the contact, the points between the baseline window
        # and the detected point of contact are kept (baseline near the
        # contact for the fit of the point of contact).
        n_pts = len(delta)
        contact = numpy.nonzero(delta >= 0)[0]
        poc_k = int(contact[0]) if len(contact) else n_pts
        poc_k = max(poc_k, stop)
        if fit_params["fit_range_type"] == 0:
            # Force range, percentages of the setpoint
            start_f = o_params["force_start"] / 100 * setpoint
            stop_f = o_params["force_stop"] / 100 * setpoint
            values = data_to_fit
            from_contact = o_params["force_start"] <= 0
        else:
            # Indentation range [m], no stop: up to the maximum indentation
            start_f = o_params["indentation_start"] * 1e-09
            stop_f = o_params["indentation_stop"] * 1e-09 or numpy.inf
            values = delta
            from_contact = o_params["indentation_start"] <= 0
        after = numpy.nonzero(values[poc_k:] >= start_f)[0]
        s_k = poc_k + int(after[0]) if len(after) else n_pts
        beyond = numpy.nonzero(values[s_k:] > stop_f)[0]
        e_k = s_k + int(beyond[0]) if len(beyond) else n_pts
        filtered_df = pd.DataFrame({"data_to_fit": data_to_fit[s_k:e_k],
                                    "delta": delta[s_k:e_k]})
        near_baseline = slice(stop, s_k) if from_contact else slice(0, 0)

        # Indentation reached, in nm (delta is in m): up to the selected stop
        # in indentation range mode, else the end of the fitted range
        if fit_params["fit_range_type"] == 1 and \
                o_params["indentation_stop"] != 0:
            indentation_val = min(o_params["indentation_stop"],
                                  numpy.max(delta) * 1e09)
        else:
            indentation_val = filtered_df["delta"].max() * 1e09

        masked_data_to_fit = numpy.r_[data_to_fit[0:stop],
                                      data_to_fit[near_baseline],
                                      filtered_df["data_to_fit"].values]
        masked_delta = numpy.r_[delta[0:stop], delta[near_baseline],
                                filtered_df["delta"].values]

        # The fit can move the point of contact (delta0_fit is relative to
        # the detected point of contact): later by a part of the indentation
        # at most, earlier down to the end of the baseline window. The
        # detected point is late by construction (threshold above the
        # noise of the baseline), by much more than this part for the
        # models whose force rises slowly (cone, pyramid): the earlier
        # point of contact of the fit has to be accepted
        window = o_params["fit_poc_window"] / 100.0 * max_indentation
        earliest = delta[min(stop, len(delta) - 1)] if len(delta) else 0.0

        if too_stiff:
            # No free fit of the point of contact (it can not be found by
            # the fit without indentation)
            delta0_fit, e0_fit = numpy.nan, numpy.nan
        else:
            # Global least squares in this range (E linear for each point
            # of contact)
            if exponent is None:
                unit = sm.unit
            else:
                unit = lambda d: sm.coeff * numpy.power(  # noqa: E731
                    numpy.clip(d, 0, None), exponent)
            delta0_fit, e0_fit = fit_tools.fit_poc_profile(
                unit, masked_delta, masked_data_to_fit, earliest, window)

        # A modulus not positive is not a fit of a contact (no contact in
        # the curve)
        if numpy.isnan(e0_fit) or not e0_fit > 0 or delta0_fit > window or \
                delta0_fit < earliest:
            # Too stiff sample, failed fit, or point of contact moved too
            # far: fit of the elasticity only, with the detected point of
            # contact
            poc_index_fit = poc_index_i
            # Point of contact of the fit: the detected one, on the
            # baseline (zero of the force of the fit)
            real_poc = [poc_pos_i[0],
                        final_fit[0] * poc_pos_i[0] + final_fit[1]]
            if exponent is None:
                e_fixed = fit_tools.fit_e_fixed_poc_unit(
                    masked_data_to_fit, sm.unit(masked_delta))
            else:
                e_fixed = fit_tools.fit_e_fixed_poc(
                    masked_delta, masked_data_to_fit, sm.coeff, exponent)
            stiffness = [e_fixed] if e_fixed > 0 else []
            fit_delta0, fit_e, fixed_poc = 0.0, e_fixed, 1.0
        else:
            poc_index_fit = (numpy.abs(delta - delta0_fit)).argmin() + first
            # Point of contact of the fit: delta is z - z_poc - (d -
            # baseline(z)) from the detected point (see
            # prepare_data_hertzfit), so delta0_fit is its shift along z,
            # with the deflection on the baseline (not a point of the
            # noisy curve)
            z_fit = poc_pos_i[0] + delta0_fit * 1e09
            real_poc = [z_fit, final_fit[0] * z_fit + final_fit[1]]
            stiffness = [e0_fit]
            fit_delta0, fit_e, fixed_poc = delta0_fit, e0_fit, 0.0

        # Quality of the fit: R^2 on the fitted points after the contact,
        # and 1 if the point of contact was fixed (fallback), 0 if fitted
        # and shift of the point of contact of the fit from the detected
        # one along z [nm] (the detected point is then known for the plots)
        fit_quality = [numpy.nan, numpy.nan, numpy.nan]
        if stiffness:
            fit_quality = [r_squared(
                filtered_df["data_to_fit"].values,
                fit_model(filtered_df["delta"].values, fit_delta0, fit_e)),
                fixed_poc, fit_delta0 * 1e09]

    # Apply corrections to get the POC corrected value
    piez = dt["piezo_image"]
    if poc_found:
        # app was converted to m by prepare_data_hertzfit, real_poc is in nm
        poc_corrected = piez + app[0][len(app[0]) - 1] * 1e09 - real_poc[0]
    else:
        # No POC: no topography value
        poc_corrected = numpy.nan

    return i, j, poc_index_fit, real_poc, final_fit, stiffness, poc_corrected, smoothing_error, indentation_val, \
        fit_quality


def r_squared(data, model):
    """Coefficient of determination of a fit (nan with less than 3 points
    or without variation of the data)."""
    data = numpy.asarray(data, dtype=float)
    if len(data) < 3:
        return numpy.nan
    ss_tot = numpy.sum((data - numpy.mean(data)) ** 2)
    if not ss_tot > 0:
        return numpy.nan
    return float(1.0 - numpy.sum((data - numpy.asarray(model)) ** 2) / ss_tot)

