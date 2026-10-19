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
from ...tools import search_poc
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
    exponent = {0: 1.5, 1: 2.0, 2: 2.0, 4: 1.0}.get(model_selected)

    if model_selected == 0:
        # Hertz with sphere
        fit_model = sm.hertz_sphere
        fit_model_lmfit = sm.hertz_sphere_lmfit
        div = numpy.sqrt(fit_params["tip_radius"] * 1e-9)
        sm.coeff = (4.0 / 3.0) * div / (1.0 - fit_params["poisson_ratio"] ** 2)

    elif model_selected == 1:
        # Sneddon (cone)
        fit_model = sm.sneddon_cone
        fit_model_lmfit = sm.sneddon_cone_lmfit
        div = numpy.tan(fit_params["tip_angle"] * numpy.pi / 180)
        sm.coeff = (2.0 / numpy.pi) * div / (1.0 - fit_params["poisson_ratio"] ** 2)

    elif model_selected == 2:
        # Barber & Billings (4-sided pyramid)
        fit_model = sm.barber_billings_pyramid
        fit_model_lmfit = sm.barber_billings_pyramid_lmfit
        div = numpy.tan(fit_params["tip_angle"] * numpy.pi / 180)
        sm.coeff = (1.0 / numpy.sqrt(2)) * div / (1.0 - fit_params["poisson_ratio"] ** 2)

    elif model_selected == 4:
        # Flat punch
        fit_model = sm.flat_punch
        fit_model_lmfit = sm.flat_punch_lmfit
        div = fit_params["tip_radius"] * 1e-9
        sm.coeff = 2 * div / (1.0 - fit_params["poisson_ratio"] ** 2)

    else:
        # Sphere (exact), blunted cone and pyramid, needle and nanowire:
        # force per unit modulus of contact_models
        fit_model = sm.tabulated
        fit_model_lmfit = sm.tabulated_lmfit
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
        fit_quality = [numpy.nan, numpy.nan]

    else:

        # A setpoint that later will be used to select the data for the fit is defined
        setpoint = fit_params["trig_threshold"] * dt["deflection_sensitivity"] * dt["spring_constant"]

        data = {"data_to_fit": data_to_fit[stop::], "delta": delta[stop::]}
        df = pd.DataFrame(data)

        if fit_params["fit_range_type"] == 0:

            # Variables used for selecting the data for the fit based on force range
            start_f = o_params["force_start"]/100 * setpoint
            stop_f = o_params["force_stop"]/100 * setpoint

            filtered_df = df[(start_f <= df.data_to_fit) & (df.data_to_fit <= stop_f)]

        elif fit_params["fit_range_type"] == 1:

            # Variables used for selecting the data for the fit based on indentation range
            start_f = o_params["indentation_start"] * 1e-09   # 0 is the default value.
            stop_f = o_params["indentation_stop"] * 1e-09

            if stop_f == 0:
                # No stop: up to the maximum indentation
                filtered_df = df[start_f <= df.delta]
            else:
                filtered_df = df[(start_f <= df.delta) & (df.delta <= stop_f)]

        # Indentation reached, in nm (delta is in m): up to the selected stop
        # in indentation range mode, else the end of the fitted range
        if fit_params["fit_range_type"] == 1 and \
                o_params["indentation_stop"] != 0:
            indentation_val = min(o_params["indentation_stop"],
                                  numpy.max(delta) * 1e09)
        else:
            indentation_val = filtered_df["delta"].max() * 1e09

        masked_data_to_fit = numpy.r_[data_to_fit[0:stop], filtered_df["data_to_fit"].values]
        masked_delta = numpy.r_[delta[0:stop], filtered_df["delta"].values]

        if too_stiff:
            # No free fit of the point of contact (it can not be found by
            # the fit without indentation)
            delta0_fit, e0_fit = numpy.nan, numpy.nan
        else:
            delta0_fit, e0_fit = fit_tools.do_hertz_fit(
                fit_model, fit_model_lmfit, masked_delta, masked_data_to_fit)

        # The fit can move the point of contact by a part of the indentation
        # (delta0_fit is relative to the detected point of contact)
        window = o_params["fit_poc_window"] / 100.0 * max_indentation

        if numpy.isnan(e0_fit) or abs(delta0_fit) > window:
            # Too stiff sample, failed fit, or point of contact moved too
            # far: fit of the elasticity only, with the detected point of
            # contact
            poc_index_fit = poc_index_i
            real_poc = poc_pos_i
            if exponent is None:
                e_fixed = fit_tools.fit_e_fixed_poc_unit(
                    masked_data_to_fit, sm.unit(masked_delta))
            else:
                e_fixed = fit_tools.fit_e_fixed_poc(
                    masked_delta, masked_data_to_fit, sm.coeff, exponent)
            stiffness = [] if numpy.isnan(e_fixed) else [e_fixed]
            fit_delta0, fit_e, fixed_poc = 0.0, e_fixed, 1.0
        else:
            poc_index_fit = (numpy.abs(delta - delta0_fit)).argmin() + first
            real_poc = [app[0][poc_index_fit] * 1e09, app[1][poc_index_fit] * 1e09]
            stiffness = [e0_fit]
            fit_delta0, fit_e, fixed_poc = delta0_fit, e0_fit, 0.0

        # Quality of the fit: R^2 on the fitted points after the contact,
        # and 1 if the point of contact was fixed (fallback), 0 if fitted
        fit_quality = [numpy.nan, numpy.nan]
        if stiffness:
            fit_quality = [r_squared(
                filtered_df["data_to_fit"].values,
                fit_model(filtered_df["delta"].values, fit_delta0, fit_e)),
                fixed_poc]

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

