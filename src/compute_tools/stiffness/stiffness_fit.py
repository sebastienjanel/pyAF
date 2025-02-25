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


def main_fit_stiffness(i, j, app, ret, dt, fit_params, o_params):
    """Gets the points of contact and computes the stiffness fitting
    Hertz model using Levenberg-Marquard(least squares) algorithm."""

    model_selected = o_params["model_selected"]

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
        # Bilodeau (Pyramid)
        #print("fit_params:", fit_params)
        fit_model = sm.bilodeau_pyramid
        fit_model_lmfit = sm.bilodeau_pyramid_lmfit
        div = numpy.tan(fit_params["tip_angle"] * numpy.pi / 180)
        sm.coeff = (1.0 / numpy.sqrt(2)) * div / (1.0 - fit_params["poisson_ratio"] ** 2)

    # Get approach curve.
    app, _, _, _, _, smoothing_error = curve_tools.get_curve(None, [i, j], app, ret,
                                                             mode="compute", dt=dt, get_time=False)

    # Get data for the fit.
    delta, data_to_fit, final_fit, first, stop = fit_tools.prepare_data_hertzfit(app, dt, fit_params)

    if len(delta) <= 1 or len(data_to_fit) <= 1:
        # Check if there is no data to fit data. Skip if no data is found.
        return

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

            filtered_df = df[(start_f <= df.delta) & (df.delta <= stop_f)]
        
        indentation_val = 0
        if o_params["indentation_stop"] != 0:
            indentation_val = o_params["indentation_stop"] * 1e-09
        
        else:
            indentation_val = filtered_df["delta"].max()

        masked_data_to_fit = numpy.r_[data_to_fit[0:stop], filtered_df["data_to_fit"].values]
        masked_delta = numpy.r_[delta[0:stop], filtered_df["delta"].values]

        delta0_fit, e0_fit = fit_tools.do_hertz_fit(fit_model, fit_model_lmfit, masked_delta, masked_data_to_fit)

        poc_index_fit = (numpy.abs(delta - delta0_fit)).argmin() + first
        real_poc = [app[0][poc_index_fit] * 1e09, app[1][poc_index_fit] * 1e09]

    # Apply corrections to get the POC corrected value
    piez = dt["piezo_image"]
    poc_corrected = piez + app[0][len(app[0]) - 1] - real_poc[0]

    return i, j, poc_index_fit, real_poc, final_fit, [e0_fit], poc_corrected, smoothing_error, indentation_val

