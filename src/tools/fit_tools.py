"""
Python application for fitting Hertz model to different force curves
"""
import numpy as np
from lmfit import Parameters, minimize
from scipy.optimize import curve_fit
from . import curve_tools
from . import search_poc

def get_contact_point_RoV():
    pass


def prepare_data_hertzfit(app, dt, fit_params):
    """
    Prepare data for hertz fit.
    """

    # Calculate initial guess POC
    poc_index_before, real_poc, final_fit, _ = search_poc.get_POC(
        app, dt["approach_position"], fit_params)

    # Initial values for the POC index and POC pos
    poc_index_i = poc_index_before
    poc_pos_i = real_poc

    # Calculate force curve
    k = dt["spring_constant"] * 1e9
    fcapp = curve_tools.get_force_curves(app[0], app[1], poc_pos_i, k)

    app[0] = app[0] * 1e-09
    app[1] = app[1] * 1e-09
    fcapp[0] = fcapp[0] * 1e-09
    fcapp[1] = fcapp[1] * 1e-09

    # Get the index where the non corrupted data of the approach curve begins.
    first = int(dt["approach_position"][0])

    # Select data for fit
    start_m = app[0][first] + fit_params["poc_skip_start"] * 1e-09
    fit_len_m = fit_params["poc_fit_length"] * 1e-09
    stop_m = start_m + fit_len_m

    # Initial position 0 if there is no corrupted data.
    start = 0
    
    for ii in range(len(app[0])):
        if (first + ii) < len(app[0]) and app[0][first + ii] >= start_m:
            start = ii - 1
            break
    # Skip corrupted data
    if first > start:
        stop_m = app[0][first] + fit_len_m

    for ii in range(len(app[0])):
        # If first + ii exceed the length of the approach curve, use the inital poc position as stop
        if first + ii >= len(app[0]):
            stop = poc_index_i
            break
        elif app[0][first + ii] >= stop_m:
            stop = ii - 1
            break
        else:
            stop = poc_index_i

    stop = stop - first

    # If first and stop are the same value or stop is smaller than first.
    # In that case use the initial poc position as stop.
    if stop <= first:
        stop = poc_index_i

    delta = fcapp[0][first::]
    data_to_fit = fcapp[1][first::]

    # Calculate tilt and substract it from the curve. Careful!! Is this done somewhere else?
    # if dt["substract_tilt"]:
    offset_no_tilt = np.median(data_to_fit[0:stop])
    data_to_fit = data_to_fit - offset_no_tilt

    return delta, data_to_fit, final_fit, first, stop


def do_hertz_fit(fit_model, fit_model_lmfit, masked_delta, masked_data_to_fit):
    """
    Function to obtain delta0 (POC) and E0 by goodness of fit.
    """

    # Initial guesses
    e0 = 10000
    delta0 = 0

    # List containing the initial parameters for the fit.
    initial_params = [delta0, e0]

    try:
        # Perform fit using numpy.optimize.curve_fit function.
        # Levenberg-Marquardt algorithm is used for the least squares fit.
        fit_result = curve_fit(fit_model, masked_delta, masked_data_to_fit,
                               method="lm", p0=initial_params,
                               absolute_sigma=True)
        ans, _ = fit_result
        delta0_fit, e0_fit = ans

    except RuntimeError:
        params = Parameters()
        params.add('delta0_fit', value=0)
        params.add('e0_fit', value=e0)
        fit_result = minimize(fit_model_lmfit, params, args=(masked_delta, masked_data_to_fit))
        e0_fit = fit_result.params['e0_fit'].value
        delta0_fit = fit_result.params['delta0_fit'].value

    return delta0_fit, e0_fit