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


"""
Python application for fitting Hertz model to different force curves
"""
import warnings
import numpy as np
from lmfit import Parameters, minimize
from scipy.optimize import curve_fit, OptimizeWarning
from . import curve_tools
from . import search_poc

def prepare_data_hertzfit(app, dt, fit_params):
    """
    Prepare data for hertz fit.
    """

    # Calculate initial guess POC
    poc_index_before, real_poc, final_fit, fit_segment = search_poc.get_POC(
        app, dt["approach_position"], fit_params)
    poc_found = search_poc.poc_found(fit_segment)

    # Initial values for the POC index and POC pos
    poc_index_i = poc_index_before
    poc_pos_i = real_poc

    # Calculate force curve
    k = dt["spring_constant"] * 1e9
    fcapp = curve_tools.get_force_curves(app[0], app[1], poc_pos_i, k)

    if poc_found:
        # Baseline of the search of the point of contact (deflection vs
        # piezo, fitted just before the contact): subtracted from the
        # deflection, so that the force is 0 before the contact as for the
        # point of contact (the tilt correction leaves a slope)
        baseline = final_fit[0] * np.asarray(app[0]) + final_fit[1]
        baseline -= final_fit[0] * poc_pos_i[0] + final_fit[1]
        fcapp[1] = fcapp[1] - baseline * k
        fcapp[0] = fcapp[0] + baseline

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

    if not poc_found:
        # No baseline: the median of the beginning of the curve is used
        offset_no_tilt = np.median(data_to_fit[0:stop])
        data_to_fit = data_to_fit - offset_no_tilt

    return delta, data_to_fit, final_fit, first, stop, poc_index_i, poc_pos_i, \
        poc_found


def fit_e_fixed_poc(delta, force, coeff, exponent):
    """Elasticity fitted with the point of contact fixed at delta = 0.

    The model F = coeff * E * delta ** exponent (0 before the contact) is
    linear in E: least squares solution, E = sum(F g) / sum(g g) with
    g = coeff * delta ** exponent. Returns nan without indentation.
    """
    if exponent is None:
        return np.nan
    g = coeff * np.power(np.clip(np.asarray(delta, float), 0, None), exponent)
    norm = np.sum(g * g)
    if not np.isfinite(norm) or norm <= 0:
        return np.nan
    return float(np.sum(np.asarray(force, float) * g) / norm)


def fit_e_fixed_poc_unit(force, unit):
    """Elasticity fitted with the point of contact fixed, for a model
    F = E * unit (unit: force per unit modulus at each indentation, see
    contact_models.unit_force). Returns nan without indentation."""
    g = np.asarray(unit, float)
    norm = np.sum(g * g)
    if not np.isfinite(norm) or norm <= 0:
        return np.nan
    return float(np.sum(np.asarray(force, float) * g) / norm)


def do_hertz_fit(fit_model, fit_model_lmfit, masked_delta, masked_data_to_fit):
    """
    Function to obtain delta0 (POC) and E0 by goodness of fit.

    Returns nan for both values if the fit failed, i.e. if the covariance of
    the parameters could not be estimated. This happens when the parameters
    are not constrained by the data (no contact, too few points, ...): the
    returned values would then be meaningless (often the initial guess).
    """

    # Initial guesses
    e0 = 10000
    delta0 = 0

    # List containing the initial parameters for the fit.
    initial_params = [delta0, e0]

    try:
        # Perform fit using numpy.optimize.curve_fit function.
        # Levenberg-Marquardt algorithm is used for the least squares fit.
        with warnings.catch_warnings():
            # A failed fit is detected below with the covariance matrix
            warnings.simplefilter("ignore", OptimizeWarning)
            fit_result = curve_fit(fit_model, masked_delta,
                                   masked_data_to_fit, method="lm",
                                   p0=initial_params, absolute_sigma=True)
        ans, pcov = fit_result
        delta0_fit, e0_fit = ans
        fit_ok = np.all(np.isfinite(pcov))

    except RuntimeError:
        params = Parameters()
        params.add('delta0_fit', value=0)
        params.add('e0_fit', value=e0)
        fit_result = minimize(fit_model_lmfit, params, args=(masked_delta, masked_data_to_fit))
        e0_fit = fit_result.params['e0_fit'].value
        delta0_fit = fit_result.params['delta0_fit'].value
        fit_ok = fit_result.success and fit_result.errorbars

    if not fit_ok:
        return np.nan, np.nan

    return delta0_fit, e0_fit