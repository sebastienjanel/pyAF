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

"""Computes the point of contact with a linear fit."""

import numpy
from . import math_tools


def get_POC(curve, cor_pos, params):
    """Get the point of contact.

    If the start fitting parameter is in the part of the curve where the curve
    is corrupted (see approach_positions array), the start parameter is shifted
    and set to the first value of data which is not corrupted.
    """
    curve_x = curve[0]
    curve_y = curve[1]
    first = cor_pos[0]

    poc_skip_start = params["poc_skip_start"]
    poc_fit_length = params["poc_fit_length"]
    poc_refit_option = params["poc_refit_option"]
    poc_noise_multiplicator = params["poc_noise_multiplicator"]
    poc_refit_times = params["poc_refit_times"]

    # Default value if no fit is found
    default_real = [curve_x[len(curve_x) - 2], curve_y[len(curve_x) - 2]]
    default = (len(curve_x) - 2), default_real, [0, 0, 0], [0, 0]

    # Get distance between to points on the curve
    xdist = math_tools.get_x_dist(curve_x, first)

    if xdist is None:
        # Some corrupted curves have bad data with no distance between points.
        return default

    # Parameters for the fit
    skip_start = int(poc_skip_start / xdist)
    fit_length = int(poc_fit_length / xdist)
    refit_option = int(poc_refit_option / xdist)

    # Skip corrupted data
    if first > skip_start:
        skip_start = first

    if skip_start + fit_length >= len(curve_x):
        # Fitting params too long, return default
        return default

    poc_index_before = 0
    start = skip_start
    stop = start + fit_length
    seg = slice(start, stop)

    coeffs, r_squared = math_tools.fit_linear(curve_x[seg], curve_y[seg])
    startfit = start
    endfit = stop

    if r_squared is None:
        # In rare cases no fit can be found, return default
        return default

    # Get the "noise" of the fit
    error = numpy.sqrt(r_squared / len(curve_x[seg])) * poc_noise_multiplicator

    # Make a fitting line, shifted by the noise
    y1 = numpy.polyval([coeffs[0], coeffs[1] + error], curve_x)

    for i in range(len(curve_y) - 1, 0, -1):
        if curve_y[i] <= y1[i]:
            poc_index_before = i
            break

    # Refit option
    if refit_option != 0:
        for _ in range(poc_refit_times):
            newstart = poc_index_before - fit_length - refit_option
            if first > newstart:
                newstart = first
            newstop = newstart + fit_length
            seg = slice(newstart, newstop)

            coeffs, r_squared = \
                math_tools.fit_linear(curve_x[seg], curve_y[seg])

            if r_squared is None:
                # In rare cases no fit can be found, return default
                return default

            error = numpy.sqrt(
                r_squared / len(curve_x[seg])) * poc_noise_multiplicator
            y1 = numpy.polyval([coeffs[0], coeffs[1] + error], curve_x)

            for i in range(len(curve_y) - 1, 0, -1):
                if curve_y[i] <= y1[i]:
                    poc_index_before = i
                    break

            startfit = newstart
            endfit = newstop

    # Security to always have a POC
    if poc_index_before >= len(curve_x) - 1:
        poc_index_before = len(curve_x) - 2

    final_fit = [coeffs[0], coeffs[1], error]

    # Search for real poc in nm (nicer)
    a, b = math_tools.find_line_coefs(
        [curve_x[poc_index_before], curve_y[poc_index_before]],
        [curve_x[poc_index_before + 1], curve_y[poc_index_before + 1]])

    if a is not None and b is not None:
        real_pocx, real_pocy = math_tools.find_intersection(
            [coeffs[0], coeffs[1] + error], [a, b])
    else:
        real_pocx, real_pocy = None, None

    # Intersection found too far away or not found
    if (real_pocx is None and real_pocy is None) or \
            real_pocx > curve_x[poc_index_before + 1] or \
            real_pocx < curve_x[poc_index_before]:
        real_pocx = curve_x[poc_index_before]
        real_pocy = curve_y[poc_index_before]

    return poc_index_before, [real_pocx, real_pocy], \
        final_fit, [startfit, endfit]
