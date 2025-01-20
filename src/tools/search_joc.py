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

"""Get jump of contact on a curve."""

import numpy
from . import math_tools


def get_JOC(curve, cor_pos, params):
    """Function used to search for a jump of contact on a curve."""
    curve_x = curve[0]
    curve_y = curve[1]
    first = cor_pos[0]

    joc_skip_start = params["joc_skip_start"]
    joc_fit_length = params["joc_fit_length"]
    joc_refit_option = params["joc_refit_option"]
    joc_noise_multiplicator = params["joc_noise_multiplicator"]
    joc_refit_times = params["joc_refit_times"]

    # Default value if no fit is found
    default = len(curve_x) - 3, len(curve_x) - 1, \
        [curve_x[len(curve_x) - 3], curve_y[len(curve_x) - 3]], \
        [curve_x[len(curve_x) - 1], curve_y[len(curve_x) - 1]], \
        [0, 0, 0], [0, len(curve_x)]

    # Get distance between to points on the curve
    xdist = math_tools.get_x_dist(curve_x, first)

    if xdist is None:
        # Some corrupted curves have bad data with no distance between points.
        return default

    # Parameters for the fit
    joc_skip_start = int(joc_skip_start / xdist)
    joc_fit_length = int(joc_fit_length / xdist)
    joc_refit_option = int(joc_refit_option / xdist)

    # Skip corrupted data
    if first > joc_skip_start:
        joc_skip_start = first

    # Return directly if the curve is too short
    if joc_skip_start + joc_fit_length >= len(curve_x):
        return default

    start = joc_skip_start
    stop = start + joc_fit_length
    seg = slice(start, stop)

    coeffs, r_squared = math_tools.fit_linear(curve_x[seg], curve_y[seg])
    startfit = start
    endfit = stop

    # No fit was found
    if r_squared is None:
        return default

    # Get the "noise" of the fit
    error = numpy.sqrt(r_squared / len(curve_x[seg])) * joc_noise_multiplicator

    # Make fitting lines, shifted by the noise
    y2 = numpy.polyval([coeffs[0], coeffs[1] - error], curve_x)
    y = numpy.polyval([coeffs[0], coeffs[1]], curve_x)

    joc1_indice = len(curve_x) - 3
    # Find intersection with -yerror
    for i in range(start, len(curve_y) - 1):
        if curve_y[i] <= y2[i]:
            joc1_indice = i
            break
    # Go back and find intersection with y
    for i in range(joc1_indice, 0, -1):
        if curve_y[i] >= y[i]:
            joc1_indice = i
            break

    # Refit option
    if joc_refit_option != 0:
        for _ in range(joc_refit_times):
            newstart = joc1_indice - joc_fit_length - joc_refit_option
            if first > newstart:
                newstart = first
            newstop = newstart + joc_fit_length
            seg = slice(newstart, newstop)

            # Get the "noise" of the fit
            coeffs, r_squared = math_tools.fit_linear(curve_x[seg],
                                                      curve_y[seg])

            # No fit was found
            if r_squared is None:
                return default

            error = numpy.sqrt(r_squared / len(curve_x[seg])) *\
                joc_noise_multiplicator

            # Make fitting lines, shifted by the noise
            y2 = numpy.polyval([coeffs[0], coeffs[1] - error], curve_x)
            y = numpy.polyval([coeffs[0], coeffs[1]], curve_x)

            # Find intersection with -yerror
            for i in range(newstart, len(curve_y) - 1):
                if curve_y[i] <= y2[i]:
                    joc1_indice = i
                    break

            # Go back and find intersection with y
            for i in range(joc1_indice, 0, -1):
                if curve_y[i] >= y[i]:
                    joc1_indice = i
                    break
            startfit = newstart
            endfit = newstop

    final_fit = [coeffs[0], coeffs[1], - error]

    # Finish Joc2 ------------------------------------------------------------
    joc2_indice = len(curve_y) - 1
    for i in range(len(curve_y) - 1, 0, -1):
        if curve_y[i] <= y[i]:
            joc2_indice = i
            break

    if joc2_indice <= joc1_indice:
        joc1_indice = len(curve_y) - 2
        joc2_indice = len(curve_y) - 1
        real_joc1x, real_joc1y = curve_x[len(curve_x) - 3], \
            curve_y[len(curve_x) - 3]
        real_joc2x, real_joc2y = curve_x[len(curve_x) - 2], \
            curve_y[len(curve_x) - 2]
        final_fit = [0, 0, 0]
    else:
        if joc2_indice + 1 == len(curve_x):
            joc2_indice = joc2_indice - 1

        # Search for real joc1 in nm (nicer)
        a, b = math_tools.find_line_coefs([curve_x[joc1_indice],
                                           curve_y[joc1_indice]],
                                          [curve_x[joc1_indice + 1],
                                           curve_y[joc1_indice + 1]])

        if a is not None and b is not None:
            real_joc1x, real_joc1y = \
                math_tools.find_intersection([coeffs[0], coeffs[1]], [a, b])
        else:
            real_joc1x, real_joc1y = None, None

        # Search for real joc2 in nm (nicer)
        a, b = math_tools.find_line_coefs([curve_x[joc2_indice],
                                           curve_y[joc2_indice]],
                                          [curve_x[joc2_indice + 1],
                                           curve_y[joc2_indice + 1]])

        if a is not None and b is not None:
            real_joc2x, real_joc2y = \
                math_tools.find_intersection([coeffs[0], coeffs[1]], [a, b])
        else:
            real_joc2x, real_joc2y = None, None

        # Intersection found too far away or not found
        not_found = False

        if (real_joc1x is None and real_joc1y is None) or \
                real_joc1x > curve_x[joc1_indice + 1] or \
                real_joc1x < curve_x[joc1_indice]:
            real_joc1x = curve_x[joc1_indice - 3]
            real_joc1y = curve_y[joc1_indice - 3]
            joc1_indice = len(curve_y) - 2
            not_found = True

        if (real_joc2x is None and real_joc2y is None) or \
                real_joc2x > curve_x[joc2_indice + 1] or \
                real_joc2x < curve_x[joc2_indice] or not_found:
            real_joc2x = curve_x[joc2_indice - 1]
            real_joc2y = curve_y[joc2_indice - 1]
            joc2_indice = len(curve_y) - 1

    # Return the results
    return joc1_indice, joc2_indice, \
        [real_joc1x, real_joc1y], [real_joc2x, real_joc2y], \
        final_fit, [startfit, endfit]
