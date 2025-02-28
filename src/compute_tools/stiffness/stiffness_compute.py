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

"""Computes the stiffness."""

import numpy
from ...tools import search_poc
from ...tools import curve_tools
from ...tools import math_tools


def main_compute_stiff(i, j, app, ret, dt, fit_params, o_params):
    """Gets the points of contact and computes the stiffness."""
    # Get the curve
    app, _, _, _, _, smoothing_error = curve_tools.get_curve(
        None, [i, j], app, ret, mode="compute", dt=dt)

    # Get POC
    poc_index_before, real_poc, final_fit, _ = search_poc.get_POC(
        app, dt["approach_position"], fit_params)

    # Get force curve
    curve_x = app[0]
    resforce = curve_tools.get_force_curves(
        curve_x,
        app[1],
        real_poc,
        dt["spring_constant"])

    piez = dt["piezo_image"]
    poc_corrected = piez + curve_x[len(curve_x) - 1] - real_poc[0]

    # Get stiffness
    segment = slice(poc_index_before, len(resforce[1]))
    curve_x = resforce[0][segment]
    curve_y = resforce[1][segment]
    curve_x[0] = 0
    curve_y[0] = 0

    indentation_val = 0
    if o_params["indentation_stop"] != 0:
        indentation_val = o_params["indentation_stop"]
    
    else:
        indentation_val = numpy.max(curve_x)

    parts = segment_curve(
        curve_x,
        curve_y,
        o_params["indentation_start"],
        o_params["indentation_stop"],
        o_params["indentation_step"])

    if o_params["model_selected"] == 3:
        # For the linear fit (slope model), pass all the curve
        curve = [curve_x, curve_y]
    else:
        curve = None

    stiffness_array = compute_stiffness(
        parts,
        o_params["coeff"],
        o_params["model_selected"],
        o_params["tomography"],
        curve)

    return i, j, poc_index_before, real_poc, final_fit, stiffness_array, \
           poc_corrected, smoothing_error, indentation_val

def get_x_parts(start, stop, step, max_x):
    """Segments the curve in small parts (x axis)

    Once the contact point determined, the curve is segmented in parts. The
    function returns a list of parts in one dimension (x).
    """
    parts_x = []
    nbr_segments = 0
    newstop = 0

    # No step, no stop, take whole indentation curve
    if step == 0 and stop == 0:
        parts_x.append(start)
        parts_x.append(max_x)
        if start >= max_x:
            # Very rare case were the poc is set as the last point of the curve
            # Then we return no segment
            return []
        return parts_x

    # All existing segments from start to end
    if step != 0 and stop == 0:
        count = 0
        while count + start + step <= max_x:
            newstop = count + start + step
            count += step

    # Stop
    if step != 0 and stop != 0:
        count = 0
        while count + start + step <= stop and count + start + step <= max_x:
            newstop = count + start + step
            count += step

    if newstop != 0:
        nbr_segments = (newstop - start) // step
        for seg in range(nbr_segments + 1):
            parts_x.append(seg * step + start)
        return parts_x
    else:
        return parts_x


def segment_curve(curve_x, curve_y, start, stop, step):
    """Returns a list of x, y positions delimiting the segments on a curve."""
    # poc_indice is only the indice of the value before the real poc,
    # so we have to replace the first value of the curve with the real poc
    # value (0,0)
    curve_x[0] = 0.0
    curve_y[0] = 0.0

    parts_x = get_x_parts(start, stop, step, max(curve_x))
    parts_y = []
    newstart = 0

    if parts_x == []:
        return [[], []]
    else:
        for x in parts_x:
            if x == 0:
                parts_y.append(0)
            else:
                for i in range(newstart, len(curve_x)):
                    if curve_x[i] >= x:
                        coefficients = numpy.polyfit(
                            curve_x[i - 1:i + 1], curve_y[i - 1:i + 1], 1)
                        polynomial = numpy.poly1d(coefficients)
                        y = polynomial(x)
                        parts_y.append(y)
                        newstart = i
                        break

    return [parts_x, parts_y]


def compute_stiffness(parts, coeff, model_selected, tomography, curve):
    """Computes the actual stiffness for each segment.

    Returns an array of Young moduli.
    """
    parts_x = parts[0]
    parts_y = parts[1]
    stiffness_array = []

    if parts_x != []:
        if model_selected != 3:
            for i in range(len(parts_x) - 1):
                if tomography:
                    delta_F = (parts_y[i + 1] - parts_y[i])
                    delta_x = (parts_x[i + 1] - parts_x[i]) * 1e-9
                else:
                    delta_F = (parts_y[i + 1] - parts_y[0])
                    delta_x = (parts_x[i + 1] - parts_x[0]) * 1e-9

                if model_selected == 0:
                    # Hertz with sphere
                    young_modulus = coeff * (delta_F / pow(delta_x, 1.5))
                elif model_selected == 1:
                    # Hertz with cone
                    young_modulus = coeff * (delta_F / (delta_x ** 2))
                elif model_selected == 2:
                    # Pyramid
                    young_modulus = coeff * (delta_F / (delta_x ** 2))

                stiffness_array.append(young_modulus)

        elif model_selected == 3:
            curve_x = curve[0]
            curve_y = curve[1]

            count = 0
            x = parts_x[0]
            while curve_x[count] < x:
                count += 1

            for i in range(len(parts_x) - 1):
                seg_x = []
                seg_y = []

                # Go from parts_x[i] to the last point before parts_x[i+1]
                # For example when there is only one segment, the while loop
                # would go out of bounds for curve_x[count]. So we have to
                # use < (and not <=)
                while x < parts_x[i + 1]:
                    x = curve_x[count]
                    y = curve_y[count]

                    seg_x.append(x * 1e-9)
                    seg_y.append(y)

                    count += 1

                seg_x.append(parts_x[i + 1] * 1e-9)
                seg_y.append(parts_y[i + 1])

                coeffs, _ = math_tools.fit_linear(seg_x, seg_y)

                # Slope
                stiffness_array.append(coeffs[0])

    return stiffness_array
