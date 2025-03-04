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

"""
Tools to detect the events on a force curve.

Two methods are proposed. The first one applies a convolution on the curve,
the other one uses the MSF algorithm.
In the two cases, the force detection is made with linear fits once the
event is detected.

"""

import numpy
from . import math_tools


def apply_kernel(kernel, curve_y, mode=None):
    """Convolves the curves.

    Different modes are aviable, depending on the border condition you want to
    have.
    """
    if mode is None:
        # Normal calculation
        active_mode = "same"
    elif mode == "special":
        active_mode = "same"
    else:
        # Force the mode, for example for display, it's only the valid part
        active_mode = mode

    convolution = numpy.convolve(curve_y, kernel, mode=active_mode)

    if mode == "special":
        diff = len(kernel) + 1
        fillval1 = convolution[diff // 2 + 1]
        for i in range(0, diff // 2):
            convolution[i] = fillval1
        fillval2 = convolution[len(convolution) - diff // 2 - 1]
        for i in range(len(convolution) - diff // 2, len(convolution)):
            convolution[i] = fillval2

    return convolution


def find_peaks(curve_x, curve_y, threshold, stop):
    """Finds the peaks one one convolution.

    Makes a linear fit in the flat part of the curve and checks if convolution
    peaks are above the defined threshold.

    The inputed curve_x and curve_y should not contain the corrupted part of
    the curve. This has to be removed before calling find_peaks, or this will
    result in a bad peak detection.
    """
    # Default value to return in case no fit is found due to a corrupted
    # curve
    default = [], None

    # Get distance between two points on the curve
    xdist = math_tools.get_x_dist(curve_x, first=None)

    if xdist is None:
        # Some corrupted curves have bad data with no distance between points.
        return default

    values = []
    count = 0
    lastfound = 0

    # Don't search peaks after the joc in the curve (stop)
    # Also, check if the stop limit is not longer than the curve
    if stop >= len(curve_y):
        stop = len(curve_y)

    # Parse the curve to find the peaks
    # Skip 10 values, there can be noise at the very beggining of the curve,
    # leading to wrongly detected events
    for i in range(10, stop):
        # Pass under the fitting limit
        if curve_y[i] <= threshold and count != 0:
            startindex = i - count
            endindex = i
            middleindex = int(startindex + abs(endindex - startindex) / 2.0)
            if middleindex < 1.0:
                middleindex = endindex

            lastfound = endindex
            values.append([startindex, middleindex, endindex])
            count = 0

        # Pass above the fitting limit
        if curve_y[i] > threshold and i > lastfound + 1:
            count += 1

    # Make sure the first position is at least one point on the left of the
    # middle position
    for val in values:
        if val[0] == val[1]:
            val[0] -= 1

    return values


def get_events_positions_by_kernel(curve_x, curve_y, stop, first, parameters):
    r"""
    Gets the positions of the peaks on the defl-ext curve.

    The curve is convoluted with a kernel.

    The curves are stripped from their corrupted part to be able to find the
    positions. The returned indices are shifted back to correspond to the
    right indicie position on the full curve.

    Adaptive method replaces a linear threshold by a threshold based on local
    average of the convolution. A savitzky golay smoothing of the convolution
    curve is removed from the convolution itself.

    Linear threshold:
    ::    ..           .
    ::---.--.---..----. .-------------.---
    ::  ..  ..... . .    .             .
    ::              .     .          .
    ::                     .        .
    ::                      .      .
    ::                       .    .
    ::                        .   .
    ::                         . .
    ::                          .

    Adaptive threshold:
    ::    _.._    __     .
    ::---/.  .\--/..\---. .-------------.---
    ::  ..  ..... . .     .             .
    ::              .     .           .
    ::                     .        .
    ::                      .      .
    ::                       .    .
    ::                        .   .
    ::                         . .
    ::                          .

    """

    # Remove corrupted part of the curve
    curve_x = curve_x[first:-1]
    curve_y = curve_y[first:-1]

    conv = get_convolution_curve(curve_y, parameters["kernel_size"])

    # Adaptive convolution correction
    if parameters["adaptive_threshold_option"]:
        curve_sg = math_tools.sg_filter(
            curve_x,
            conv,
            parameters["adaptive_smoothing_order"],
            parameters["adaptive_smoothing_window"],
            True)

        conv = conv - curve_sg

    # Find the peaks
    positions = find_peaks(curve_x, conv, parameters["threshold"], stop)

    # The algorithm has worked on the non corrupted part of the curve
    # so that now we have to shift the position of the indice.
    pos = []
    for i in range(len(positions)):
        pos.append(first + positions[i])

    return pos


def find_fits_around_middle_pos(curve, positions, fit_size):
    """Find the positions and the coefficients of the fits around the events.

    seg_left_fit_base, seg_right_fit_base are the boundaries (indexes) of the
    points used by the fit.
    seg_left_fit, seg_right_fit are the boundaries (indexes) of the fit line
    itself, to be displayed to the user and used for the determination of the
    event's force. It is bigger than the base fits by the value of step.
    """
    curve_x = curve[0]
    curve_y = curve[1]

    # Default length of the fit, could be smaller at the end if the fit is
    # corrected by the algorithm.
    length = fit_size
    # Move away 3 points from the middle event position
    step = 3
    fits = []

    for k in range(len(positions)):
        # Middle index
        pos = positions[k][1]

        # LEFT FIT ############################################################

        # Default value
        effective_length = length
        # Check if we go out of the boundaries of the curve (on the left)
        if pos - length < 0:
            effective_length = pos
        if k > 0:
            # Check for another event on the left
            if pos - length < positions[k - 1][1]:
                effective_length = pos - positions[k - 1][1]

        # Small check if we trim too much the fit, go back a little bit.
        if effective_length <= step:
            effective_length = step + 1

        # The fit is done one the indices defined in seg_left_fit_base
        # (used in the slice defined as segment_left)
        # The fit itself is "interpreted" on the segment seg_left_fit
        # which is a little is bigger by a value of "step" indices.
        seg_left_fit_base = [pos - effective_length, pos - step]
        seg_left_fit = [pos - effective_length, pos]
        segment_left = slice(seg_left_fit_base[0], seg_left_fit_base[1])

        # Get the left fit's coefficients
        seg1_x = curve_x[segment_left]
        seg1_y = curve_y[segment_left]
        coeffs_left, _ = math_tools.fit_linear(seg1_x, seg1_y)

        # Sometimes the left fit is too long, and if two events are very close
        # it can happen that the slope of this fit is negative. Check for this,
        # and reduce the fit size. (Do not reduce too much by limiting the
        # effective length to be greater than the step size)
        while coeffs_left[0] < -0.0001 and effective_length > step + 1:
            effective_length -= 1
            seg_left_fit_base = [pos - effective_length, pos - step]
            seg_left_fit = [pos - effective_length, pos]
            segment_left = slice(seg_left_fit_base[0], seg_left_fit_base[1])
            seg1_x = curve_x[segment_left]
            seg1_y = curve_y[segment_left]
            coeffs_left, _ = math_tools.fit_linear(seg1_x, seg1_y)

        # RIGHT FIT ###########################################################

        # Default value
        effective_length = length
        if k < len(positions) - 1:
            # Check for another event on the right
            if pos + length > positions[k + 1][1]:
                # In this case, the new effective lenght is redefined.
                # +1 because it can happen that effective_length == length
                # We need at least 2 values to make a fitting segment.
                effective_length = positions[k + 1][1] - pos + 1

        seg_right_fit_base = [pos + step, pos + effective_length]
        seg_right_fit = [pos, pos + effective_length]
        segment_right = slice(seg_right_fit_base[0], seg_right_fit_base[1])

        # RIGHT FIT CORRECTION ################################################

        # Check if the fit goes above the maximum force of the event itself
        # or if the fit goes below 2 * the force delta of the event.
        #
        # ******* curve
        # - - - - fit
        # X boundaries of right fit
        #
        # Case up, right fit is limited by left fit
        #
        # ******* - - - X - -
        #       *     *
        #       *   *
        #       * *
        #       X
        #
        # Case down, right fit is limited 2 * delta
        #
        # ******* - - - - - -
        #       *         |
        #       *         |  delta
        #       X         |
        # -------------------
        #         *       |
        #           *     | 2 * delta
        #             *   |
        #               X |
        # -------------------

        correct = False
        # To have an estimation of the force delta, use the start and end
        # positions of the event found by the peak detection algorithm
        pos_start = positions[k][0]
        pos_end = positions[k][2]

        # Check if pos_start == pos_middle. Depening on the threshold,
        # the peak detection algorithm (find_peaks method) can find
        # start and middle indices which are very close.
        # In this case move pos_start two indices to the left to have bigger
        # force estimation (else, delta = 0). This hinders that the right fit
        # correction is triggered in wrong cases. This is like "stretching"
        # the event's size.
        if pos_start == positions[k][1]:
            pos_start -= 2

        # Loop through the curve positions and check if we meet one of the
        # two conditions. In the case we meet one of the conditions, set
        # correct to false, and define a new effective_length
        for a in range(int(pos), int(pos + effective_length), 1):
            ref = 0
            # print(f"Out of bounds: a={a}, pos_start={pos_start}")
            current = curve_y[a] - curve_y[pos_start]
            delta2 = abs(2 * (curve_y[pos_start] - curve_y[pos_end]))
            # Case up
            if (ref - current) < 0:
                effective_length = a - pos
                correct = True
                break
            # Case down
            elif (ref - current) > delta2:
                effective_length = a - pos
                correct = True
                break

        if correct:
            # If the correction is too strong, use a default fit size of
            # step + 1
            if effective_length <= step:
                effective_length = step + 1
            # Define new right fits in case we need to correct the fit.
            seg_right_fit_base = [pos + step, pos + effective_length]
            seg_right_fit = [pos, pos + effective_length]
            segment_right = slice(seg_right_fit_base[0], seg_right_fit_base[1])

        # Get the right fit's coefficients
        seg2_x = curve_x[segment_right]
        seg2_y = curve_y[segment_right]

        coeffs_right, _ = math_tools.fit_linear(seg2_x, seg2_y)

        # Check for a fit which is negative on the right size. This happens
        # sometimes when two fits are too close. In this case reduce the
        # fit size to try to improve a little bit the fit. Do not reduce
        # less than the step size.
        while coeffs_right[0] < -0.0001 and effective_length > step + 1:
            effective_length -= 1
            seg_right_fit_base = [pos + step, pos + effective_length]
            seg_right_fit = [pos, pos + effective_length]
            segment_right = slice(seg_right_fit_base[0], seg_right_fit_base[1])
            seg1_x = curve_x[segment_right]
            seg1_y = curve_y[segment_right]
            coeffs_right, _ = math_tools.fit_linear(seg1_x, seg1_y)

        # SAVE FIT ############################################################

        a = {
            "coeffs_left": coeffs_left,
            "coeffs_right": coeffs_right,
            "seg_left_fit_base": seg_left_fit_base,
            "seg_right_fit_base": seg_right_fit_base,
            "seg_left_fit": seg_left_fit,
            "seg_right_fit": seg_right_fit}

        fits.append(a)

    return fits


def find_force_from_fits(curve, fits, positions):
    """Computes the forces for all the events of the curve.

    Returns the exact position of the chosen force segment [[x1, y1], [x2, y2]]
    and the value of the force == abs(y2 - y1)

    positions contains for each event the estimation of the start, middle and
    end position of the event (in indices)
    """
    curve_x = curve[0]
    curve_y = curve[1]

    forces = []
    fit_pos = []
    fit_pos_vert = []

    # Compute the force between Y1 and Y2
    # The position of Y1 and Y2 is determined as the intersection of the
    # horizontal fits, given as input parameters, and a vertical fit,
    # which follows the step (this fit is done with the y1, y2 estimation
    # done by one of the algorithms (kernel or msf))
    #
    # -*-*-*-*-*-*-*--Y1
    #                *
    #                 X              *
    #                   *       *
    #                 Y2--*-*-*-----

    for i in range(len(fits)):
        fit = fits[i]

        coef_left = fit["coeffs_left"]
        coef_right = fit["coeffs_right"]

        segment = slice(positions[i][0], positions[i][2] + 1)
        seg_x = curve_x[segment]
        seg_y = curve_y[segment]
        coeffs, _ = math_tools.fit_linear(seg_x, seg_y)

        # Position of Y1, Y2
        xleft, yleft = math_tools.find_intersection(coeffs, coef_left)
        xright, yright = math_tools.find_intersection(coeffs, coef_right)

        # Get the force
        found_force = abs(yleft - yright)
        forces.append(found_force)
        fit_pos.append([
            [curve_x[positions[i][1]], yleft],
            [curve_x[positions[i][1]], yright]])

        fit_pos_vert.append([[xleft, xright], [yleft, yright]])

    return fit_pos, fit_pos_vert, forces


def get_events_positions_by_msf(curve, threshold, window_size, stop, first):
    """Detect the position of the events with the MSF algorithm."""
    curve_x = curve[0]
    curve_y = curve[1]

    # Remove corrupted part of the curve
    curve_x = curve_x[first:-1]
    curve_y = curve_y[first:-1]

    # Apply the kernel
    msf = apply_msf(curve_x, curve_y, window_size)

    curve_x = curve_x[int(window_size / 2.0):-int(window_size / 2.0)]

    # Find the peaks
    stop = stop - int(window_size / 2.0)  # Has to be shifted
    positions = find_peaks(curve_x, msf, threshold, stop)

    # The algorithm has worked on the non corrupted part of the curve
    # so that now we have to shift the position of the indice.
    pos = []
    for i in range(len(positions)):
        pos.append(first + positions[i] + int(window_size / 2.0))

    return pos


def apply_msf(curve_x, curve_y, window_size):
    """Fitting routine performed on the window function for each datapoint.
    Linear fitting 2 times first using the step function (algebric least square
    calculation of the angular coefficient m and of the 2 intercepts csx and
    cdx) and then the linear fit.
    Then storage of the fit functions in the 2 variables f (concatenation
    of the left and right fits), and g (the whole fit).
    d is the difference between the differences between the 2 fits and the data
    points.
    """
    window_size = int(numpy.fix(window_size / 2) * 2)

    l1 = int(window_size / 2.0)
    d = []

    i1 = list(range(0, l1))
    i2 = list(range(l1, window_size))

    curve_x = numpy.array(curve_x)
    curve_y = numpy.array(curve_y)
    coef = 2 / float(window_size)
    norm = (float(window_size) * pow((0.01), 3))

    for c in range(len(curve_x) - 2 * l1):
        xc = curve_x[c:(c + window_size)]
        yc = curve_y[c:(c + window_size)]

        # Least square calculation for the step function "fi"
        div = l1 * numpy.sum(xc ** 2) - numpy.sum(xc[i1]) ** 2 - \
            numpy.sum(xc[i2]) ** 2
        m = (
            l1 * numpy.sum(xc * yc) - numpy.sum(xc[i1]) * numpy.sum(yc[i1]) -
            numpy.sum(xc[i2]) * numpy.sum(yc[i2])) / div

        csx = coef * (numpy.sum(yc[i1]) - m * numpy.sum(xc[i1]))
        cdx = coef * (numpy.sum(yc[i2]) - m * numpy.sum(xc[i2]))

        # Automatic least square (polyfit order 1), for the simple linear fit
        p = numpy.polyfit(xc, yc, 1)

        fsx = m * xc[i1] + csx
        fdx = m * xc[i2] + cdx
        f = numpy.concatenate((fsx, fdx))
        g = p[0] * xc + p[1]

        # Calculation of the residues used for discriminating the steps.
        resg = numpy.sum((g - yc) ** 2)
        resf = numpy.sum((f - yc) ** 2)

        d.append((resg - resf) * (csx - cdx) / norm)
        # The normalization factor has been estimated on calculating the
        # residues considering that on a perfect step the distance between
        # the step function and the curve is 0, while between the straight
        # line and the curve is proportional to Window size/step^2. Then,
        # considering in the calculation of d there is a further factor
        # step size, I divised for a factor of 10pN as step height at the
        # third power, and by window size.)

    return numpy.array(d)


def normalize_convolution(curve):
    """You can do 2 experiments on different days and have two different noise
    levels on you curves ... or use a different sampling in z (nbr points on
    the curve). You should NEVER do this but still, some people try ...
    Still, the events should be the same. To prevent problems, the convolutions
    are normalized on their flat part to take into account the noise.
    The solution is of course not perfect, you should ensure to have the same
    noise level and sampling on you curves to be able to compare your data.
    This is only a security making it a little bit more clean when you have big
    differences in noise or sampling.

    This method was proposed by Simone Bovio and Michka Popoff.
    """
    # Normalize the convolution
    start = 100
    stop = int(len(curve) / 2.0)

    if start > stop:
        # This can happen for Peak force files with only 128 points ...
        start = 0
        stop = int(len(curve) / 4.0)

    length = stop - start + 1

    val = numpy.sum(abs(curve[start:stop]))

    if val != 0:
        div = (val / length)
        curve = curve / div

    return curve


def create_kernel(kernel_size):
    """Create kernel with the kernel size chosen by the user.

    Example for kernel_size = 1          or 2          or 3

    ::
                                                             *
                                             *                *
                                 *            *                *
         zero :  ---------------  * --------   *  -----------   *
                                   *            *                *
                                                 *                *
                                                                   *

    Only the symetry is important, K value as no effect
    """
    n = 2 * int(kernel_size)

    k = -100.0
    r = 2 * abs(k)
    kernel = []
    kernel.append(k)
    for _ in range(n):
        k = k + r // n
        kernel.append(k)

    return kernel


def get_convolution_curve(curve_y, kernel_size):
    """Apply convolution to curve for kernel detection method.

    The curve are also normalized here.
    """
    kernel = create_kernel(kernel_size)

    # Apply the convolution
    curve = apply_kernel(kernel, curve_y, mode="special")

    # Normalize the convolutions
    normalized_curve = normalize_convolution(curve)

    return normalized_curve
