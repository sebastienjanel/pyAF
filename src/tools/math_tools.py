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

"""Some math tools."""

import math
import numpy
import itertools
from ..tools.utils import module_exists
if module_exists("scipy"):
    from scipy import stats


def in_list(liste, val):
    """Basic function to check if a value is present in a list.

    Very slow if the list is long ...
    """
    if liste == []:
        return False
    for item in liste:
        if item == val:
            return True
    return False


def get_x_dist(curve_x, first):
    """Computes the distance between two points.

    This is needed for each curve because JPK curves have non uniform
    distances.

    Do not take into account the corrupted values of the curve for this
    computation (see bug #357). (The first value gives the position of the
    first non corrupted point on the curve)

    Needs abs value because in rare cases it can happen that there are only
    negative values in the curves (JPK bug).

    If the "first" argument is set to None, the curve_x is already the
    non-corrupted curve, so we can get the distance on the whole curve.
    """
    if first is None:
        start = 0
    else:
        start = int(first)

    dist_nm = abs(curve_x[-1] - curve_x[start])
    dist_ind = len(curve_x) - start

    if dist_nm == 0 or dist_ind == 0:
        return None
    else:
        return abs(dist_nm / dist_ind)


def compute_noises(curve_x, curve_y):
    """Function used to compute 4 types of noise."""
    if len(curve_y) < 2:
        return None, None, None, None

    # Root mean square
    # http://en.wikipedia.org/wiki/Root_mean_square
    rms = math.sqrt(sum(n * n for n in curve_y) / len(curve_y))

    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.std.html
    std = numpy.std(curve_y)

    # Coefficient of determination (R squared)
    # http://en.wikipedia.org/wiki/Coefficient_of_determination
    # http://stackoverflow.com/questions/893657/
    coeffs = numpy.polyfit(curve_x, curve_y, 1)
    p = numpy.poly1d(coeffs)
    yhat = p(curve_x)
    ybar = numpy.sum(curve_y) / len(curve_y)
    ssreg = numpy.sum((yhat - ybar) ** 2)
    sstot = numpy.sum((curve_y - ybar) ** 2)
    determination = ssreg / sstot

    # SEM
    # docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.stats.sem.html
    # http://en.wikipedia.org/wiki/Standard_error
    if module_exists("scipy"):
        sem = stats.sem(curve_y, axis=None, ddof=0)
    else:
        sem = None

    return rms, std, determination, sem


def flatten(data, order=1):
    """Flatten a meshgrid.

    You can define 3 orders, and there is only one flattening direction for the
    moment (Y).
    """
    vect_x = numpy.arange(data.shape[1])
    m_slope1 = 0
    m_slope2 = 0
    m_slope3 = 0
    # Determine the slope of the drift.
    if order == 1:
        slope1 = [numpy.polyfit(vect_x, line, 1)[0] for line in data]
        m_slope1 = numpy.mean(slope1)
        # Correct the image with the average slope.
        new_data = [line - m_slope1 * vect_x for line in data]
    elif order == 2:
        slope1 = [numpy.polyfit(vect_x, line, 2)[0] for line in data]
        m_slope1 = numpy.mean(slope1)
        slope2 = [numpy.polyfit(vect_x, line, 2)[1] for line in data]
        m_slope2 = numpy.mean(slope2)
        # Correct the image with the average slope(s).
        new_data = [line - m_slope1 * vect_x ** 2 - m_slope2 * vect_x
                    for line in data]
    elif order == 3:
        slope1 = [numpy.polyfit(vect_x, line, 3)[0] for line in data]
        m_slope1 = numpy.mean(slope1)
        slope2 = [numpy.polyfit(vect_x, line, 3)[1] for line in data]
        m_slope2 = numpy.mean(slope2)
        slope3 = [numpy.polyfit(vect_x, line, 3)[2] for line in data]
        m_slope3 = numpy.mean(slope3)
        # Correct the image with the average slope(s).
        new_data = [line - m_slope1 * vect_x ** 3 - m_slope2 * vect_x ** 2
                    - m_slope3 * vect_x for line in data]

    # Go to zero after flatten
    minimum = numpy.amin(new_data)
    new_data = new_data - minimum

    return numpy.array(new_data)


def correct_tilt(approach, retraction, limit_1, limit_2, mode,
                 approach_pos, retraction_pos):
    """
    This function corrects the tilt of the force curves.

    The same tilt is applied on the trace and the retrace, but you can chose
    with the mode option if the first fit is done on the trace or the retrace.

    The corrupted parts of the curve are not taken into account.

    """

    # Some informations on the tilt correction
    info = {
        "is_corrected": False,
        "length_of_fit_nbr_points": 0,
        "detected_slope": 0}

    # The default, to be returned in some special cases were no tilt correction
    # can be done
    default = (approach, retraction, info)

    if mode == "trace":
        curve_first = approach
        curve_second = retraction
        start_pos = approach_pos[0]
    elif mode == "retrace":
        curve_first = retraction
        curve_second = approach
        start_pos = retraction_pos[0]
    else:
        return default

    xlimit1 = 0
    xlimit2 = 0

    for k in range(start_pos, len(curve_first[0]), 1):
        if curve_first[0][k] >= limit_1:
            xlimit1 = k
            break
    nbr_points = 0
    for k in range(xlimit1, len(curve_first[0])):
        if curve_first[0][k] >= limit_2:
            xlimit2 = k
            break
        nbr_points += 1

    if xlimit1 == xlimit2:
        # No fit was found because there is no data
        # In this case, do not apply the fit
        return default

    # Get a fit between the two limits
    segment = slice(xlimit1, xlimit2)
    coeffs, _ = fit_linear(curve_first[0][segment],
                           curve_first[1][segment])

    # Get the values to substract. Curve_first can have a different
    # length than curve_second so we get the values for the two curves.
    subs_y1 = numpy.polyval([coeffs[0], coeffs[1]], curve_first[0])
    subs_y2 = numpy.polyval([coeffs[0], coeffs[1]], curve_second[0])

    # Substract
    if mode == "trace":
        result_curve_first = [curve_first[0], curve_first[1] - subs_y1]
        result_curve_second = [curve_second[0], curve_second[1] - subs_y2]
    elif mode == "retrace":
        result_curve_second = [curve_first[0], curve_first[1] - subs_y1]
        result_curve_first = [curve_second[0], curve_second[1] - subs_y2]

    info["is_corrected"] = True
    info["length_of_fit_nbr_points"] = nbr_points
    info["detected_slope"] = coeffs[0]

    return result_curve_first, result_curve_second, info


def fit_linear(curve_x, curve_y):
    """Makes a linear fit on a curve.

    Returns the coefficients of the fit and the squared residual.
    """
    # Fit, get the coefficients and some infos
    fit = numpy.polyfit(curve_x, curve_y, 1, full=True)
    coefficients = fit[0]

    # Sometimes with corrupted curves, you can get (for example) :
    # curve_x = [100, 100, 100, 100, 100, 100, 100, 100, 100]
    # curve_y = [50, 50, 50, 50, 50, 50, 50, 50, 50]
    # Info is then = [], so I return None as r_squared.
    if len(fit[1]) == 0:
        r_squared = None
    else:
        r_squared = fit[1][0]

    return [coefficients, r_squared]


def find_non_flat(curve):
    """Find non flat part of a curve."""
    indice = 0
    prev_val = curve[indice]
    next_val = curve[indice + 1]
    while prev_val == next_val and indice + 1 <= len(curve):
        indice += 1
        next_val, prev_val = curve[indice], next_val
    return indice


def find_intersection(line1, line2):
    """Find intersection between two lines, given their coefficients."""
    # y1 = ax + b and y2 = cx + d, if y1 = y2 => x(a-c) = d-b
    # x = (d-b)/(a-c)
    if line2[0] - line1[0] != 0:
        x = (line1[1] - line2[1]) / (line2[0] - line1[0])
        y = line1[0] * x + line1[1]
    else:
        # No intersection found, slopes are the same
        x, y = None, None
    return x, y


def find_line_coefs(point1, point2):
    """Find the coefficients of a line, given the position of 2 points.

    Point1 = [x1, y1]
    Point2 = [x2, y2]
    """
    div = point1[0] - point2[0]
    if div == 0:
        return None, None
    else:
        # y = ax + b
        a = (point1[1] - point2[1]) / div
        b = point1[1] - a * point1[0]
        return a, b


def polyfit2d(x, y, z, order=3, linear=False):
    """Two-dimensional polynomial fit.

    Based upon code provided by Joe Kington.

    Reference:
        http://stackoverflow.com/questions/7997152/
        python-3d-polynomial-surface-fit-order-dependent/7997925#7997925
    """
    ncols = (order + 1) ** 2
    G = numpy.zeros((x.size, ncols))
    ij = itertools.product(list(range(order + 1)), list(range(order + 1)))
    for k, (i, j) in enumerate(ij):
        G[:, k] = x ** i * y ** j
        if linear & (i != 0.) & (j != 0.):
            G[:, k] = 0
    m, _, _, _ = numpy.linalg.lstsq(G, z)

    return m


def polyval2d(x, y, m):
    """Values to two-dimensional polynomial fit.

    Based upon code provided by Joe Kington.
    """
    order = int(numpy.sqrt(len(m))) - 1
    ij = itertools.product(list(range(order + 1)), list(range(order + 1)))
    z = numpy.zeros_like(x)
    for a, (i, j) in zip(m, ij):
        z += a * x ** i * y ** j
    return z


def sg_filter(x, y, order, w, spacing):
    """The function SG_filter is an algorithm performing a Savitzky-Golay
    smoothing of a set of data y.

    Originally written by Simone Bovio in MatLab, adapted to python by
    Michka Popoff (with the help of the SG filter proposal for scipy).
    """
    # The algorithm calculate a polynomial fit of nth order on a window of data
    # of lenght w centered on one point of the curve: the value of the point is
    # then substituted with the value taken by the polynomial p at the same
    # position. The procedure is repeated for each point within the dataset.
    # (Actually, the SG filter does not perform a polyfit for every interval,
    # but actually a convolution of the date belonging to the window,
    # using a kernel calculated considering only the number of point in the
    # window and the order of the polynomial (in the case of equally spaced
    # points).
    # The output, ys, is a colums vector containing the smoothed data.
    #
    # The intervale lenght w must be an integer odd number.
    # The order of the polynomial must be n<w-1, or the fit would be
    # ill-conditioned.
    # spaceing is used to indicate to the software if to consider the points as
    # equally spaced (spaceing='unif', to insert using the '' fonts)
    # (just one single kernel is calculated for  all the intervals), or not
    # (spaceing='nunif'), in which case a new kernel is calculated for every
    # window, considering the real spaceing between subsequent points.
    #
    # Note1: x,y inputs can be either column or line vectors, since inside the
    # code they will be reshaped.
    # Note2: the multiplication of the matrix A by a constant, does not affect
    # the values of the first line of H, i.e. of the convolution kernel. If an
    # are the coefficient of the polinomial calculated as ai=Sum_j(H_ij*y_j),
    # ai', the coefficients obtained multiplying the A matrix by c will be
    # related the ai as:
    #    ai'=(c^i)*ai,
    # so a0', the coefficient we calculate with the convolution, that is, the
    # value of the polynomial evaluated at the central point the interval, is:
    #    a0'=a0.
    #

    #################################
    # New x and y vectors: at beginning and end of x and y, w/2 points are
    # mirrored in order to have a simmetric interval for the first and the last
    # point. This part can be improved allowing for a variation of the interval
    # at the extremes, or by performing a cut of the bad points from ys at the
    # end of the computations.
    firstvals = x[0] - numpy.abs(x[1:int(w / 2) + 1][::-1] - x[0])
    lastvals = x[-1] + numpy.abs(x[-int(w / 2) - 1:-1][::-1] - x[-1])
    new_x = numpy.concatenate((firstvals, x, lastvals))

    firstvals = y[0] - numpy.abs(y[1:int(w / 2) + 1][::-1] - y[0])
    lastvals = y[-1] + numpy.abs(y[-int(w / 2) - 1:-1][::-1] - y[-1])
    new_y = numpy.concatenate((firstvals, y, lastvals))

    #################################

    if spacing:
        A = numpy.zeros([w, order + 1])

        wind = numpy.arange(-int(w / 2), int(w / 2) + 1, 1)

        for ind in range(w):
            for o in range(order + 1):
                A[ind][o] = wind[ind] ** o

        ######################################

        # Part of the code used to calculate the matrix H, whos first line
        # contains the elements of the convolution kernel.
        # The calculation is now performed directly using the function pinv,
        # that should make the same matrix procedures commented below, but
        # avoids error messages concernig 0 on inf values obtained in the case
        # of the explicit calculation. After some proofs the results seem to be
        # the same.
        # To further check.

        #     AT=transpose(A);
        #     B=AT*A;
        #     C=inv(B);
        #     H=C*AT;

        H = numpy.linalg.pinv(A)

        h0 = H[0]

        ######################################

        ys = []

        c = 0

        for _ in range(len(x)):
            yc = new_y[c:(c + w)]
            y0 = numpy.dot(h0, yc)
            ys.append(y0)
            c += 1

    else:
        # Calculation of the smoothing for non-equally spaced points:
        # the kernel is calculated inside the for loop, one time for each
        # interval, so numel(x) times.

        ######################################

        # Mean step size: used to divide the x data used for calculating the
        # matrix in the case of non-equally spaced points (u==0) to avoid to
        # have to small numbers
        dx = numpy.mean(x)

        ys = []

        c = 0

        for _ in range(len(x)):
            A = numpy.zeros([w, order + 1])

            xc = new_x[c:(c + w)]
            xc = xc - xc[int(w / 2)]
            xc = xc / dx

            for ind in range(w):
                for o in range(order + 1):
                    A[ind][o] = xc[ind] ** o

            ######################################

            # Part of the code used to calculate the matrix H, whos first line
            # contains the elements of the convolution kernel.
            # The calculation is now performed directly using the function
            # pinv, that should make the same matrix procedures commented
            # below, but avoids error messages concernig 0 on inf values
            # obtained in the case of the explicit calculation. After some
            # proofs the results seem to be the same.
            # To further check.

            #     AT=transpose(A);
            #     B=AT*A;
            #     C=inv(B);
            #     H=C*AT;

            H = numpy.linalg.pinv(A)

            h0 = H[0]

            yc = new_y[c:(c + w)]
            y0 = numpy.dot(h0, yc)
            ys.append(y0)
            c += 1

    return numpy.array(ys)
