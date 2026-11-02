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

"""Filters of the image stacks (8 bits, (frames, rows, cols)), as the ITK
filters used before: BinaryThresholdImageFilter (same result), and
DiscreteGaussianImageFilter with its default parameters (at most 2 gray levels
of difference, from the rounding of the sums; isosurfaces of 10 microscopy
stacks, sigma from 0.5 to 8: at most 0.025 % of the voxels of the region and
1.1 % of the volume)."""

import math

import numpy

# Default parameters of itk::DiscreteGaussianImageFilter
GAUSSIAN_MAX_ERROR = 0.01
GAUSSIAN_MAX_KERNEL_WIDTH = 32
# Accuracy of the downward recurrence of the Bessel functions
BESSEL_ACCURACY = 40


def binary_threshold(image, lower, upper, outside, inside):
    """inside for the values from lower to upper, else outside."""
    if lower > upper:
        raise ValueError("Lower threshold cannot be greater than upper "
                         "threshold.")
    # Not wrapped around (300 would become 44)
    if not all(0 <= value <= 255 for value in (lower, upper, outside, inside)):
        raise ValueError("The thresholds and the values must be from 0 to "
                         "255.")
    return numpy.where((image >= lower) & (image <= upper), inside,
                       outside).astype(numpy.uint8)


def bessel_i0(y):
    """Modified Bessel function of order 0 (polynomial approximation, as
    itk::GaussianOperator; not scipy.special.iv: the rounding of the kernel
    decides the truncated values, see discrete_gaussian)."""
    if y < 3.75:
        m = (y / 3.75) ** 2
        return 1.0 + m * (3.5156229 + m * (3.0899424 + m * (1.2067492 + m * (0.2659732 + m * (0.360768e-1 + m * 0.45813e-2)))))
    m = 3.75 / y
    return (math.exp(y) / math.sqrt(y)) * (0.39894228 + m * (0.1328592e-1 + m * (0.225319e-2 + m * (-0.157565e-2 + m * (0.916281e-2 + m * (-0.2057706e-1 + m * (0.2635537e-1 + m * (-0.1647633e-1 + m * 0.392377e-2))))))))


def bessel_i1(y):
    """Modified Bessel function of order 1 (see bessel_i0)."""
    if y < 3.75:
        m = (y / 3.75) ** 2
        return y * (0.5 + m * (0.87890594 + m * (0.51498869 + m * (0.15084934 + m * (0.2658733e-1 + m * (0.301532e-2 + m * 0.32411e-3))))))
    m = 3.75 / y
    a = 0.2282967e-1 + m * (-0.2895312e-1 + m * (0.1787654e-1 - m * 0.420059e-2))
    a = 0.39894228 + m * (-0.3988024e-1 + m * (-0.362018e-2 + m * (0.163801e-2 + m * (-0.1031555e-1 + m * a))))
    return a * (math.exp(y) / math.sqrt(y))


def bessel_in(order, y):
    """Modified Bessel function of order >= 2: downward recurrence (see
    bessel_i0)."""
    to_y = 2.0 / y
    previous = result = 0.0
    current = 1.0
    for j in range(2 * (order + int(math.sqrt(BESSEL_ACCURACY * order))), 0, -1):
        previous, current = current, previous + j * to_y * current
        if abs(current) > 1e10:
            result *= 1e-10
            current *= 1e-10
            previous *= 1e-10
        if j == order:
            result = previous
    return result * bessel_i0(y) / current


def gaussian_kernel(variance):
    """Discrete Gaussian kernel (Lindeberg; not a sampled Gaussian), up to
    the maximal error or the maximal width, normalized."""
    scale = math.exp(-variance)
    half = [scale * bessel_i0(variance), scale * bessel_i1(variance)]
    total = half[0] + 2.0 * half[1]
    while total < 1.0 - GAUSSIAN_MAX_ERROR:
        half.append(scale * bessel_in(len(half), variance))
        total += 2.0 * half[-1]
        if len(half) > GAUSSIAN_MAX_KERNEL_WIDTH:
            break
    half = [value / total for value in half]
    return numpy.array(half[:0:-1] + half)


def discrete_gaussian(image, sigma):
    """Gaussian smoothing, sigma in pixels, one axis after the other, the
    last pixels repeated at the edges. As in ITK, the values are truncated to
    integers after each axis: the isosurfaces are drawn at the level 1,
    which depends on it (+14 % of volume with a plain
    scipy.ndimage.gaussian_filter)."""
    kernel = gaussian_kernel(sigma ** 2)
    radius = len(kernel) // 2
    values = image.astype(float)
    for axis in range(values.ndim):
        width = [(0, 0)] * values.ndim
        width[axis] = (radius, radius)
        padded = numpy.pad(values, width, mode="edge")
        result = numpy.zeros(values.shape)
        part = [slice(None)] * values.ndim
        # Sum in the order of the kernel (rounding, see bessel_i0)
        for start, coefficient in enumerate(kernel):
            part[axis] = slice(start, start + values.shape[axis])
            result += coefficient * padded[tuple(part)]
        values = numpy.floor(result)
    return numpy.clip(values, 0, 255).astype(numpy.uint8)
