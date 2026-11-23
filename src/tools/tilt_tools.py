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

"""Automatic tilt correction: the baseline of each curve found without
limits, and one slope for the file.

The baseline of a curve is its first part (far from the surface), up to
the first sustained deviation from its line (contact, adhesion, events):
fitted on a first part of the curve, then extended or shortened to the
first points which stay further than a few times its noise from the line,
until it does not change. Its length follows each curve (cells of
different heights, short ramps).

The tilt comes from the instrument (optics, interferences on the sample):
its slope is about the same for all the curves of a file. The slope of the
file is the weighted median of the slopes of the baselines (weight: number
of points times the square of the height range, the precision of a slope),
precise even with short baselines; each curve keeps its own offset."""

import numpy

# Points of a baseline at least
MIN_POINTS = 8
# Deviation from the line [times the noise] starting the end of the baseline
DEVIATION = 4.0
# Part of the points of a window which deviate: sustained deviation
SUSTAINED = 0.8
# Part of the baseline left before its end (start of the deviation)
MARGIN = 0.05
# Baselines shorter than this are reported (too few points)
SHORT_BASELINE = 20
# Slopes further than this from the slope of the file are reported [times
# the spread of the slopes]
OUTLIER = 5.0
# A baseline shorter than this part of the curve is searched again after
# its end (noise of the turnaround of the piezo at the start of the curves
# of some instruments)
RETRY_BELOW = 0.3
RETRIES = 3
# Level of a baseline found again: not further towards the contact than
# this from the first one [part of the deflection range of the curve]
SAME_LEVEL = 0.04


def line_fit(z, d):
    """Slope and intercept of a line, refitted without the points further
    than 3 times the noise (spikes)."""
    slope, intercept = numpy.polyfit(z, d, 1)
    residuals = d - (slope * z + intercept)
    noise = 1.4826 * numpy.median(numpy.abs(residuals))
    if noise > 0:
        kept = numpy.abs(residuals) <= 3 * noise
        if kept.sum() >= MIN_POINTS:
            slope, intercept = numpy.polyfit(z[kept], d[kept], 1)
    return slope, intercept


def window_noise(d, window):
    """Noise of the mean of window consecutive points of a curve, from the
    second differences of the means of its successive windows: a straight
    part (tilt, contact) does not count, nor a few curved parts (median).
    At the scale of the window: also the slow noise (drift, correlated
    noise) and the steps of the digitizer, unlike the noise of the
    points."""
    blocks = len(d) // window
    if blocks < 4:
        return 0.0
    means = d[:blocks * window].reshape(blocks, window).mean(axis=1)
    second = numpy.diff(means, 2)
    return 1.4826 * numpy.median(numpy.abs(second - numpy.median(second))) \
        / numpy.sqrt(6.0)


def baseline(z, d, slope=None):
    """Baseline of a curve: z and deflection from the start of the valid
    points (far from the surface first). With slope None, the slope of the
    baseline is fitted; otherwise only its offset. Returns (slope, offset,
    first point, end (excluded) of the baseline in the valid points, noise
    of the mean of a window of points), or None (too few points).

    A short baseline is searched again after its end (noise of the
    turnaround of the piezo at the start of the curve): the longest one is
    kept, if not shifted towards the contact (not a part of the contact).
    With a large tilt over the curve, the shift of the tilt hides it: the
    first baseline is kept (short, reported)."""
    z = numpy.asarray(z, dtype=float)
    d = numpy.asarray(d, dtype=float)
    valid = numpy.isfinite(z) & numpy.isfinite(d)
    z, d = z[valid], d[valid]
    found = baseline_from_start(z, d, slope)
    if found is None:
        return None
    best_start, best = 0, found
    level = numpy.median(d[:found[2]])
    tolerance = SAME_LEVEL * (numpy.max(d) - numpy.min(d))
    # Direction of the contact (end of the curve) from the baseline
    towards = numpy.sign(numpy.median(d[-max(MIN_POINTS, len(d) // 20):])
                         - level)
    start = 0
    for _ in range(RETRIES):
        if best[2] >= RETRY_BELOW * len(z):
            break
        start += found[2]
        found = baseline_from_start(z[start:], d[start:], slope)
        if found is None:
            break
        part = slice(start, start + found[2])
        # Not shifted towards the contact (not a part of the contact)
        shift = (numpy.median(d[part]) - level) * towards
        if found[2] > best[2] and shift <= tolerance:
            best_start, best = start, found
    s, offset, number, noise = best
    return s, offset, best_start, best_start + number, noise


def baseline_from_start(z, d, slope=None):
    """Baseline from the first point of the curve (finite values): returns
    (slope, offset, number of points, noise), or None.

    The line of each first part of the curve (all the lengths at once, with
    cumulative sums) is compared with its last points: the baseline ends
    where they depart from it more than their noise allows (contact,
    adhesion, events)."""
    number = len(z)
    if number < MIN_POINTS:
        return None
    window = max(5, number // 100)
    noise = window_noise(d, window)
    if noise <= 0:
        noise = numpy.std(d) / numpy.sqrt(window) or 1e-12

    # Sums of the first k points (k = 0 ... number)
    def sums(values):
        return numpy.concatenate([[0.0], numpy.cumsum(values)])
    # Centered heights: precision of the sums
    z0 = z - z[0]
    s1 = numpy.arange(number + 1, dtype=float)
    sz, sd = sums(z0), sums(d)
    szz, szd = sums(z0 * z0), sums(z0 * d)
    k = numpy.arange(max(MIN_POINTS, 2 * window), number + 1)
    if len(k) == 0:
        return None
    if slope is None:
        det = s1[k] * szz[k] - sz[k] ** 2
        det[det == 0] = numpy.nan
        a = (s1[k] * szd[k] - sz[k] * sd[k]) / det
    else:
        a = numpy.full(len(k), float(slope))
    b = (sd[k] - a * sz[k]) / s1[k]
    # Mean residual of the last window of points of each first part
    last_z = sz[k] - sz[k - window]
    last_d = sd[k] - sd[k - window]
    residual = (last_d - a * last_z - b * window) / window
    departs = numpy.abs(residual) > DEVIATION * noise
    departs |= ~numpy.isfinite(residual)
    # First sustained departure (over a window of lengths)
    run = numpy.convolve(departs, numpy.ones(window) / window, "valid")
    sustained = numpy.nonzero(run >= SUSTAINED)[0]
    if len(sustained):
        end = int(k[sustained[0]]) - window
    else:
        end = number
    end = min(number, max(MIN_POINTS, int(end * (1 - MARGIN))))
    if slope is None:
        s, intercept = line_fit(z0[:end], d[:end])
    else:
        s = float(slope)
        intercept = numpy.median(d[:end] - s * z0[:end])
    # Offset of the line in the heights of the curve
    return s, intercept - s * z[0], end, noise


def correct(curve_first, curve_second, start, slope):
    """Curves corrected with the slope of the file and the offset of the
    baseline of curve_first (from the index start). None if the curve has
    no baseline."""
    found = baseline(curve_first[0][start:], curve_first[1][start:], slope)
    if found is None:
        return None
    offset = found[1]
    first = [curve_first[0],
             curve_first[1] - (slope * curve_first[0] + offset)]
    second = [curve_second[0],
              curve_second[1] - (slope * curve_second[0] + offset)]
    return first, second


def weighted_median(values, weights):
    """Median of values with weights."""
    order = numpy.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = numpy.cumsum(weights)
    return values[numpy.searchsorted(cumulative, cumulative[-1] / 2.0)]


def file_slope(baselines):
    """Slope of a file from the baselines of its curves: [(slope, number of
    points, height range)]. Returns (slope, [indexes of the baselines too
    short], [indexes of the slopes far from the others]), slope None
    without baselines."""
    if not baselines:
        return None, [], []
    slopes = numpy.array([b[0] for b in baselines], dtype=float)
    points = numpy.array([b[1] for b in baselines], dtype=float)
    ranges = numpy.array([b[2] for b in baselines], dtype=float)
    weights = points * ranges ** 2
    if not numpy.any(weights > 0):
        weights = numpy.ones(len(slopes))
    slope = weighted_median(slopes, weights)
    spread = 1.4826 * numpy.median(numpy.abs(slopes - slope))
    short = [k for k in range(len(slopes)) if points[k] < SHORT_BASELINE]
    far = [] if spread <= 0 else \
        [k for k in range(len(slopes))
         if abs(slopes[k] - slope) > OUTLIER * spread]
    return float(slope), short, far


def data_baselines(data, mode, progress=None):
    """Baselines of all the curves of a file (raw deflection, as recorded):
    mode "trace" (approach) or "retrace" (retraction). Returns
    [(i, j, slope, number of points, height range)] of the curves with a
    baseline (not discarded, in the region of the map)."""
    if mode == "trace":
        curves, positions = data.curves_approach, data.approach_positions
    else:
        curves, positions = data.curves_retraction, data.retraction_positions
    result = []
    for i in range(data.nbr_pixels_x):
        row = curves[i]
        for j in range(data.nbr_pixels_y):
            if data.discarded_curves[i][j] or data.is_outside(i, j):
                continue
            start, stop = int(positions[i][j][0]), int(positions[i][j][1])
            z = numpy.asarray(row[j][0][start:stop], dtype=float)
            d = numpy.asarray(row[j][1][start:stop], dtype=float)
            if len(z) < MIN_POINTS or numpy.ptp(z) == 0:
                continue
            found = baseline(z, d)
            if found is None:
                continue
            slope, _, first, end, _ = found
            result.append((i, j, slope, end - first,
                           abs(z[end - 1] - z[first])))
        if progress is not None:
            progress()
    return result
