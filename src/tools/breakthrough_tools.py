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

"""Breakthrough of supported lipid bilayers on the approach curves (see
Gumi-Audenis & Giannotti, in Biomembrane Simulations: Computational
Studies of Biological Membranes, 2019).

The tip compresses the bilayer elastically, then punctures it: the
separation (tip-sample distance, z - deflection) jumps by the thickness of
the bilayer while the force drops slightly, and the tip is on the support.
The breakthrough force Fb is the force at the rupture, the thickness th the
jump of the separation (several breakthroughs for stacked bilayers).

A breakthrough is detected on the contact part of the approach when, within
a few points, the separation advances more than its local trend by at least
the minimum thickness and by n_sigma times the noise of its point-to-point
steps (robust standard deviation), and the force drops (by more than twice
its noise). Fb is the maximum of the force just before the jump, th the
median separation of the hard contact after the jump minus the separation
at Fb.
"""

import numpy

# Points over which a jump is searched, points of the hard contact after it
# (thickness), points of the local trend of the separation before it
JUMP_POINTS = 3
AFTER_POINTS = 8
TREND_POINTS = 10

# Minimum force drop at the rupture [nN] (besides twice the noise)
MIN_DROP = 0.05

# Values of an event in the results: force [N], thickness [nm], indices of
# the force peak and of the end of the jump on the approach curve
EVENT_VALUES = ("force", "thickness", "index_peak", "index_end")


def robust_sigma(values):
    """Robust standard deviation (1.4826 x median absolute deviation)."""
    values = values[numpy.isfinite(values)]
    if values.size == 0:
        return numpy.nan
    return 1.4826 * float(numpy.median(numpy.abs(values -
                                                 numpy.median(values))))


def force_distance(z, deflection, spring_constant):
    """Force - distance curve of the approach without point of contact:
    force [nN] from the baseline (first third of the approach), distance
    [nm] = separation from the support (median separation of the last
    points of the approach). Same length as the given arrays (nan for the
    corrupted points)."""
    z = numpy.asarray(z, dtype=float)
    deflection = numpy.asarray(deflection, dtype=float)
    valid = numpy.isfinite(z) & numpy.isfinite(deflection) & \
        (numpy.abs(deflection) < 1e5) & (numpy.abs(z) < 1e6)
    force = numpy.full(z.shape, numpy.nan)
    distance = numpy.full(z.shape, numpy.nan)
    indices = numpy.flatnonzero(valid)
    if indices.size < 2 * AFTER_POINTS:
        return distance, force
    d, separation = deflection[indices], z[indices] - deflection[indices]
    baseline = numpy.median(d[: len(d) // 3])
    force[indices] = (d - baseline) * spring_constant
    distance[indices] = separation - numpy.median(separation[-AFTER_POINTS:])
    return distance, force


def detect(z, deflection, spring_constant, n_sigma=6.0, min_thickness=1.0,
           min_force=0.2, max_events=5):
    """Breakthroughs of an approach curve.

    z: piezo position [nm], increasing toward the sample; deflection [nm];
    spring_constant [N/m]; min_thickness [nm]; min_force [nN].
    Returns the events, in the order of the approach: list of (force [nN],
    thickness [nm], index of the force peak, index of the end of the jump),
    the indices on the given arrays.
    """
    z = numpy.asarray(z, dtype=float)
    deflection = numpy.asarray(deflection, dtype=float)
    # Valid points (no corrupted values)
    valid = numpy.isfinite(z) & numpy.isfinite(deflection) & \
        (numpy.abs(deflection) < 1e5) & (numpy.abs(z) < 1e6)
    indices = numpy.flatnonzero(valid)
    if indices.size < 50:
        return []
    z, deflection = z[indices], deflection[indices]

    # Force from the baseline (first third of the approach), separation
    baseline = deflection[: len(deflection) // 3]
    noise = robust_sigma(baseline)
    force = (deflection - numpy.median(baseline)) * spring_constant
    separation = z - deflection
    noise_force = noise * spring_constant

    # Contact part of the approach
    above = force > max(min_force, 5 * noise_force)
    if not above.any():
        return []
    start = int(numpy.argmax(above))
    sep, f = separation[start:], force[start:]
    steps = numpy.diff(sep)
    sigma = robust_sigma(steps)
    threshold = max(min_thickness, n_sigma * sigma)
    min_drop = max(2 * noise_force, MIN_DROP)

    events = []
    i = 0
    while i < len(sep) - JUMP_POINTS - 1 and len(events) < max_events:
        # Advance of the separation beyond its local trend (it advances
        # with the piezo before a firm contact)
        trend = float(numpy.median(steps[max(0, i - TREND_POINTS):i])) \
            if i >= 3 else 0.0
        advance = sep[i + 1:i + 1 + JUMP_POINTS] - sep[i] - \
            trend * numpy.arange(1, JUMP_POINTS + 1)
        j = int(numpy.argmax(advance))
        drop = f[i] - numpy.min(f[i + 1:i + 2 + j])
        if advance[j] >= threshold and drop >= min_drop and \
                f[i] >= min_force:
            # Force peak just before the jump, end of the jump
            low = max(0, i - 3)
            peak = low + int(numpy.argmax(f[low:i + 1]))
            end = i + 1 + j
            while end + 1 < len(sep) and sep[end + 1] - sep[end] - trend > \
                    max(min_thickness / 3.0, 3 * sigma):
                end += 1
            thickness = float(numpy.median(sep[end:end + AFTER_POINTS]) -
                              sep[peak])
            if thickness >= min_thickness:
                events.append((float(f[peak]), thickness,
                               int(indices[start + peak]),
                               int(indices[start + end])))
            i = end + 1
            continue
        i += 1
    return events
