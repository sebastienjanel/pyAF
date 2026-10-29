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

"""Viscous drag coefficient of the cantilever at the surface, b(0), for the
microrheology (Alcaraz et al., Langmuir 2002 and Biophys. J. 2003).

The cantilever is oscillated at several heights h above the surface,
without contact (JPK modulation segments separated by retract segments).
At each height, the force is the drag of the liquid only:
b(h) = Im(F / delta) / (2 pi f), delta = z - d (PyFMRheo, ComputeBh,
routine ViscousDragSteps). b(h) is then extrapolated to the surface with
b(h) = 6 pi eta a_eff^2 / (h + h_eff), fitted by least squares, and
b(0) = 6 pi eta a_eff^2 / h_eff.

Units: SI (m, N, s).
"""

import math

import numpy
from scipy.optimize import curve_fit

from .. import shared
from . import microrheology
from ..tools import curve_tools


def available(data):
    """True if the file has modulations at two heights above the surface or
    more (viscous drag measurement)."""
    heights = data.modulation_heights or []
    return len({h for h in heights if h > 0}) >= 2


def drag_files():
    """Loaded files with a viscous drag measurement."""
    return [data for data in shared.exp.list if available(data)]


def curve(data, i, j):
    """b(h) of a curve: list of (height [m], frequency [Hz], b [N s/m]) for
    the modulations above the surface, None without modulation."""
    segments = curve_tools.get_modulation_segments(data, i, j)
    if segments is None:
        return None
    k = data.spring_constant * 1e9
    results = []
    for (frequency, z, d, t), height in zip(segments, data.modulation_heights):
        if height <= 0 or len(z) <= 8 or not numpy.ptp(t) > 0:
            continue
        zd, dd, td = microrheology.detrend(frequency, z, d, t)
        if len(zd) <= 8:
            continue
        h, _ = microrheology.transfer_function(frequency, zd, dd, td, "fft",
                                               k)
        results.append((height * 1e-9, frequency,
                        h.imag / (2.0 * math.pi * frequency)))
    return results


def compute(data):
    """b(h) of all the curves of the file (median and standard deviation at
    each height) and its extrapolation to the surface. Returns a dict:
    heights [m], b, b_std [N s/m], nbr (curves), b0 [N s/m], h_eff [m]
    and r2 of the fit (nan if the fit is not possible)."""
    values = {}
    nbr = 0
    for i in range(data.nbr_pixels_x):
        for j in range(data.nbr_pixels_y):
            if data.discarded_curves[i][j]:
                continue
            result = curve(data, i, j)
            if not result:
                continue
            nbr += 1
            for height, _, b in result:
                values.setdefault(round(height, 12), []).append(b)
    heights = numpy.array(sorted(values))
    b = numpy.array([numpy.nanmedian(values[h]) for h in heights])
    b_std = numpy.array([numpy.nanstd(values[h]) for h in heights])
    b0, h_eff, r2 = fit(heights, b)
    return {"heights": heights, "b": b, "b_std": b_std, "nbr": nbr,
            "b0": b0, "h_eff": h_eff, "r2": r2}


def fit(heights, b):
    """Least squares fit of b(h) = A / (h + h_eff) (Alcaraz et al. 2002):
    b(0) = A / h_eff, h_eff and R2 (nan if not possible). Initial values:
    1 / b linear in h."""
    heights = numpy.asarray(heights, float)
    b = numpy.asarray(b, float)
    ok = numpy.isfinite(b) & (b > 0)
    nan = (numpy.nan, numpy.nan, numpy.nan)
    if numpy.sum(ok) < 2:
        return nan
    # Fitted in nm and uN s/m (values of order 1)
    h, values = heights[ok] * 1e9, b[ok] * 1e6
    slope, intercept = numpy.polyfit(h, 1.0 / values, 1)
    if slope > 0 and intercept > 0:
        init = [1.0 / slope, intercept / slope]
    else:
        init = [values[0] * h[0], h[0]]
    try:
        params, _ = curve_fit(lambda x, a, he: a / (x + he), h, values,
                              p0=init, bounds=([0, 0], [numpy.inf,
                                                        numpy.inf]))
    except (RuntimeError, ValueError):
        return nan
    a, h_eff = params
    if not h_eff > 0:
        return nan
    residuals = values - a / (h + h_eff)
    ss_tot = numpy.sum((values - values.mean()) ** 2)
    r2 = 1.0 - numpy.sum(residuals ** 2) / ss_tot if ss_tot > 0 \
        else numpy.nan
    return a / h_eff * 1e-6, h_eff * 1e-9, r2


def model(heights, b0, h_eff):
    """b(h) of the fit [N s/m]."""
    return b0 * h_eff / (numpy.asarray(heights) + h_eff)
