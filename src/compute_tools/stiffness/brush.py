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

"""Brush model of the cells (Sokolov et al., Appl. Phys. Lett. 2007; Iyer
et al., Nat. Nanotechnol. 2009): the cell body is indented through a layer
of pericellular brush (glycocalyx, microvilli, microridges), with a
spherical probe.

1. The brush is squeezed near the maximum load: the cell body is fitted
   above a threshold, a percentage of the maximum force of the curve
   (paraboloid or sphere), giving its modulus E and the position of its
   surface.
2. Below the threshold, the distance h between the probe and the cell body
   is the indentation of the surface of the body for the force (from the
   fit) minus the indentation reached; the force of the brush is fitted
   with F = 100 kT R N^(3/2) L exp(-2 pi h / L) for 0.1 < h / L < 0.8
   (Makarova et al., ACS Appl. Mater. Interfaces 2023; linear fit of ln F
   vs h, iterated as L is fitted): length L of the brush and grafting
   density N (effective number of chains per unit area).

The curves are the force - indentation curves of the elasticity fit
(point of contact and baseline of its search, see bec.approach_points),
with the zero of the force on the baseline line (the point of contact is
detected within the brush). The refit of the baseline of the point of
contact moves it into the brush: to be turned off for the brush model.
The brush points are the contiguous ones before the threshold, above
FLOOR_SIGMA x noise.
"""

import os

import numpy
from scipy import optimize

from ... import shared
from ... import widgets_list
from ...tools import misc_tools
from ...widgets.progressbar import Progressbar
from . import bec
from . import contact_models

# Models of the cell body (the brush force is the one of a sphere)
BRUSH_MODELS = (contact_models.PARABOLOID, contact_models.SPHERE)

# Values of a curve in the results: modulus of the body [Pa], length of the
# brush [nm], grafting density [1/um^2], R^2 of the brush fit, surface of
# the body [nm, from the point of contact]
VALUES = ("modulus", "length", "density", "r2", "surface")

# Thermal energy (25 degrees C) [J]
KT = 1.380649e-23 * 298.15

# Lowest force of the brush points: FLOOR_SIGMA x noise of the baseline
FLOOR_SIGMA = 5.0

# Prefactor of the brush force F = PREFACTOR kT R N^(3/2) L exp(-2 pi h /
# L) (Makarova et al., ACS Appl. Mater. Interfaces 2023, eq. 6; 50 in Guz
# et al., Biophys. J. 2014), valid for H_MIN < h / L < H_MAX
PREFACTOR = 100.0
H_MIN = 0.1
H_MAX = 0.8

# Iterations of the fit of the brush in its window of validity
MAX_ITERATIONS = 20

# Minimum number of points of the fits
MIN_POINTS = 10

# Minimum R^2 of the brush fit (ln F vs h): below, no exponential decay of
# the force (no brush, or brush and body not separated), no length and
# density
MIN_R2 = 0.8


def available(data, used=False):
    """True if the brush model can be used with the model of the
    elasticity."""
    model = data.used_stiffness_model_selected if used else \
        data.stiffness_model_selected
    return model in BRUSH_MODELS


def body_unit_force(model, radius, poisson_ratio, delta):
    """Force per unit modulus of the body [N/Pa], delta [m]."""
    return contact_models.unit_force(model, radius, 0.0, poisson_ratio,
                                     delta)


def fit_body(model, radius, poisson_ratio, delta, force):
    """Modulus [Pa] and surface [m] of the body, fitted on the points above
    the threshold (paraboloid: linear fit of F^(2/3); sphere: least squares
    from the paraboloid)."""
    coeff = 4.0 / 3.0 * numpy.sqrt(radius) / (1.0 - poisson_ratio ** 2)
    slope, intercept = numpy.polyfit(delta, force ** (2.0 / 3.0), 1)
    if not slope > 0:
        return numpy.nan, numpy.nan
    surface = -intercept / slope
    modulus = slope ** 1.5 / coeff
    if model == contact_models.PARABOLOID:
        return modulus, surface

    scale = numpy.max(force)

    def residuals(p):
        return (p[1] * 1e3 * body_unit_force(
            model, radius, poisson_ratio, delta - p[0]) - force) / scale

    try:
        fit = optimize.least_squares(residuals, [surface, modulus / 1e3],
                                     x_scale=[1e-7, modulus / 1e3])
    except ValueError:
        return modulus, surface
    if fit.success and fit.x[1] > 0:
        return fit.x[1] * 1e3, fit.x[0]
    return modulus, surface


def body_indentation(model, radius, poisson_ratio, modulus, force):
    """Indentation of the body [m] for the forces [N]."""
    force = numpy.clip(force, 0.0, None)
    if model == contact_models.PARABOLOID:
        coeff = 4.0 / 3.0 * numpy.sqrt(radius) / (1.0 - poisson_ratio ** 2)
        return (force / (coeff * modulus)) ** (2.0 / 3.0)
    # Inverse of the tabulated force (increasing)
    top = 2.0 * radius
    grid = numpy.r_[0.0, numpy.logspace(numpy.log10(top) - 6,
                                        numpy.log10(top), 400)]
    table = modulus * body_unit_force(model, radius, poisson_ratio, grid)
    return numpy.interp(force, table, grid)


def fit_brush(h, force, brush):
    """Linear fit of ln F vs h on the brush points, restricted to the
    validity of the brush force, H_MIN < h / L < H_MAX (iterated, L being
    fitted). Returns (slope, intercept, indices of the fitted points), or
    None."""
    points = brush
    fit = None
    for _ in range(MAX_ITERATIONS):
        if points.size < MIN_POINTS:
            return None
        slope, intercept = numpy.polyfit(h[points],
                                         numpy.log(force[points]), 1)
        if not slope < 0:
            return None
        fit = (slope, intercept, points)
        length = -2.0 * numpy.pi / slope
        window = brush[(h[brush] > H_MIN * length) &
                       (h[brush] < H_MAX * length)]
        if numpy.array_equal(window, points):
            break
        points = window
    return fit


def analyse(delta, force, noise, model, radius, poisson_ratio, fraction,
            full=False):
    """Brush model of a curve: delta [m] and force [N] from the point of
    contact, noise of the force [N], radius [m], fraction of the maximum
    force of the curve above which the cell body is fitted.

    Returns a dict with VALUES (nan when not fitted) and the force of the
    threshold [N], None if the body can not be fitted. With full, also the
    curve ("delta", "force"), the distance to the body h [m], the indices
    of the brush points before the threshold ("candidates") and of the
    fitted ones, within 0.1 < h / L < 0.8 ("brush"), and the fit ("slope",
    "intercept" of ln F vs h).
    """
    ok = numpy.isfinite(delta) & numpy.isfinite(force)
    delta, force = delta[ok], force[ok]
    if force.size < MIN_POINTS:
        return None
    threshold = fraction * numpy.max(force)
    below = numpy.flatnonzero(force <= threshold)
    start = below[-1] + 1 if below.size else 0
    if len(force) - start < MIN_POINTS:
        return None
    modulus, surface = fit_body(model, radius, poisson_ratio,
                                delta[start:], force[start:])
    if not (numpy.isfinite(modulus) and modulus > 0):
        return None
    result = {"modulus": modulus, "surface": surface * 1e9,
              "length": numpy.nan, "density": numpy.nan, "r2": numpy.nan,
              "threshold": threshold}

    # Distance between the probe and the body
    h = surface + body_indentation(model, radius, poisson_ratio, modulus,
                                   force) - delta
    # Brush: contiguous points before the threshold, above the noise
    floor = FLOOR_SIGMA * noise
    low = numpy.flatnonzero(force[:start] <= floor)
    first = low[-1] + 1 if low.size else 0
    brush = numpy.arange(first, start)
    brush = brush[(force[brush] > floor) & (h[brush] > 0)]
    fit = fit_brush(h, force, brush)
    if full:
        result.update({"delta": delta, "force": force, "h": h,
                       "candidates": brush,
                       "brush": fit[2] if fit else brush[:0]})
    if fit is None:
        return result

    slope, intercept, points = fit
    length = -2.0 * numpy.pi / slope
    amplitude = numpy.exp(intercept)
    density = (amplitude / (PREFACTOR * KT * radius * length)) ** \
        (2.0 / 3.0)
    log_force = numpy.log(force[points])
    total = numpy.sum((log_force - log_force.mean()) ** 2)
    r2 = 1.0 - numpy.sum((log_force - (slope * h[points] + intercept)) **
                         2) / total if total > 0 else numpy.nan
    result["r2"] = r2
    if full:
        result.update({"slope": slope, "intercept": intercept})
    if not r2 >= MIN_R2:
        return result
    result.update({"length": length * 1e9, "density": density * 1e-12})
    return result


def settings(data, used=True):
    """Model, radius [m], Poisson ratio and threshold (fraction of the
    maximum force)."""
    if used:
        return (data.used_stiffness_model_selected,
                data.used_tip_radius * 1e-9, data.used_poisson_ratio,
                data.used_brush_threshold_pct / 100.0)
    return (data.stiffness_model_selected, data.tip_radius * 1e-9,
            data.poisson_ratio, data.brush_threshold_pct / 100.0)


def curve(data, i, j, full=False):
    """Brush model of the curve of a pixel (used settings), None without
    point of contact or body fit."""
    points = bec.approach_points(data, i, j)
    if points is None:
        return None
    delta, force, z = points
    # Zero of the force on the baseline line of the fit of the point of
    # contact (approach_points), not at the detected point of contact: it
    # is above the noise band of the baseline, within the brush (its weak
    # force would be lost)
    # Noise of the force: window of the baseline fit of the point of
    # contact (as bec.fit_points)
    skip = data.used_fitparam_poc_skip_start or 0
    length = data.used_fitparam_poc_fit_length or 0
    window = (delta <= 0) & (z >= z[0] + skip) & \
        (z <= z[0] + skip + length)
    baseline = force[window] if numpy.count_nonzero(window) > 2 else \
        force[delta <= 0]
    noise = 1.4826 * numpy.median(numpy.abs(baseline -
                                            numpy.median(baseline))) \
        if baseline.size > 2 else 0.0
    model, radius, poisson_ratio, fraction = settings(data)
    return analyse(delta, force, noise, model, radius, poisson_ratio,
                   fraction, full=full)


def compute(calc_id):
    """Brush model of the curves of a file, after its elasticity (points of
    contact). Returns the number of curves with a brush."""
    data = shared.exp.list[calc_id]
    misc_tools.save_used_params(data, "brush")
    nx, ny = data.nbr_pixels_x, data.nbr_pixels_y
    results = numpy.full((nx, ny, len(VALUES)), numpy.nan,
                         dtype=numpy.float32)

    Progressbar(os.path.basename(data.filename))
    widgets_list.widget_progressbar.set_label("Brush model")
    widgets_list.widget_progressbar.set_range(0, nx * ny)

    found = 0
    for i in range(nx):
        for j in range(ny):
            if data.discarded_curves[i][j] == 0:
                result = curve(data, i, j)
                if result is not None:
                    results[i, j] = [result[name] for name in VALUES]
                    found += bool(numpy.isfinite(result["length"]))
            widgets_list.widget_progressbar.update()

    data.brush_results = results
    data.brush_calculated = True
    data.update()
    widgets_list.widget_progressbar.close()
    return found
