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


"""Dynamic force spectroscopy: Bell-Evans fit of the rupture forces against
the loading rates.

The most probable rupture force depends on the loading rate r as
    F = kBT / x * ln(r * x / (kBT * koff))
so a linear fit F = a * ln(r) + b gives the distance to the barrier
x = kBT / a and the dissociation rate koff = exp(-b / a) / a.

Two ways to get the points of the fit:
- "medians": one point per file (or group) displayed in the plots, at the
  median loading rate and the median force of its events;
- "binned": the events of the displayed files (or groups) are put in bins of
  loading rate with the same number of events (quantiles of ln(r)); one
  point per bin, at its median loading rate and median force.
The logarithm needs positive loading rates: the events with r <= 0 are not
used (their number is given with the fit).
"""

import numpy
from .. import shared

# Boltzmann constant [J/K]
KB = 1.380649e-23

METHODS = {"medians": "Medians per file", "binned": "Binned events"}


def get_items(mode):
    """The files (mode "single") or groups (mode "groups") ticked in their
    table, with their forces [pN] and loading rates [pN/s]:
    list of (name, forces, loading rates, color, temperatures [K])."""
    from ..experiment import get_short_name

    if mode == "groups":
        all_data = shared.groups_data
        items = shared.exp.groups_list
        first = 1
    else:
        all_data = shared.single_data
        items = shared.exp.results_list
        first = 0

    result = []
    if all_data is None:
        return result

    for i in range(first, min(len(all_data), len(items))):
        item = items[i]
        data = all_data[i]
        if not item.display or len(data) != 2:
            continue
        forces = numpy.asarray(data[0], dtype=float).ravel()
        lrs = numpy.asarray(data[1], dtype=float).ravel()
        if forces.size == 0 or forces.size != lrs.size:
            continue

        if mode == "groups":
            name = item.name
            temps = [shared.exp.list[res.data_id].temperature
                     for res in shared.exp.results_list
                     if res.group == item.group_id]
        else:
            name = get_short_name(item)
            temps = [shared.exp.list[item.data_id].temperature]

        result.append((name, forces, lrs, item.color, temps))

    return result


def median_point(forces, lrs):
    """Median loading rate and force of the events, with their quartiles:
    (r, F, [r25, r75], [F25, F75])."""
    return (numpy.median(lrs), numpy.median(forces),
            numpy.percentile(lrs, [25, 75]),
            numpy.percentile(forces, [25, 75]))


def get_points(mode, method, nbr_bins=10):
    """Points of the fit [pN/s, pN], the number of events used and the number
    of events with a loading rate <= 0 (not used), the temperatures."""
    items = get_items(mode)
    temps = [t for item in items for t in item[4]]
    points_r, points_f = [], []
    used = 0
    negative = 0

    if method == "medians":
        for _, forces, lrs, _, _ in items:
            r, f, _, _ = median_point(forces, lrs)
            if r > 0:
                points_r.append(r)
                points_f.append(f)
                used += forces.size
            else:
                negative += forces.size

    elif method == "binned":
        forces = numpy.concatenate([item[1] for item in items]) \
            if items else numpy.array([])
        lrs = numpy.concatenate([item[2] for item in items]) \
            if items else numpy.array([])
        positive = lrs > 0
        negative = int(numpy.count_nonzero(~positive))
        forces, lrs = forces[positive], lrs[positive]
        used = forces.size
        if used >= nbr_bins:
            # Bins with the same number of events
            order = numpy.argsort(lrs)
            for indexes in numpy.array_split(order, nbr_bins):
                points_r.append(numpy.median(lrs[indexes]))
                points_f.append(numpy.median(forces[indexes]))

    return (numpy.array(points_r), numpy.array(points_f), used, negative,
            temps)


def fit(points_r, points_f, temperature):
    """Linear fit F = a ln(r) + b (SI units) and the Bell-Evans parameters
    with their standard errors (NaN with less than 3 points)."""
    ln_r = numpy.log(points_r)
    n = len(ln_r)
    a, b = numpy.polyfit(ln_r, points_f, 1)

    residuals = points_f - (a * ln_r + b)
    sxx = numpy.sum((ln_r - ln_r.mean()) ** 2)
    ss_tot = numpy.sum((points_f - points_f.mean()) ** 2)
    r2 = 1 - numpy.sum(residuals ** 2) / ss_tot if ss_tot > 0 else numpy.nan
    if n > 2:
        s2 = numpy.sum(residuals ** 2) / (n - 2)
        var_a = s2 / sxx
        var_b = s2 * (1.0 / n + ln_r.mean() ** 2 / sxx)
        cov_ab = -ln_r.mean() * s2 / sxx
    else:
        var_a = var_b = cov_ab = numpy.nan

    kbt = KB * temperature
    x = kbt / a
    koff = numpy.exp(-b / a) / a
    sx = kbt * numpy.sqrt(var_a) / a ** 2
    # ln(koff) = -b/a - ln(a)
    d_a = b / a ** 2 - 1 / a
    d_b = -1 / a
    var_ln = d_a ** 2 * var_a + d_b ** 2 * var_b + 2 * d_a * d_b * cov_ab
    skoff = koff * numpy.sqrt(var_ln) if var_ln >= 0 else numpy.nan
    if a <= 0:
        # The force does not increase with the loading rate: no Bell-Evans
        # parameters
        x = sx = koff = skoff = numpy.nan

    return {"a": a, "b": b, "sa": numpy.sqrt(var_a), "sb": numpy.sqrt(var_b),
            "x": x, "sx": sx, "koff": koff, "skoff": skoff, "r2": r2,
            "temperature": temperature}


def new_fit(mode, method, nbr_bins=10):
    """Computes a fit with the displayed files (or groups). Returns the fit,
    or an error message."""
    points_r, points_f, used, negative, temps = get_points(
        mode, method, nbr_bins)
    if len(points_r) < 2:
        if method == "medians":
            what = "speeds" if mode == "groups" else "files"
            return ("At least 2 " + what + " with a positive median "
                    "loading rate must be displayed (ticked) for the fit.")
        return ("Not enough events with a positive loading rate for " +
                str(nbr_bins) + " bins.")
    if numpy.ptp(numpy.log(points_r)) == 0:
        return "The points of the fit have the same loading rate."

    # Values in N and N/s for the fit
    result = fit(points_r * 1e-12, points_f * 1e-12, float(numpy.mean(temps)))
    result.update({
        "method": method, "nbr_bins": nbr_bins, "nbr_points": len(points_r),
        "nbr_events": used, "nbr_negative": negative,
        "temperatures": [min(temps), max(temps)],
        "points_r": points_r.tolist(), "points_f": points_f.tolist(),
        "display": True})
    return result


def get_fits(mode):
    """The fits of the plots tab (mode "single") or of the grouped plots tab
    (mode "groups"). The fits of the old files (lists) are converted: their
    line was not the Bell-Evans fit and is not drawn."""
    exp = shared.exp
    name = "list_dfs_" + mode + "_fits"
    fits = getattr(exp, name)
    olds = getattr(exp, "list_dfs_" + mode + "_results", [])
    for i, item in enumerate(fits):
        if isinstance(item, dict):
            continue
        x, koff = olds[i] if i < len(olds) else (numpy.nan, numpy.nan)
        fits[i] = {"method": "old", "x": x, "koff": koff, "display": False,
                   "nbr_points": numpy.nan, "r2": numpy.nan,
                   "sx": numpy.nan, "skoff": numpy.nan}
    return fits
