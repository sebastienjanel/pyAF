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


"""Dynamic force spectroscopy: models of the rupture forces against the
loading rates (see Mesbah et al., in Atomic Force Microscopy: Methods and
Protocols, 2nd ed., chap. 5).

Bell-Evans (Evans and Ritchie 1997): the most probable rupture force
depends on the loading rate r as
    F* = kBT / x * ln(r * x / (kBT * koff))
so a linear fit F = a * ln(r) + b gives the distance to the barrier
x = kBT / a and the dissociation rate koff = exp(-b / a) / a. The points
are the median forces, whose exact value for the Bell-Evans distribution of
the rupture forces is
    F_median = kBT / x * ln(1 + ln 2 * r * x / (kBT * koff))
x and koff are fitted with this expression (started from the linear fit):
at high forces it is the line with koff multiplied by ln 2, and it stays
right at low loading rates, where the line is biased.

Dudko-Hummer-Szabo (DHS, Phys Rev Lett 2006): the rate of rupture
    k(F) = koff (1 - nu F x / dG) ** (1 / nu - 1)
           * exp(dG / kBT * (1 - (1 - nu F x / dG) ** (1 / nu)))
with the height of the barrier dG, nu = 1/2 (cusp barrier) or 2/3
(linear-cubic). Curved spectrum (x decreases with the force); Bell-Evans
for nu = 1. Needs a wide range of loading rates. The median of its rupture
forces is exact:
    F_median = dG / (nu x) * (1 - [1 - kBT / dG
                                   * ln(1 + ln 2 * r x / (kBT koff))] ** nu)
(the mean of the book, Eq. 1.5, is an approximation, with exp(-g) instead
of ln 2, g: Euler's constant).

Friddle-Noy-De Yoreo (PNAS 2012): rebinding near equilibrium,
    F = feq + kBT / x * ln(1 + c * r x / (k(feq) kBT))
with k(feq) = koff exp(feq x / kBT): plateau at the equilibrium force feq
at low loading rates, Bell-Evans at high rates. c = exp(-g) for the mean
force (their approximation), ln 2 here for the median (exact in the
Bell-Evans limit). feq = sqrt(2 kc dGbu), kc: effective spring constant
(cantilever and linker), dGbu: free energy of the bond.

All the models are fitted on the median forces: the mean is pulled up by
the large forces of the multiple bonds and the nonspecific adhesions, and
all the models are then fitted on the same points.

Two ways to get the points of the fit:
- "medians": one point per file (or group) displayed in the plots, at the
  median loading rate and the median (or mean) force of its events;
- "binned": the events of the displayed files (or groups) are put in bins of
  loading rate with the same number of events (quantiles of ln(r)); one
  point per bin, at its median loading rate and median (or mean) force.
The logarithm needs positive loading rates: the events with r <= 0 are not
used (their number is given with the fit).
"""

import warnings

import numpy
from scipy import optimize
from .. import shared


KB = 1.380649e-23

METHODS = {"medians": "Medians per file", "binned": "Binned events"}

# Models: name, statistic of the forces of the points
MODELS = {"bell": ("Bell-Evans", "median"),
          "dhs_cusp": ("DHS (\u03bd = 1/2)", "median"),
          "dhs_cubic": ("DHS (\u03bd = 2/3)", "median"),
          "fnd": ("Friddle-Noy-De Yoreo", "median")}
NU = {"dhs_cusp": 0.5, "dhs_cubic": 2.0 / 3.0}


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


def get_points(mode, method, nbr_bins=10, statistic="median"):
    """Points of the fit [pN/s, pN] (median loading rate, median or mean
    force), the number of events used and the number of events with a
    loading rate <= 0 (not used), the temperatures, and the quartiles of
    the loading rates and forces of each point (error bars)."""
    average = numpy.mean if statistic == "mean" else numpy.median
    items = get_items(mode)
    temps = [t for item in items for t in item[4]]
    points_r, points_f = [], []
    quartiles_r, quartiles_f = [], []
    used = 0
    negative = 0

    if method == "medians":
        for _, forces, lrs, _, _ in items:
            r = numpy.median(lrs)
            if r > 0:
                points_r.append(r)
                points_f.append(average(forces))
                quartiles_r.append(numpy.percentile(lrs, [25, 75]))
                quartiles_f.append(numpy.percentile(forces, [25, 75]))
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
                points_f.append(average(forces[indexes]))
                quartiles_r.append(numpy.percentile(lrs[indexes], [25, 75]))
                quartiles_f.append(
                    numpy.percentile(forces[indexes], [25, 75]))

    return (numpy.array(points_r), numpy.array(points_f), used, negative,
            temps, numpy.array(quartiles_r), numpy.array(quartiles_f))


def fit(points_r, points_f, temperature):
    """Bell-Evans fit of the median forces (SI units): the linear fit
    F = a ln(r) + b, then x and koff with their standard errors from the
    exact median (see the top; NaN errors with less than 3 points)."""
    ln_r = numpy.log(points_r)
    n = len(ln_r)
    a, b = numpy.polyfit(ln_r, points_f, 1)

    residuals = points_f - (a * ln_r + b)
    sxx = numpy.sum((ln_r - ln_r.mean()) ** 2)
    ss_tot = numpy.sum((points_f - points_f.mean()) ** 2)
    if n > 2:
        s2 = numpy.sum(residuals ** 2) / (n - 2)
        var_a = s2 / sxx
        var_b = s2 * (1.0 / n + ln_r.mean() ** 2 / sxx)
    else:
        var_a = var_b = numpy.nan

    kbt = KB * temperature
    result = {"a": a, "b": b, "sa": numpy.sqrt(var_a),
              "sb": numpy.sqrt(var_b), "x": numpy.nan, "sx": numpy.nan,
              "koff": numpy.nan, "skoff": numpy.nan, "r2": numpy.nan,
              "temperature": temperature, "median_exact": True}
    if a <= 0:
        # The force does not increase with the loading rate: no Bell-Evans
        # parameters
        return result

    # Exact median, in nm and log10(koff), from the line (koff x ln 2)
    def function(r, x_nm, log_koff):
        return median_force(r, x_nm * 1e-9, 10 ** log_koff, kbt)

    start = [kbt / a * 1e9,
             numpy.log10(numpy.log(2) * numpy.exp(-b / a) / a)]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            params, cov = optimize.curve_fit(function, points_r, points_f,
                                             p0=start, maxfev=20000)
    except (RuntimeError, ValueError):
        params, cov = numpy.array(start), numpy.full((2, 2), numpy.nan)
    if params[0] <= 0:
        return result
    errors = numpy.sqrt(numpy.diag(cov)) if n > 2 and \
        numpy.all(numpy.isfinite(cov)) else numpy.full(2, numpy.nan)
    koff = 10 ** params[1]
    sse = numpy.sum((points_f - function(points_r, *params)) ** 2)
    result.update({
        "x": params[0] * 1e-9, "sx": errors[0] * 1e-9, "koff": koff,
        "skoff": koff * numpy.log(10) * errors[1],
        "r2": 1 - sse / ss_tot if ss_tot > 0 else numpy.nan})
    return result


def median_force(r, x, koff, kbt):
    """Median rupture force of the Bell-Evans distribution (SI units)."""
    return kbt / x * numpy.log1p(numpy.log(2) * r * x / (kbt * koff))


def new_fit(mode, method, nbr_bins=10, model="bell"):
    """Computes a fit with the displayed files (or groups). Returns the fit,
    or an error message."""
    statistic = MODELS[model][1]
    points_r, points_f, used, negative, temps, quartiles_r, quartiles_f = \
        get_points(mode, method, nbr_bins, statistic)
    needed = 2 if model == "bell" else 4
    if len(points_r) < needed:
        what = ("speeds" if mode == "groups" else "files") \
            if method == "medians" else "bins"
        text = "At least " + str(needed) + " points are needed"
        if model != "bell":
            text += " (3 parameters)"
        if method == "medians":
            return (text + ": display (tick) " + str(needed) + " " + what +
                    " with a positive median loading rate.")
        return (text + ": not enough events with a positive loading rate "
                "for " + str(nbr_bins) + " bins.")
    if numpy.ptp(numpy.log(points_r)) == 0:
        return "The points of the fit have the same loading rate."

    temperature = float(numpy.mean(temps))
    if model == "bell":
        # Values in N and N/s for the fit
        result = fit(points_r * 1e-12, points_f * 1e-12, temperature)
    else:
        result = fit_model(model, points_r, points_f, temperature)
        if isinstance(result, str):
            return result
        if model == "fnd":
            # Free energy of the bond with the spring constant of the
            # cantilever (softer with the linker: larger dGbu)
            kc = numpy.mean(spring_constants(mode)) * 1e9  # N/m
            result["kc"] = kc
            result["dGbu"] = (result["feq"] * 1e-12) ** 2 / (2 * kc) / \
                (KB * temperature)
    result.update({
        "model": model, "statistic": statistic,
        "method": method, "nbr_bins": nbr_bins, "nbr_points": len(points_r),
        "nbr_events": used, "nbr_negative": negative,
        "temperatures": [min(temps), max(temps)],
        "points_r": points_r.tolist(), "points_f": points_f.tolist(),
        "points_r_q": quartiles_r.tolist(), "points_f_q": quartiles_f.tolist(),
        "display": True})
    return result


def spring_constants(mode):
    """Spring constants [N/nm] of the files of the displayed items."""
    if mode == "groups":
        ids = {g.group_id for g in shared.exp.groups_list[1:] if g.display}
        results = [r for r in shared.exp.results_list if r.group in ids]
    else:
        results = [r for r in shared.exp.results_list if r.display]
    values = [shared.exp.list[r.data_id].spring_constant for r in results]
    return values or [numpy.nan]


def dhs_force(r, x, koff, dg, nu, kbt):
    """Median rupture force of the DHS model (SI units), dg in J."""
    with numpy.errstate(divide="ignore", invalid="ignore"):
        inner = 1.0 - kbt / dg * numpy.log1p(
            numpy.log(2) * r * x / (kbt * koff))
    inner = numpy.clip(inner, 0.0, None)
    return dg / (nu * x) * (1.0 - inner ** nu)


def fnd_force(r, x, koff, feq, kbt):
    """Median rupture force of the Friddle-Noy-De Yoreo model (SI
    units)."""
    k_eq = koff * numpy.exp(feq * x / kbt)
    return feq + kbt / x * numpy.log1p(r * x * numpy.log(2) /
                                       (k_eq * kbt))


def model_force(fit_result, r):
    """Force [pN] of a fit at the loading rates r [pN/s]."""
    model = fit_result.get("model", "bell")
    r = numpy.asarray(r, dtype=float)
    kbt = KB * fit_result.get("temperature", 298.15)
    if model == "bell":
        if fit_result.get("median_exact") and numpy.isfinite(fit_result["x"]):
            return median_force(r * 1e-12, fit_result["x"],
                                fit_result["koff"], kbt) * 1e12
        # Line of the older fits
        return (fit_result["a"] * numpy.log(r * 1e-12) +
                fit_result["b"]) * 1e12
    if model in NU:
        force = dhs_force(r * 1e-12, fit_result["x"], fit_result["koff"],
                          fit_result["dG"] * kbt, NU[model], kbt)
    else:
        force = fnd_force(r * 1e-12, fit_result["x"], fit_result["koff"],
                          fit_result["feq"] * 1e-12, kbt)
    return force * 1e12


def fit_model(model, points_r, points_f, temperature):
    """Least-squares fit of DHS or Friddle-Noy-De Yoreo to the median
    points [pN/s, pN]. Parameters x [nm], log10 koff and dG [kBT] (DHS) or
    feq [pN] (Friddle-Noy-De Yoreo), started from the Bell-Evans fit and
    from several values of the third parameter. Returns the parameters with
    their standard errors, or an error message."""
    kbt = KB * temperature
    r_si = points_r * 1e-12

    if model in NU:
        nu = NU[model]

        def function(r, x_nm, log_koff, dg_kt):
            return dhs_force(r * 1e-12, x_nm * 1e-9, 10 ** log_koff,
                             dg_kt * kbt, nu, kbt) * 1e12
        thirds = (5.0, 10.0, 20.0, 40.0, 80.0)
        bounds = ([1e-3, -15.0, 1.0], [100.0, 8.0, 500.0])
    else:
        def function(r, x_nm, log_koff, feq):
            return fnd_force(r * 1e-12, x_nm * 1e-9, 10 ** log_koff,
                             feq * 1e-12, kbt) * 1e12
        low = max(float(numpy.min(points_f)), 1e-3)
        thirds = (0.2 * low, 0.5 * low, 0.8 * low, 0.95 * low)
        bounds = ([1e-3, -15.0, 0.0], [100.0, 8.0,
                                       float(numpy.max(points_f))])

    # Start: Bell-Evans fit of the median forces
    bell = fit(r_si, points_f * 1e-12, temperature)
    if numpy.isfinite(bell["x"]) and numpy.isfinite(bell["koff"]):
        x0 = bell["x"] * 1e9
        log_koff0 = numpy.log10(bell["koff"])
    else:
        x0, log_koff0 = 0.5, 0.0
    x0 = float(numpy.clip(x0, 2e-3, 50.0))
    log_koff0 = float(numpy.clip(log_koff0, -14.0, 7.0))

    best = None
    for third in thirds:
        start = [x0, log_koff0, float(numpy.clip(
            third, bounds[0][2] + 1e-6, bounds[1][2] - 1e-6))]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                params, cov = optimize.curve_fit(
                    function, points_r, points_f, p0=start, bounds=bounds,
                    maxfev=20000)
        except (RuntimeError, ValueError):
            continue
        sse = float(numpy.sum((points_f - function(points_r, *params)) ** 2))
        if best is None or sse < best[0]:
            best = (sse, params, cov)
    if best is None:
        return "The fit of the model did not converge."

    sse, params, cov = best
    errors = numpy.sqrt(numpy.diag(cov)) if numpy.all(numpy.isfinite(cov)) \
        else numpy.full(3, numpy.nan)
    if len(points_r) <= 3:
        errors = numpy.full(3, numpy.nan)
    ss_tot = float(numpy.sum((points_f - points_f.mean()) ** 2))
    koff = 10 ** params[1]
    result = {"x": params[0] * 1e-9, "sx": errors[0] * 1e-9,
              "koff": koff, "skoff": koff * numpy.log(10) * errors[1],
              "r2": 1 - sse / ss_tot if ss_tot > 0 else numpy.nan,
              "temperature": temperature}
    at_bound = [name for name, value, lo, hi in zip(
        ("x", "koff", "third"), params, bounds[0], bounds[1])
        if numpy.isclose(value, lo, rtol=1e-3, atol=1e-9) or
        numpy.isclose(value, hi, rtol=1e-3, atol=1e-9)]
    # Limits where the model is Bell-Evans: no curvature (DHS, dG at its
    # maximum) or no plateau (Friddle-Noy-De Yoreo, feq at 0)
    if model in NU:
        result.update({"dG": params[2], "sdG": errors[2], "nu": NU[model]})
        names = {"third": "\u0394G\u2021"}
        bell_limit = params[2] >= 0.99 * bounds[1][2]
    else:
        result.update({"feq": params[2], "sfeq": errors[2]})
        names = {"third": "feq"}
        bell_limit = params[2] <= 0.01 * float(numpy.median(points_f))
    if bell_limit:
        at_bound = [n for n in at_bound if n != "third"]
    result["bell_limit"] = bool(bell_limit)
    result["at_bound"] = [names.get(n, n) for n in at_bound]
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
            if "model" not in item and item.get("method") != "old":
                # Bell-Evans fit of an older version: koff of the median
                # forces, corrected (x ln 2, see the top)
                item["model"] = "bell"
                item["koff"] *= numpy.log(2)
                item["skoff"] *= numpy.log(2)
                item["koff_corrected"] = True
            continue
        x, koff = olds[i] if i < len(olds) else (numpy.nan, numpy.nan)
        fits[i] = {"method": "old", "x": x, "koff": koff, "display": False,
                   "nbr_points": numpy.nan, "r2": numpy.nan,
                   "sx": numpy.nan, "skoff": numpy.nan}
    return fits
