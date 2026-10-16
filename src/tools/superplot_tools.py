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


"""SuperPlots (Lord et al., J Cell Biol 2020): the values of the maps
(one value per map: the median of its curves, or the mean for the number of
events per curve), colored by experiment, the means of the replicates, and
the statistics done on the replicates (n = number of replicates).

Hierarchy: condition > replicate (the groups) > map (file) > curves.
"""

import numpy
from scipy import stats

from .. import shared

DEFAULT_OPTIONS = {
    "control": 1,           # Condition used as control for the tests
    "paired": False,        # Replicates paired by experiment number
    "log": False,           # Log scale (plot and tests)
    "error": "SD",          # Error bars of the conditions: SD or 95% CI
    "curves": False,        # Violin of all the curves in the background
    "background": "Default",
}


def get_options():
    """Options of the SuperPlots (saved with the experiment)."""
    options = getattr(shared.exp, "superplot_options", None)
    if not isinstance(options, dict):
        options = {}
        shared.exp.superplot_options = options
    for key, value in DEFAULT_OPTIONS.items():
        options.setdefault(key, value)
    return options


def map_value(values):
    """Value of a map: median of its curves (mean for the number of events
    per curve, which is mostly 0 or 1)."""
    if shared.exp.results_type == "events_per_curve":
        return float(numpy.mean(values))
    return float(numpy.median(values))


def experiment_number(group, order):
    """Experiment of a replicate (pairing and colors): the one given in the
    replicates table, else its order in its condition."""
    number = getattr(group, "experiment", None)
    return number if number else order


def get_data():
    """Displayed conditions with their replicates and maps.

    Returns a list of conditions:
    {"id", "name", "replicates": [{"id", "name", "experiment", "value",
    "maps": [{"name", "value", "nbr_curves"}]}], "curves": all the values}
    The replicates without values are not listed.
    """
    from ..experiment import get_short_name

    if shared.single_data is None or \
            shared.exp.results_type == "loading_rates":
        return []

    conditions = []
    for condition in shared.exp.conditions_list[1:]:
        if not condition.display:
            continue
        replicates = []
        curves = []
        groups = [g for g in shared.exp.groups_list[1:]
                  if g.condition == condition.condition_id]
        for order, group in enumerate(groups, start=1):
            if not group.display:
                continue
            maps = []
            for i, result in enumerate(shared.exp.results_list):
                if result.group != group.group_id or \
                        i >= len(shared.single_data):
                    continue
                values = numpy.asarray(shared.single_data[i], dtype=float)
                values = values[numpy.isfinite(values)]
                if values.size == 0:
                    continue
                maps.append({"name": get_short_name(result),
                             "value": map_value(values),
                             "nbr_curves": int(values.size)})
                curves.append(values)
            if maps:
                replicates.append({
                    "id": group.group_id, "name": group.name,
                    "experiment": experiment_number(group, order),
                    "value": float(numpy.mean([m["value"] for m in maps])),
                    "maps": maps})
        if replicates:
            conditions.append({
                "id": condition.condition_id, "name": condition.name,
                "replicates": replicates,
                "curves": numpy.concatenate(curves)})
    return conditions


def summary(condition):
    """Mean, SD and 95% confidence interval of the replicate values."""
    values = numpy.array([r["value"] for r in condition["replicates"]])
    n = values.size
    mean = float(values.mean())
    sd = float(values.std(ddof=1)) if n > 1 else numpy.nan
    ci = float(stats.t.ppf(0.975, n - 1) * sd / numpy.sqrt(n)) \
        if n > 1 else numpy.nan
    return {"n": n, "mean": mean, "sd": sd, "ci": ci,
            "maps": sum(len(r["maps"]) for r in condition["replicates"]),
            "curves": sum(m["nbr_curves"] for r in condition["replicates"]
                          for m in r["maps"])}


def stars(p):
    """Significance of a p-value."""
    if p is None or not numpy.isfinite(p):
        return "–"
    for limit, text in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"),
                        (0.05, "*")):
        if p < limit:
            return text
    return "ns"


def compare(conditions, options):
    """Tests of each condition against the control, on the replicate values
    (log10 of the values with the log scale).

    Unpaired: Welch t-test for one comparison, Dunnett's test for several.
    Paired (replicates with the same experiment number): paired t-tests,
    Holm correction for several comparisons.
    Returns a list of {"name", "test", "n", "p", "note"}.
    """
    control = next((c for c in conditions if c["id"] == options["control"]),
                   None)
    if control is None:
        return []
    others = [c for c in conditions if c is not control]
    if not others:
        return []

    def values(condition):
        result = {r["experiment"]: r["value"] for r in condition["replicates"]}
        if options["log"]:
            result = {k: (numpy.log10(v) if v > 0 else numpy.nan)
                      for k, v in result.items()}
        return result

    ctrl = values(control)
    results = []
    if options["paired"]:
        pvalues = []
        for other in others:
            vals = values(other)
            common = sorted(k for k in ctrl if k in vals and
                            numpy.isfinite(ctrl[k]) and numpy.isfinite(vals[k]))
            entry = {"name": other["name"] + " vs " + control["name"],
                     "test": "Paired t-test", "n": len(common), "p": numpy.nan,
                     "note": ""}
            if len(common) >= 2:
                entry["p"] = float(stats.ttest_rel(
                    [vals[k] for k in common], [ctrl[k] for k in common]).pvalue)
            else:
                entry["note"] = "Less than 2 paired experiments"
            results.append(entry)
            pvalues.append(entry["p"])
        if len(results) > 1:
            # Holm correction
            finite = [i for i, p in enumerate(pvalues) if numpy.isfinite(p)]
            order = sorted(finite, key=lambda i: pvalues[i])
            m = len(order)
            running = 0.0
            for rank, i in enumerate(order):
                running = max(running, min(1.0, (m - rank) * pvalues[i]))
                results[i]["p"] = running
                results[i]["test"] = "Paired t-test, Holm"
        return results

    ctrl_values = [v for v in ctrl.values() if numpy.isfinite(v)]
    groups = []
    for other in others:
        vals = [v for v in values(other).values() if numpy.isfinite(v)]
        groups.append(vals)
        results.append({"name": other["name"] + " vs " + control["name"],
                        "test": "Welch t-test" if len(others) == 1
                        else "Dunnett", "n": len(vals), "p": numpy.nan,
                        "note": ""})
    valid = len(ctrl_values) >= 2 and all(len(g) >= 2 for g in groups)
    if not valid:
        for entry, vals in zip(results, groups):
            if len(vals) < 2 or len(ctrl_values) < 2:
                entry["note"] = "Less than 2 replicates"
        if len(others) > 1:
            return results
    if len(others) == 1:
        if valid:
            results[0]["p"] = float(stats.ttest_ind(
                groups[0], ctrl_values, equal_var=False).pvalue)
    else:
        res = stats.dunnett(*groups, control=ctrl_values)
        for entry, p in zip(results, res.pvalue):
            entry["p"] = float(p)
    return results


def export_rows(conditions):
    """Rows of the csv export: one per map, with its replicate and condition."""
    rows = []
    for condition in conditions:
        for rep in condition["replicates"]:
            for item in rep["maps"]:
                rows.append({"condition": condition["name"],
                             "replicate": rep["name"],
                             "experiment": rep["experiment"],
                             "map": item["name"],
                             "map_value": item["value"],
                             "nbr_curves": item["nbr_curves"],
                             "replicate_mean": rep["value"]})
    return rows
