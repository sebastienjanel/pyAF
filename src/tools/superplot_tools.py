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


"""SuperPlots (Lord et al., J Cell Biol 2020): the values of the files
(force maps, QI images or single curves; one value per file: the median of
its curves, or the mean for the number of
events per curve), colored by experiment, the means of the replicates, and
the statistics done on the replicates (n = number of replicates).

Hierarchy: condition > replicate (the groups) > file > curves (the "maps"
of the data: one per file, a single curve being a file of one curve).
"""

import numpy
from scipy import stats

from .. import shared

# Control option: all the pairs of conditions are compared (no control; the
# ids of the conditions start at 1, 0 is "no condition")
ALL_PAIRS = -1

DEFAULT_OPTIONS = {
    "control": 1,           # Condition used as control for the tests
                            # (ALL_PAIRS: all pairwise comparisons)
    "paired": False,        # Replicates paired by experiment number
    "log": False,           # Log scale (plot and tests)
    "error": "SD",          # Error bars of the conditions: SD or 95% CI
    "curves": False,        # Violin of all the curves in the background
    "background": "Default",
}


# Tooltips of the columns of the Plots and Replicates tabs (hierarchy of
# the SuperPlots: condition > replicate > file > curves)
REPLICATE_TIP = (
    "Replicate: independent repetition of a condition (e.g. one dish or "
    "one gel, prepared and measured separately),\ngroup of files. It is "
    "the statistical unit of the SuperPlots: n = number of replicates.\n"
    "Several files of the same dish are not independent: put them in the "
    "same replicate.")
EXPERIMENT_TIP = (
    "Experiment: replicates of different conditions done together (same "
    "day, cells, tip, calibration)\nhave the same number: colors of the "
    "SuperPlots and pairing of the paired tests.\nEmpty: order of the "
    "replicate in its condition.")
CONDITION_TIP = (
    "Condition: treatment or sample compared in the SuperPlots (control, "
    "drug...), set of replicates.")
MODE_TIP = (
    "Mode: most frequent value, maximum of the density of the values "
    "(kernel density estimate,\nsame smoothing as the violins).")
FILE_TIP = (
    "File: force map, QI image or single curve. Its value is the median "
    "of its curves (mean for the number of events per curve).")
COLUMN_TOOLTIPS = {
    "plots": {
        "Replicate": "Replicate of the file. " + REPLICATE_TIP,
        "Mode": MODE_TIP + "\nOf the curves of the file."},
    "replicates": {
        "Replicate": REPLICATE_TIP,
        "Files": "Number of files in the replicate (column Replicate of "
                 "the Plots tab)",
        "Experiment": EXPERIMENT_TIP,
        "Condition": CONDITION_TIP,
        "Mode": MODE_TIP + "\nOf all the curves of the replicate."},
    "superplots": {
        "Condition": CONDITION_TIP,
        "Replicates": "Number of replicates of the condition: n of the "
                      "tests. " + REPLICATE_TIP,
        "Files": "Number of files of the condition (small dots). " +
                 FILE_TIP,
        "Curves": "Number of curves of the condition",
        "95% CI": "95 % confidence interval of the mean of the replicate "
                  "values (half-width)"},
}


def set_header_tooltips(table, tab):
    """Tooltips of the headers of a table of the Plots ("plots"),
    Replicates ("replicates") or SuperPlots ("superplots") tab, from the
    texts of the headers."""
    tips = COLUMN_TOOLTIPS[tab]
    for column in range(table.columnCount()):
        item = table.horizontalHeaderItem(column)
        if item is not None and item.text() in tips:
            item.setToolTip(tips[item.text()])


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


def fieller(m1, m0, v11, v00, v10, crit):
    """Confidence interval of the ratio m1 / m0 of two means (Fieller):
    v11, v00 variances of the means, v10 their covariance, crit the
    critical value of the test (the interval excludes 1 exactly when the
    test of the difference is significant). nan if unbounded (the control
    mean is not significantly different from 0)."""
    a = m0 ** 2 - crit ** 2 * v00
    b = -2.0 * (m1 * m0 - crit ** 2 * v10)
    c = m1 ** 2 - crit ** 2 * v11
    delta = b ** 2 - 4.0 * a * c
    if not a > 0 or not delta >= 0:
        return numpy.nan, numpy.nan
    root = numpy.sqrt(delta)
    return float((-b - root) / (2 * a)), float((-b + root) / (2 * a))


def fold_change(a, b, v11, v00, v10, crit, log):
    """Fold change of the values a vs the control values b, and its
    confidence interval, from the variances of the means and the critical
    value of the test. With the log scale (log10 values): ratio of the
    geometric means. Else ratio of the means (Fieller interval)."""
    if log:
        d = float(numpy.mean(a) - numpy.mean(b))
        se = numpy.sqrt(v11 + v00 - 2 * v10)
        return 10 ** d, 10 ** (d - crit * se), 10 ** (d + crit * se)
    low, high = fieller(float(numpy.mean(a)), float(numpy.mean(b)), v11, v00,
                        v10, crit)
    return float(numpy.mean(a) / numpy.mean(b)), low, high


def compare(conditions, options):
    """Tests of each condition against the control, on the replicate values
    (log10 of the values with the log scale), with the fold change vs the
    control and its 95 % confidence interval (same test: the interval
    excludes 1 when the test is significant).

    Unpaired: Welch t-test for one comparison, Dunnett's test (simultaneous
    intervals) for several.
    Paired (replicates with the same experiment number): paired t-tests,
    Holm correction for several comparisons (Bonferroni intervals, Holm
    gives no interval).
    Fold change: ratio of the means, ratio of the geometric means with the
    log scale.
    Returns a list of {"name", "test", "n", "p", "note", "fc", "fc_low",
    "fc_high", "fc_note"}.
    """
    if options["control"] == ALL_PAIRS:
        return compare_all_pairs(conditions, options)
    control = next((c for c in conditions if c["id"] == options["control"]),
                   None)
    if control is None:
        return []
    others = [c for c in conditions if c is not control]
    if not others:
        return []
    log = options["log"]

    def values(condition):
        result = {r["experiment"]: r["value"] for r in condition["replicates"]}
        if log:
            result = {k: (numpy.log10(v) if v > 0 else numpy.nan)
                      for k, v in result.items()}
        return result

    def entry(other, test, n):
        return {"name": other["name"] + " vs " + control["name"],
                "test": test, "n": n, "p": numpy.nan, "note": "",
                "fc": numpy.nan, "fc_low": numpy.nan, "fc_high": numpy.nan,
                "fc_note": ""}

    ctrl = values(control)
    results = []
    if options["paired"]:
        pvalues = []
        m = len(others)
        for other in others:
            vals = values(other)
            common = sorted(k for k in ctrl if k in vals and
                            numpy.isfinite(ctrl[k]) and numpy.isfinite(vals[k]))
            item = entry(other, "Paired t-test", len(common))
            if len(common) >= 2:
                a = numpy.array([vals[k] for k in common])
                b = numpy.array([ctrl[k] for k in common])
                n = len(common)
                item["p"] = float(stats.ttest_rel(a, b).pvalue)
                # Bonferroni critical value for several comparisons
                crit = stats.t.ppf(1 - 0.025 / m, n - 1)
                cov = numpy.cov(a, b, ddof=1)
                item["fc"], item["fc_low"], item["fc_high"] = fold_change(
                    a, b, cov[0, 0] / n, cov[1, 1] / n, cov[0, 1] / n, crit,
                    log)
                if m > 1:
                    item["fc_note"] = "Simultaneous 95 % CI (Bonferroni)"
            else:
                item["note"] = "Less than 2 paired experiments"
            results.append(item)
            pvalues.append(item["p"])
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
        return finish_fold_changes(results)

    ctrl_values = numpy.array([v for v in ctrl.values() if numpy.isfinite(v)])
    groups = []
    for other in others:
        vals = numpy.array([v for v in values(other).values()
                            if numpy.isfinite(v)])
        groups.append(vals)
        results.append(entry(other, "Welch t-test" if len(others) == 1
                             else "Dunnett", len(vals)))
    valid = len(ctrl_values) >= 2 and all(len(g) >= 2 for g in groups)
    if not valid:
        for item, vals in zip(results, groups):
            if len(vals) < 2 or len(ctrl_values) < 2:
                item["note"] = "Less than 2 replicates"
        if len(others) > 1:
            return results
    n0 = len(ctrl_values)
    if len(others) == 1:
        if valid:
            a, b = groups[0], ctrl_values
            results[0]["p"] = float(stats.ttest_ind(
                a, b, equal_var=False).pvalue)
            v1, v0 = a.var(ddof=1) / len(a), b.var(ddof=1) / n0
            df = (v1 + v0) ** 2 / (v1 ** 2 / (len(a) - 1) +
                                   v0 ** 2 / (n0 - 1))
            crit = stats.t.ppf(0.975, df)
            item = results[0]
            item["fc"], item["fc_low"], item["fc_high"] = fold_change(
                a, b, v1, v0, 0.0, crit, log)
    else:
        res = stats.dunnett(*groups, control=ctrl_values)
        ci = res.confidence_interval(0.95)
        # Pooled variance of the test, critical value of each comparison
        all_groups = groups + [ctrl_values]
        dof = sum(len(g) for g in all_groups) - len(all_groups)
        mse = sum((len(g) - 1) * g.var(ddof=1) for g in all_groups) / dof
        for i, (item, a) in enumerate(zip(results, groups)):
            item["p"] = float(res.pvalue[i])
            se = numpy.sqrt(mse * (1 / len(a) + 1 / n0))
            crit = (ci.high[i] - ci.low[i]) / 2 / se
            item["fc"], item["fc_low"], item["fc_high"] = fold_change(
                a, ctrl_values, mse / len(a), mse / n0, 0.0, crit, log)
            item["fc_note"] = "Simultaneous 95 % CI (Dunnett)"
    return finish_fold_changes(results)


def compare_all_pairs(conditions, options):
    """Tests of all the pairs of conditions (no natural control), on the
    replicate values (log10 with the log scale), with the fold change of
    each pair and its 95 % confidence interval (same test).

    Unpaired: Tukey's HSD test (simultaneous intervals). Paired
    (replicates with the same experiment number): paired t-tests, Holm
    correction (Bonferroni intervals). Each pair is "B vs A", A being the
    first of the two conditions in the plot.
    Returns the same entries as compare, with "pair": (index A, index B).
    """
    log = options["log"]
    if len(conditions) < 2:
        return []

    def values(condition):
        result = {r["experiment"]: r["value"] for r in condition["replicates"]}
        if log:
            result = {k: (numpy.log10(v) if v > 0 else numpy.nan)
                      for k, v in result.items()}
        return result

    pairs = [(i, j) for i in range(len(conditions))
             for j in range(i + 1, len(conditions))]
    results = []
    for i, j in pairs:
        results.append({"name": conditions[j]["name"] + " vs " +
                        conditions[i]["name"], "pair": (i, j),
                        "test": "", "n": 0, "p": numpy.nan, "note": "",
                        "fc": numpy.nan, "fc_low": numpy.nan,
                        "fc_high": numpy.nan, "fc_note": ""})
    vals = [values(c) for c in conditions]

    if options["paired"]:
        m = len(pairs)
        pvalues = []
        for item, (i, j) in zip(results, pairs):
            common = sorted(k for k in vals[i] if k in vals[j] and
                            numpy.isfinite(vals[i][k]) and
                            numpy.isfinite(vals[j][k]))
            item["n"] = len(common)
            item["test"] = "Paired t-test, Holm" if m > 1 else \
                "Paired t-test"
            if len(common) >= 2:
                a = numpy.array([vals[j][k] for k in common])
                b = numpy.array([vals[i][k] for k in common])
                n = len(common)
                item["p"] = float(stats.ttest_rel(a, b).pvalue)
                crit = stats.t.ppf(1 - 0.025 / m, n - 1)
                cov = numpy.cov(a, b, ddof=1)
                item["fc"], item["fc_low"], item["fc_high"] = fold_change(
                    a, b, cov[0, 0] / n, cov[1, 1] / n, cov[0, 1] / n, crit,
                    log)
                if m > 1:
                    item["fc_note"] = "Simultaneous 95 % CI (Bonferroni)"
            else:
                item["note"] = "Less than 2 paired experiments"
            pvalues.append(item["p"])
        # Holm correction
        finite = [k for k, p in enumerate(pvalues) if numpy.isfinite(p)]
        order = sorted(finite, key=lambda k: pvalues[k])
        running = 0.0
        for rank, k in enumerate(order):
            running = max(running, min(1.0, (len(order) - rank) * pvalues[k]))
            results[k]["p"] = running
        return finish_fold_changes(results)

    groups = [numpy.array([v for v in g.values() if numpy.isfinite(v)])
              for g in vals]
    for item, (i, j) in zip(results, pairs):
        item["n"] = len(groups[j])
        item["test"] = "Tukey"
    if any(len(g) < 2 for g in groups):
        for item, (i, j) in zip(results, pairs):
            if len(groups[i]) < 2 or len(groups[j]) < 2:
                item["note"] = "Less than 2 replicates"
        return results
    res = stats.tukey_hsd(*groups)
    ci = res.confidence_interval(0.95)
    # Pooled variance of the test, critical value of each comparison
    dof = sum(len(g) for g in groups) - len(groups)
    mse = sum((len(g) - 1) * g.var(ddof=1) for g in groups) / dof
    for item, (i, j) in zip(results, pairs):
        a, b = groups[j], groups[i]
        item["p"] = float(res.pvalue[j, i])
        se = numpy.sqrt(mse * (1 / len(a) + 1 / len(b)))
        crit = (ci.high[j, i] - ci.low[j, i]) / 2 / se
        item["fc"], item["fc_low"], item["fc_high"] = fold_change(
            a, b, mse / len(a), mse / len(b), 0.0, crit, log)
        item["fc_note"] = "Simultaneous 95 % CI (Tukey)"
    return finish_fold_changes(results)


def finish_fold_changes(results):
    """Fold changes are positive (ratio of elasticities, forces...): lower
    bound at 0; note for the unbounded intervals."""
    for item in results:
        if numpy.isfinite(item["fc_low"]) and item["fc_low"] < 0:
            item["fc_low"] = 0.0
        if numpy.isfinite(item["fc"]) and not numpy.isfinite(item["fc_low"]):
            item["fc_note"] = "Unbounded interval: the mean of the control " \
                "is too uncertain (use the log scale for paired " \
                "experiments with a day-to-day variation)"
    return results


def export_rows(conditions):
    """Rows of the csv export: one per file, with its replicate and
    condition."""
    rows = []
    for condition in conditions:
        for rep in condition["replicates"]:
            for item in rep["maps"]:
                rows.append({"condition": condition["name"],
                             "replicate": rep["name"],
                             "experiment": rep["experiment"],
                             "file": item["name"],
                             "file_value": item["value"],
                             "nbr_curves": item["nbr_curves"],
                             "replicate_mean": rep["value"]})
    return rows


def test_rows(conditions, options):
    """Rows of the csv export of the tests: one per comparison, with the
    p-value and the fold change vs the control with its 95 % confidence
    interval."""
    rows = []
    for test in compare(conditions, options):
        rows.append({"comparison": test["name"], "test": test["test"],
                     "n": test["n"], "p_value": test["p"],
                     "significance": stars(test["p"]),
                     "fold_change": test["fc"],
                     "fold_change_ci95_low": test["fc_low"],
                     "fold_change_ci95_high": test["fc_high"],
                     "scale": "log (geometric means)" if options["log"]
                     else "linear (means)",
                     "note": "; ".join(n for n in (test["note"],
                                                   test["fc_note"]) if n)})
    return rows
