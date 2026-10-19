# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
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

"""Exports the results as text files (txt files, or a csv table for R)."""

import os
import shutil

import pandas as pd
import numpy as np

from .. import shared
from .. import widgets_list as wl
from ..tools import results_sorting
from ..widgets.progressbar import Progressbar


def write_file(path, header, lines):
    """Writes a tab separated text file."""
    with open(path, "wb") as newfile:
        newfile.write((header + "\r\n").encode(encoding="UTF-8"))
        for line in lines:
            newfile.write((line + "\r\n").encode(encoding="UTF-8"))


def exported_rows():
    """The checked rows of the plots tab: (row index, result, data)."""
    rows = []
    for row, result in enumerate(shared.exp.results_list):
        if result.display:
            rows.append((row, result, shared.exp.list[result.data_id]))
    return rows


def group_and_condition(result):
    """The group and the condition (objects) of a row."""
    group = shared.exp.groups_list[result.group]
    condition = shared.exp.conditions_list[group.condition]
    return group, condition


def elasticity_title(data):
    """Title of the elasticity column (the slope model gives a stiffness)."""
    if data.used_stiffness_model_selected == 3:
        return "Stiffness (N/m)"
    return "Elasticity (Pa)"


def slice_label(data, i, j, k):
    """Indentation segment [nm] of the slice k of the curve [i, j]."""
    start, step, tomography = data.get_used_indentation()
    if step == 0:
        # A single segment, up to the indentation reached by the curve
        return str(start) + " - " + str(round(float(data.indentation[i][j]), 1))
    low = start + k * step if tomography else start
    return str(low) + " - " + str(start + (k + 1) * step)


def export_results(folder):
    """Exports the results as text files.

    The results are saved as filtered (results displayed in the plots) and
    unfiltered (raw results). The results of the replicates (the groups) are
    also saved, with a file listing the replicates, their condition and their
    files.
    """
    rows = exported_rows()
    Progressbar()
    wl.widget_progressbar.set_label("Exporting")
    wl.widget_progressbar.set_range(0, 2 * len(rows) + 8)

    # Create folders
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)
    os.mkdir(folder + "/filtered")
    os.mkdir(folder + "/unfiltered")
    os.mkdir(folder + "/replicates")

    # Unfiltered data (all the curves, one file per checked row)
    for _, result, data in rows:
        write_unfiltered(folder + "/unfiltered/" + data.filename, data)
        wl.widget_progressbar.update()

    # Events per scan (unfiltered)
    if any(data.events_calculated for _, _, data in rows):
        lines = []
        for index, (_, result, data) in enumerate(rows):
            if data.events_calculated:
                val = sum(len(data.events_forces[i][j])
                          for i in range(data.nbr_pixels_x)
                          for j in range(data.nbr_pixels_y))
                lines.append(str(index + 1) + "\t" + result.name + "\t" +
                             str(val))
        write_file(folder + "/unfiltered/Events_per_scan.txt",
                   "Index\tName\tNumber of events", lines)

    # Filtered data (as in the plots: ROI, slice, filters of the row)
    per_scan = []
    for index, (row, result, data) in enumerate(rows):
        pre = folder + "/filtered/" + result.name
        nbr = write_filtered(pre, row, data)
        if nbr is not None:
            per_scan.append(str(index + 1) + "\t" + result.name + "\t" +
                            str(nbr))
        wl.widget_progressbar.update()

    if per_scan:
        write_file(folder + "/filtered/Events_per_scan.txt",
                   "Index\tName\tNumber of events", per_scan)

    # Replicates
    for res_type in ("stiffness", "work", "rupture_force", "events_forces",
                     "events_per_curve", "events_rupture_force",
                     "loading_rates", "events_distance") + \
            tuple(results_sorting.VISCO_TYPES):
        write_group_file(res_type, folder)
        wl.widget_progressbar.update()

    # Infos
    write_info_file(folder)

    wl.widget_progressbar.close()

    return True


def write_unfiltered(pre, data):
    """Writes the raw results of a file (all the curves)."""
    nx, ny = data.nbr_pixels_x, data.nbr_pixels_y

    def pixels():
        index = 0
        for i in range(nx):
            for j in range(ny):
                yield index, i, j
                index += 1

    def start(index, i, j):
        return str(index + 1) + "\t" + str(i + 1) + "\t" + str(j + 1) + "\t"

    if data.stiffness_calculated:
        lines = [start(index, i, j) + slice_label(data, i, j, k) + "\t" +
                 str(data.stiffness_array[i][j][k])
                 for index, i, j in pixels()
                 for k in range(len(data.stiffness_array[i][j]))]
        write_file(pre + "_elasticity.txt",
                   "Index\tXpos\tYpos\tSlice (nm)\t" + elasticity_title(data),
                   lines)

    if data.work_and_rupture_force1_calculated:
        write_file(pre + "_work.txt", "Index\tXpos\tYpos\tWork (J)",
                   [start(index, i, j) + str(data.work[i][j])
                    for index, i, j in pixels()])
        write_file(pre + "_rupture_force1.txt",
                   "Index\tXpos\tYpos\tRupture force (N)",
                   [start(index, i, j) + str(data.rupture_force1[i][j])
                    for index, i, j in pixels()])

    if data.events_calculated:
        write_file(pre + "_events_forces.txt",
                   "Index\tXpos\tYpos\tEvent Force (N)",
                   [start(index, i, j) + str(force)
                    for index, i, j in pixels()
                    for force in data.events_forces[i][j]])
        write_file(pre + "_event_distance.txt",
                   "Index\tXpos\tYpos\tDistance (m)",
                   [start(index, i, j) + str(dist)
                    for index, i, j in pixels()
                    for dist in data.events_forces_distance[i][j]])

    if data.loading_rates_calculated:
        write_file(pre + "_loading_rates.txt",
                   "Index\tXpos\tYpos\tLoading rate (N/s)",
                   [start(index, i, j) + str(lr)
                    for index, i, j in pixels()
                    for lr in data.loading_rates[i][j]])

    if data.visco_calculated:
        # All the parameters of the Ting fit (nan: not fitted)
        maps = [data.visco_e0, data.visco_fluidity, data.visco_tan_delta,
                data.visco_storage, data.visco_loss,
                data.get_visco_result("tc"), data.get_visco_result("f0"),
                data.get_visco_result("r2"), data.get_visco_result("height")]
        freq = "%g Hz" % data.visco_frequency
        write_file(pre + "_viscoelasticity.txt",
                   "Index\tXpos\tYpos\tE0 (Pa)\tFluidity\t"
                   "Loss tangent (tan delta)\tStorage modulus E' at " +
                   freq + " (Pa)\tLoss modulus E'' at " + freq + " (Pa)\t"
                   "Time of contact (s)\t"
                   "F0 (N)\tR2\tB.E.C. height (nm)",
                   [start(index, i, j) +
                    "\t".join(str(float(m[i][j])) for m in maps)
                    for index, i, j in pixels()])


def write_filtered(pre, row, data):
    """Writes the results of a row as in the plots. Returns the number of
    kept events (None without events)."""
    def get(force_type):
        util = results_sorting.GetResults(export="single")
        return util.get_results_single(result_id=row, force_type=force_type)

    def lines(array, value_col=3):
        return [str(line[0] + 1) + "\t" + str(line[1] + 1) + "\t" +
                str(line[2] + 1) + "\t" + str(line[value_col])
                for line in array]

    if data.stiffness_calculated:
        array = get("stiffness")
        write_file(pre + "_elasticity.txt",
                   "Index\tXpos\tYpos\tSlice (nm)\t" + elasticity_title(data),
                   [str(line[0] + 1) + "\t" + str(line[1] + 1) + "\t" +
                    str(line[2] + 1) + "\t" +
                    slice_label(data, line[1], line[2], line[3]) + "\t" +
                    str(line[4])
                    for line in array])

    if data.work_and_rupture_force1_calculated:
        write_file(pre + "_work.txt", "Index\tXpos\tYpos\tWork (J)",
                   lines(get("work")))
        write_file(pre + "_rupture_force1.txt",
                   "Index\tXpos\tYpos\tRupture force (N)",
                   lines(get("rupture_force")))

    nbr_events = None
    if data.events_calculated:
        write_file(pre + "_events_forces.txt",
                   "Index\tXpos\tYpos\tEvent Force (N)",
                   lines(get("events_forces")))
        per_curve = get("events_per_curve")
        nbr_events = int(sum(line[3] for line in per_curve))
        write_file(pre + "_events_per_curve.txt",
                   "Index\tXpos\tYpos\tEvents per curve", lines(per_curve))
        write_file(pre + "_event_max_force.txt",
                   "Index\tXpos\tYpos\tEvent max. force (N)",
                   lines(get("events_rupture_force")))
        write_file(pre + "_event_distance.txt",
                   "Index\tXpos\tYpos\tDistance (m)",
                   lines(get("events_distance")))

    if data.loading_rates_calculated:
        write_file(pre + "_loading_rates.txt",
                   "Index\tXpos\tYpos\tLoading rate (N/s)",
                   lines(get("loading_rates")))

    if data.visco_calculated:
        for res_type, values in results_sorting.VISCO_TYPES.items():
            write_file(pre + "_" + values["column"] + ".txt",
                       "Index\tXpos\tYpos\t" + values["export"],
                       lines(get(res_type)))

    return nbr_events


def export_results_r(path, stats_flag=True):
    """Export the results of the checked rows as a csv table for R.

    One line per curve, or per event for the curves with events (the values
    of the curve are repeated on the lines of its events). The values are
    filtered as in the plots and are in SI units (elasticity in Pa, or the
    stiffness in N/m for the slope model, work in J, forces in N, loading
    rates in N/s, distances in m, indentation in nm). The descriptive
    statistics of each row are saved in a second file (_stats).
    """
    if not path.lower().endswith(".csv"):
        path += ".csv"

    tables = []
    stats = []
    for row, result, data in exported_rows():
        curves, events, row_stats = get_r_tables(row, data)

        if curves is None and events is None:
            continue
        if curves is None:
            table = events
        elif events is None:
            table = curves
        else:
            table = pd.merge(curves, events, how="outer", on=CURVE_KEYS)
        keys = [key for key in CURVE_KEYS + ["event"] if key in table.columns]
        if "event" in table.columns:
            # Integers, empty for the curves without events
            table["event"] = table["event"].astype("Int64")
        table = table.sort_values(keys)
        table = table[keys + [c for c in table.columns if c not in keys]]

        info = file_info(result, data)
        for pos, (name, value) in enumerate(info.items()):
            table.insert(pos, name, value)
            row_stats.insert(pos, name, value)
        tables.append(table)
        stats.append(row_stats)

    data_df = pd.concat(tables, ignore_index=True) if tables \
        else pd.DataFrame()
    data_df.to_csv(path, index=False)

    if stats_flag:
        stats_df = pd.concat(stats, ignore_index=True) if stats \
            else pd.DataFrame()
        stats_df.to_csv(path[:-4] + "_stats.csv", index=False)


CURVE_KEYS = ["curve_index", "x_pos", "y_pos"]


def file_info(result, data):
    """Columns describing the row: condition, replicate and file."""
    group, condition = group_and_condition(result)
    return {"condition": condition.name,
            "replicate": group.name,
            "name": result.name,
            "file_name": data.filename,
            "condition_id": condition.condition_id,
            "replicate_id": group.group_id,
            "file_id": result.data_id}


def get_r_tables(row, data):
    """Tables of a row for the R export: values per curve, values per event
    and descriptive statistics."""
    def get(force_type):
        util = results_sorting.GetResults(export="r_format")
        return util.get_results_single(result_id=row, force_type=force_type)

    def table(array, name, event=False):
        columns = CURVE_KEYS + [name] + (["event"] if event else [])
        df = pd.DataFrame([line[:len(columns)] for line in array],
                          columns=columns)
        # Indexes from 1, as in the txt export
        for key in CURVE_KEYS + (["event"] if event else []):
            df[key] = df[key].astype(int) + 1
        return df

    curve_tables = []
    event_tables = []
    stats = []

    if data.stiffness_calculated:
        array = get("stiffness")
        df = pd.DataFrame([line[:5] for line in array],
                          columns=CURVE_KEYS + ["elasticity", "slice"])
        stats.append(describe(df, "elasticity"))
        for key in CURVE_KEYS:
            df[key] = df[key].astype(int) + 1
        if len(df) and df["slice"].max() > 0:
            # Tomography: one column per slice
            df = df.pivot_table(index=CURVE_KEYS, columns="slice",
                                values="elasticity", aggfunc="first")
            df.columns = ["elasticity_slice_" + str(int(k) + 1)
                          for k in df.columns]
            df = df.reset_index()
        else:
            df = df.drop(columns="slice")
        # Indentation reached by the curves
        df["indentation"] = [float(data.indentation[x - 1][y - 1])
                             for x, y in zip(df["x_pos"], df["y_pos"])]
        curve_tables.append(df)

    if data.work_and_rupture_force1_calculated:
        for force_type, name in (("work", "work"),
                                 ("rupture_force", "rupture_force")):
            df = table(get(force_type), name)
            stats.append(describe(df, name))
            curve_tables.append(df)

    if data.visco_calculated:
        for force_type, values in results_sorting.VISCO_TYPES.items():
            df = table(get(force_type), values["column"])
            stats.append(describe(df, values["column"]))
            curve_tables.append(df)

    if data.events_calculated:
        for force_type, name in (("events_per_curve", "events_per_curve"),
                                 ("events_rupture_force", "event_max_force")):
            df = table(get(force_type), name)
            stats.append(describe(df, name))
            curve_tables.append(df)

        types = [("events_forces", "events_forces"),
                 ("events_distance", "events_distance")]
        if data.loading_rates_calculated:
            types.insert(1, ("loading_rates", "loading_rates"))
        for force_type, name in types:
            df = table(get(force_type), name, event=True)
            stats.append(describe(df, name))
            event_tables.append(df)

    curves = merge_tables(curve_tables, CURVE_KEYS)
    events = merge_tables(event_tables, CURVE_KEYS + ["event"])
    row_stats = pd.concat(stats, axis=1) if stats else pd.DataFrame([{}])
    return curves, events, row_stats


def merge_tables(tables, keys):
    """Outer merge of the tables on the keys (None without table)."""
    if not tables:
        return None
    merged = tables[0]
    for df in tables[1:]:
        merged = pd.merge(merged, df, how="outer", on=keys)
    return merged


def describe(df, name):
    """Descriptive statistics of a column, as one line (count_x, mean_x, std_x,
    min_x, 25%_x, 50%_x, 75%_x, max_x)."""
    values = df[name].describe()
    return pd.DataFrame([{stat + "_" + name: values[stat]
                          for stat in values.index}])


def write_group_file(res_type, folder):
    """Write file containing the results of the replicates (the groups).

    Each replicate is written in a column. The columns are separated by
    semi-colons; a replicate with less values has empty cells.

    NOTE: The results_sorting.py file provides the values used here.
    All the values have been scaled to be displayed directly in the sotware.
    We need to scale them back for exporting in m, Pa, N, N/s or J.
    """
    if res_type == "stiffness":
        model3 = any(shared.exp.list[res.data_id].used_stiffness_model_selected
                     == 3 for res in shared.exp.results_list)
        if shared.exp.hist_log:
            # log10(E [Pa]) in the plots
            convert = lambda v: 10 ** v  # noqa: E731
        elif model3:
            convert = lambda v: v  # noqa: E731
        else:
            convert = lambda v: v * 1e3  # noqa: E731
        title = "Stiffness (N/m)" if model3 else "Elasticity (Pa)"
    else:
        factor, title = {
            "work": (1e-15, "Detachment work (J)"),
            "rupture_force": (1e-12, "Detachment force (N)"),
            "events_forces": (1e-12, "Events forces (N)"),
            "events_per_curve": (1, "Events per curve"),
            "events_rupture_force": (1e-12, "Event max. force (N)"),
            "loading_rates": (1e-12, "Loading rates (N/s)"),
            "events_distance": (1e-9, "Event's distance (m)"),
            **{name: (1 / values["factor"], values["export"])
               for name, values in results_sorting.VISCO_TYPES.items()},
        }[res_type]
        convert = lambda v: v * factor  # noqa: E731

    results_util = results_sorting.GetResults(export="groups")
    [groups, _] = results_util.get_results_groups(force_type=res_type)

    # Values of the groups (group 0 is the None group)
    columns = []
    for group_id in range(1, len(shared.exp.groups_list)):
        values = groups[group_id]
        if res_type == "loading_rates" and len(values) > 0:
            # Forces and loading rates
            values = values[1]
        if len(values) > 0:
            name = shared.exp.groups_list[group_id].name
            columns.append((name, [str(convert(v)) for v in values]))

    if not columns:
        return

    max_len = max(len(values) for _, values in columns)
    lines = [";".join(values[line_id] if line_id < len(values) else ""
                      for _, values in columns)
             for line_id in range(max_len)]
    write_file(folder + "/replicates/Replicates_" + res_type + ".txt",
               title + "\r\n" + ";".join(name for name, _ in columns), lines)


def write_info_file(folder):
    """Writes a file containing the list of replicates (with their condition)
    and of their rows."""
    lines = []
    for group in shared.exp.groups_list[1:]:
        condition = shared.exp.conditions_list[group.condition]
        lines.append(group.name + " (condition: " + condition.name + ")")
        for result in shared.exp.results_list:
            if result.group == group.group_id:
                lines.append(result.name)
        lines.append("------------------------------------------")
    write_file(folder + "/info.txt", "Replicates and their files", lines)
