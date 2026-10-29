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

"""Exports of the results (see the Export menu): CSV tables (one row per
curve or event, summary per row), SuperPlot values and tests, DFS fits,
maps (32-bit TIFF of the values for ImageJ / Fiji, and figures) and the
parameters (JSON).

Every data export writes its parameters next to it (name.params.json, see
tools.params_export). The values are in SI units (elasticity in Pa, or the
stiffness in N/m for the slope model, work in J, forces in N, loading rates
in N/s, distances in m, indentation in nm).
"""

import csv
import os

import matplotlib
import numpy as np
import pandas as pd
import tifffile

from .. import shared
from ..tools import dfs_tools
from ..tools import map_values
from ..compute_tools import microrheology
from ..tools import params_export
from ..tools import results_sorting
from ..tools import superplot_tools
from ..tools.breakthrough_tools import EVENT_VALUES
from ..tools.colortables import ColorTables


def csv_path(path):
    """Path with the .csv extension."""
    return path if path.lower().endswith(".csv") else path + ".csv"


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


def data_tables(rows=None):
    """Tables of the rows of the Plots tab (the checked ones by default):
    values (one line per curve, or per event for the curves with events,
    the values of the curve being repeated on the lines of its events;
    filtered as in the plots) and descriptive statistics of each row."""
    tables = []
    stats = []
    if rows is None:
        selected = exported_rows()
    else:
        selected = [(row, shared.exp.results_list[row],
                     shared.exp.list[shared.exp.results_list[row].data_id])
                    for row in rows]
    for row, result, data in selected:
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
    stats_df = pd.concat(stats, ignore_index=True) if stats \
        else pd.DataFrame()
    return data_df, stats_df


def checked_rows():
    """Indices of the checked rows of the Plots tab."""
    return [row for row, _, _ in exported_rows()]


def export_data_table(path, rows=None):
    """Values of the rows (the checked ones by default), one line per curve
    (or event). Returns the written files."""
    path = csv_path(path)
    rows = checked_rows() if rows is None else rows
    data_df, _ = data_tables(rows)
    data_df.to_csv(path, index=False)
    json_path = params_export.write_sidecar(
        path, "Data table (one row per curve or event)", rows=rows)
    return [path, json_path]


def export_summary(path, rows=None):
    """Descriptive statistics of each row (the checked ones by default).
    Returns the written files."""
    path = csv_path(path)
    rows = checked_rows() if rows is None else rows
    _, stats_df = data_tables(rows)
    stats_df.to_csv(path, index=False)
    json_path = params_export.write_sidecar(
        path, "Summary per row", rows=rows)
    return [path, json_path]


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

    for force_type, values in results_sorting.PIXEL_TYPES.items():
        if force_type in results_sorting.MICRORHEO_TYPES:
            # All the frequencies, below
            continue
        if getattr(data, values["calculated"]):
            df = table(get(force_type), values["column"])
            stats.append(describe(df, values["column"]))
            curve_tables.append(df)

    if data.microrheo_calculated:
        # G' and G'' of the curves at each frequency (Pa)
        curves = get("microrheo_storage")
        for frequency in microrheology.frequencies(data):
            for name, column in (("storage", "G_storage"),
                                 ("loss", "G_loss")):
                column += "_" + format(frequency, "g") + "Hz"
                values = data.get_microrheo_result(name, frequency)
                df = table([[index, i, j, float(values[i][j])]
                            for index, i, j, _ in curves], column)
                stats.append(describe(df, column))
                curve_tables.append(df)

    if data.breakthrough_calculated:
        # Number of breakthroughs per curve, and all the breakthroughs
        count = data.breakthrough_count
        curves = get("breakthrough_force")
        df = table([[index, i, j, count[i][j]] for index, i, j, _ in curves],
                   "breakthroughs_per_curve")
        stats.append(describe(df, "breakthroughs_per_curve"))
        curve_tables.append(df)
        force = EVENT_VALUES.index("force")
        thickness = EVENT_VALUES.index("thickness")
        lines = [[index, i, j, float(event[force]), float(event[thickness]),
                  k] for index, i, j, _ in curves
                 for k, event in enumerate(data.breakthrough_results[i, j])
                 if np.isfinite(event[force])]
        columns = CURVE_KEYS + ["breakthrough_forces", "breakthrough_jumps",
                                "event"]
        df = pd.DataFrame(lines, columns=columns)
        for key in CURVE_KEYS + ["event"]:
            df[key] = df[key].astype(int) + 1
        event_tables.append(df)

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


def superplot_rows():
    """Rows of the Plots table in the SuperPlots: in a displayed replicate
    of a displayed condition."""
    groups = {g.group_id for g in shared.exp.groups_list[1:]
              if g.display and g.condition and
              shared.exp.conditions_list[g.condition].display}
    return [row for row, result in enumerate(shared.exp.results_list)
            if result.group in groups]


def write_rows(path, rows):
    """Write a list of dictionaries as a CSV table."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export_superplot(path):
    """Values of the files of the SuperPlot (with their replicate and
    condition) and the tests (_tests.csv). Returns the written files, empty
    without SuperPlot."""
    conditions = superplot_tools.get_data()
    rows = superplot_tools.export_rows(conditions)
    if not rows:
        return []
    path = csv_path(path)
    files = [path]
    write_rows(path, rows)
    options = superplot_tools.get_options()
    tests = superplot_tools.test_rows(conditions, options)
    if tests:
        tests_path = path[:-4] + "_tests.csv"
        write_rows(tests_path, tests)
        files.append(tests_path)
    control = next((c.name for c in shared.exp.conditions_list[1:]
                    if c.condition_id == options["control"]), None)
    files.append(params_export.write_sidecar(
        path, "SuperPlot values and tests", rows=superplot_rows(),
        options={"superplot": {
            "control": "all pairs"
            if options["control"] == superplot_tools.ALL_PAIRS else control,
            "replicates": "paired" if options["paired"] else "unpaired",
            "test": "mixed model" if options["test"] == "mixed"
            else "replicate means",
            "log_scale": bool(options["log"])}}))
    return files


def dfs_rows(mode):
    """Rows of the Plots table of the DFS fits: the displayed files, or the
    files of the displayed replicates (speeds)."""
    if mode == "groups":
        groups = {g.group_id for g in shared.exp.groups_list[1:]
                  if g.display}
        return [row for row, result in enumerate(shared.exp.results_list)
                if result.group in groups]
    return checked_rows()


def export_dfs_fits(path, mode=None):
    """DFS fits of the Plots tab (mode "single") or of the Replicates tab
    ("groups"): model, parameters with their standard errors, R2, and their
    points (_points.csv). Returns the written files, empty without fit."""
    mode = mode or shared.exp.dfs_mode
    fits = [f for f in dfs_tools.get_fits(mode) if f.get("method") != "old"]
    if not fits:
        return []
    path = csv_path(path)
    lines, points = [], []
    for number, fit in enumerate(fits, start=1):
        model = fit.get("model", "bell")
        lines.append({
            "fit": number, "model": dfs_tools.MODELS[model][0],
            "points_from": "files" if mode == "single" else "speeds",
            "method": fit["method"], "bins": fit.get("nbr_bins")
            if fit["method"] == "binned" else "",
            "statistic": fit.get("statistic", "median"),
            "nbr_points": fit["nbr_points"], "nbr_events": fit["nbr_events"],
            "x_m": fit["x"], "x_se_m": fit["sx"],
            "koff_s-1": fit["koff"], "koff_se_s-1": fit["skoff"],
            "dG_kBT": fit.get("dG", ""), "dG_se_kBT": fit.get("sdG", ""),
            "feq_N": fit["feq"] * 1e-12 if "feq" in fit else "",
            "feq_se_N": fit["sfeq"] * 1e-12 if "sfeq" in fit else "",
            "reduces_to_bell_evans": bool(fit.get("bell_limit", False)),
            "r2": fit["r2"], "temperature_K": fit["temperature"]})
        quartiles = fit.get("points_f_q", [])
        for k, (r, f) in enumerate(zip(fit["points_r"], fit["points_f"])):
            point = {"fit": number, "loading_rate_N_s": r * 1e-12,
                     "force_N": f * 1e-12}
            if k < len(quartiles):
                point["force_q25_N"] = quartiles[k][0] * 1e-12
                point["force_q75_N"] = quartiles[k][1] * 1e-12
            points.append(point)
    write_rows(path, lines)
    points_path = path[:-4] + "_points.csv"
    keys = max(points, key=len).keys()
    write_rows(points_path, [{k: p.get(k, "") for k in keys} for p in points])
    json_path = params_export.write_sidecar(
        path, "DFS fits", rows=dfs_rows(mode),
        options={"dfs": {"points_from": "files" if mode == "single"
                         else "speeds"}})
    return [path, points_path, json_path]


def map_type_name(map_type):
    """Short name of a map for the file names."""
    return map_type.replace("stiffness_corr", "elasticity_bec").replace(
        "stiffness", "elasticity")


def map_lut(data, map_type):
    """Colour table of a map (as in pyAF) for ImageJ: (3, 256) uint8."""
    opts = data.get_color_opts(map_type)
    tables = ColorTables()
    name = tables.colortables_mpl_list[opts[0]]
    if name in ("Rainbow 1", "Cool-Warm"):
        factor = map_values.factor(data, map_type)
        cmap, _ = ColorTables(opts[0], None, None, opts[2] * factor,
                              opts[3] * factor).get_color_as_cdict()
    else:
        cmap = matplotlib.colormaps[name]
    colors = cmap(np.linspace(0.0, 1.0, 256))[:, :3]
    return np.round(colors * 255).astype(np.uint8).T


def export_map(path, data_id, map_type):
    """Values of a map of a file as a 32-bit TIFF (NaN where there is no
    value; in the units of the colour bar, see map_values) for ImageJ /
    Fiji: pixel size in um, the colour table and the range of the colour
    scale of pyAF (the values are not changed), and the parameters in the
    information of the image (Image > Show Info). Returns the written
    file."""
    data = shared.exp.list[data_id]
    factor = map_values.factor(data, map_type)
    values = map_values.values(data, map_type) * factor
    # Rows of the image: y, from the top
    image = np.ascontiguousarray(np.flipud(values.T), dtype=np.float32)
    opts = data.get_color_opts(map_type)
    low, high = float(opts[1] * factor), float(opts[2] * factor)
    if not high > low:
        finite = image[np.isfinite(image)]
        low, high = (float(finite.min()), float(finite.max())) \
            if finite.size else (0.0, 1.0)
    pixel_x = data.scan_size_x / data.nbr_pixels_x / 1000.0
    pixel_y = data.scan_size_y / data.nbr_pixels_y / 1000.0
    info = params_export.tiff_description(
        "Map " + map_values.label(data, map_type), path, [data_id],
        options={"map": map_type, "unit": map_values.label(data, map_type),
                 "pixel_size_um": [pixel_x, pixel_y]})
    tifffile.imwrite(
        path, image, imagej=True,
        resolution=(1.0 / pixel_x, 1.0 / pixel_y),
        metadata={"unit": "um", "min": low, "max": high,
                  "Ranges": (low, high), "LUTs": [map_lut(data, map_type)],
                  "Info": info})
    return [path]


def export_map_figure(path, data_id, dpi=300):
    """Figure of the map of a file as displayed in the Results tab (axes,
    colour bar), PNG, PDF or SVG (by the extension), with its parameters
    next to it. Returns the written files."""
    from ..plots.PYAFPlot import PYAFPlot

    previous = shared.exp.id_selected
    red_square = shared.exp.meshgrid_display_red_square
    shared.exp.id_selected = data_id
    # Without the square of the selected curve (on the screen only)
    shared.exp.meshgrid_display_red_square = False
    try:
        with matplotlib.rc_context({"savefig.dpi": dpi,
                                    "savefig.bbox": "tight"}):
            PYAFPlot(None, "meshgrid", save_path=path)
    finally:
        shared.exp.id_selected = previous
        shared.exp.meshgrid_display_red_square = red_square
    data = shared.exp.list[data_id]
    json_path = params_export.write_sidecar(
        path, "Map figure", data_ids=[data_id],
        options={"map": data.meshgrid_type})
    return [path, json_path]


def export_parameters(path):
    """Parameters of all the files. Returns the written file."""
    if not path.lower().endswith(".json"):
        path += ".json"
    params_export.write(path, params_export.build(
        "Parameters", path, data_ids=list(range(len(shared.exp.list))),
        rows=list(range(len(shared.exp.results_list)))))
    return [path]


# Result types of the SuperPlots, with the computation they need
RESULT_TYPES = (("stiffness", "stiffness_calculated"),
                ("work", "work_and_rupture_force1_calculated"),
                ("rupture_force", "work_and_rupture_force1_calculated"),
                ("events_forces", "events_calculated"),
                ("events_per_curve", "events_calculated"),
                ("events_rupture_force", "events_calculated"),
                ("events_distance", "events_calculated"),
                ("visco_e0", "visco_calculated"),
                ("visco_fluidity", "visco_calculated"),
                ("visco_tan_delta", "visco_calculated"),
                ("microrheo_storage", "microrheo_calculated"),
                ("microrheo_loss", "microrheo_calculated"),
                ("microrheo_tan_delta", "microrheo_calculated"),
                ("breakthrough_force", "breakthrough_calculated"),
                ("breakthrough_thickness", "breakthrough_calculated"),
                ("brush_modulus", "brush_calculated"),
                ("brush_length", "brush_calculated"),
                ("brush_density", "brush_calculated"))


def computed_result_types():
    """Result types computed for at least one file."""
    return [name for name, flag in RESULT_TYPES
            if any(getattr(data, flag, False) for data in shared.exp.list)]


def export_superplots(folder):
    """SuperPlot of each computed result type (superplot_<type>.csv), the
    values of the rows being computed for each type, then for the
    displayed one again. Returns the written files."""
    from . import stat_tools

    if not any(g.condition for g in shared.exp.groups_list[1:]):
        return []
    files = []
    displayed = shared.exp.results_type
    try:
        for name in computed_result_types():
            shared.exp.results_type = name
            stat_tools.get_values()
            files += export_superplot(os.path.join(
                folder, "superplot_" + map_type_name(name) + ".csv"))
    finally:
        shared.exp.results_type = displayed
        stat_tools.get_values()
    return files


def export_all(folder):
    """Everything computed in a folder: data table and summary of all the
    rows of the Plots tab (all the computed results of each row), SuperPlot
    of each computed result type, DFS fits, all the computed maps of the
    files (values in maps/, figures in figures/) and the parameters.
    Returns the written files."""
    os.makedirs(folder, exist_ok=True)
    files = []
    rows = list(range(len(shared.exp.results_list)))
    if rows:
        files += export_data_table(os.path.join(folder, "data_table.csv"),
                                   rows)
        files += export_summary(os.path.join(folder, "summary.csv"), rows)
    files += export_superplots(folder)
    for mode, name in (("single", "dfs_fits_files.csv"),
                       ("groups", "dfs_fits_speeds.csv")):
        files += export_dfs_fits(os.path.join(folder, name), mode)
    maps = os.path.join(folder, "maps")
    figures = os.path.join(folder, "figures")
    os.makedirs(maps, exist_ok=True)
    os.makedirs(figures, exist_ok=True)
    for data_id, data in enumerate(shared.exp.list):
        base = os.path.splitext(os.path.basename(str(data.filename)))[0]
        displayed = data.meshgrid_type
        try:
            for map_type in map_values.available_maps(data):
                name = base + "_" + map_type_name(map_type)
                files += export_map(os.path.join(maps, name + ".tif"),
                                    data_id, map_type)
                # The figure of the map (drawn for the map of the file)
                data.meshgrid_type = map_type
                files += export_map_figure(
                    os.path.join(figures, name + ".png"), data_id)
        finally:
            data.meshgrid_type = displayed
    files += export_parameters(os.path.join(folder, "parameters.json"))
    return files
