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


"""Plant cells (see tools.plant_tools): apparent stiffness k_as of the
curves (approach and retraction), and, for the curves recorded in a cell
(ROI of a map, the map of the file or the map where a single curve was
recorded), the modulus E and the turgor pressure P of the cell wall
(Yamasaki et al. 2025), and the pressure from the elasticity of the Hertz
fit (Beauzamy et al. 2015).

The values of the curves are computed once (data.plant_results [i, j,
CURVE_VALUES]); the values of the cells are computed from them when they
are used (cell_values), with the current ROIs and wall thickness: drawing a
cell or changing the thickness needs no new computation."""

import os

import numpy

from .. import shared
from .. import widgets_list
from ..tools import misc_tools
from ..tools import plant_tools
from ..widgets.progressbar import Progressbar
from . import multiproc

# Values of each curve, stored (nan when not found): k_as [N/m] (mean of
# the approach and retraction), fit of each one, elasticity of the Hertz
# fit E Hertz [Pa]
CURVE_VALUES = ("k_as", "k_approach", "k_retraction", "r2_approach",
                "r2_retraction", "modulus_hertz")
# Values of the curves in a cell, computed when used: E [Pa], P [Pa] and
# tau (Yamasaki), P [Pa] and tau from E Hertz (Beauzamy), structural factor
# f of the cell, number of its ROI (from 1)
CELL_VALUES = ("modulus", "pressure", "tau", "pressure_hertz", "tau_hertz",
               "f", "roi")
VALUES = CURVE_VALUES + CELL_VALUES

# Settings of the computation (attributes plant_* of the files); the
# thickness of the wall (plant_wall_thickness) is applied when the values
# are used
SETTINGS = ("method", "smoothing", "force_start", "force_stop")

# Values of the cells of the files: {id of the file: (key, values)} (see
# cell_values), not saved
_CELLS = {}


def settings(data, used=False):
    """Settings of the computation (the ones used for the maps with
    used=True)."""
    prefix = "used_plant_" if used else "plant_"
    return {name: getattr(data, prefix + name) for name in SETTINGS}


def check_file(data):
    """None if the file can be computed, else the reason."""
    if not data.stiffness_calculated:
        return "the elasticity has not been computed (point of contact " \
            "and E Hertz)"
    return None


def cell_rois(data):
    """ROIs of a map drawn on cells: [(number from 0, ROI)], all but the
    substrate of the B.E.C."""
    return [(number, roi) for number, roi in enumerate(data.roi_list)
            if number + 1 != data.roi_glass_id and len(roi.values) >= 6]


def in_si(geometry, point):
    """Geometry of a cell at a point (see plant_tools.cell_geometry) in m
    and 1/m."""
    values = {"L_major": geometry["L_major"] * 1e-6,
              "L_minor": geometry["L_minor"] * 1e-6, "f": point["f"],
              "dome": point["dome"]}
    for name in ("k_major", "k_minor", "kM"):
        values[name] = point[name] * 1e6
    return values


def cells(data, pixel=None):
    """Cell of the curves of a file (of one curve with pixel=(i, j)):
    {(i, j): (geometry at the curve (see in_si), map data, number of the
    ROI from 0)}. A single curve is in the cell (ROI) of the map where it
    was recorded, the curves of a map in the cells of the map."""
    found = {}
    if data.is_single:
        on_map = plant_tools.map_of_curve(data)
        if on_map is None:
            return found
        _, map_data, position = on_map
        for number, roi in cell_rois(map_data):
            if [position["i"], position["j"]] in \
                    [[int(i), int(j)] for i, j in roi.values]:
                geometry = plant_tools.cell_geometry(
                    map_data, roi.values, [[position["x"], position["y"]]])
                if geometry is not None:
                    found[(0, 0)] = (in_si(geometry, geometry["points"][1]),
                                     map_data, number)
                break
        return found
    size_x = data.scan_size_x / data.nbr_pixels_x
    size_y = data.scan_size_y / data.nbr_pixels_y
    for number, roi in cell_rois(data):
        pixels = [(int(i), int(j)) for i, j in roi.values]
        if pixel is not None:
            if tuple(pixel) not in pixels:
                continue
            pixels = [tuple(pixel)]
        points = [[(i + 0.5) * size_x, (j + 0.5) * size_y]
                  for i, j in pixels]
        geometry = plant_tools.cell_geometry(data, roi.values, points)
        if geometry is None:
            continue
        for key, point in zip(pixels, geometry["points"][1:]):
            found.setdefault(key, (in_si(geometry, point), data, number))
        if pixel is not None:
            break
    return found


def elasticity(data, i, j):
    """Elasticity of the Hertz fit of a curve [Pa] (first segment), nan
    without value."""
    array = data.stiffness_array
    if array is None or not len(array[i][j]):
        return float("nan")
    return float(array[i][j][0])


def cell_mechanics(data, values, cell):
    """Values of a curve in a cell (see CELL_VALUES), from its k_as and E
    Hertz (values, completed), the cell (see cells) and the thickness of
    the wall of the file."""
    geometry, _, number = cell
    thickness = data.plant_wall_thickness * 1e-9
    poisson_ratio = data.used_poisson_ratio
    values["f"] = geometry["f"]
    values["roi"] = number + 1
    values["modulus"], values["pressure"], values["tau"] = \
        plant_tools.pressure_yamasaki(values["k_as"], geometry, thickness,
                                      poisson_ratio)
    values["pressure_hertz"], values["tau_hertz"] = \
        plant_tools.pressure_beauzamy(values["k_as"],
                                      values["modulus_hertz"], geometry,
                                      thickness, poisson_ratio)


def curve(data, i, j, used=True):
    """Apparent stiffness of a curve (see plant_tools.apparent_stiffness)
    with its values (dict, see VALUES) in result["values"], and its cell in
    result["cell"] (see cells, None outside of a cell). None without point
    of contact."""
    current = settings(data, used)
    result = plant_tools.apparent_stiffness(
        data, i, j, current["method"], current["smoothing"] / 100.0,
        current["force_start"], current["force_stop"])
    if result is None:
        return None
    values = dict.fromkeys(VALUES, float("nan"))
    values["k_as"] = result["k_as"]
    for name in ("approach", "retraction"):
        if result[name] is not None:
            values["k_" + name] = result[name]["k_as"]
            values["r2_" + name] = result[name]["r2"]
    values["modulus_hertz"] = elasticity(data, i, j)
    cell = cells(data, (i, j)).get((i, j))
    if cell is not None:
        cell_mechanics(data, values, cell)
    result["values"] = values
    result["cell"] = cell
    return result


def rois_key(data):
    """What the cells of a map depend on: its ROIs, its substrate and its
    heights (see plant_tools.map_heights)."""
    return (tuple(tuple(tuple(int(k) for k in pixel) for pixel in roi.values)
                  for roi in data.roi_list),
            data.roi_glass_id, data.stiffness_calculated,
            id(data.stiffness_array),
            data.flatten_applied and data.applied_flatten_motor)


def cell_values(data):
    """Values of the cells (see CELL_VALUES) of the curves of a file:
    {name: array [i][j]} (nan outside of a cell). Computed again when the
    ROIs, the wall thickness or the curves change."""
    nx, ny = data.nbr_pixels_x, data.nbr_pixels_y
    if not data.plant_calculated or data.plant_results is None:
        return {name: numpy.full((nx, ny), numpy.nan)
                for name in CELL_VALUES}
    if data.is_single:
        on_map = plant_tools.map_of_curve(data)
        maps = None if on_map is None else \
            (on_map[0], on_map[2]["i"], on_map[2]["j"], rois_key(on_map[1]))
    else:
        maps = rois_key(data)
    key = (maps, data.plant_wall_thickness, data.used_poisson_ratio,
           id(data.plant_results),
           numpy.asarray(data.discarded_curves).tobytes())
    cached = _CELLS.get(id(data))
    if cached is not None and cached[0] == key:
        return cached[1]

    values = {name: numpy.full((nx, ny), numpy.nan) for name in CELL_VALUES}
    for (i, j), cell in cells(data).items():
        if data.discarded_curves[i][j] != 0:
            continue
        curve_values = dict(zip(CURVE_VALUES,
                                (float(v) for v in data.plant_results[i, j])))
        if not numpy.isfinite(curve_values["k_as"]):
            continue
        cell_mechanics(data, curve_values, cell)
        for name in CELL_VALUES:
            values[name][i, j] = curve_values[name]
    _CELLS[id(data)] = (key, values)
    return values


def compute(calc_id):
    """Plant cells for the curves of a file. Returns the numbers of curves
    {"curves": computed, "k_as": with k_as, "cell": in a cell (with k_as),
    "pressure": with the pressure (Yamasaki), "saddle": in a cell but not on
    a dome}."""
    data = shared.exp.list[calc_id]
    misc_tools.save_used_params(data, "plant")
    nx, ny = data.nbr_pixels_x, data.nbr_pixels_y
    results = numpy.full((nx, ny, len(CURVE_VALUES)), numpy.nan,
                         dtype=numpy.float32)
    counts = dict.fromkeys(("curves", "k_as", "cell", "pressure",
                            "saddle"), 0)

    Progressbar(os.path.basename(data.filename))
    widgets_list.widget_progressbar.set_label("Plant cells")
    widgets_list.widget_progressbar.set_range(0, nx * ny)
    current = settings(data, used=True)

    def curves():
        # Force vs indentation of the curves (read in this process), fitted
        # in the processes of the pool
        for i in range(nx):
            for j in range(ny):
                if data.discarded_curves[i][j] != 0:
                    widgets_list.widget_progressbar.update()
                    continue
                counts["curves"] += 1
                parts = plant_tools.force_indentation(data, i, j)
                if parts is None:
                    widgets_list.widget_progressbar.update()
                    continue
                yield i, j, parts

    def store(found):
        for i, j, values in found:
            results[i, j] = values + [elasticity(data, i, j)]
            counts["k_as"] += bool(numpy.isfinite(values[0]))
        widgets_list.widget_progressbar.update(len(found))

    options = {"method": current["method"],
               "window": current["smoothing"] / 100.0,
               "start": current["force_start"], "stop": current["force_stop"]}
    multiproc.compute_curves(plant_tools.stiffness_task, options, curves(),
                             store, nx * ny)

    data.plant_results = results
    data.plant_calculated = True
    found = cells(data)
    pressure = cell_values(data)["pressure"]
    for (i, j), cell in found.items():
        if data.discarded_curves[i][j] == 0 and \
                numpy.isfinite(results[i, j, 0]):
            counts["cell"] += 1
            counts["pressure"] += bool(numpy.isfinite(pressure[i, j]))
            counts["saddle"] += not cell[0]["dome"]
    data.update()
    widgets_list.widget_progressbar.close()
    return counts


def refresh_results():
    """Maps, curve and tables of the results after a change of the cells
    (ROIs) or of the wall thickness: the values of the cells change without
    a new computation."""
    results = widgets_list.widget_results
    if results is None or not any(data.plant_calculated
                                  for data in shared.exp.list):
        return
    data = shared.exp.current_data
    if data.meshgrid_type.startswith("plant_"):
        results.update_MPL("MPL_canvas1")
    if data.display_curve_type == "plant":
        results.update_MPL("MPL_canvas2")
    # Values of the plots tabs: computed again when they are displayed
    main = widgets_list.widget_main
    if main is not None:
        if main.tabs.currentIndex() in (3, 4):
            widgets_list.widget_results_single.schedule_refresh()
        else:
            main.update_results = True
