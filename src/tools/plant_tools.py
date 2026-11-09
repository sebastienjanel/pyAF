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

"""Plant cells: position of the force curves on a map, and geometry of the
cells (ROIs drawn on the topography).

Positions (JPK files). The microscope records the position of the tip
scanner, of the sample scanner and of the motorized stage. The position of a
curve on the sample is

    position = tip scanner - stage - sample scanner

(missing terms are 0, e.g. without motorized stage), so that single curves
recorded after a map can be placed on it, as long as the stage did not move.
Only files of the same microscope, recorded less than SAME_SESSION hours
apart, are compared (unrelated files without motorized stage would be at
the same positions). A sample moved by hand can not be detected.

Cell geometry (Yamasaki et al. 2025, Planta 261:131; Beauzamy et al. 2015,
Biophys J 108:2448): the axes of the cell are the principal axes of its ROI,
and the half-lengths L are the distances from the centre to the edge of the
ROI along these axes (Yamasaki Fig. 2b). The height of the ROI is fitted with
a quadratic surface z = a x² + b y² + c xy + d x + e y + g, and the curvatures
along the axes, at a point, are the normal curvatures of this surface
(positive for a dome). The mean and Gaussian curvatures are
kM = (k_major + k_minor) / 2 and kG = k_major k_minor (Yamasaki Eq. 13-14),
and f is the structural factor (Eq. 15, 1 for a sphere).
"""

import datetime
import math

import numpy

from .. import shared

# Distance between the stages of two files to consider them at the same
# place [m]
SAME_STAGE = 1e-7
# Maximum time between a map and a single curve recorded on it [h]
SAME_SESSION = 12.0


def frame_from_header(header):
    """Positions of the stage and of the scanners in a JPK header (dict of
    strings): {"tip", "sample", "stage": [x, y] in m or None, "z_motor": m
    or None, "coordinates": name of the frame of the grid or None}."""

    def pair(suffix_x, suffix_y):
        x = y = None
        for key, value in header.items():
            if x is None and key.endswith(suffix_x):
                x = value
            if y is None and key.endswith(suffix_y):
                y = value
        if x in (None, "") or y in (None, ""):
            return None
        return [float(x), float(y)]

    frame = {
        "tip": pair("tip-scanner.scanner-displacement.x",
                    "tip-scanner.scanner-displacement.y"),
        "sample": pair("sample-scanner.scanner-displacement.x",
                       "sample-scanner.scanner-displacement.y"),
        "stage": pair("motor-stage-moved.parameters.transform.translatex",
                      "motor-stage-moved.parameters.transform.translatey"),
        "z_motor": None, "coordinates": None}
    for key, value in header.items():
        if frame["z_motor"] is None and \
                key.endswith("z-scanner-map.z-motors.position"):
            frame["z_motor"] = float(value)
        if frame["coordinates"] is None and \
                key.endswith("position-pattern.coordinates.name"):
            frame["coordinates"] = value
    return frame


def sample_position(tip, frame):
    """Position [x, y] (m) on the sample of a point at the position tip of
    the tip scanner (see the module)."""
    x, y = tip
    for name in ("stage", "sample"):
        if frame.get(name) is not None:
            x -= frame[name][0]
            y -= frame[name][1]
    return [x, y]


def map_grid(grid, frame):
    """Grid of a map (xcenter, ycenter, ulength, vlength, theta in m and
    rad) on the sample (see sample_position)."""
    grid = dict(grid)
    if frame.get("coordinates") != "surface":
        grid["xcenter"], grid["ycenter"] = sample_position(
            [grid["xcenter"], grid["ycenter"]], frame)
    return grid


def same_stage(frame_1, frame_2):
    """True if the stages of the two files are at the same position (or if
    both files have no stage)."""
    stage_1, stage_2 = frame_1.get("stage"), frame_2.get("stage")
    if stage_1 is None or stage_2 is None:
        return stage_1 is None and stage_2 is None
    return math.hypot(stage_1[0] - stage_2[0],
                      stage_1[1] - stage_2[1]) < SAME_STAGE


def same_session(data_1, data_2):
    """True if two files were recorded with the same microscope, less than
    SAME_SESSION hours apart (the positions of files without motorized stage
    can not be compared otherwise)."""
    if not data_1.microscope_name or \
            data_1.microscope_name != data_2.microscope_name:
        return False
    try:
        delay = abs((datetime.datetime.fromisoformat(str(data_1.date)) -
                     datetime.datetime.fromisoformat(str(data_2.date)))
                    .total_seconds())
    except ValueError:
        return False
    return delay <= SAME_SESSION * 3600.0


def locate(curve_data, map_data):
    """Position of a single curve on a map: {"x", "y": position on the map
    in nm (as on the plots, from the corner), "i", "j": pixel}, or None if
    the curve is not on the map (other microscope or session, other stage
    position, outside of the map, rotated map, or positions not known)."""
    position = getattr(curve_data, "sample_position", None)
    grid = getattr(map_data, "sample_grid", None)
    if position is None or grid is None or curve_data is map_data:
        return None
    if not same_session(curve_data, map_data):
        return None
    if not same_stage(curve_data.frame or {}, map_data.frame or {}):
        return None
    if abs(grid.get("theta", 0.0)) > 1e-9:
        return None
    x = (position[0] - grid["xcenter"] + grid["ulength"] / 2) * 1e9
    y = (position[1] - grid["ycenter"] + grid["vlength"] / 2) * 1e9
    if not (0 <= x < map_data.scan_size_x and 0 <= y < map_data.scan_size_y):
        return None
    i = int(x / map_data.scan_size_x * map_data.nbr_pixels_x)
    j = int(y / map_data.scan_size_y * map_data.nbr_pixels_y)
    return {"x": x, "y": y, "i": i, "j": j}


def curves_on_map(map_data):
    """Single curves loaded on a map: [(index of the file, data, position
    (see locate))]."""
    located = []
    for index, data in enumerate(shared.exp.list):
        if data.is_single:
            position = locate(data, map_data)
            if position is not None:
                located.append((index, data, position))
    return located


def map_of_curve(curve_data):
    """First loaded map on which a single curve is: (index of the file,
    data, position (see locate)), or None."""
    for index, data in enumerate(shared.exp.list):
        if not data.is_single:
            position = locate(curve_data, data)
            if position is not None:
                return index, data, position
    return None


def map_heights(data):
    """Heights of a map [nm] for the cell geometry: the zero-force
    topography when the elasticity is computed, else the piezo height (at
    the setpoint), without flatten. (array [x][y], name)."""
    if data.stiffness_calculated and data.topography is not None:
        return numpy.array(data.topography, dtype=float), \
            "zero-force topography"
    return numpy.array(data.piezo_image, dtype=float), "piezo height"


def structural_factor(kappa_m, kappa_g):
    """Structural factor f of an ellipsoidal shell (Yamasaki Eq. 15,
    Beauzamy Eq. 1), 1 for a sphere."""
    ratio = kappa_g / kappa_m ** 2
    if abs(ratio) < 1e-12:
        # Limit of a cylinder (kG = 0)
        return 0.75
    root = math.sqrt(max(1.0 - ratio, 0.0))
    return (2.0 + (root - 1.0) / ratio) / (root + 1.0)


def cell_geometry(data, pixels, points=None):
    """Geometry of a cell from the pixels [[i, j], ...] of its ROI.

    points: positions [[x, y], ...] (nm, as on the plots) where the
    curvatures are computed, after the centre of the ROI.

    Returns a dict (lengths in µm, curvatures in 1/µm): pixels, area,
    centre, angle (of the major axis, degrees from x), L_major, L_minor
    (half-lengths), rms (of the fit, nm), heights (name of the map), and
    points: for each point, x, y, k_major, k_minor, kM, kG, f, k1, k2
    (principal curvatures), dome (both principal curvatures positive). None
    if there are less than 6 pixels."""
    pixels = numpy.array(pixels, dtype=int).reshape(-1, 2)
    if len(pixels) < 6:
        return None
    heights, name = map_heights(data)
    size_x = data.scan_size_x / data.nbr_pixels_x / 1000.0
    size_y = data.scan_size_y / data.nbr_pixels_y / 1000.0

    # Centres of the pixels [µm]
    x = (pixels[:, 0] + 0.5) * size_x
    y = (pixels[:, 1] + 0.5) * size_y
    z = heights[pixels[:, 0], pixels[:, 1]] / 1000.0
    inside = data.inside_region()
    used = numpy.isfinite(z)
    if inside is not None:
        used &= inside[pixels[:, 0], pixels[:, 1]]
    x, y, z = x[used], y[used], z[used]
    if len(x) < 6:
        return None
    centre = numpy.array([x.mean(), y.mean()])

    # Axes of the cell: principal axes of the pixels of the ROI
    values, vectors = numpy.linalg.eigh(numpy.cov(numpy.vstack(
        [x - centre[0], y - centre[1]])))
    major = vectors[:, numpy.argmax(values)]
    if major[0] < 0:
        major = -major
    minor = numpy.array([-major[1], major[0]])

    # Half-lengths: from the centre to the edge of the ROI, mean of the two
    # directions
    roi = numpy.zeros((data.nbr_pixels_x, data.nbr_pixels_y), bool)
    roi[pixels[:, 0], pixels[:, 1]] = True
    step = min(size_x, size_y) / 4.0

    def half_length(direction):
        lengths = []
        for sign in (1, -1):
            distance = 0.0
            while True:
                point = centre + sign * (distance + step) * direction
                i = int(math.floor(point[0] / size_x))
                j = int(math.floor(point[1] / size_y))
                if not (0 <= i < roi.shape[0] and 0 <= j < roi.shape[1]) \
                        or not roi[i, j]:
                    break
                distance += step
            lengths.append(distance + step / 2.0)
        return sum(lengths) / 2.0

    # Quadratic surface, in µm, around the centre
    dx, dy = x - centre[0], y - centre[1]
    matrix = numpy.column_stack([dx ** 2, dy ** 2, dx * dy, dx, dy,
                                 numpy.ones_like(dx)])
    coeffs, *_ = numpy.linalg.lstsq(matrix, z, rcond=None)
    a, b, c, d, e, _ = coeffs
    rms = float(numpy.sqrt(numpy.mean((matrix @ coeffs - z) ** 2))) * 1000.0
    hessian = numpy.array([[2 * a, c], [c, 2 * b]])

    points = [[centre[0] * 1000.0, centre[1] * 1000.0]] + list(points or [])
    results = []
    for point in points:
        px, py = point[0] / 1000.0 - centre[0], point[1] / 1000.0 - centre[1]
        gradient = numpy.array([2 * a * px + c * py + d,
                                c * px + 2 * b * py + e])
        norm = 1.0 + gradient @ gradient

        def normal_curvature(direction):
            # Normal curvature along the direction (unit, in the plane of
            # the map), positive for a dome
            slope = gradient @ direction
            return float(-(direction @ hessian @ direction) /
                         (math.sqrt(norm) * (1.0 + slope ** 2)))

        k_major = normal_curvature(major)
        k_minor = normal_curvature(minor)
        kappa_m = (k_major + k_minor) / 2.0
        kappa_g = k_major * k_minor
        # Principal curvatures (shape operator of the surface)
        metric = numpy.eye(2) + numpy.outer(gradient, gradient)
        shape = numpy.linalg.solve(metric, -hessian / math.sqrt(norm))
        k1, k2 = sorted(numpy.real(numpy.linalg.eigvals(shape)),
                        reverse=True)
        results.append({
            "x": float(point[0]), "y": float(point[1]),
            "k_major": k_major, "k_minor": k_minor, "kM": kappa_m,
            "kG": kappa_g,
            "f": structural_factor(kappa_m, kappa_g)
            if kappa_m != 0 else float("nan"),
            "k1": float(k1), "k2": float(k2),
            "dome": bool(k1 > 0 and k2 > 0)})

    return {"pixels": int(len(x)),
            "area": float(len(pixels) * size_x * size_y),
            "centre": [float(centre[0]), float(centre[1])],
            "angle": float(math.degrees(math.atan2(major[1], major[0]))),
            "L_major": half_length(major), "L_minor": half_length(minor),
            "rms": rms, "heights": name, "points": results}
