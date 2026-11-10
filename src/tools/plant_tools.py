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


# Apparent stiffness of the cell (slope of the force vs the indentation at
# large indentation): Yamasaki et al. 2025 (Eq. 16), the derivative of the
# force is fitted with k = k_as (1 - exp(-alpha (delta - beta))); Beauzamy
# et al. 2015, linear fit of the force between 75 and 100 % of the maximum.
STIFFNESS_METHODS = ("Derivative (Yamasaki)", "Force range (Beauzamy)")


def force_indentation(data, i, j):
    """Indentation [m] and force [N] of the loading (approach) and of the
    unloading (retraction) of a curve, from its point of contact
    (elasticity computed): ((delta, force) of the loading, (delta, force)
    of the unloading), each in the order of the indentation, or None."""
    from ..compute_tools.viscoelasticity.main import prepare_curve
    if not data.stiffness_calculated:
        return None
    curve = prepare_curve(data, i, j, {"drag_correction": False},
                          whole=True)
    if curve is None:
        return None
    delta, force, k_tm = curve
    parts = []
    for part in (slice(0, k_tm + 1), slice(k_tm, None)):
        d, f = numpy.asarray(delta[part]), numpy.asarray(force[part])
        inside = numpy.isfinite(d) & numpy.isfinite(f) & (d > 0)
        d, f = d[inside], f[inside]
        order = numpy.argsort(d, kind="stable")
        parts.append((d[order], f[order]))
    if min(len(part[0]) for part in parts) < 10:
        return None
    return parts


def derivative(delta, force, nbr_points=200, window=0.05):
    """Derivative dF/ddelta [N/m] on a regular grid of the indentation
    (mean of the force in nbr_points bins, then a Savitzky-Golay filter of
    the given width, fraction of the indentation range): (delta, k)."""
    from scipy.signal import savgol_filter
    edges = numpy.linspace(delta[0], delta[-1], nbr_points + 1)
    index = numpy.clip(numpy.searchsorted(edges, delta, side="right") - 1,
                       0, nbr_points - 1)
    counts = numpy.bincount(index, minlength=nbr_points)
    sums = numpy.bincount(index, weights=force, minlength=nbr_points)
    filled = counts > 0
    centres = (edges[:-1] + edges[1:]) / 2.0
    grid, values = centres[filled], sums[filled] / counts[filled]
    if len(grid) < 7:
        return None
    length = max(5, int(round(window * len(grid))) | 1)
    length = min(length, len(grid) - (1 - len(grid) % 2))
    step = (grid[-1] - grid[0]) / (len(grid) - 1)
    slope = savgol_filter(numpy.interp(
        numpy.linspace(grid[0], grid[-1], len(grid)), grid, values),
        length, 2, deriv=1, delta=step)
    return numpy.linspace(grid[0], grid[-1], len(grid)), slope


# End of the indentation not fitted (fraction of the range): the derivative
# is not reliable where the piezo turns back (edge of the filter)
END_EXCLUDED = 0.03


def stiffness_derivative(delta, force, window=0.05):
    """Apparent stiffness from the derivative of the force (Yamasaki Eq.
    16, with k = 0 before beta: the point of contact can be found before
    the force rises): {"k_as" [N/m], "alpha" [1/m], "beta" [m], "r2",
    "delta", "k" (derivative), "fit" (fitted derivative)}, or None."""
    from scipy.optimize import curve_fit
    curve = derivative(delta, force, window=window)
    if curve is None:
        return None
    grid, slope = curve
    kept = grid <= grid[0] + (grid[-1] - grid[0]) * (1.0 - END_EXCLUDED)
    grid, slope = grid[kept], slope[kept]
    if len(grid) < 7:
        return None

    def model(d, k_as, alpha, beta):
        return k_as * (1.0 - numpy.exp(-alpha * numpy.maximum(d - beta,
                                                               0.0)))

    top = slope[int(len(slope) * 0.75):]
    k0 = float(numpy.median(top)) if len(top) else float(slope.max())
    if not k0 > 0:
        return None
    rising = numpy.nonzero(slope > 0.1 * k0)[0]
    beta0 = float(grid[rising[0]]) if len(rising) else 0.0
    guess = [k0, 3.0 / (grid[-1] - beta0), beta0]
    try:
        params, _ = curve_fit(
            model, grid, slope, p0=guess,
            bounds=([0.0, 0.0, 0.0], [numpy.inf, numpy.inf, grid[-1]]),
            maxfev=5000)
    except (RuntimeError, ValueError):
        return None
    fit = model(grid, *params)
    residual = numpy.sum((slope - fit) ** 2)
    total = numpy.sum((slope - slope.mean()) ** 2)
    return {"k_as": float(params[0]), "alpha": float(params[1]),
            "beta": float(params[2]),
            "r2": float(1.0 - residual / total) if total > 0 else
            float("nan"), "delta": grid, "k": slope, "fit": fit}


def stiffness_force_range(delta, force, start=75.0, stop=100.0):
    """Apparent stiffness: slope of a linear fit of the force between start
    and stop % of its maximum (Beauzamy et al. 2015): {"k_as" [N/m], "r2",
    "delta", "force" (fitted part), "fit"}, or None."""
    maximum = force.max()
    used = (force >= maximum * start / 100.0) & \
        (force <= maximum * stop / 100.0)
    if used.sum() < 5:
        return None
    d, f = delta[used], force[used]
    slope, intercept = numpy.polyfit(d, f, 1)
    fit = slope * d + intercept
    total = numpy.sum((f - f.mean()) ** 2)
    return {"k_as": float(slope),
            "r2": float(1.0 - numpy.sum((f - fit) ** 2) / total)
            if total > 0 else float("nan"),
            "delta": d, "force": f, "fit": fit}


def stiffness_of_curves(curves, method=0, window=0.05, start=75.0,
                        stop=100.0):
    """Apparent stiffness of the loading and of the unloading (curves, see
    force_indentation): {"approach", "retraction", "k_as"} (see
    apparent_stiffness)."""
    results = {}
    for name, (delta, force) in zip(("approach", "retraction"), curves):
        if method == 0:
            results[name] = stiffness_derivative(delta, force, window)
        else:
            results[name] = stiffness_force_range(delta, force, start, stop)
    found = [results[name]["k_as"] for name in ("approach", "retraction")
             if results[name] is not None]
    results["k_as"] = float(numpy.mean(found)) if found else float("nan")
    return results


def stiffness_task(task):
    """Apparent stiffness of curves, in a process of the computation (see
    multiproc.compute_curves). task: settings of stiffness_of_curves
    ("common") and list of (i, j, curves) ("curves"). Returns the list of
    (i, j, values): k_as, k and R^2 of the approach and of the
    retraction (nan when not found)."""
    results = []
    for i, j, curves in task["curves"]:
        result = stiffness_of_curves(curves, **task["common"])
        values = [result["k_as"]]
        for key in ("k_as", "r2"):
            for name in ("approach", "retraction"):
                part = result[name]
                values.append(part[key] if part else float("nan"))
        results.append((i, j, values))
    return results


def apparent_stiffness(data, i, j, method=0, window=0.05, start=75.0,
                       stop=100.0):
    """Apparent stiffness of the curve (i, j), loading and unloading (see
    STIFFNESS_METHODS): {"approach": result or None, "retraction": result
    or None, "k_as": mean of the two (or the one found), "curves":
    (loading, unloading) from force_indentation}, or None without point of
    contact."""
    curves = force_indentation(data, i, j)
    if curves is None:
        return None
    results = stiffness_of_curves(curves, method, window, start, stop)
    results["curves"] = curves
    return results


# Turgor pressure and modulus of the cell wall, shell of thickness t and
# mean curvature kM indented at a point (Vella et al. 2012; Beauzamy et al.
# 2015; Yamasaki et al. 2025, Eq. 1):
#
#     k_as = pi P f / (kM log(P f c / (E t^2 kM^2))),  c = sqrt(3 (1 - nu^2))
#
# valid for a large dimensionless pressure tau = P f c / (2 E t^2 kM^2)
# (Eq. 2, tau >> 1).
# Yamasaki: the strain of the wall, from its curvature (Eq. 9), gives the
# second equation (Eq. 6, P = E t kM (e_major + e_minor) / f), and E and P
# are explicit (Eq. 10-11).
# Beauzamy: E is the elasticity of the Hertz fit (shallow indentation, below
# the thickness of the wall), and P is the root of Eq. 1.


def wall_strain(half_length, curvature):
    """Strain of the cell wall along an axis (Yamasaki Eq. 9), from its
    half-length L [m] and its curvature [1/m]: 2 (asin(L k) / (L k) - 1),
    nan if L k is not in ]0, 1]."""
    product = half_length * curvature
    if not 0 < product <= 1:
        return float("nan")
    return 2.0 * (math.asin(product) / product - 1.0)


def pressure_yamasaki(k_as, geometry, thickness, poisson_ratio):
    """Modulus E [Pa], pressure P [Pa] and dimensionless pressure tau of the
    wall (Yamasaki Eq. 10, 11 and 2) from the apparent stiffness [N/m], the
    geometry of the cell at the curve (lengths [m], curvatures [1/m]:
    L_major, L_minor, k_major, k_minor, kM, f) and the thickness of the
    wall [m]. nan without solution (saddle, flat or too curved wall)."""
    nan = (float("nan"),) * 3
    kappa_m, f = geometry["kM"], geometry["f"]
    strain = wall_strain(geometry["L_major"], geometry["k_major"]) + \
        wall_strain(geometry["L_minor"], geometry["k_minor"])
    if not (kappa_m > 0 and f > 0 and numpy.isfinite(strain) and
            strain > 0 and k_as > 0):
        return nan
    c = math.sqrt(3.0 * (1.0 - poisson_ratio ** 2))
    argument = strain * c / (thickness * kappa_m)
    if argument <= 1:
        return nan
    log = math.log(argument)
    modulus = k_as * log / (math.pi * thickness * strain)
    pressure = k_as * kappa_m * log / (math.pi * f)
    return modulus, pressure, argument / 2.0


def pressure_beauzamy(k_as, modulus, geometry, thickness, poisson_ratio):
    """Pressure P [Pa] and dimensionless pressure tau, root of Eq. 1 with
    the modulus E [Pa] of the wall (Hertz fit). With u = P f c / (E t^2
    kM^2) = 2 tau, Eq. 1 is u / log(u) = k_as c / (pi E t^2 kM), solved on
    the branch u > e (large pressures). nan without solution."""
    from scipy.special import lambertw
    nan = (float("nan"),) * 2
    kappa_m, f = geometry["kM"], geometry["f"]
    if not (kappa_m > 0 and f > 0 and k_as > 0 and modulus > 0):
        return nan
    c = math.sqrt(3.0 * (1.0 - poisson_ratio ** 2))
    ratio = k_as * c / (math.pi * modulus * thickness ** 2 * kappa_m)
    if ratio <= math.e:
        # u / log(u) >= e: k_as too small for this modulus
        return nan
    u = float(numpy.real(-ratio * lambertw(-1.0 / ratio, -1)))
    pressure = u * modulus * thickness ** 2 * kappa_m ** 2 / (f * c)
    return pressure, u / 2.0
