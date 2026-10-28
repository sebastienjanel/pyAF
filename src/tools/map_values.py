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

"""Values of the maps of the results tab as 2D arrays ([x][y], raw units,
nan without value), their labels and display factors, and the maps computed
for a file (3D view)."""

import numpy

from . import planefit

# Name, color bar label and factor (displayed value = raw value * factor) of
# the maps that can be drawn as a surface (not the slices of the elasticity)
MAPS = {
    "piezo": ("Piezo height", "Piezo height [nm]", 1.0),
    "topo": ("Topography (zero force)", "Zero-force topography [nm]", 1.0),
    "stiffness": ("Elasticity", "E [kPa]", 1e-3),
    "stiffness_corr": ("Elasticity (BEC)", "E [kPa]", 1e-3),
    "work": ("Detachment work", "Detachment work [fJ]", 1e15),
    "rupture_force": ("Detachment force", "Detachment force [pN]", 1e12),
    "events_per_curve": ("Events per curve", "Events per curve", 1.0),
    "events_rupture_force": ("Event max. force", "Force [pN]", 1e12),
    "visco_e0": ("E₀", "E₀ [kPa]", 1e-3),
    "visco_fluidity": ("Fluidity", "Fluidity", 1.0),
    "visco_tan_delta": ("Loss tangent", "tan δ", 1.0),
    "visco_storage": ("Storage modulus E′", "E′ [kPa]", 1e-3),
    "visco_loss": ("Loss modulus E″", "E″ [kPa]", 1e-3),
    "breakthrough_force": ("Breakthrough force", "Fb [nN]", 1e9),
    "breakthrough_thickness": ("Breakthrough jump", "Jump [nm]", 1.0),
    "breakthrough_count": ("Breakthroughs per curve",
                           "Breakthroughs per curve", 1.0),
    "fit_r2": ("Fit R²", "R² of the fit", 1.0),
    "fit_poc_fixed": ("Fit: fixed point of contact",
                      "Point of contact fixed (1)", 1.0),
}

# Elasticity maps (one value per indentation segment with the Tomography)
ELASTICITY_MAPS = ("stiffness", "stiffness_corr")


def available_maps(data):
    """Maps computed for a file, in the order of MAPS."""
    computed = {
        "piezo": True,
        "topo": data.stiffness_calculated,
        "stiffness": data.stiffness_calculated,
        "stiffness_corr": data.stiffness_corrected,
        "work": data.work_and_rupture_force1_calculated,
        "rupture_force": data.work_and_rupture_force1_calculated,
        "events_per_curve": data.events_calculated,
        "events_rupture_force": data.events_calculated,
        "breakthrough_force": data.breakthrough_calculated,
        "breakthrough_thickness": data.breakthrough_calculated,
        "breakthrough_count": data.breakthrough_calculated,
        "fit_r2": data.has_fit_quality(),
        "fit_poc_fixed": data.has_fit_quality()}
    return [mt for mt in MAPS
            if computed.get(mt, mt.startswith("visco_") and
                            data.visco_calculated)]


def label(data, map_type):
    """Color bar label of a map (stiffness in N/m with the slope model)."""
    if map_type in ELASTICITY_MAPS and \
            data.used_stiffness_model_selected == 3:
        return "Apparent stiffness [N/m]"
    return MAPS[map_type][1]


def factor(data, map_type):
    """Factor of the displayed values (raw value * factor)."""
    if map_type in ELASTICITY_MAPS and \
            data.used_stiffness_model_selected == 3:
        return 1.0
    return MAPS[map_type][2]


def values(data, map_type):
    """Raw values of a map ([x][y], nan without value). The elasticity at the
    indentation segment of the 3D view (stiffness_depth_view)."""
    nx, ny = data.nbr_pixels_x, data.nbr_pixels_y
    if map_type == "piezo":
        return numpy.array(planefit.piezo_map(data), dtype=float)
    if map_type == "topo":
        return planefit.topo_map(data)
    if map_type in ELASTICITY_MAPS:
        array = data.stiffness_corrected_array \
            if map_type == "stiffness_corr" else data.stiffness_array
        depth = data.stiffness_depth_view
        result = numpy.full((nx, ny), numpy.nan)
        for i in range(nx):
            for j in range(ny):
                if len(array[i][j]) > depth:
                    result[i, j] = array[i][j][depth]
        return result
    sources = {"work": "work", "rupture_force": "rupture_force1",
               "events_per_curve": "events_per_curve",
               "events_rupture_force": "rupture_force2"}
    result = getattr(data, sources.get(map_type, map_type))
    return numpy.array(numpy.ma.filled(numpy.ma.asarray(result, dtype=float),
                                       numpy.nan), dtype=float)


def height(data):
    """Height of the 3D surface [nm]: piezo height or zero-force
    topography (opengl_height), with the planefit. Pixels without height
    (no point of contact) at the lowest level, as VTK needs finite
    coordinates."""
    use_topo = getattr(data, "opengl_height", "topo") == "topo" and \
        data.stiffness_calculated
    array = values(data, "topo" if use_topo else "piezo")
    finite = numpy.isfinite(array)
    array[~finite] = array[finite].min() if finite.any() else 0.0
    return array
