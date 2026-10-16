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


"""Planefit of the height maps (piezo height and zero-force topography).

The substrate ROI of a file (ROI manager, "Define as substrate") is used as
the zero plane: a plane (offset and tilt) is fitted on the values of the map
in the ROI and subtracted from the whole map. Only the displayed maps are
corrected, the data used for the computations are not changed.
"""

import numpy
from . import math_tools


def get_planefit_roi(data):
    """Substrate ROI used for the planefit of a file (None if the planefit
    is not applied, or if the ROI is no longer a substrate)."""
    if not getattr(data, "planefit_applied", False):
        return None
    index = getattr(data, "planefit_roi", None)
    if index is None or not 0 <= index < len(data.roi_list):
        return None
    roi = data.roi_list[index]
    if roi.glass_coeffs is None or not roi.values:
        return None
    return roi


def get_coeffs(data, array, roi):
    """Coefficients of the plane fitted on the values of array in the ROI
    (None if there are less than 3 values)."""
    x_size = data.scan_size_x / data.nbr_pixels_x
    y_size = data.scan_size_y / data.nbr_pixels_y
    pixels = numpy.array(roi.values, dtype=int)
    z = numpy.asarray(array, dtype=float)[pixels[:, 0], pixels[:, 1]]
    finite = numpy.isfinite(z)
    if numpy.count_nonzero(finite) < 3:
        return None
    x = pixels[finite, 0] * x_size
    y = pixels[finite, 1] * y_size
    return math_tools.polyfit2d(x, y, z[finite], order=1, linear=True)


def get_plane(data, array, roi):
    """Plane fitted on the values of array in the ROI, over the whole map
    (None if there are less than 3 values)."""
    coeffs = get_coeffs(data, array, roi)
    if coeffs is None:
        return None
    x_size = data.scan_size_x / data.nbr_pixels_x
    y_size = data.scan_size_y / data.nbr_pixels_y
    grid_x, grid_y = numpy.meshgrid(
        numpy.arange(array.shape[0]) * x_size,
        numpy.arange(array.shape[1]) * y_size, indexing="ij")
    return math_tools.polyval2d(grid_x, grid_y, coeffs)


def topo_coeffs(data):
    """Coefficients of the plane removed from the topography by the
    planefit (None without planefit)."""
    roi = get_planefit_roi(data)
    if roi is None:
        return None
    return get_coeffs(data, data.topography, roi)


def apply_planefit(data, array):
    """The map with the planefit of the file (a copy), or the map itself
    when no planefit is applied."""
    roi = get_planefit_roi(data)
    if roi is None:
        return array
    array = numpy.array(array, dtype=float)
    plane = get_plane(data, array, roi)
    if plane is None:
        return array
    return array - plane


def piezo_map(data, flatten=True):
    """Displayed piezo height (with the flatten of the compute tab when
    flatten is True, and the planefit)."""
    if flatten and data.flatten_applied:
        array = math_tools.flatten(data.piezo_image,
                                   data.applied_flatten_order)
    else:
        array = data.piezo_image
    return apply_planefit(data, numpy.array(array, dtype=float))


def topo_map(data):
    """Displayed zero-force topography (with the planefit)."""
    return apply_planefit(data, numpy.array(data.topography, dtype=float))


def reset_color_scales(data):
    """Color scales of the height maps after a change of the planefit: the
    range of the planefitted values, or the default scales without
    planefit."""
    if get_planefit_roi(data) is not None:
        for array, opts in ((piezo_map(data), data.color_opts_piezo),
                            (topo_map(data) if data.stiffness_calculated
                             else None, data.color_opts_topo)):
            if array is None or not numpy.isfinite(array).any():
                continue
            low, high = numpy.nanmin(array), numpy.nanmax(array)
            opts[1], opts[2], opts[3] = low, high, (low + high) / 2.0
    else:
        piezo = numpy.asarray(data.piezo_image, dtype=float)
        high = numpy.median(piezo) + 3 * numpy.std(piezo)
        opts = data.color_opts_piezo
        opts[1], opts[2], opts[3] = 0, high, high / 2.0
        if data.stiffness_calculated and data.max_topo is not None:
            opts = data.color_opts_topo
            opts[1], opts[2], opts[3] = 0, data.max_topo, data.max_topo / 2.0


def get_volume(data, roi):
    """Volume between the zero-force topography and the plane of a substrate
    ROI (fitted when the ROI was defined as substrate), over the whole map.

    Returns a dict (volumes in nm^3): the net volume (heights above and
    below the plane), the volumes above and below the plane, the number of
    pixels used and the number of pixels without topography (no contact
    point: not counted).
    """
    x_size = data.scan_size_x / data.nbr_pixels_x
    y_size = data.scan_size_y / data.nbr_pixels_y
    topo = numpy.array(data.topography, dtype=float)
    grid_x, grid_y = numpy.meshgrid(
        numpy.arange(topo.shape[0]) * x_size,
        numpy.arange(topo.shape[1]) * y_size, indexing="ij")
    plane = math_tools.polyval2d(grid_x, grid_y,
                                 numpy.asarray(roi.glass_coeffs, dtype=float))
    heights = topo - plane
    finite = numpy.isfinite(heights)
    heights = heights[finite]
    area = x_size * y_size
    return {"net": float(heights.sum() * area),
            "above": float(heights[heights > 0].sum() * area),
            "below": float(heights[heights < 0].sum() * area),
            "pixels": int(heights.size),
            "no_topography": int((~finite).sum()),
            "max_height": float(heights.max()) if heights.size else numpy.nan,
            "pixel_size": (x_size, y_size)}
