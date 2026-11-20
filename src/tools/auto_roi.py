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

"""ROI assistance: the pixels of a map sorted into the ROIs from a few
pixels marked by the user in each ROI (pixel classification, as
ilastik).

Each pixel is described by the maps computed for the file (topography,
elasticity, fluidity...), normalized on the map, and by their values
smoothed around it at a few scales in pixels (context, edges), so that the
same marks work on a cell and on a small scan of molecules. A random forest
learns the ROIs from the marked pixels and gives each pixel a probability
for each ROI; the pixels below the confidence threshold (borders, debris)
are left uncertain, in no ROI. Pixels without values (nan) are handled by
the forest, they do not block the sorting.

Cells that look alike (plant tissue) are not told apart by the forest: they
are sorted by the domes of the height instead (watershed: one region per
dome, the boundaries in the grooves between them), and each ROI takes the
whole domes its marks touch."""

import heapq

import numpy
from scipy import ndimage

from . import map_values

# Smoothing scales of the context [pixels]
SCALES = (1.0, 2.0, 4.0)
# Maps of moduli, spread over decades: sorted on their logarithm
LOG_MAPS = ("stiffness", "stiffness_corr", "visco_e0", "visco_storage",
            "visco_loss", "microrheo_storage", "microrheo_loss",
            "microrheo_crossover", "brush_modulus", "plant_k_as",
            "plant_modulus", "plant_pressure", "plant_pressure_hertz")
# Pixels marked at least, in at least two ROIs
MIN_MARKS = 2


class AutoRoiError(Exception):
    """The proposal can not be made (message for the user)."""


def available():
    """True if scikit-learn can be imported."""
    try:
        import sklearn.ensemble  # noqa: F401
    except ImportError:
        return False
    return True


def normalized(values, log=False):
    """Map centered on its median and divided by its spread (median
    absolute deviation, standard deviation if most pixels have the same
    value); None without values or without contrast."""
    values = numpy.array(values, dtype=float)
    if log:
        with numpy.errstate(divide="ignore", invalid="ignore"):
            values = numpy.where(values > 0, numpy.log10(values), numpy.nan)
    finite = values[numpy.isfinite(values)]
    if finite.size < 2:
        return None
    median = numpy.median(finite)
    spread = 1.4826 * numpy.median(numpy.abs(finite - median))
    if spread <= 0:
        spread = finite.std()
    if spread <= 0:
        return None
    return (values - median) / spread


def smoothed(values, sigma):
    """Gaussian smoothing ignoring the pixels without value (nan where no
    pixel with a value is near)."""
    finite = numpy.isfinite(values)
    weights = ndimage.gaussian_filter(finite.astype(float), sigma,
                                      mode="nearest")
    total = ndimage.gaussian_filter(numpy.where(finite, values, 0.0), sigma,
                                    mode="nearest")
    with numpy.errstate(divide="ignore", invalid="ignore"):
        return numpy.where(weights > 1e-3, total / weights, numpy.nan)


def channel_features(values):
    """Features of one normalized map: the values, smoothed at each scale,
    the slope (edges) and the local spread (texture)."""
    result = [values]
    for sigma in SCALES:
        result.append(smoothed(values, sigma))
    near = result[1]
    grad_x, grad_y = numpy.gradient(near)
    result.append(numpy.hypot(grad_x, grad_y))
    square = smoothed(values ** 2, SCALES[1])
    result.append(numpy.sqrt(numpy.clip(square - result[2] ** 2, 0, None)))
    return result


def features(data):
    """Features of the pixels of a file: array [x][y][feature], from all the
    maps computed for it."""
    shape = (data.nbr_pixels_x, data.nbr_pixels_y)
    channels = []
    for map_type in map_values.available_maps(data):
        values = normalized(map_values.values(data, map_type),
                            map_type in LOG_MAPS)
        if values is not None and values.shape == shape:
            channels.extend(channel_features(values))
    if not channels:
        raise AutoRoiError("No map with values for this file.")
    return numpy.stack(channels, axis=-1)


def marks_labels(shape, marks):
    """Labels of the marked pixels: array [x][y], index of the ROI in marks
    (list of sets of pixels (i, j)), -1 not marked. A pixel marked in two
    ROIs keeps the last one."""
    labels = numpy.full(shape, -1, dtype=int)
    for index, pixels in enumerate(marks):
        for i, j in pixels:
            if 0 <= i < shape[0] and 0 <= j < shape[1]:
                labels[i, j] = index
    return labels


def probabilities(feats, labels, inside=None):
    """Probability of each ROI for each pixel: array [x][y][ROI], learned on
    the marked pixels (labels, see marks_labels); nan outside of the region
    of the map (inside: boolean array, None for the whole map)."""
    from sklearn.ensemble import RandomForestClassifier

    nbr_classes = labels.max() + 1
    counts = numpy.bincount(labels[labels >= 0], minlength=nbr_classes)
    if numpy.count_nonzero(counts >= MIN_MARKS) < 2:
        raise AutoRoiError("Mark pixels in at least two ROIs (Draw column) "
                           "to get a proposal.")
    marked = labels >= 0
    forest = RandomForestClassifier(n_estimators=50, min_samples_leaf=2,
                                    class_weight="balanced", n_jobs=-1,
                                    random_state=0)
    forest.fit(feats[marked], labels[marked])
    shape = labels.shape
    valid = numpy.ones(shape, dtype=bool) if inside is None else \
        numpy.asarray(inside, dtype=bool)
    result = numpy.full(shape + (nbr_classes,), numpy.nan)
    if valid.any():
        proba = forest.predict_proba(feats[valid])
        # Columns of the ROIs marked (a ROI without marks: probability 0)
        values = numpy.zeros((proba.shape[0], nbr_classes))
        values[:, forest.classes_] = proba
        result[valid] = values
    return result


def classify(proba, labels, threshold, min_size):
    """ROI of each pixel: array [x][y], index of the ROI, -1 uncertain
    (probability below the threshold, or in a group of less than min_size
    pixels of the ROI), -2 outside of the region of the map. The marked
    pixels stay in their ROI."""
    outside = numpy.isnan(proba).all(axis=-1)
    filled = numpy.nan_to_num(proba, nan=0.0)
    result = numpy.argmax(filled, axis=-1)
    result[filled.max(axis=-1) < threshold] = -1
    if min_size > 1:
        for index in range(proba.shape[-1]):
            groups, number = ndimage.label(result == index)
            if number:
                sizes = numpy.bincount(groups.ravel())
                small = sizes < min_size
                small[0] = False
                result[small[groups]] = -1
    marked = labels >= 0
    result[marked] = labels[marked]
    result[outside] = -2
    return result


def height(data):
    """Height of the domes [nm]: zero-force topography once computed, else
    piezo height (with the planefit), nan without value."""
    return map_values.values(
        data, "topo" if data.stiffness_calculated else "piezo")


def reconstruction(seed, mask):
    """Grayscale reconstruction by dilation of seed under mask."""
    result = numpy.minimum(seed, mask)
    while True:
        grown = numpy.minimum(ndimage.grey_dilation(result, size=(3, 3)),
                              mask)
        if numpy.array_equal(grown, result):
            return result
        result = grown


def domes(values, min_height, seeds=()):
    """Domes of a height map: array [x][y], number of the dome of each pixel
    (from 1), and the number of domes. The tops of the domes higher than
    min_height (fraction of the range of the heights) above the grooves
    around them are flooded downhill (watershed): the boundaries are in the
    grooves. Lower domes (noise) belong to their neighbors. Each set of
    pixels of seeds is flooded too, as the region number (number of domes +
    1 + its index): a region without top (flat substrate) stops at the
    grooves around it."""
    z = smoothed(numpy.array(values, dtype=float), SCALES[0])
    finite = numpy.isfinite(z)
    if not finite.any():
        raise AutoRoiError("No height for this file.")
    z = numpy.where(finite, z, z[finite].min())
    depth = max(min_height * (z.max() - z.min()), 1e-12)
    # Tops: where the domes are at least depth high (h-maxima)
    dome = z - reconstruction(z - depth, z)
    tops, number = ndimage.label(dome >= depth * (1 - 1e-9))
    labels = tops.astype(int)
    nbr_x, nbr_y = z.shape
    for index, pixels in enumerate(seeds):
        for i, j in pixels:
            if 0 <= i < nbr_x and 0 <= j < nbr_y:
                labels[i, j] = number + 1 + index
    heap = [(-z[i, j], i, j) for i, j in numpy.argwhere(labels > 0)]
    heapq.heapify(heap)
    while heap:
        _, i, j = heapq.heappop(heap)
        for a, b in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
            if 0 <= a < nbr_x and 0 <= b < nbr_y and labels[a, b] == 0:
                labels[a, b] = labels[i, j]
                heapq.heappush(heap, (-z[a, b], a, b))
    return labels, number


def dome_classes(dome_labels, labels, inside=None, banned=None):
    """ROI of each pixel (as classify) from the domes: each ROI takes the
    domes under its marks (labels, see marks_labels), a dome marked in two
    ROIs is uncertain; -1 the domes not marked, -2 outside of the region of
    the map (inside: boolean array, None for the whole map); banned: the
    domes a ROI does not take ({index of the ROI: domes}). And the number of
    domes marked in several ROIs."""
    result = numpy.full(labels.shape, -1, dtype=int)
    banned = banned or {}
    owners = {}
    for index in range(labels.max() + 1):
        for dome in numpy.unique(dome_labels[labels == index]):
            if dome not in banned.get(index, ()):
                owners.setdefault(dome, set()).add(index)
    shared = 0
    for dome, rois in owners.items():
        shared += len(rois) > 1
        result[dome_labels == dome] = rois.pop() if len(rois) == 1 else -1
    marked = labels >= 0
    result[marked] = labels[marked]
    if inside is not None:
        result[~numpy.asarray(inside, dtype=bool)] = -2
    return result, shared
