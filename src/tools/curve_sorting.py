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

"""Curve sorting: the good and bad curves of an experiment learned from a
few curves marked by the user (in any of its files, maps or single curves),
and the other curves sorted.

Bad curves differ from one experiment to the other (indentation of cells,
unbinding of molecules...), so nothing is assumed about them: the shape
of each curve, approach and retraction resampled to a few points over the
recorded points (no point of contact needed), relative to its baseline and
to the typical force range of the experiment, is described by the typical
shapes of the experiment it contains and how far it is from them
(principal components of all the curves); with a few numbers (force range,
adhesion, noise of the baseline, height range, largest jumps of the
approach and retraction, quality of the elasticity fit when computed).
Extremely randomized trees learn the marks, and give each curve a
probability to be bad. The raw points of the shapes are not used: with a
few marks, the noise of some of them would separate the marks by chance.
About 5 to 10 marks of each kind of curve (good, and each kind of bad
curve) are needed."""

import numpy

# Points of each part of the curve (approach, retraction)
POINTS = 32
# Principal components (typical shapes)
COMPONENTS = 10
# Curves used to find the typical shapes (at most)
PCA_CURVES = 20000
# Marks of each kind at least
MIN_MARKS = 2


class CurveSortingError(Exception):
    """The curves can not be sorted (message for the user)."""


def available():
    """True if scikit-learn can be imported."""
    try:
        import sklearn.ensemble  # noqa: F401
        import sklearn.decomposition  # noqa: F401
    except ImportError:
        return False
    return True


def resampled(values, points=POINTS):
    """Curves [curve][point] resampled to points values: mean of each of
    points consecutive parts (interpolated with fewer points)."""
    number = values.shape[1]
    if number < points:
        positions = numpy.linspace(0, number - 1, points)
        return numpy.array([numpy.interp(positions, numpy.arange(number), row)
                            for row in values])
    used = number // points * points
    start = (number - used) // 2
    return values[:, start:start + used].reshape(
        len(values), points, -1).mean(axis=2)


def measures(approach, retraction):
    """Measures of the curves of a row of a map (or of a single curve):
    approach and retraction [curve][z, deflection (, time)][point]. Shape
    (deflection of the approach and retraction relative to the baseline of
    the approach, resampled: [curve][2 * POINTS]), and [curve][force
    range, adhesion, noise of the baseline, largest jump of the approach,
    largest jump of the retraction, height range]; nan for the missing
    curves (no height range), and for all the curves without at least two
    points on the approach and on the retraction."""
    approach = numpy.asarray(approach, dtype=float)
    retraction = numpy.asarray(retraction, dtype=float)
    if approach.ndim != 3 or retraction.ndim != 3 or \
            approach.shape[-1] < 2 or retraction.shape[-1] < 2:
        return numpy.full((len(approach), 2 * POINTS), numpy.nan), \
            numpy.full((len(approach), 6), numpy.nan)
    z_app, defl_app = approach[:, 0], approach[:, 1]
    defl_ret = retraction[:, 1]
    number = defl_app.shape[1]
    baseline = numpy.median(defl_app[:, :max(number // 10, 1)], axis=1)
    shape = numpy.concatenate(
        [resampled(defl_app - baseline[:, None]),
         resampled(defl_ret - baseline[:, None])], axis=1)
    low, high = numpy.percentile(defl_app, [1, 99], axis=1)
    start = defl_app[:, :max(number // 3, 2)]
    numbers = numpy.stack(
        [high - low, baseline - defl_ret.min(axis=1),
         numpy.diff(start, axis=1).std(axis=1),
         numpy.abs(numpy.diff(defl_app, axis=1)).max(axis=1),
         numpy.abs(numpy.diff(defl_ret, axis=1)).max(axis=1),
         numpy.ptp(z_app, axis=1)], axis=1)
    missing = numpy.ptp(z_app, axis=1) == 0
    shape[missing] = numpy.nan
    numbers[missing] = numpy.nan
    return shape, numbers


def file_measures(data):
    """Measures of all the curves of a file (see measures): [x][y][...]."""
    nbr_x, nbr_y = data.nbr_pixels_x, data.nbr_pixels_y
    shape = numpy.full((nbr_x, nbr_y, 2 * POINTS), numpy.nan)
    numbers = numpy.full((nbr_x, nbr_y, 6), numpy.nan)
    for i in range(nbr_x):
        # A whole row of the map at once
        shape[i], numbers[i] = measures(data.curves_approach[i],
                                        data.curves_retraction[i])
    return shape, numbers


def features(all_measures, fit_quality=None):
    """Features of the curves of the files sorted together: for each file,
    [x][y][feature]. all_measures: (shape, numbers) of each file (see
    file_measures); fit_quality: R^2 map of each file (nan without), or
    None. The deflections are relative to the typical force range of all
    the files, the heights to their typical height range. The first
    feature is nan for the missing curves."""
    from sklearn.decomposition import PCA

    shapes = [shape.reshape(-1, shape.shape[-1]) for shape, _ in all_measures]
    numbers = [values.reshape(-1, values.shape[-1])
               for _, values in all_measures]
    every = numpy.concatenate(numbers)
    valid = numpy.isfinite(every).all(axis=1)
    if not valid.any():
        raise CurveSortingError("No curve to sort.")
    force = numpy.median(every[valid, 0])
    force = force if force > 0 else 1.0
    height = numpy.median(every[valid, -1])
    height = height if height > 0 else 1.0
    scale = numpy.array([force] * 5 + [height])
    shapes = [shape / force for shape in shapes]
    numbers = [values / scale for values in numbers]

    # Typical shapes of the experiment
    every = numpy.concatenate(shapes)
    finite = numpy.isfinite(every).all(axis=1)
    rows = every[finite]
    if len(rows) > PCA_CURVES:
        rows = rows[numpy.random.default_rng(0).choice(
            len(rows), PCA_CURVES, replace=False)]
    pca = PCA(n_components=min(COMPONENTS, len(rows), rows.shape[1]))
    pca.fit(rows)

    result = []
    for index, (shape, values) in enumerate(zip(shapes, numbers)):
        nbr_x, nbr_y = all_measures[index][0].shape[:2]
        finite = numpy.isfinite(shape).all(axis=1)
        components = numpy.full((len(shape), pca.n_components_), numpy.nan)
        distance = numpy.full((len(shape), 1), numpy.nan)
        if finite.any():
            components[finite] = pca.transform(shape[finite])
            rebuilt = pca.inverse_transform(components[finite])
            distance[finite, 0] = numpy.sqrt(
                ((shape[finite] - rebuilt) ** 2).mean(axis=1))
        quality = numpy.full((len(shape), 1), numpy.nan)
        if fit_quality is not None and fit_quality[index] is not None:
            quality[:, 0] = numpy.asarray(fit_quality[index],
                                          dtype=float).ravel()
        result.append(numpy.concatenate(
            [components, distance, values, quality],
            axis=1).reshape(nbr_x, nbr_y, -1))
    return result


def bad_probability(all_features, marks):
    """Probability of each curve to be bad: for each file, [x][y] (nan for
    the missing curves). marks: (index of the file, i, j, True if bad) of
    the curves marked."""
    from sklearn.ensemble import ExtraTreesClassifier

    rows, labels = [], []
    for index, i, j, bad in marks:
        values = all_features[index][i, j]
        if numpy.isfinite(values[0]):
            rows.append(values)
            labels.append(int(bad))
    labels = numpy.array(labels, dtype=int)
    if numpy.count_nonzero(labels == 1) < MIN_MARKS or \
            numpy.count_nonzero(labels == 0) < MIN_MARKS:
        raise CurveSortingError(
            "Mark at least %d good and %d bad curves (Good and Bad)."
            % (MIN_MARKS, MIN_MARKS))
    # All the features at each split (the informative ones are found even
    # with a few marks)
    forest = ExtraTreesClassifier(n_estimators=200, max_features=None,
                                  class_weight="balanced", n_jobs=-1,
                                  random_state=0)
    forest.fit(numpy.array(rows), labels)
    result = []
    for values in all_features:
        nbr_x, nbr_y = values.shape[:2]
        flat = values.reshape(nbr_x * nbr_y, -1)
        valid = numpy.isfinite(flat[:, 0])
        proba = numpy.full(nbr_x * nbr_y, numpy.nan)
        if valid.any():
            proba[valid] = forest.predict_proba(flat[valid])[:, 1]
        result.append(proba.reshape(nbr_x, nbr_y))
    return result
