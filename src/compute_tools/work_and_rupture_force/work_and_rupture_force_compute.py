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

"""Computes the work and rupture force."""

import numpy
from ...tools import search_joc
from ...tools import curve_tools
from ...tools import math_tools


def main_compute_work(i, j, app, ret, dt, fit_params):
    """Utilitary function to compute the work and rupture force for one curve."""
    _, ret, _, _, _, _ = curve_tools.get_curve(
        None, [i, j], app, ret, mode="compute", dt=dt)

    # Get JOC
    joc1_indice, joc2_indice, joc1_real, joc2_real, fit_joc, _ = \
        search_joc.get_JOC(ret, dt["retraction_position"], fit_params)

    work, rupture_force1 = get_surface_and_force(ret[0],
                                                 ret[1],
                                                 joc1_indice,
                                                 joc2_indice,
                                                 joc1_real,
                                                 joc2_real,
                                                 dt["spring_constant"])

    return i, j, joc1_indice, joc2_indice, \
        joc1_real, joc2_real, fit_joc, work, rupture_force1


def get_surface_and_force(x, y, ind1, ind2, joc1_real, joc2_real, spring_ct):
    """Computes the work and the rupture force, both relative to the line
    between the two jumps of contact.

    The surface is computed with the trapeze's method.
    """
    curve_x, curve_y = x, y
    joc1_indice, joc2_indice = ind1, ind2

    new_curve_x = curve_x.copy()
    new_curve_y = curve_y.copy()

    # Segment curve
    segment = slice(joc1_indice, joc2_indice)
    new_curve_x = new_curve_x[segment]
    new_curve_y = new_curve_y[segment]

    # Get force curve
    new_curve_x = new_curve_x - new_curve_y
    new_curve_y = new_curve_y * spring_ct

    joc1_real_force_x = joc1_real[0] - joc1_real[1]
    joc1_real_force_y = joc1_real[1] * spring_ct
    joc2_real_force_x = joc2_real[0] - joc2_real[1]
    joc2_real_force_y = joc2_real[1] * spring_ct

    # Insert joc 1 at x=0,y=0
    new_curve_x[0] = joc1_real_force_x
    new_curve_y[0] = joc1_real_force_y
    new_curve_x[len(new_curve_x) - 1] = joc2_real_force_x
    new_curve_y[len(new_curve_y) - 1] = joc2_real_force_y

    # Go to 0
    new_curve_x = new_curve_x - new_curve_x[0]
    new_curve_y = new_curve_y - new_curve_y[0]

    # Distance between joc1 and joc2 (signed: the line goes through both)
    distx = joc2_real_force_x - joc1_real_force_x
    disty = joc2_real_force_y - joc1_real_force_y

    # The fit was the line between the two jocs
    a, b = math_tools.find_line_coefs([0.0, 0.0], [distx, disty])
    if a is not None and b is not None:
        y = numpy.polyval([a, b], new_curve_x)
    else:
        return 0, 0

    # Work: surface between the line and the curve (trapezes); rupture
    # force: largest distance between the line and the curve, so that both
    # use the same reference when the baseline drifts
    diff = y - new_curve_y
    # dx is in nm, convert it to m so that the work is in J
    surface_trapezes = numpy.sum(
        numpy.diff(new_curve_x) * 1e-9 * (diff[:-1] + diff[1:]) / 2.0)
    rupture_force1 = max(numpy.max(diff), 0.0) if len(diff) else 0.0

    return abs(surface_trapezes), rupture_force1
