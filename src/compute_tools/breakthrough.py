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

"""Breakthroughs of supported lipid bilayers on the approach curves of a
file (see tools.breakthrough_tools): force [N], jump of the separation
[nm] and positions of each event, stored in data.breakthrough_results."""

import os

import numpy

from .. import shared
from .. import widgets_list
from ..tools import breakthrough_tools
from ..tools import curve_tools
from ..tools import misc_tools
from ..widgets.progressbar import Progressbar


def compute(calc_id):
    """Breakthroughs of the curves of a file. Returns the number of curves
    with at least one breakthrough."""
    data = shared.exp.list[calc_id]
    misc_tools.save_used_params(data, "breakthrough")

    nx, ny = data.nbr_pixels_x, data.nbr_pixels_y
    max_events = max(1, int(data.breakthrough_max_events))
    results = numpy.full((nx, ny, max_events,
                          len(breakthrough_tools.EVENT_VALUES)),
                         numpy.nan, dtype=numpy.float32)
    spring_constant = data.spring_constant * 1e9  # N/m

    Progressbar(os.path.basename(data.filename))
    widgets_list.widget_progressbar.set_label("Breakthroughs")
    widgets_list.widget_progressbar.set_range(0, nx * ny)

    found = 0
    for i in range(nx):
        for j in range(ny):
            if data.discarded_curves[i][j] == 0:
                approach, _, _, _, _, _ = curve_tools.get_curve(
                    data, [i, j], mode="curve_results")
                events = breakthrough_tools.detect(
                    approach[0], approach[1], spring_constant,
                    n_sigma=data.breakthrough_n_sigma,
                    min_thickness=data.breakthrough_min_thickness,
                    min_force=data.breakthrough_min_force,
                    max_events=max_events)
                for k, (force, thickness, peak, end) in enumerate(events):
                    # Force in N
                    results[i, j, k] = (force * 1e-9, thickness, peak, end)
                found += bool(events)
            widgets_list.widget_progressbar.update()

    data.breakthrough_results = results
    data.breakthrough_calculated = True
    data.update()
    widgets_list.widget_progressbar.close()
    return found
