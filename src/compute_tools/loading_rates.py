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

"""This module computes the loading rates for each event."""

import os
import numpy
from ..tools import math_tools
from ..tools.curve_tools import get_force_curves, get_curve
from ..widgets.progressbar import Progressbar
from .. import widgets_list
from .. import shared


class GetLoadingRate:
    """Computes the loading rates for each event, on each curve.

    The events have to be computed before.
    """

    def __init__(self):
        self.temp_file = shared.exp.temp_file

    def compute(self, calc_id):
        """Does the actual computation.

        Goes through each curve adn cals the get_loading_rate_for_curve
        method to get the values. It stores them in the current dataset.
        """
        data = shared.exp.list[calc_id]
        speed = data.retraction_velocity * 1e-9  # Is in nm/s -> m/s
        lr_coef = data.lr_coef
        data.used_lr_coef = data.lr_coef
        # Use the same spring constant as the one for the events computation
        # Convert spring constant to N/m
        spring_constant = data.used_spring_constant * 1e9

        # If we recalculate the loading rates, it is better to remove the data
        # and recreate the tables.
        tmp = self.temp_file
        if data.loading_rates_calculated:
            tmp.delete_tables_for_results(str(calc_id), "loading_rates")
            tmp.create_tables_for_results(str(calc_id), "loading_rates")

        # Create a progressbar which is displayed during the computations
        Progressbar(os.path.basename(data.filename))
        widgets_list.widget_progressbar.set_label("Getting loading rates")
        nbr_curves = data.nbr_pixels_x * data.nbr_pixels_y
        widgets_list.widget_progressbar.set_range(0, nbr_curves)

        # Get the tables
        st = "/data/_" + str(calc_id)
        loading_rates = tmp.file.get_node(st + "/results/loading_rates")

        # Temporary table to store the data
        res_loading_rate = numpy.zeros([data.nbr_pixels_y], list)

        for i in range(data.nbr_pixels_x):
            for j in range(data.nbr_pixels_y):
                _, ret, _, _, _, _ = get_curve(
                    data, [i, j], mode="compute_no_shared")

                # X and Y axis in meters
                ret[0] = ret[0] * 1e9
                ret[1] = ret[1] * 1e9

                # The computation is done on the force curve
                force_curve = get_force_curves(
                    ret[0],
                    ret[1],
                    data.events_jocs2_real[i][j],
                    spring_constant)

                # Get the values for the loading rates of all the events of
                # the curve
                res_loading_rate[j], _ = \
                    get_loading_rate_for_curve(
                        force_curve,
                        data.events_positions_start[i][j],
                        data.events_positions_stop[i][j],
                        speed, lr_coef)

                widgets_list.widget_progressbar.update()

            for k in range(len(res_loading_rate)):
                loading_rates.append(res_loading_rate[k])

        # Close the progressbar
        widgets_list.widget_progressbar.close()

        # Flush the file and update the values in the experiment class
        self.temp_file.flush_file()
        data.loading_rates_calculated = True
        data.update()


def get_loading_rate_for_curve(curve, pos1, pos2, speed, lr_coef):
    """Computes the loading rate."""
    curve_x = curve[0]
    curve_z = curve[1]

    # Make sure pos1 and pos2 are integers.
    pos1 = pos1.astype("uint32")
    pos2 = pos2.astype("uint32")

    loading_rates = []
    fit_boundaries = []

    # Get loading rates
    nbr_events = len(pos1)

    for event in range(nbr_events):
        indextop = pos1[event]
        indexbottom = pos2[event]

        total_force = abs(curve_z[indextop] - curve_z[indexbottom])
        force_step = lr_coef * total_force

        check_next_event = False
        if event < nbr_events - 1:
            check_next_event = True

        if check_next_event:
            max_search_index = pos1[event + 1]
        else:
            max_search_index = len(curve_z)

        ztop = 0

        for z in range(indexbottom, max_search_index, 1):
            if curve_z[z] >= curve_z[indexbottom] + force_step:
                ztop = z
                break

        if ztop == 0:
            if event + 1 == nbr_events:
                # No next event, so ztop is directly the next point
                ztop = indexbottom + 1
            else:
                # Max z of the loading rate is first point of the next event
                ztop = pos1[event + 1]
            a, b = math_tools.find_line_coefs(
                [curve_x[indexbottom], curve_z[indexbottom]],
                [curve_x[ztop], curve_z[ztop]])
            if a is None:
                x, y = curve_x[ztop], curve_z[ztop]
            else:
                x, y = math_tools.find_intersection([0, curve_z[ztop]], [a, b])
        else:
            a, b = math_tools.find_line_coefs(
                [curve_x[ztop - 1], curve_z[ztop - 1]],
                [curve_x[ztop], curve_z[ztop]])
            if a is None:
                x, y = curve_x[ztop], curve_z[ztop]
            else:
                pos = curve_z[indexbottom] + force_step
                x, y = math_tools.find_intersection([0, pos], [a, b])

        if x is not None or y is not None:
            a, b = math_tools.find_line_coefs(
                [curve_x[indexbottom], curve_z[indexbottom]], [x, y])
        else:
            a = None

        if a is None:
            loading_rates.append(0)
        else:
            loading_rates.append(a * speed)

        fit_boundaries.append(
            [[curve_x[indexbottom], x], [curve_z[indexbottom], y]])

    return loading_rates, fit_boundaries
