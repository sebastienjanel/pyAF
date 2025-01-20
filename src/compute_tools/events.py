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

"""Computes the events on the curves."""

import os
import numpy
from ..tools import misc_tools
from ..tools import events_tools
from ..tools import search_joc2_events
from ..tools.curve_tools import get_curve, get_force_curves
from ..tools import events_refresher
from ..widgets.progressbar import Progressbar
from .. import widgets_list as wl
from .. import shared


class GetEvents:
    """Computes the events for each file."""

    def __init__(self):
        self.temp_file = shared.exp.temp_file

    def compute(self, calc_id):
        """Gets the events for each event on each curve.

        Computes also the distance between each event and the Joc2. The Joc2 is
        the intersection of the curve with a fit on the flat part of the curve.
        """
        data = shared.exp.list[calc_id]

        misc_tools.save_used_params(data, "events")

        # If we recalculate the events, it is better to remove the data
        # and recreate the tables.
        if data.events_calculated:
            self.temp_file.delete_tables_for_results(str(calc_id), "events")
            self.temp_file.create_tables_for_results(str(calc_id), "events")

        # Create a progressbar which is displayed during the computations
        Progressbar(os.path.basename(data.filename))
        wl.widget_progressbar.set_label("Getting events")
        wl.widget_progressbar.set_range(
            0, data.nbr_pixels_x * data.nbr_pixels_y)

        # Get the tables
        st = "/data/_" + str(calc_id)
        tf = self.temp_file.file
        events_positions_middle = \
            tf.get_node(st + "/results/events_positions_middle")
        events_positions_start = \
            tf.get_node(st + "/results/events_positions_start")
        events_positions_stop = \
            tf.get_node(st + "/results/events_positions_stop")
        events_forces = tf.get_node(st + "/results/events_forces")
        events_slopes = tf.get_node(st + "/results/events_slopes")
        events_distance = tf.get_node(st + "/results/events_forces_distance")
        jocs2_indices = tf.get_node(st + "/results/events_jocs2_indices")
        jocs2_real = tf.get_node(st + "/results/events_jocs2_real")
        events_fits_joc = tf.get_node(st + "/results/events_fits_joc")

        # Temporary table to store the data
        resevents_positions_middle = numpy.zeros([data.nbr_pixels_y], list)
        resevents_positions_start = numpy.zeros([data.nbr_pixels_y], list)
        resevents_positions_stop = numpy.zeros([data.nbr_pixels_y], list)
        resevents_forces = numpy.zeros([data.nbr_pixels_y], list)
        resevents_slopes = numpy.zeros([data.nbr_pixels_y], list)
        resevents_jocs2_indices = numpy.zeros([data.nbr_pixels_y], int)
        resevents_jocs2_real = numpy.zeros([data.nbr_pixels_y, 2], float)
        resevents_events_distance = numpy.zeros([data.nbr_pixels_y], list)
        resevents_fits_joc = numpy.zeros([data.nbr_pixels_y, 3], float)

        fit_params = {
            "joc_skip_start": data.fitparam_events_joc_skip_start,
            "joc_fit_length": data.fitparam_events_joc_fit_length,
            "joc_refit_option": data.fitparam_events_joc_refit_option,
            "joc_noise_multiplicator":
            data.fitparam_events_joc_noise_multiplicator,
            "joc_refit_times": data.fitparam_events_joc_refit_times}

        if data.events_algorithm == "kernel":
            # Get threshold
            if data.adaptive_threshold_option:
                threshold = \
                    data.events_kernel_adaptive_threshold
            else:
                threshold = data.events_kernel_detection_threshold

            convolution_params = {
                "kernel_size":
                data.kernel_size,
                "threshold": threshold,
                "adaptive_threshold_option":
                data.adaptive_threshold_option,
                "adaptive_smoothing_window":
                data.adaptive_smoothing_window,
                "adaptive_smoothing_order":
                data.adaptive_smoothing_order}

        for i in range(data.nbr_pixels_x):
            for j in range(data.nbr_pixels_y):
                app, ret, _, _, _, _ = get_curve(data, [i, j],
                                           mode="compute_no_shared")

                # Postion of the jump of contact
                # Needed for the stop value. Beyond the joc no event will be
                # searched. Can also be used to filter the events by distance.
                retpos = data.retraction_positions[i][j][0]
                joc2_i, joc2_r, fit_joc, _ = \
                    search_joc2_events.get_JOC_events(ret, retpos, fit_params)

                # The computation is done on the force curve
                force_curve = get_force_curves(
                    ret[0], ret[1], joc2_r, data.spring_constant * 1e9)

                # Force curve units: nN

                if app is None and ret is None:
                    # If no curve found (discarded curve for example),
                    # do not save any event
                    positions = []

                else:
                    # Get the events positions
                    if data.events_algorithm == "kernel":
                        # Get events position
                        positions = \
                            events_tools.get_events_positions_by_kernel(
                                force_curve[0],
                                force_curve[1],
                                joc2_i,
                                retpos,
                                convolution_params)

                    elif data.events_algorithm == "msf":
                        positions = events_tools.get_events_positions_by_msf(
                            force_curve,
                            data.events_msf_detection_threshold,
                            data.events_msf_window_size,
                            joc2_i,
                            retpos)

                # Get forces and slopes
                pos_start = []
                pos_stop = []
                pos_middle_indice = []
                dist = []
                ang = []

                fit_size = data.used_events_fit_size
                fits = events_tools.find_fits_around_middle_pos(
                    force_curve, positions, fit_size)
                _, _, forces = events_tools.find_force_from_fits(
                    force_curve, fits, positions)

                # Store forces in N
                forces = numpy.array(forces) * 1e-9

                # pyplot.plot(force_curve[0], force_curve[1])

                slopes = []
                for fit in fits:
                    # a = slice(
                    # fit["seg_left_fit"][0], fit["seg_left_fit"][1] + 1)
                    # seg_x = force_curve[0][a]
                    # y = numpy.polyval(
                    # [fit["coeffs_left"][0], fit["coeffs_left"][1]], seg_x)
                    # pyplot.plot(seg_x, y, color="r", linewidth=1)

                    slopes.append(fit["coeffs_left"][0])

                    # a = slice(
                    # fit["seg_right_fit"][0], fit["seg_right_fit"][1] + 1)
                    # seg_x = force_curve[0][a]
                    # y = numpy.polyval(
                    # [fit["coeffs_right"][0], fit["coeffs_right"][1]], seg_x)
                    # pyplot.plot(seg_x, y, color="r", linewidth=1)

                # pyplot.show()

                for k in range(len(positions)):
                    # Middle position of the event (indice)
                    pos_middle_indice.append(positions[k][1])

                    # Estimation of the start and end positon (in indices)
                    pos_start.append(positions[k][0])
                    pos_stop.append(positions[k][2])

                    # Distance between Joc2 and event on force curve
                    # The joc2 is at 0, 0. The value is given in m
                    dist.append(abs(force_curve[0][positions[k][1]]) * 1e-9)

                    # Save the angle
                    ang.append(0)

                resevents_positions_middle[j] = pos_middle_indice
                resevents_positions_start[j] = pos_start
                resevents_positions_stop[j] = pos_stop
                resevents_slopes[j] = slopes
                resevents_forces[j] = forces
                resevents_events_distance[j] = dist
                resevents_jocs2_indices[j] = joc2_i
                resevents_jocs2_real[j] = joc2_r
                resevents_fits_joc[j] = fit_joc

                wl.widget_progressbar.update()

            # Arrays are VLarrays with only one dimension
            # Store the data in the hdf5 file
            for k in range(len(resevents_forces)):
                events_positions_middle.append(resevents_positions_middle[k])
                events_positions_start.append(resevents_positions_start[k])
                events_positions_stop.append(resevents_positions_stop[k])
                events_slopes.append(resevents_slopes[k])
                events_forces.append(resevents_forces[k])
                events_distance.append(resevents_events_distance[k])

            jocs2_indices[i, :] = resevents_jocs2_indices
            jocs2_real[i, :] = resevents_jocs2_real
            events_fits_joc[i, :] = resevents_fits_joc

        # Close the progressbar
        wl.widget_progressbar.close()

        # Flush the file and update the values in the experiment class
        self.temp_file.flush_file()
        data.events_calculated = True

        # Reset datadict
        data.datadict = {}

        # Update the events array (with a progressbar)
        events_refresher.update_events(data_id=data.unique_id)

        data.update()
