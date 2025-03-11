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

"""Tools used to fetch and sort the results."""

import numpy
import pandas as pd
from .. import shared
from ..tools import math_tools
from ..tools import curve_tools
from .. import widgets_list


class GetResults:
    """Results fetching class."""

    def __init__(self, export=None, force_type=None):
        # Can be set to "single" or "groups"
        self.export = export

        if force_type is not None:
            if not export:
                # Meshgrids
                self.type = force_type
        else:
            # Results (histograms)
            if shared.exp.results_type is not None:
                self.type = shared.exp.results_type
            else:
                self.type = None

        self.data = None

        self.roi_id = None
        self.slice_nbr = None
        self.corr = None
        self.select = None
        self.log = None

        self.filter_dist_left = None
        self.filter_dist_right = None
        self.filter_dist_keep_middle = None

        self.filter_slope_min = None
        self.filter_slope_max = None

    def load_data(self, data_id=None, result_id=None):
        """Define some variables."""
        if data_id is not None and result_id is None:
            self.data = shared.exp.list[data_id]

            self.filter_dist_left = self.data.events_results_filter_dist_left
            self.filter_dist_right = self.data.events_results_filter_dist_right
            self.filter_dist_keep_middle = \
                self.data.events_results_filter_dist_keep_middle

            self.filter_slope_min = self.data.events_results_filter_slope_min
            self.filter_slope_max = self.data.events_results_filter_slope_max

            self.roi_id = 0  # No Roi
            self.slice_nbr = 0  # No slice

        elif data_id is None and result_id is not None:
            result = shared.exp.results_list[result_id]
            self.data = shared.exp.list[result.data_id]

            self.filter_dist_left = result.filter_dist_left
            self.filter_dist_right = result.filter_dist_right
            self.filter_dist_keep_middle = result.filter_dist_keep_middle

            self.filter_slope_min = result.filter_slope_min
            self.filter_slope_max = result.filter_slope_max

            self.slice_nbr = result.slice
            self.roi_id = result.roi
            self.corr = result.sub_type
            self.select = result.dist_to_joc

    def get_results_single(
            self, data_id=None, result_id=None, force_type=None):
        """
        Get the results for a single file.

        """

        # Used by DFS widget and exporting
        if force_type is not None:
            self.type = force_type

        self.load_data(data_id, result_id)

        self.log = shared.exp.hist_log

        result = []
        forces = []
        lr = []

        if self.type == "stiffness":
            result = self.prepare_stiffness()
        elif self.type == "work":
            result = self.prepare_work()
        elif self.type == "rupture_force":
            result = self.prepare_rupture_force1()
        elif self.type == "events_forces":
            result = self.prepare_event_forces()
        elif self.type == "events_per_curve":
            result = self.prepare_events_per_curve()
        elif self.type == "events_rupture_force":
            result = self.prepare_rupture_force2()
        elif self.type == "loading_rates":
            forces, lr = self.prepare_loading_rates()
            result = [forces, lr]
        elif self.type == "events_distance":
            result = self.prepare_events_distance()

        result = self.remove_zeros(result)

        if self.export == "r_format":
            return result

        if self.export == "single":
            # Return for single export function
            return result

        else:
            if self.type == "loading_rates":
                forces, lr = result[0], result[1]

        # Normal data fetching
        if self.type != "loading_rates":
            if len(result) == 0:
                val = numpy.nan
                return [[], [val, val, val, 0]]

            else:
                values = [numpy.mean(result),
                          numpy.median(result),
                          numpy.std(result), 0]

                return [numpy.array(result), values]

        elif self.type == "loading_rates":
            if len(forces) == 0:
                val = numpy.nan
                return [[[], []], [val, val, val, 0]]

            else:
                values = [numpy.mean(lr), numpy.median(lr), numpy.std(lr), 0]

                return [[forces, lr], values]

    def get_results_groups(self, force_type=None):
        """Gets the results for a group."""
        # Used by DFS widget and exporting
        if force_type is not None:
            self.type = force_type

        values = []
        groups = []
        groups_forces = []
        groups_lr = []
        for i in range(len(shared.exp.groups_list)):
            groups.append([])
            # For loading rates
            groups_forces.append([])
            groups_lr.append([])

        for i in range(len(shared.exp.results_list)):
            res = shared.exp.results_list[i]

            if res.group != 0:
                if self.export != "groups":
                    # Get the values from single_data which has already been
                    # prepared. Is really faster
                    if shared.single_data is not None:
                        result = shared.single_data[i]
                    else:
                        # When the refresh at load checkbox is unclicked,
                        # single data is still set to None. Return nothing
                        # in this case.
                        result = []
                else:
                    # Refetch the data for exporting
                    result = self.get_results_single(result_id=i)[0]

                if self.type != "loading_rates":
                    groups[res.group].extend(result)
                else:
                    groups_forces[res.group].extend(result[0])
                    groups_lr[res.group].extend(result[1])

                if widgets_list.widget_progressbar is not None:
                    widgets_list.widget_progressbar.update()

        if self.type != "loading_rates":
            for i in range(len(groups)):
                if groups[i]:
                    values.append([numpy.mean(groups[i]),
                                   numpy.median(groups[i]),
                                   numpy.std(groups[i]), 0])
                else:
                    values.append([0, 0, 0, 0])

        elif self.type == "loading_rates":
            for i in range(len(groups)):
                if groups_forces[i]:
                    values.append([numpy.mean(groups_lr[i]),
                                   numpy.median(groups_lr[i]),
                                   numpy.std(groups_lr[i]), 0])

                    groups[i] = numpy.array([groups_forces[i], groups_lr[i]])
                else:
                    values.append([0, 0, 0, 0])
        # print(f"Groups:{[groups, values]}")
        return [groups, values]

    def get_results_conditions(self, force_type=None):
        """"Gets the results for a experimental condition."""
        # Used by DFS widget and exporting
        if force_type is not None:
            self.type = force_type

        values = []
        conditions = []
        conditions_forces = []
        conditions_lr = []
        for i in range(len(shared.exp.conditions_list)):
            conditions.append([])
            # For loading rates
            conditions_forces.append([])
            conditions_lr.append([])

        for i in range(len(shared.exp.groups_list)):
            grp = shared.exp.groups_list[i]

            if grp.condition != 0:
                if self.export != "conditions":  # How to manage this
                    # Get the values from single_data which has already been
                    # prepared. Is really faster
                    if shared.groups_data is not None:
                        result = shared.groups_data[i]
                    else:
                        # When the refresh at load checkbox is unclicked,
                        # single data is still set to None. Return nothing
                        # in this case.
                        result = []
                else:
                    # Refetch the data for exporting
                    result = self.get_results_groups(force_type)[0][i]

                if self.type != "loading_rates":
                    conditions[grp.condition].extend(result)
                else:
                    conditions_forces[grp.condition].extend(result[0])
                    conditions_lr[grp.condition].extend(result[1])

                if widgets_list.widget_progressbar is not None:
                    widgets_list.widget_progressbar.update()

        if self.type != "loading_rates":
            for i in range(len(conditions)):
                if conditions[i]:
                    values.append([len(conditions[i]),
                                   numpy.mean(conditions[i]),
                                   numpy.median(conditions[i]),
                                   numpy.std(conditions[i]), 0])
                else:
                    values.append([0, 0, 0, 0, 0])

        elif self.type == "loading_rates":
            for i in range(len(conditions)):
                if conditions_forces[i]:
                    values.append([len(conditions[i]),
                                   numpy.mean(conditions_lr[i]),
                                   numpy.median(conditions_lr[i]),
                                   numpy.std(conditions_lr[i]), 0])

                    conditions[i] = numpy.array([conditions_forces[i], conditions_lr[i]])
                else:
                    values.append([0, 0, 0, 0, 0])
        # print(f"Conditions:{[conditions, values]}")
        return [conditions, values]

    def get_results_statistical(self, force_type=None):
        """"Gets the results for a experimental condition."""
        # Used by DFS widget and exporting
        if force_type is not None:
            self.type = force_type

        groups = []
        conditions = []
        data = []

        for i in range(1, len(shared.exp.groups_list)):
            grp = shared.exp.groups_list[i]
            if shared.groups_data is not None:
                dt = numpy.asarray(shared.groups_data[i])
                grps = numpy.ones(dt.shape, dtype=int) * grp.group_id
                conds = numpy.ones(dt.shape, dtype=int) * grp.condition

                groups.extend(grps)
                conditions.extend(conds)
                data.extend(dt)

        # Construct data frame used for plots
        shared.statistics_data = pd.DataFrame({'Replicate': groups, 'Condition': conditions, 'Data': data})

    def prepare_stiffness(self):
        """Get the stiffness results."""
        dt = self.data

        if (not self.corr and dt.stiffness_calculated is False) or \
                (self.corr and dt.stiffness_corrected is False):
            return []

        if self.slice_nbr is not None:
            if self.slice_nbr == 0:
                slice_nbr = None
            else:
                slice_nbr = self.slice_nbr - 1

        if not self.log:
            if self.data.used_stiffness_model_selected == 3:
                factor = 1.0  # In N/m
            else:
                factor = 1e-3  # In kPa
        else:
            # Log(E [Pa])
            factor = 1.0

        array = []
        export_array = []
        index = 0
        for i in range(dt.nbr_pixels_x):
            for j in range(dt.nbr_pixels_y):
                # Filter out discarded curves
                if dt.discarded_curves[i][j] == 0:
                    # Filter out by ROI
                    if self.roi_id != 0:
                        roi_vals = dt.roi_list[self.roi_id - 1].values
                        inlist = math_tools.in_list(roi_vals, [i, j])

                    if self.roi_id == 0 or inlist:
                        if not self.corr:
                            array_vals = dt.stiffness_array[i][j]
                        else:
                            array_vals = dt.stiffness_corrected_array[i][j]

                        # Filter by slice
                        if slice_nbr is None:
                            # Take all the slices
                            therange = list(range(0, len(array_vals)))

                        else:
                            # Take only one slice (if present)
                            if len(array_vals) > slice_nbr:
                                therange = list(range(slice_nbr, slice_nbr + 1))
                            else:
                                # Do nothing
                                therange = list(range(0))

                        # Take all the values
                        for k in therange:
                            # Fetch the data
                            found = False
                            if not self.log:
                                array.append(array_vals[k])
                                found = True
                            elif self.log:
                                if array_vals[k] > 0:
                                    array.append(array_vals[k])
                                    found = True

                            # To export the data to a file, prepare the array
                            if self.export == "single" and found:
                                z = k * dt.indentation_step
                                if not dt.tomography:
                                    sl = "0 - " + str(z + dt.indentation_step)
                                else:
                                    sl = str(z) + " - " + \
                                        str(z + dt.indentation_step)

                                export_array.append(
                                    [index, i, j, sl, array[-1]])

                            if self.export == "r_format" and found:
                                export_array.append(
                                    [index, i, j, array[-1]])

                    index += 1

            if widgets_list.widget_progressbar is not None:
                widgets_list.widget_progressbar.update()

        # Returning for export function
        if self.export == "single" or self.export == "r_format":
            return export_array

        else:
            # Normal return
            if self.log:
                return numpy.log10(numpy.array(array) * factor)
            else:
                return numpy.array(array) * factor

    def prepare_work(self):
        """Get the work."""
        dt = self.data

        if dt.work_and_rupture_force1_calculated is False:
            return []

        array = []
        export_array = []
        index = 0
        for i in range(dt.nbr_pixels_x):
            for j in range(dt.nbr_pixels_y):
                if dt.discarded_curves[i][j] == 0:
                    if self.roi_id != 0:
                        roi_vals = dt.roi_list[self.roi_id - 1].values
                        inlist = math_tools.in_list(roi_vals, [i, j])

                    if self.roi_id == 0 or inlist:
                        array.append(dt.work[i][j])

                        if self.export == "single" or self.export == "r_format":
                            export_array.append([index, i, j, dt.work[i][j]])

                index += 1

            if widgets_list.widget_progressbar is not None:
                widgets_list.widget_progressbar.update()

        if self.export == "single" or self.export == "r_format":
            # Returning for export function
            return export_array

        else:
            return numpy.array(array) * 1e15  # in fJ

    def prepare_rupture_force1(self):
        """Get the rupture forces."""
        dt = self.data

        if dt.work_and_rupture_force1_calculated is False:
            return []

        array = []
        export_array = []
        index = 0
        for i in range(dt.nbr_pixels_x):
            for j in range(dt.nbr_pixels_y):
                if dt.discarded_curves[i][j] == 0:
                    if self.roi_id != 0:
                        roi_vals = dt.roi_list[self.roi_id - 1].values
                        inlist = math_tools.in_list(roi_vals, [i, j])

                    if self.roi_id == 0 or inlist:
                        array.append(dt.rupture_force1[i][j])

                        if self.export == "single" or self.export == "r_format":
                            val = dt.rupture_force1[i][j]
                            export_array.append([index, i, j, val])

                index += 1

            if widgets_list.widget_progressbar is not None:
                widgets_list.widget_progressbar.update()

        if self.export == "single" or self.export == "r_format":
            # Returning for export function
            return export_array

        else:
            return numpy.array(array) * 1e12  # in pN

    def prepare_event_forces(self):
        """Get the event forces."""
        dt = self.data

        if dt.events_calculated is False:
            return []

        array = []
        export_array = []
        index = 0
        for i in range(dt.nbr_pixels_x):
            for j in range(dt.nbr_pixels_y):
                if dt.discarded_curves[i][j] == 0:
                    if self.roi_id != 0:
                        roi_vals = dt.roi_list[self.roi_id - 1].values
                        inlist = math_tools.in_list(roi_vals, [i, j])

                    if self.roi_id == 0 or inlist:
                        forces = self.get_data_from_curve(
                            i, j, "forces_per_curve")

                        if forces:
                            array.extend(forces)

                        if self.export == "single" or self.export == "r_format":
                            for k in range(len(forces)):
                                export_array.append([index, i, j, forces[k]])

                index += 1

            if widgets_list.widget_progressbar is not None:
                widgets_list.widget_progressbar.update()

        if self.export == "single" or self.export == "r_format":
            # Returning for export function
            return export_array

        else:
            return numpy.array(array) * 1e12  # in pN

    def prepare_events_per_curve(self):
        """Get the number of events per curve."""
        dt = self.data

        if dt.events_calculated is False:
            return []

        tp = "events_per_curve"
        array = []
        for i in range(dt.nbr_pixels_x):
            for j in range(dt.nbr_pixels_y):
                if dt.discarded_curves[i][j] == 0:
                    if self.roi_id != 0:
                        roi_vals = dt.roi_list[self.roi_id - 1].values
                        inlist = math_tools.in_list(roi_vals, [i, j])

                    if self.roi_id == 0 or inlist:
                        nbr = self.get_data_from_curve(i, j, tp)

                        array.append(nbr)

            if widgets_list.widget_progressbar is not None:
                widgets_list.widget_progressbar.update()

        return numpy.array(array)

    def prepare_rupture_force2(self):
        """Get the rupture force from the events."""
        dt = self.data

        if dt.events_calculated is False:
            return []

        tp = "events_rupture_force"
        array = []
        for i in range(dt.nbr_pixels_x):
            for j in range(dt.nbr_pixels_y):
                if dt.discarded_curves[i][j] == 0:
                    if self.roi_id != 0:
                        roi_vals = dt.roi_list[self.roi_id - 1].values
                        inlist = math_tools.in_list(roi_vals, [i, j])

                    if self.roi_id == 0 or inlist:
                        force = self.get_data_from_curve(i, j, tp)

                        array.append(force)

            if widgets_list.widget_progressbar is not None:
                widgets_list.widget_progressbar.update()

        return numpy.array(array) * 1e12  # in pN

    def prepare_events_distance(self):
        """Get the distances to joc for the events."""
        dt = self.data

        if dt.events_calculated is False:
            return []

        f_dist = dt.events_forces_distance

        array = []
        export_array = []
        index = 0
        for i in range(dt.nbr_pixels_x):
            for j in range(dt.nbr_pixels_y):
                if dt.discarded_curves[i][j] == 0:
                    if self.roi_id != 0:
                        roi_vals = dt.roi_list[self.roi_id - 1].values
                        inlist = math_tools.in_list(roi_vals, [i, j])

                    if self.roi_id == 0 or inlist:
                        res = self.get_data_from_curve(
                            i, j, "events_distance")

                        if res:
                            array.extend(res)

                        if self.export == "single" or self.export == "r_format":
                            if res:
                                for k in range(len(res)):
                                    export_array.append(
                                        [index, i, j, f_dist[i][j][k]])

                index += 1

            if widgets_list.widget_progressbar is not None:
                widgets_list.widget_progressbar.update()

        if self.export == "single" or self.export == "r_format":
            # Returning for export function
            return export_array

        else:
            return numpy.array(array) * 1e9  # in nm

    def count_events_per_scan(self):
        """This function counts the number of events per scan.

        It parses each curve with the count_events_per_curve function and
        returns the value.
        """
        val = 0
        for i in range(self.data.nbr_pixels_x):
            for j in range(self.data.nbr_pixels_y):
                nbr = self.get_data_from_curve(i, j, "events_per_curve")
                val = val + nbr

        return val

    def prepare_loading_rates(self):
        """Get the loading rates."""
        dt = self.data

        if dt.loading_rates_calculated is False:
            return [], []

        array_forces = []
        array_lr = []
        export_array = []
        index = 0
        for i in range(dt.nbr_pixels_x):
            for j in range(dt.nbr_pixels_y):
                if dt.discarded_curves[i][j] == 0:
                    if self.roi_id != 0:
                        roi_vals = dt.roi_list[self.roi_id - 1].values
                        inlist = math_tools.in_list(roi_vals, [i, j])

                    if self.roi_id == 0 or inlist:
                        forces, lr = self.get_data_from_curve(i, j,
                                                              "loading_rates")

                        if lr:
                            array_forces.extend(forces)
                            array_lr.extend(lr)

                        if self.export == "single" or self.export == "r_format":
                            for loading_rate in lr:
                                export_array.append(
                                    [index, i, j, loading_rate])

                index += 1

            if widgets_list.widget_progressbar is not None:
                widgets_list.widget_progressbar.update()

        if self.export == "single" or self.export == "r_format":
            # Returning for export function
            return export_array, []

        else:
            return numpy.array(
                array_forces) * 1e12, numpy.array(array_lr) * 1e12

    def get_data_from_curve(self, i, j, type_of_data):
        """Get the wanted data from a curve.

        Filtering can be done by slope or distance limit.
        """
        _, retraction, _, _,  _, _ = \
            curve_tools.get_curve(self.data, [i, j], mode="results_sorting")

        if retraction is None:
            if type_of_data == "events_per_curve":
                return 0
            elif type_of_data == "forces_per_curve":
                return []
            elif type_of_data == "events_rupture_force":
                return 0
            elif type_of_data == "events_distance":
                return []
            elif type_of_data == "loading_rates":
                return [], []

        joc2_real = self.data.events_jocs2_real[i][j]

        re = curve_tools.get_force_curves(retraction[0], retraction[1],
                                          joc2_real, self.data.spring_constant)

        ret_pos = self.data.retraction_positions[i][j][0]

        # Get distance between two points on the curve
        xdist = math_tools.get_x_dist(re[0], ret_pos)

        if xdist is None:
            # Some corrupted curves have bad data with no distance between
            # points. Return default (empty) values.
            tp = type_of_data
            if tp == "forces_per_curve" or tp == "loading_rates":
                return [], []
            elif tp == "events_per_curve" or tp == "events_rupture_force":
                return 0
            elif tp == "events_distance":
                return []

        forces = self.data.events_forces[i][j]
        event_pos = self.data.events_positions_middle[i][j]
        dist_to_joc = self.data.events_forces_distance[i][j]
        slopes = self.data.events_slopes[i][j]

        # Make sure positions are integers.
        event_pos = event_pos.astype("uint32")

        if type_of_data == "events_per_curve":
            nbr_events = 0

            for z in range(len(forces)):
                if self.filter_event(re[0][event_pos[z]], slopes[z]):
                    nbr_events += 1

            return nbr_events

        elif type_of_data == "forces_per_curve":
            resforces = []

            for z in range(len(forces)):
                if self.filter_event(re[0][event_pos[z]], slopes[z]):
                    resforces.append(forces[z])

            return resforces

        elif type_of_data == "events_rupture_force":
            force = 0
            for z in range(len(forces)):
                if forces[z] > force:
                    if self.filter_event(re[0][event_pos[z]], slopes[z]):
                        force = forces[z]

            return force

        elif type_of_data == "events_distance":
            result = []

            if self.select == 0:
                # All
                for z in range(len(dist_to_joc)):
                    if self.filter_event(re[0][event_pos[z]], slopes[z]):
                        result.append(dist_to_joc[z])

            elif self.select == 1:
                # Last one
                for z in range(len(dist_to_joc)):
                    if self.filter_event(re[0][event_pos[z]], slopes[z]):
                        result.append(dist_to_joc[z])
                if result:
                    result = [result[0]]

            return result

        elif type_of_data == "loading_rates":
            result_forces = []
            result_lr = []

            for z in range(len(self.data.events_forces[i][j])):
                if self.filter_event(re[0][event_pos[z]], slopes[z]):
                    result_forces.append(self.data.events_forces[i][j][z])
                    result_lr.append(self.data.loading_rates[i][j][z])

            return result_forces, result_lr

    def remove_zeros(self, array):
        """Removes zeros and negatives from a results array."""
        if self.export != "single" and self.export != "r_format":
            if self.type == "loading_rates" and numpy.array(array).size == 0:
                return [[], []]
            elif self.type != "loading_rates" and numpy.array(array).size == 0:
                return []

        if shared.exp.hist_remove_zeros:
            if self.export == "single" or self.export == "r_format":
                if self.type == "loading_rates":
                    array = array[0]  # In this case the wanted data is in 0

                new_array = []

                for line in array:
                    # Stiffness result arrays incorporate an additional parameter sl (segment length)
                    if line[-1] > 0:
                        new_array.append(line)

                array = new_array

            else:
                # For normal usage or group exports

                if self.type != "loading_rates":
                    array = array[array > 0]
                else:
                    new_f = []
                    new_lr = []

                    # Sort on array[0] == forces
                    for i in range(len(array[0])):
                        if array[0][i] > 0:
                            new_f.append(array[0][i])
                            new_lr.append(array[1][i])
                    array[0] = new_f
                    array[1] = new_lr

        return array

    def filter_event(self, pos, slope):
        """Cheks if an event on the curve can be kept or not.

        Some filters can remove events from the curve so that it is not
        counted.

        pos is a negative value (the events are situaded before the joc2 == 0)

        slope is the value of the slope of the fit on the right of the event.
        The user can ask to filter by keeping only slopes between two values.
        """
        if self.filter_slope_max != 0 and self.filter_slope_min != 0:
            if slope > self.filter_slope_max or slope < self.filter_slope_min:
                return False

        if self.filter_dist_keep_middle:
            if self.filter_dist_left < pos and pos < self.filter_dist_right:
                return True
        else:
            if self.filter_dist_left > pos or pos > self.filter_dist_right:
                return True

        return False
