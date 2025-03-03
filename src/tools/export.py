# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
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

"""Exports the results as text files."""

import os
import shutil

import pandas as pd
import numpy as np

from .. import shared
from .. import widgets_list as wl
from ..tools import math_tools as mt
from ..tools import results_sorting
from ..widgets.progressbar import Progressbar


def export_results(folder):
    """Exports the results as text files.

    The results are saved as filtered (results displayed in the histograms)
    and unfiltered (raw results). The groups histograms are also saved.

    A file is also saved with the list of which file belongs to which group.
    """
    Progressbar()
    wl.widget_progressbar.set_label("Exporting")

    progressbar_range = 0
    for row in range(len(shared.exp.results_list)):
        data_id = shared.exp.results_list[row].data_id
        if shared.exp.results_list[row].display:
            data = shared.exp.list[data_id]
            # Filtered and unfiltered single
            if data.stiffness_calculated:
                progressbar_range += 2
            if data.stiffness_corrected:
                progressbar_range += 2
            if data.work_and_rupture_force1_calculated:
                progressbar_range += 4
            if data.events_calculated:
                progressbar_range += 2
            if data.loading_rates_calculated:
                progressbar_range += 2
            # Filtered and unfiltered events per scan
            if data.events_calculated:
                progressbar_range += 2
            # Filtered and unfilterd distances of the events
            if data.events_calculated:
                progressbar_range += 2

        # Groups
        groups = []
        if shared.exp.results_list[row].group != 0:
            grp = shared.exp.results_list[row].group
            if mt.in_list(groups, grp) is False:
                groups.append(shared.exp.results_list[row].group)
                data = shared.exp.list[data_id]
                if data.stiffness_calculated:
                    progressbar_range += 1
                if data.stiffness_corrected:
                    progressbar_range += 1
                if data.work_and_rupture_force1_calculated:
                    progressbar_range += 2
                if data.events_calculated:
                    progressbar_range += 2
                if data.loading_rates_calculated:
                    progressbar_range += 1
                # Distance of the events
                if data.events_calculated:
                    progressbar_range += 1

    wl.widget_progressbar.set_range(0, progressbar_range)

    # Create folders
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)
    os.mkdir(folder + "/filtered")
    os.mkdir(folder + "/unfiltered")
    os.mkdir(folder + "/groups")

    # Save unfiltered data
    for row in range(len(shared.exp.results_list)):
        data_id = shared.exp.results_list[row].data_id
        if shared.exp.results_list[row].display:
            data = shared.exp.list[data_id]

            pre = folder + "/unfiltered/" + data.filename

            # Stiffness
            if data.stiffness_calculated:
                newfile = open(pre + "_elasticity.txt", "wb")
                newfile.write("Index\tXpos\tYpos\tSlice\tElasticity (Pa)\r\n".encode(encoding='UTF-8',errors='strict'))
                index = 0
                for i in range(len(data.stiffness_array[0])):
                    for j in range(len(data.stiffness_array[0])):
                        for k in range(0, len(data.stiffness_array[i][j])):
                            E = data.stiffness_array[i][j][k]
                            z = k * data.indentation_step
                            if not data.tomography:
                                sl = "0 - " + str(z + data.indentation_step)
                            else:
                                sl = str(z) + " - " + \
                                    str(z + data.indentation_step)

                            text = str(index + 1) + "\t" + \
                                str(i + 1) + "\t" + \
                                str(j + 1) + "\t" + \
                                sl + "\t" + \
                                str(E) + "\r\n"

                            newfile.write(text.encode(encoding='UTF-8',errors='strict'))
                        index += 1
                newfile.close()
                wl.widget_progressbar.update()

            # Work
            if data.work_and_rupture_force1_calculated:
                newfile = open(pre + "_work.txt", "wb")
                newfile.write("Index\tXpos\tYpos\tWork (J)\r\n".encode(encoding='UTF-8',errors='strict')
)
                index = 0
                for i in range(len(data.work[0])):
                    for j in range(len(data.work[0])):
                        text = str(index + 1) + "\t" + \
                            str(i + 1) + "\t" + \
                            str(j + 1) + "\t" + \
                            str(data.work[i][j]) + "\r\n"

                        newfile.write(text.encode(encoding='UTF-8',errors='strict'))
                        index += 1
                newfile.close()
                wl.widget_progressbar.update()

            # Rupture force1
            if data.work_and_rupture_force1_calculated:
                newfile = open(pre + "_rupture_force1.txt", "wb")
                newfile.write("Index\tXpos\tYpos\tRupture force (N)\r\n".encode(encoding='UTF-8',errors='strict'))
                index = 0
                for i in range(len(data.rupture_force1[0])):
                    for j in range(len(data.rupture_force1[0])):
                        text = str(index + 1) + "\t" + \
                            str(i + 1) + "\t" + \
                            str(j + 1) + "\t" + \
                            str(data.rupture_force1[i][j]) + "\r\n"

                        newfile.write(text.encode(encoding='UTF-8',errors='strict'))

                        index += 1
                newfile.close()
                wl.widget_progressbar.update()

            # Events force
            if data.events_calculated:
                newfile = open(pre + "_events_forces.txt", "wb")
                newfile.write("Index\tXpos\tYpos\tEvent Force (N)\r\n".encode(encoding='UTF-8',errors='strict'))
                index = 0
                for i in range(len(data.events_forces[0])):
                    for j in range(len(data.events_forces[0])):
                        for k in range(len(data.events_forces[i][j])):
                            text = str(index + 1) + "\t" + \
                                str(i + 1) + "\t" + \
                                str(j + 1) + "\t" + \
                                str(data.events_forces[i][j][k]) + "\r\n"

                            newfile.write(text.encode(encoding='UTF-8',errors='strict'))

                        index += 1
                newfile.close()
                wl.widget_progressbar.update()

            # Loading rates
            if data.loading_rates_calculated:
                newfile = open(pre + "_loading_rates.txt", "wb")
                newfile.write("Index\tXpos\tYpos\tLoading rate (N/s)\r\n".encode(encoding='UTF-8',errors='strict'))
                index = 0
                for i in range(len(data.events_forces[0])):
                    for j in range(len(data.events_forces[0])):
                        for k in range(len(data.events_forces[i][j])):
                            text = str(index + 1) + "\t" + \
                                str(i + 1) + "\t" + \
                                str(j + 1) + "\t" + \
                                str(data.loading_rates[i][j][k]) + "\r\n"

                            newfile.write(text.encode(encoding='UTF-8',errors='strict'))

                        index += 1
                newfile.close()
                wl.widget_progressbar.update()

            # Event's distance
            if data.events_calculated:
                newfile = open(pre + "_event_distance.txt", "wb")
                newfile.write("Index\tXpos\tYpos\tDistance (m)\r\n".encode(encoding='UTF-8',errors='strict'))
                index = 0
                for i in range(len(data.events_forces_distance[0])):
                    for j in range(len(data.events_forces_distance[0])):
                        for k in range(len(data.events_forces_distance[i][j])):
                            text = str(index + 1) + "\t" + \
                                str(i + 1) + "\t" + \
                                str(j + 1) + "\t" + \
                                str(data.events_forces_distance[i][j][k]) + \
                                "\r\n"

                            newfile.write(text.encode(encoding='UTF-8',errors='strict'))

                        index += 1
                newfile.close()
                wl.widget_progressbar.update()

    # Events per scan (unfiltered)
    found_events = False
    for i in range(len(shared.exp.results_list)):
        if shared.exp.results_list[i].display:
            data_id = shared.exp.results_list[i].data_id
            if shared.exp.list[data_id].events_calculated:
                found_events = True

    if found_events:
        index = 0
        newfile = open(folder + "/unfiltered/Events_per_scan.txt", "wb")
        newfile.write("Index\tFilename\tNumber of events\r\n".encode(encoding='UTF-8',errors='strict'))
        for a in range(len(shared.exp.results_list)):
            if shared.exp.results_list[a].display:
                val = 0
                data_id = shared.exp.results_list[a].data_id
                data = shared.exp.list[data_id]
                name = shared.exp.results_list[a].name
                if data.events_calculated:
                    for i in range(len(data.events_forces[0])):
                        for j in range(len(data.events_forces[0])):
                            for k in range(len(data.events_forces[i][j])):
                                val += 1
                text = str(index) + "\t" + name + "\t" + str(val) + "\r\n"
                newfile.write(text.encode(encoding='UTF-8',errors='strict'))
                wl.widget_progressbar.update()
            index += 1
        newfile.close()

    # Save filtered data
    for row in range(len(shared.exp.results_list)):
        data_id = shared.exp.results_list[row].data_id
        if shared.exp.results_list[row].display:
            data = shared.exp.list[data_id]
            name = shared.exp.results_list[row].name
            pre = folder + "/filtered/" + name
            # Stiffness
            if data.stiffness_calculated:
                newfile = open(pre + "_elasticity.txt", "wb")
                results_util = results_sorting.GetResults(export="single")
                newfile.write("Index\tXpos\tYpos\tSlice\tElasticity (Pa)\r\n".encode(encoding='UTF-8',errors='strict'))
                array = results_util.get_results_single(
                    result_id=row, force_type="stiffness")
                write_single_filtered(array, newfile, "stiffness")
                newfile.close()

            # Work
            if data.work_and_rupture_force1_calculated:
                newfile = open(pre + "_work.txt", "wb")
                results_util = results_sorting.GetResults(export="single")
                newfile.write("Index\tXpos\tYpos\tWork (J)\r\n".encode(encoding='UTF-8',errors='strict'))
                array = results_util.get_results_single(
                    result_id=row, force_type="work")
                write_single_filtered(array, newfile, "work")
                newfile.close()

            # Rupture force1
            if data.work_and_rupture_force1_calculated:
                newfile = open(pre + "_rupture_force1.txt", "wb")
                results_util = results_sorting.GetResults(export="single")
                newfile.write("Index\tXpos\tYpos\tRupture force (N)\r\n".encode(encoding='UTF-8',errors='strict'))
                array = results_util.get_results_single(
                    result_id=row, force_type="rupture_force")
                write_single_filtered(array, newfile, "rupture_force")
                newfile.close()

            # Events force
            if data.events_calculated:
                newfile = open(pre + "_events_forces.txt", "wb")
                results_util = results_sorting.GetResults(export="single")
                newfile.write("Index\tXpos\tYpos\tEvent Force (N)\tSlope\r\n".encode(encoding='UTF-8',errors='strict'))
                array = results_util.get_results_single(
                    result_id=row, force_type="events_forces")
                write_single_filtered(array, newfile, "events_forces")
                newfile.close()

            # Loading rates
            if data.loading_rates_calculated:
                newfile = open(pre + "_loading_rates.txt", "wb")
                results_util = results_sorting.GetResults(export="single")
                newfile.write("Index\tXpos\tYpos\tLoading rate (N/s)\r\n".encode(encoding='UTF-8',errors='strict'))
                array = results_util.get_results_single(
                    result_id=row, force_type="loading_rates")
                write_single_filtered(array, newfile, "loading_rates")
                newfile.close()

            # Event's distance
            if data.events_calculated:
                newfile = open(pre + "_event_distance.txt", "wb")
                results_util = results_sorting.GetResults(export="single")
                newfile.write("Index\tXpos\tYpos\tDistance (m)\r\n".encode(encoding='UTF-8',errors='strict'))
                array = results_util.get_results_single(
                    result_id=row, force_type="events_distance")
                write_single_filtered(array, newfile, "events_distance")
                newfile.close()

            wl.widget_progressbar.update()

    # Groups
    write_group_file("stiffness", folder)
    write_group_file("work", folder)
    write_group_file("rupture_force", folder)
    write_group_file("events_forces", folder)
    write_group_file("events_per_curve", folder)
    write_group_file("loading_rates", folder)
    write_group_file("events_distance", folder)

    # Infos
    write_info_file(folder)

    wl.widget_progressbar.close()

    return True


def export_results_r(folder, stats_flag=True):
    """
    Export results in a friendly R format.
    """

    data_df = pd.DataFrame()

    if stats_flag:
        data_df_stats = pd.DataFrame()

    for row in range(len(shared.exp.results_list)):
        data_id = shared.exp.results_list[row].data_id
        if shared.exp.results_list[row].display:
            data = shared.exp.list[data_id]
            row_df = pd.DataFrame()

            if stats_flag:
                row_df_stats = pd.DataFrame()

            results_util = results_sorting.GetResults(export="r_format")

            if data.stiffness_calculated:
                print("Shape of data.indentation:", data.indentation.shape)

                # Convert CArray to a NumPy array and flatten it
                flattened_indentation = np.array(data.indentation).flatten()
                print("Length of flattened indentation:", flattened_indentation.shape)

                # Get the filtered stiffness data (ROI only)
                stiffness_array = results_util.get_results_single(result_id=row, force_type="stiffness")
                stiffness = pd.DataFrame(stiffness_array, columns=["curve_index", "x_pos", "y_pos", "elasticity"])

                # Filter the indentation data to match the ROI
                # Assuming stiffness_array contains the ROI-filtered data, we need to filter indentation accordingly
                # Create a mask for the ROI
                roi_mask = np.zeros(data.indentation.shape, dtype=bool)
                for curve in stiffness_array:
                    i, j = int(curve[1]), int(curve[2])  # x_pos and y_pos
                    roi_mask[i, j] = True

                # Flatten the mask and filter the indentation data
                flattened_mask = roi_mask.flatten()
                roi_indentation = flattened_indentation[flattened_mask]

                # Assign the filtered indentation data to the stiffness DataFrame
                stiffness["indentation"] = roi_indentation
                print("Length of elasticity indentation:", len(stiffness["indentation"]))

                if stats_flag:
                    stiffness_stats = stiffness.describe()
                    ft_stats = export_descriptive_stats(stiffness_stats, "elasticity")
                    row_df_stats = pd.concat([row_df_stats, ft_stats], axis=1)

                row_df = merge_dfs(row_df, stiffness)

            # Repeat the same logic for other data types (work, rupture_force, etc.)
            if data.work_and_rupture_force1_calculated:
                work_array = results_util.get_results_single(result_id=row, force_type="work")
                work = pd.DataFrame(work_array, columns=["curve_index", "x_pos", "y_pos", "work"])

                if stats_flag:
                    work_stats = work.describe()
                    ft_stats = export_descriptive_stats(work_stats, "work")
                    row_df_stats = pd.concat([row_df_stats, ft_stats], axis=1)

                row_df = merge_dfs(row_df, work)

            if data.work_and_rupture_force1_calculated:
                rupture_force_array = results_util.get_results_single(result_id=row, force_type="rupture_force")
                rupture_force = pd.DataFrame(rupture_force_array, columns=["curve_index", "x_pos", "y_pos", "rupture_force"])

                if stats_flag:
                    rupture_force_stats = rupture_force.describe()
                    ft_stats = export_descriptive_stats(rupture_force_stats, "rupture_force")
                    row_df_stats = pd.concat([row_df_stats, ft_stats], axis=1)

                row_df = merge_dfs(row_df, rupture_force)

            if data.loading_rates_calculated:
                events_forces_array = results_util.get_results_single(result_id=row, force_type="events_forces")
                events_forces = pd.DataFrame(events_forces_array, columns=["curve_index", "x_pos", "y_pos", "events_forces"])

                if stats_flag:
                    events_forces_stats = events_forces.describe()
                    ft_stats = export_descriptive_stats(events_forces_stats, "events_forces")
                    row_df_stats = pd.concat([row_df_stats, ft_stats], axis=1)

                row_df = merge_dfs(row_df, events_forces)

            if data.events_calculated:
                loading_rates_array = results_util.get_results_single(result_id=row, force_type="loading_rates")
                loading_rates = pd.DataFrame(loading_rates_array, columns=["curve_index", "x_pos", "y_pos", "loading_rates"])

                if stats_flag:
                    loading_rates_stats = loading_rates.describe()
                    ft_stats = export_descriptive_stats(loading_rates_stats, "loading_rates")
                    row_df_stats = pd.concat([row_df_stats, ft_stats], axis=1)

                row_df = merge_dfs(row_df, loading_rates)

            if data.events_calculated:
                events_distance_array = results_util.get_results_single(result_id=row, force_type="events_distance")
                events_distance = pd.DataFrame(events_distance_array, columns=["curve_index", "x_pos", "y_pos", "events_distance"])

                if stats_flag:
                    events_distance_stats = events_distance.describe()
                    ft_stats = export_descriptive_stats(events_distance_stats, "events_distance")
                    row_df_stats = pd.concat([row_df_stats, ft_stats], axis=1)

                row_df = merge_dfs(row_df, events_distance)

            # Add information about the file
            res = shared.exp.results_list[row]
            group_id = res.group

            grp = shared.exp.groups_list[group_id]
            condition_id = grp.condition

            row_df.insert(0, 'condition', condition_id)
            row_df.insert(1, 'group', group_id)
            row_df.insert(2, 'file_id', data_id)
            row_df.insert(3, 'file_name', data.filename)

            data_df = merge_dfs(data_df, row_df, on=None)

            if stats_flag:
                row_df_stats.insert(0, 'condition', condition_id)
                row_df_stats.insert(1, 'group', group_id)
                row_df_stats.insert(2, 'file_id', data_id)
                row_df_stats.insert(3, 'file_name', data.filename)

                data_df_stats = merge_dfs(data_df_stats, row_df_stats, on=None)

    data_df.to_csv(folder)

    if stats_flag:
        data_df_stats.to_csv(folder + "_stats")


def export_descriptive_stats(variable_df_stats, variable_name):
    """
    Export descriptive statistics (Mean, Median, Mod, Min, Max) in a friendly R format.
    """
    out_df = pd.DataFrame()
    for stat_field in variable_df_stats.index.values.tolist():
        val = variable_df_stats.loc[stat_field, variable_name]
        out_df[stat_field + "_" + variable_name] = [val]

    return out_df


def merge_dfs(df1, df2, how='outer', on=["curve_index", "x_pos", "y_pos"]):

    if df1.empty:
        mdf = df2

    elif df2.empty:
        mdf = df1

    else:
        mdf = pd.merge(df1, df2, how=how, on=on)
        mdf.drop_duplicates()

        if on is not None:
            mdf.sort_values(by=on)

    return mdf


def write_single_filtered(array, newfile, mode):
    """Writes results to the file for one line.

    Do not forget to shift the indices from 0 to 1 for the users
    """
    for line in array:
        if mode == "stiffness":
            txt = str(line[0] + 1) + "\t"\
                + str(line[1] + 1) + "\t"\
                + str(line[2] + 1) + "\t"\
                + str(line[3]) + "\t"\
                + str(line[4]) + "\r\n"

        else:
            txt = str(line[0] + 1) + "\t"\
                + str(line[1] + 1) + "\t"\
                + str(line[2] + 1) + "\t"\
                + str(line[3]) + "\r\n"\

        newfile.write(txt.encode(encoding='UTF-8',errors='strict'))


def write_group_file(res_type, folder):
    """Write file containing the grouped results.

    Each group is written in a column. The columns are separated by
    semi-colons. If there is no more data in the group, the empty values are
    written as tabs.

    NOTE: The results_sorting.py file provides the values used here.
    All the values have been scaled to be displayed directly in the sotware.
    We need to scale them back for exporting in m, Pa, N or J.
    """
    if res_type == "stiffness":
        name = "stiffness"
        factor = 1e3
        title = "Stiffness (Pa)"
    elif res_type == "work":
        name = "work"
        factor = 1e-15
        title = "Work (J)"
    elif res_type == "rupture_force":
        name = "rupture_force"
        factor = 1e-12
        title = "Rupture force (N)"
    elif res_type == "events_forces":
        name = "events_forces"
        factor = 1e-12
        title = "Events forces (N)"
    elif res_type == "events_per_curve":
        name = "events_per_curve"
        factor = 1
        title = "Events per curve"
    elif res_type == "loading_rates":
        name = "loading_rates"
        factor = 1e-12
        title = "Loading rates"
    elif res_type == "events_distance":
        name = "events_distance"
        factor = 1e9
        title = "Event's distance (m)"

    results_util = results_sorting.GetResults(export="groups")
    [groups, _] = results_util.get_results_groups(force_type=name)

    group_found = False
    for group_id in range(len(shared.exp.groups_list)):
        if groups[group_id] != []:
            group_found = True

    if group_found:
        newfile = open(folder + "/groups/Groups_" + name + ".txt", "wb")
        newfile.write((title + "\r\n").encode(encoding='UTF-8',errors='strict'))

        # Get the length of the biggest group
        max_len = 0
        # Group 0 is an empty group (None value in the list)
        for group_id in range(1, len(shared.exp.groups_list)):
            if len(groups[group_id]) > max_len:
                max_len = len(groups[group_id])

        # Add header
        header = ""
        for group_id in range(1, len(shared.exp.groups_list)):
            if groups[group_id] != []:
                is_last_group = True
                for i in range(group_id, len(shared.exp.groups_list) - 1):
                    if groups[i + 1] != []:
                        is_last_group = False
                        break

                if is_last_group:
                    # Last group of the line
                    separator = ""

                else:
                    # All the other groups
                    separator = ";"

                header += shared.exp.groups_list[group_id].name + separator

        newfile.write((header + "\r\n").encode(encoding='UTF-8',errors='strict'))

        # Write lines
        for line_id in range(max_len):
            text = ""

            for group_id in range(1, len(shared.exp.groups_list)):
                if groups[group_id] != []:
                    is_last_group = True
                    for i in range(group_id,
                                   len(shared.exp.groups_list) - 1):
                        if groups[i + 1] != []:
                            is_last_group = False
                            break

                    if is_last_group:
                        # Last group of the line
                        separator = ""
                    else:
                        # All the other groups
                        separator = ";"

                    if line_id < len(groups[group_id]):
                        if name == "loading_rates":
                            val = groups[group_id][1][line_id] * factor
                        else:
                            val = groups[group_id][line_id] * factor

                        text += str(val) + separator

                    else:
                        text += "\t" + separator

            newfile.write((text + "\r\n").encode(encoding='UTF-8',errors='strict'))

        newfile.close()


def write_info_file(folder):
    """Writes a file containing the list of groups and files."""
    newfile = open(folder + "/info.txt", "wb")

    for group in shared.exp.groups_list:
        # Skip first None groups
        if group.group_id != 0:
            newfile.write((group.name + "\r\n").encode(encoding='UTF-8',errors='strict'))

            for result in shared.exp.results_list:
                if result.group == group.group_id:
                    newfile.write((result.name + "\r\n").encode(encoding='UTF-8',errors='strict'))

        newfile.write("------------------------------------------\r\n".encode(encoding='UTF-8',errors='strict'))

    newfile.close()
