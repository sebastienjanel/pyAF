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

"""Tools used to do some statistics and further results sorting."""

import numpy
from .. import widgets_list
from .. import shared
from ..tools import results_sorting
from ..tools import misc_tools
from ..tools.utils import module_exists
from ..widgets.progressbar import Progressbar

if module_exists("scipy"):
    from scipy import stats


def get_values():
    """Fetch the results and fill the shared arrays with the values."""
    # Display a progressbar (use loading progressbar at load or
    # create a new one for the refreshing button).
    new_bar = False
    if widgets_list.widget_progressbar is None:
        new_bar = True
        Progressbar()
    widgets_list.widget_progressbar.set_label("Refreshing plots")
    if new_bar:
        # Get the range of the progressbar
        progressbar_range = misc_tools.get_hist_refresh_range(shared.exp)
        widgets_list.widget_progressbar.set_range(0, progressbar_range)

    # For the plots
    shared.single_data = []
    shared.single_values = []

    for i in range(len(shared.exp.results_list)):
        result = shared.exp.results_list[i]

        results_util = results_sorting.GetResults()

        # Will return the sorted results
        res = results_util.get_results_single(result_id=result.result_id)
        array, values = res[0], res[1]

        # For the plots
        shared.single_data.append(array)
        shared.single_values.append(list(values))

    update_single_modes()

    if new_bar:
        widgets_list.widget_progressbar.close()


def update_single_modes():
    """Computes the modes of the single results (index 3 of the values),
    with the density estimate of the violins (same smoothing)."""
    from ..experiment import get_violin_options
    bw_adjust = get_violin_options(shared.exp)["bw_adjust"]

    for i, data in enumerate(shared.single_data or []):
        if shared.exp.results_type == "loading_rates":
            # 0 = forces, 1 = loading rates
            data = data[1] if len(data) > 1 else []
        shared.single_values[i][3] = get_mode(data, bw_adjust)


def update_groups_modes():
    """Computes the modes of the groups (index 3 of the values), with the
    density estimate of their violins (same smoothing)."""
    from ..experiment import get_violin_options
    bw_adjust = get_violin_options(shared.exp, "groups")["bw_adjust"]

    for i, data in enumerate(shared.groups_data or []):
        if shared.exp.results_type == "loading_rates":
            # 0 = forces, 1 = loading rates
            data = data[1] if len(data) > 1 else []
        shared.groups_values[i][3] = get_mode(data, bw_adjust)


def get_mode(data, bw_adjust=1.0):
    """Mode of the values: maximum of their density estimate (gaussian
    kernels, Scott's rule for the bandwidth multiplied by bw_adjust).

    NaN if there are no values.
    """
    data = numpy.asarray(data, dtype=float).ravel()
    data = data[numpy.isfinite(data)]
    if data.size == 0:
        return numpy.nan
    if data.size == 1 or numpy.ptp(data) == 0:
        return data[0]
    kde = stats.gaussian_kde(data)
    kde.set_bandwidth(kde.factor * bw_adjust)
    grid = numpy.linspace(data.min(), data.max(), 1024)
    return grid[numpy.argmax(kde(grid))]


def fetch_group_data():
    """Fills the groups_data and groups_values arrays in shared."""
    results_util = results_sorting.GetResults()

    # Will return the sorted results
    [groups, values] = results_util.get_results_groups()

    shared.groups_data = groups
    shared.groups_values = values

    update_groups_modes()


