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

"""Tools used to do some statistics and further results sorting."""

import numpy
from PyQt5 import QtWidgets
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
    # Reset the pdfs, these need to be recalculated
    shared.single_pdfs_x = None
    shared.single_pdfs_y = None
    shared.groups_pdfs_x = None
    shared.groups_pdfs_y = None
    shared.conditions_pdfs_x = None
    shared.conditions_pdfs_y = None

    # Display a progressbar (use loading progressbar at load or
    # create a new one for the refreshing button).
    new_bar = False
    if widgets_list.widget_progressbar is None:
        new_bar = True
        Progressbar()
    widgets_list.widget_progressbar.set_label("Refreshing histograms")
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
        shared.single_values.append(values)

    if new_bar:
        widgets_list.widget_progressbar.close()


def get_range(mode, themax, data):
    """Define a range for the plots, depending on the option the user selected."""
    if mode == "single":
        if shared.exp.hist_single_x_mode == "auto":
            r = (0, numpy.mean(data) * 3.0)
        elif shared.exp.hist_single_x_mode == "max":
            r = (0, themax)
        elif shared.exp.hist_single_x_mode == "manual":
            r = (shared.exp.hist_single_min_x, shared.exp.hist_single_max_x)
    elif mode == "groups":
        if shared.exp.hist_groups_x_mode == "auto":
            r = (0, numpy.mean(data) * 3.0)
        elif shared.exp.hist_groups_x_mode == "max":
            r = (0, themax)
        elif shared.exp.hist_groups_x_mode == "manual":
            r = (shared.exp.hist_groups_min_x, shared.exp.hist_groups_max_x)
    elif mode == "experiment":
        if shared.exp.hist_experiment_x_mode == "auto":
            r = (0, numpy.mean(data) * 3.0)
        elif shared.exp.hist_experiment_x_mode == "max":
            r = (0, themax)
        elif shared.exp.hist_experiment_x_mode == "manual":
            r = (shared.exp.hist_experiment_min_x, shared.exp.hist_experiment_max_x)

    return r


def get_bins(mode):
    """Get the number of bins for the histograms."""
    if mode == "single":
        if shared.exp.hist_single_bins_mode == "manual":
            bins = shared.exp.hist_single_bins
        else:
            bins = 10
    elif mode == "groups":
        if shared.exp.hist_groups_bins_mode == "manual":
            bins = shared.exp.hist_groups_bins
        else:
            bins = 10
    elif mode == "experiment":
        if shared.exp.hist_experiment_bins_mode == "manual":
            bins = shared.exp.hist_experiment_bins
        else:
            bins = 10
    return int(bins)


def get_bandwidth(mode):
    """Get the bandwidth for the PDF fitting."""
    if mode == "single":
        if shared.exp.hist_single_bw_mode == "manual":
            bandwidth = shared.exp.hist_single_bw
        else:
            bandwidth = None
    elif mode == "groups":
        if shared.exp.hist_groups_bw_mode == "manual":
            bandwidth = shared.exp.hist_groups_bw
        else:
            bandwidth = None
    elif mode == "experiment":
        if shared.exp.hist_experiment_bw_mode == "manual":
            bandwidth = shared.exp.hist_experiment_bw
        else:
            bandwidth = None

    return bandwidth


def get_frequencies(mode, i):
    """Get the frequency of a value in the PDF."""
    if mode == "single":
        _, _, _, mode_index = shared.single_values[i]
        pdf = shared.single_pdfs_y[i]
    elif mode == "groups":
        _, _, _, mode_index = shared.groups_values[i]
        pdf = shared.groups_pdfs_y[i]
    elif mode == "experiment":
        _, _, _, _, mode_index = shared.conditions_values[i]
        pdf = shared.conditions_pdfs_y[i]

    return [numpy.mean(pdf), numpy.median(pdf), None, pdf[mode_index]]


def get_pdf_and_mode(array, mode):
    """Returns the probability density function and the mode."""
    if not module_exists("scipy"):
        return None, None, None

    # Get a range
    steps = 500.0
    themax = int(numpy.amax(array))
    step_x = themax / steps

    # Get bandwidth. If set to None, the Scott's rule is used
    bw = get_bandwidth(mode)

    # Get a probability density function, and the mode
    pdf = stats.kde.gaussian_kde(array, bw_method=bw)

    newpdf_x = []
    newpdf_y = []

    for i in range(int(steps)):
        newpdf_x.append(i * step_x)
        try:
            newpdf_y.append(pdf(i * step_x)[0])
        except FloatingPointError:
            newpdf_y.append(0)
        widgets_list.widget_progressbar.update()

    max_index = numpy.argmax(newpdf_y)

    return numpy.array(newpdf_x), numpy.array(newpdf_y), max_index


def update_counts(mode):
    """Update the frequencies."""
    if mode == "single" or mode == "all":
        shared.single_factors = []
        shared.single_frequencies = []

        factors = []
        frequencies = []

        bins = get_bins("single")

        for i in range(len(shared.single_data)):
            data = shared.single_data[i]
            if data != []:
                # PDFS are computed on the whole range
                hist_normed, _ = \
                    numpy.histogram(data, bins=bins, density=True)
                hist_normal, _ = \
                    numpy.histogram(data, bins=bins, density=False)

                # Get a ratio for the factor calculation.
                # The factor is the link between the normed and non normed
                # histogram. It will help multiply the PDF in the plots in case
                # of a non-normed histogram

                for j in range(len(hist_normal)):
                    if hist_normal[j] != 0 and hist_normed[j] != 0:
                        factor = hist_normal[j] / hist_normed[j]
                        break

                factors.append(factor)
                counts = get_frequencies("single", i)

                frequencies.append(counts)

            else:
                factors.append(1.0)
                frequencies.append([0, 0, 0, 0])

        shared.single_factors = factors
        shared.single_frequencies = frequencies

    if mode == "groups" or mode == "all":
        # Groups
        shared.groups_factors = []
        shared.groups_frequencies = []

        factors = []
        frequencies = []

        bins = get_bins("groups")

        for i in range(len(shared.groups_data)):
            data = shared.groups_data[i]
            if data != []:
                # PDFS are computed on the whole range
                hist_normed, _ = \
                    numpy.histogram(data, bins=bins, density=True)
                hist_normal, _ = \
                    numpy.histogram(data, bins=bins, density=False)

                # Get a ratio for the factor calculation.
                # The factor is the link between the normed and non normed
                # histogram. It will help multiply the PDF in the plots in case
                # of a non-normed histogram

                for j in range(len(hist_normal)):
                    if hist_normal[j] != 0 and hist_normed[j] != 0:
                        factor = hist_normal[j] / hist_normed[j]
                        break

                factors.append(factor)
                counts = get_frequencies("groups", i)

                frequencies.append(counts)

            else:
                factors.append(1.0)
                frequencies.append([0, 0, 0, 0])

        shared.groups_factors = factors
        shared.groups_frequencies = frequencies

    if mode == "experiment" or mode == "all":
        # Groups
        shared.conditions_factors = []
        shared.conditions_frequencies = []

        factors = []
        frequencies = []

        bins = get_bins("experiment")

        for i in range(len(shared.conditions_data)):
            data = shared.conditions_data[i]
            if data != []:
                # PDFS are computed on the whole range
                hist_normed, _ = \
                    numpy.histogram(data, bins=bins, density=True)
                hist_normal, _ = \
                    numpy.histogram(data, bins=bins, density=False)

                # Get a ratio for the factor calculation.
                # The factor is the link between the normed and non normed
                # histogram. It will help multiply the PDF in the plots in case
                # of a non-normed histogram

                for j in range(len(hist_normal)):
                    if hist_normal[j] != 0 and hist_normed[j] != 0:
                        factor = hist_normal[j] / hist_normed[j]
                        break

                factors.append(factor)
                counts = get_frequencies("experiment", i)

                frequencies.append(counts)

            else:
                factors.append(1.0)
                frequencies.append([0, 0, 0, 0])

        shared.conditions_factors = factors
        shared.conditions_frequencies = frequencies


def get_pdfs_and_modes(mode):
    """Compute the PDF and the mode."""
    if mode == "single":
        shared.single_pdfs_x = []
        shared.single_pdfs_y = []

        # Create a progressbar, this operation is slow
        prog_range = 0
        for i in range(len(shared.single_data)):
            data = shared.single_data[i]
            themin = 0
            themax = 0
            if shared.exp.results_type != "loading_rates":
                if data != []:
                    themin = numpy.amin(data)
                    themax = numpy.amax(data)
            else:
                if data != [[], []]:
                    themin = numpy.amin(data[1])
                    themax = numpy.amax(data[1])
            prog_range += int(themax - themin)

        Progressbar()
        widgets_list.widget_progressbar.set_label("Getting PDFs (single)")
        widgets_list.widget_progressbar.set_range(0, prog_range)

        for i in range(len(shared.single_data)):
            data = shared.single_data[i]

            if shared.exp.results_type != "loading_rates":
                if data != []:
                    pdf_x, pdf_y, mode_index = get_pdf_and_mode(data, "single")

                    shared.single_pdfs_x.append(pdf_x)
                    shared.single_pdfs_y.append(pdf_y)
                    shared.single_values[i][3] = mode_index

                else:
                    shared.single_pdfs_x.append([])
                    shared.single_pdfs_y.append([])

            else:
                if data != [[], []]:
                    # 0 = forces, 1 = lr
                    pdf_x, pdf_y, mode_index = get_pdf_and_mode(
                        data[1], "single")

                    shared.single_pdfs_x.append(pdf_x)
                    shared.single_pdfs_y.append(pdf_y)
                    shared.single_values[i][3] = mode_index

                else:
                    shared.single_pdfs_x.append([])
                    shared.single_pdfs_y.append([])

        widgets_list.widget_progressbar.close()

    elif mode == "groups":
        shared.groups_pdfs_x = []
        shared.groups_pdfs_y = []

        # Create a progressbar, this operation is slow
        prog_range = 0
        for i in range(len(shared.groups_data)):
            data = shared.groups_data[i]

            themin = 0
            themax = 0
            if shared.exp.results_type != "loading_rates":
                if data != []:
                    themin = numpy.amin(data)
                    themax = numpy.amax(data)
            else:
                if data != [] and data != [[], []]:
                    themin = numpy.amin(data[1])
                    themax = numpy.amax(data[1])
            prog_range += int(themax - themin)

        Progressbar()
        widgets_list.widget_progressbar.set_label("Getting PDFs (groups)")
        widgets_list.widget_progressbar.set_range(0, prog_range)

        for i in range(len(shared.groups_data)):
            data = shared.groups_data[i]

            if shared.exp.results_type != "loading_rates":
                if data != []:
                    pdf_x, pdf_y, mode_index = get_pdf_and_mode(data, "groups")

                    shared.groups_pdfs_x.append(pdf_x)
                    shared.groups_pdfs_y.append(pdf_y)
                    shared.groups_values[i][3] = mode_index

                else:
                    shared.groups_pdfs_x.append([])
                    shared.groups_pdfs_y.append([])

            else:
                if data != [] and data != [[], []]:
                    # 0 = forces, 1 = lr
                    pdf_x, pdf_y, mode_index = get_pdf_and_mode(
                        data[1], "groups")

                    shared.groups_pdfs_x.append(pdf_x)
                    shared.groups_pdfs_y.append(pdf_y)
                    shared.groups_values[i][3] = mode_index

                else:
                    shared.groups_pdfs_x.append([])
                    shared.groups_pdfs_y.append([])

        widgets_list.widget_progressbar.close()

    elif mode == "experiment":
        shared.conditions_pdfs_x = []
        shared.conditions_pdfs_y = []

        # Create a progressbar, this operation is slow
        prog_range = 0
        for i in range(len(shared.conditions_data)):
            data = shared.conditions_data[i]

            themin = 0
            themax = 0
            if shared.exp.results_type != "loading_rates":
                if data != []:
                    themin = numpy.amin(data)
                    themax = numpy.amax(data)
            else:
                if data != [] and data != [[], []]:
                    themin = numpy.amin(data[1])
                    themax = numpy.amax(data[1])
            prog_range += int(themax - themin)

        Progressbar()
        widgets_list.widget_progressbar.set_label("Getting PDFs (experiment)")
        widgets_list.widget_progressbar.set_range(0, prog_range)

        for i in range(len(shared.conditions_data)):
            data = shared.conditions_data[i]

            if shared.exp.results_type != "loading_rates":
                if data != []:
                    pdf_x, pdf_y, mode_index = get_pdf_and_mode(data, "experiment")

                    shared.conditions_pdfs_x.append(pdf_x)
                    shared.conditions_pdfs_y.append(pdf_y)
                    shared.conditions_values[i][4] = mode_index

                else:
                    shared.conditions_pdfs_x.append([])
                    shared.conditions_pdfs_y.append([])

            else:
                if data != [] and data != [[], []]:
                    # 0 = forces, 1 = lr
                    pdf_x, pdf_y, mode_index = get_pdf_and_mode(
                        data[1], "experiment")

                    shared.conditions_pdfs_x.append(pdf_x)
                    shared.conditions_pdfs_y.append(pdf_y)
                    shared.conditions_values[i][4] = mode_index

                else:
                    shared.conditions_pdfs_x.append([])
                    shared.conditions_pdfs_y.append([])

        widgets_list.widget_progressbar.close()


def do_statistical_test(test):
    """Performs selected statistical test."""

    result = []

    for grp in shared.exp.sample_groups:

        if test == "Two-sample T-test":
            ttest_ind_result = stats.ttest_ind(shared.conditions_data[grp[0]],
                                               shared.conditions_data[grp[1]],
                                               nan_policy="omit")
            # print(ttest_ind_result)
            result.append((grp, ttest_ind_result))

        elif test == "Paired T-test":
            if numpy.shape(shared.conditions_data[grp[0]]) == numpy.shape(shared.conditions_data[grp[1]]):
                ttest_rel_result = stats.ttest_rel(shared.conditions_data[grp[0]],
                                                   shared.conditions_data[grp[1]],
                                                   nan_policy="omit")
                # print(ttest_rel_result)
                result.append((grp, ttest_rel_result))
            else:
                text = """Paired T-test can only be perfomed between samples of same size."""
                QtWidgets.QMessageBox.warning(None, "Info", text)

    shared.exp.statistical_test_results[test] = result


def fetch_group_data():
    """Fills the groups_data and groups_values arrays in shared."""
    results_util = results_sorting.GetResults()

    # Will return the sorted results
    [groups, values] = results_util.get_results_groups()

    shared.groups_data = groups
    shared.groups_values = values


def fetch_conditions_data():
    """Fills the conditions_data and conditions_values arrays in shared."""
    results_util = results_sorting.GetResults()

    # Will return the sorted results
    [conditions, values] = results_util.get_results_conditions()

    shared.conditions_data = conditions
    shared.conditions_values = values

    results_util.get_results_statistical()


def get_pdf():
    """Get the probability density function for the histogram."""
    get_pdfs_and_modes("single")
    update_counts("single")

    get_pdfs_and_modes("groups")
    update_counts("groups")

    get_pdfs_and_modes("experiment")
    update_counts("experiment")
