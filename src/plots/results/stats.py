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

"""
Plots for experiment summary.
Based on: https://doi.org/10.1083/jcb.202001064
"""
import seaborn as sns
import itertools
from ..plot_main import MainPlot
from ... import shared


class PlotExperimentStats(MainPlot):
    """Plots the grouped results."""

    def __init__(self, parent):
        super().__init__(parent)

        # Plots supporting swarm plot overlay
        show_swarm_plots = ["Boxplot", "Violinplot"]

        # Plots supporting P-value display
        display_pval_plots = ["Two-sample T-test", "Paired T-test"]

        data = shared.statistics_data
        pairs = shared.exp.sample_groups

        stat_test = shared.exp.selected_statistical_test
        stat_results = shared.exp.statistical_test_results
        results = stat_results.get(stat_test)

        plot_type = shared.exp.stat_plot_type

        if len(pairs) >= 1:
            unique_grps = set(itertools.chain.from_iterable(pairs))
            filtered_data = data.loc[data['Condition'].isin(unique_grps)]

            # Seaborn stripplot scales badly, for more than 4000 points it takes a while.
            swarm_flag = len(filtered_data.index) > 4000

            show_swarm = False  # To define in experiment plots options
            show_descriptive_value = False  # To define in experiment plots options
            descriptive_value = "mean"  # This string can be: mean, median, To define in experiment plots options

            if plot_type == "Boxplot":
                sns.boxplot(x='Condition', y='Data', data=filtered_data, ax=self.axes)

            elif plot_type == "Violinplot":
                sns.violinplot(x='Condition', y='Data', data=filtered_data, ax=self.axes)

            elif plot_type == "Swarmplot":
                sns.stripplot(x='Condition', y='Data', data=filtered_data, ax=self.axes, linewidth=1)

            elif plot_type == "Superplot":

                if swarm_flag:
                    show_descriptive_value = True  # Provisional
                    sns.stripplot(x='Condition', y='Data', hue="Replicate", data=filtered_data, ax=self.axes,
                                  palette="Set2", linewidth=1, color=".25")
                else:
                    show_descriptive_value = True  # Provisional
                    sns.swarmplot(x='Condition', y='Data', hue="Replicate", data=filtered_data, ax=self.axes,
                                  palette="Set2", linewidth=1, color=".25")

            if show_swarm and plot_type in show_swarm_plots:
                if swarm_flag:
                    sns.stripplot(x='Condition', y='Data', data=filtered_data, ax=self.axes,
                                  linewidth=1, color=".25")
                else:
                    sns.swarmplot(x='Condition', y='Data', data=filtered_data, ax=self.axes,
                                  linewidth=1, color=".25")
            if show_descriptive_value:
                ReplicateDescriptiveValues = filtered_data.groupby(['Condition', 'Replicate'], as_index=False).agg(
                    {'Data': f"{descriptive_value}"})

                ax = sns.swarmplot(x='Condition', y='Data', hue='Replicate', size=15, edgecolor="k", linewidth=2,
                                   ax=self.axes, palette="Set2", data=ReplicateDescriptiveValues)
                ax.legend_. remove()

            if len(stat_results) >= 1:
                if stat_test in display_pval_plots:
                    for pair, results in results:
                        P_value = str(results[1])
                        x1, x2 = pair[0] - 1, pair[1] - 1
                        y, h, col = filtered_data['Data'].max() * 1.05, filtered_data['Data'].max() * 0.05, 'w'
                        self.axes.text((x1 + x2) * .5, y + h * 2, "P = " + P_value, ha='center', va='bottom', color=col)
                        self.axes.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)

        self.start_plot()
