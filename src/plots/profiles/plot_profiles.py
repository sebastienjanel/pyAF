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

"""Plots the profiles as a line."""

from ... import widgets_list
from ... import shared
from ..plot_main import MainPlot


class PlotProfiles(MainPlot):
    """Plots the profiles as a line."""

    def __init__(self, parent):
        super().__init__(parent)

        # Profiles of the values of the map displayed (see the profiles
        # widget)
        widget = widgets_list.widget_profiles
        self.pr_x = widget.profiles_x
        self.pr_z = widget.profiles_z
        pl = self.data.profile_list

        # Displayed profiles, with the color of each profile on the map
        colors = self.data.profiles_color_names
        profile_list = []
        for profile in range(len(pl)):
            if widgets_list.widget_profiles.checkbox_list[profile].isChecked():
                profile_list.append([self.pr_x[profile], self.pr_z[profile],
                                     colors[profile % len(colors)]])

        self.plot_type = "profiles"
        self.xlabel = "[\u03bcm]"
        self.ylabel = widget.value_label
        profiles_bar_1 = shared.exp.profiles_bar_1 / 1000.0
        profiles_bar_2 = shared.exp.profiles_bar_2 / 1000.0
        self.show_profiles_bars = shared.exp.show_profiles_bars

        if profile_list:
            # Plot the profiles
            for x, z, color in profile_list:
                self.axes.plot(x, z, color=color)

            # Plot the bars
            if shared.exp.show_profiles_bars:
                self.axes.axvline(x=profiles_bar_1, linewidth=2,
                                  color='black', linestyle='--')
                self.axes.axvline(x=profiles_bar_2, linewidth=2,
                                  color='black', linestyle='--')

        # Plot the labels
        self.plot_labels()

        # Do the actual plotting, is defined in MainPlot class
        self.start_plot()
