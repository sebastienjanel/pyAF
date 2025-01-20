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

"""Exports the figure"""

from .. import widgets_list
from ..plots.PYAFPlot import PYAFPlot


def export_figure(name):
    """Function used to open a figure in a separate window.

    (Allows zooming and saving).
    """
    # To acces the tabs, go through the parent of the data widget
    current_tab = widgets_list.widget_data.parent.tabs.currentIndex()

    if name == "meshgrid" and current_tab == 0:
        # Data tab, meshgrid
        PYAFPlot(None, "meshgrid_data")
    elif name == "curve" and current_tab == 0:
        # Data tab, curve
        PYAFPlot(None, "curve_data")

    elif name == "curve" and current_tab == 1:
        # Compute tab, curve
        PYAFPlot(None, "fit_preview")

    elif name == "meshgrid" and current_tab == 2:
        # Results tab, meshgrid
        PYAFPlot(None, "meshgrid")
    elif name == "curve" and current_tab == 2:
        # Results tab, curve
        PYAFPlot(None, "curve_results")

    elif name == "results" and current_tab == 3:
        # Plots tab, curve
        PYAFPlot(None, "results_single")

    elif name == "results" and current_tab == 4:
        # Grouped plots tab, curve
        PYAFPlot(None, "results_groups")

    elif name == "results" and current_tab == 5:
        # Grouped plots tab, curve
        PYAFPlot(None, "results_experiment")