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

"""Tools for the events"""

from .. import shared
from .. import widgets_list
from ..widgets.progressbar import Progressbar


def update_events(mode=None, data_id=None):
    """Updates the events arrays while displaying a progressbar.

    The mode can be used to force the update of all the maps,
    when set to "all".
    """
    # Display a progressbar
    Progressbar("Refresh events")

    # Count the number of pixels to be refreshed
    total = 0
    if (shared.exp.apply_to_all_data or mode == "all") and data_id is None:
        for i in range(len(shared.exp.list)):
            data = shared.exp.list[i]
            if data.events_calculated:
                total += data.nbr_pixels_x * data.nbr_pixels_y
    else:
        if data_id is None:
            data = shared.exp.current_data
        else:
            data = shared.exp.list[data_id]
        if data.events_calculated:
            total = data.nbr_pixels_x * data.nbr_pixels_y

    widgets_list.widget_progressbar.set_range(0, total)

    # Update the arrays
    if (shared.exp.apply_to_all_data or mode == "all") and data_id is None:
        for i in range(len(shared.exp.list)):
            data = shared.exp.list[i]
            if data.events_calculated:
                widgets_list.widget_progressbar.set_label(data.filename)
                data.update_events_per_curve()
                data.update_rupture_force2()
    else:
        # Update only the current data
        if data_id is None:
            data = shared.exp.current_data
        else:
            data = shared.exp.list[data_id]
        if data.events_calculated:
            widgets_list.widget_progressbar.set_label(data.filename)
            data.update_events_per_curve()
            data.update_rupture_force2()

    widgets_list.widget_progressbar.close()
