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



"""Filters of the events, by distance from the jump of contact and by slope.

Each file has the filters of the Results tab (attributes events_results_*
of the dataset) and each row of the Plots table has its own filters
(attributes filter_* of the result).

The distance filter removes the events between two distances from the jump
of contact (in nm, positive away from the surface), or keeps only them.
The slope filter keeps only the events whose slope is between two values.
A filter is off when its maximum is not larger than its minimum (0 and 0 by
default).
"""


def keep_event(distance, slope, dist_min, dist_max, dist_keep,
               slope_min, slope_max):
    """True if the event passes the filters.

    distance is the distance of the event from the jump of contact (nm).
    """
    if slope_max > slope_min:
        if slope < slope_min or slope > slope_max:
            return False

    if dist_max > dist_min:
        inside = dist_min <= distance <= dist_max
        return inside if dist_keep else not inside

    return True


def format_number(value):
    """Short text of a value (no useless decimals)."""
    return "{:g}".format(value)


def dist_filter_text(dist_min, dist_max, dist_keep):
    """Short description of a distance filter."""
    if dist_max <= dist_min:
        return "None"
    action = "Keep only" if dist_keep else "Remove"
    return action + " " + format_number(dist_min) + "–" + \
        format_number(dist_max) + " nm"


def slope_filter_text(slope_min, slope_max):
    """Short description of a slope filter."""
    if slope_max <= slope_min:
        return "None"
    return format_number(slope_min) + "–" + format_number(slope_max) + \
        " nN/nm"


def convert_old_filters(obj):
    """Convert the distance filter of an old file (dataset or result).

    The old filters had a left and a right limit on the distance axis of the
    curve, where the events are at negative values, and removed the events
    between the limits unless keep middle was checked. The new ones are
    distances from the jump of contact: the meaning is kept.
    """
    for old, new in (("events_results_filter_dist", "events_results_filter"),
                     ("filter_dist", "filter")):
        left = old + "_left"
        if left not in obj.__dict__:
            continue
        right = obj.__dict__.pop(old + "_right", 0)
        left = obj.__dict__.pop(left)
        keep = obj.__dict__.pop(old + "_keep_middle", False)
        setattr(obj, new + "_dist_min", 0 - right)
        setattr(obj, new + "_dist_max", 0 - left)
        setattr(obj, new + "_dist_keep", keep)
