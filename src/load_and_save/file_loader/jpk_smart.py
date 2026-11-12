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

"""Loads JPK Smart Mapping files (.h5-jpk).

These are HDF5 files (USID layout, written by the JPK software from version
8): the curves of a map, with a region of the map drawn by hand (the pixels
outside of the region have no curve), and the height of the head combined
with the z motor (combinedHeightMeasured, the z sensor, and combinedHeight,
the piezo), for maps higher than the range of the z piezo.

    /Measurement_000/Map: header of the map (attributes, the same keys as in
        the header.properties files of the other JPK files, with the
        multi-scan-series and shared-data prefixes), and for each curve
        (row): its position index (meta-data/valid-indices) and position
        (Position_Values, [x, y] in m)
    /Measurement_000/Map/SegmentN: style of the segment (attribute; or
        segments named Extend, Retract... in the files of the software
        before 8.1.6x), number
        of points of each curve (meta-data/num-points, the curves are
        padded to the longest one) and the channels (ChannelK/Name, one row
        per curve, the raw values are converted with the net-encoder
        attributes of ChannelK)

The curves are stored like the ones of the other JPK files (see LoaderJPK).
"""

import warnings

import numpy
import tables
from dateutil.parser import parse

from ...tools import plant_tools
from .jpk import LoaderJPK

MAP = "/Measurement_000/Map"
SERIES = "multi-scan-series."
SETTINGS = SERIES + "map.header.settings.force-settings."
GRID = SERIES + "map.header.position-pattern.grid."

# Height channels, by order of preference: the combined height (z sensor or
# piezo, with the z motor), always used when it is there
HEIGHT_CHANNELS = ("combinedHeightMeasured", "combinedHeight",
                   "measuredHeight", "height")
# Height channels from the z sensor
SENSOR_CHANNELS = ("combinedHeightMeasured", "measuredHeight")
# Height of the head without the z motor, for each combined height (the
# position of the motor is the difference)
MOTOR_CHANNELS = {"combinedHeightMeasured": "measuredHeight",
                  "combinedHeight": "height"}
# Smallest move of the z motor (nm)
MOTOR_TOLERANCE = 10.0


def read_attributes(node):
    """Attributes of a node of the file, as strings."""
    attributes = {}
    for name in node._v_attrs._f_list("user"):
        value = node._v_attrs[name]
        if isinstance(value, numpy.ndarray) and value.size == 1:
            value = value.item()
        if isinstance(value, (bytes, numpy.bytes_)):
            value = value.decode("utf-8", "replace")
        attributes[name] = str(value)
    return attributes


# Segments named by their style (files of the software before 8.1.6x:
# Extend, Retract... without style attribute), in the order of a curve
NAMED_SEGMENTS = ("extend", "pause", "modulation", "retract")


def segment_groups(afm_file):
    """Segment groups of the map ({segment number: group}): SegmentN, or the
    segments named by their style (numbered in the order of a curve)."""
    groups = {}
    named = {}
    for name, group in afm_file.get_node(MAP)._v_children.items():
        if name.startswith("Segment") and name[7:].isdigit():
            groups[int(name[7:])] = group
        elif name.lower() in NAMED_SEGMENTS:
            named[name.lower()] = group
    if not groups:
        groups = {number: named[style] for number, style in
                  enumerate(s for s in NAMED_SEGMENTS if s in named)}
    return groups


def segment_style(group):
    """Style of a segment (extend, retract, pause, modulation): its
    attribute, or its name (see NAMED_SEGMENTS)."""
    style = read_attributes(group).get("style")
    if style is None and group._v_name.lower() in NAMED_SEGMENTS:
        style = group._v_name.lower()
    return style


def channel_nodes(group):
    """Channels of a segment ({channel name: (dataset, attributes)}): the
    order of the channels is not the same in all the segments."""
    channels = {}
    for name, node in group._v_children.items():
        if not name.startswith("Channel"):
            continue
        attributes = read_attributes(node)
        datasets = [child for child in node._v_children.values()
                    if isinstance(child, tables.Leaf)]
        if "name" in attributes and datasets:
            channels[attributes["name"]] = (datasets[0], attributes)
    return channels


def check_smart_mapping(path):
    """(error, number of curves) of a Smart Mapping file, for the
    pre_check: error is "" if the file can be loaded."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", tables.exceptions.DataTypeWarning)
        return check_file(path)


def check_file(path):
    """See check_smart_mapping."""
    try:
        afm_file = tables.open_file(path, "r")
    except (OSError, tables.HDF5ExtError):
        return "Corrupted file", 0
    try:
        info = read_attributes(afm_file.root).get("file-format-info", "")
        if "Force Map" not in info or MAP not in afm_file:
            return "Not a JPK force map", 0
        styles = {segment_style(group)
                  for group in segment_groups(afm_file).values()}
        if "extend" not in styles or "retract" not in styles:
            return "No force curves in the file", 0
        nbr_curves = afm_file.get_node(MAP + "/meta-data/valid-indices") \
            .shape[0]
        return "", nbr_curves
    except (tables.NoSuchNodeError, KeyError, ValueError):
        return "Corrupted file", 0
    finally:
        afm_file.close()


class LoaderJPKSmart(LoaderJPK):
    """Class used to load JPK Smart Mapping files (see LoaderJPK, the curves
    are stored the same way)."""

    def __init__(self, file_info, unique_id):
        with warnings.catch_warnings():
            # PyTables can not read some attributes (HDF5 references), which
            # are not used
            warnings.simplefilter("ignore", tables.exceptions.DataTypeWarning)
            afm_file = tables.open_file(file_info["path"], "r")
            try:
                self.afm_file = afm_file
                self.groups = segment_groups(afm_file)
                # Channels of each segment
                self.channels_of = {}
                super().__init__(file_info, unique_id)
            finally:
                afm_file.close()

    def open_file(self):
        """The file is opened in __init__ (see LoaderJPK)."""

    def fetch_header(self):
        """Saves all the informations in the header."""
        header = read_attributes(self.afm_file.get_node(MAP))

        def value(key, default=None):
            return header.get(key, default)

        def number(key, default=None):
            text = value(key)
            return float(text) if text not in (None, "") else default

        self.data.microscope_name = value(SERIES + "description.instrument")
        self.data.version = value(SERIES + "description.source-software")
        self.data.retracted_delay = number(SETTINGS + "retracted-pause-time")
        self.data.extended_delay = number(SETTINGS + "extended-pause-time")
        self.data.pause_duration = self.data.extended_delay
        start = value(SERIES + "map.header.start-time")
        if start:
            date = parse(start, ignoretz=True)
            self.data.date = str(date).split(".")[0]

        # Grid of the map
        self.data.scan_angle = number(GRID + "theta", 0.0)
        self.data.nbr_pixels_x = int(number(GRID + "ilength"))
        self.data.nbr_pixels_y = int(number(GRID + "jlength"))
        self.data.scan_size_x = number(GRID + "ulength") * 1e9
        self.data.scan_size_y = number(GRID + "vlength") * 1e9
        for key in ("xcenter", "ycenter", "ulength", "vlength"):
            self.grid[key] = number(GRID + key)
        self.grid["theta"] = self.data.scan_angle
        self.numbering = value(SERIES + "map.header.position-pattern.numbering",
                               "back-and-forth")
        if self.numbering not in ("back-and-forth", "line-by-line"):
            raise ValueError("Unknown numbering of the curves: " +
                             self.numbering)

        closed_loop = value(SERIES + "header.force-settings.closed-loop")
        if closed_loop == "true":
            self.data.z_closed_loop = "On"
        elif closed_loop == "false":
            self.data.z_closed_loop = "Off"

        # Ramp
        z_start = number(SETTINGS + "relative-z-start")
        z_end = number(SETTINGS + "relative-z-end")
        if z_start is not None and z_end is not None:
            self.data.ramp_size = (z_start - z_end) * 1e9
        time_approach = number(SETTINGS + "extend-scan-time")
        time_retraction = number(SETTINGS + "retract-scan-time")
        self.data.trig_threshold = number(SETTINGS + "relative-setpoint")

        # Curves of the file: position index and position of each curve
        valid = self.afm_file.get_node(MAP + "/meta-data/valid-indices")
        self.position_indexes = [int(index) for index in valid[:, 0]]
        positions = self.afm_file.get_node(MAP + "/Position_Values")[:]
        self.positions = [[float(x), float(y)] for x, y in positions[:, :2]]
        self.data.real_nbr_pixels = len(self.position_indexes)

        # Segments
        for segment_nbr in sorted(self.groups):
            group = self.groups[segment_nbr]
            style = segment_style(group)
            if style is None:
                continue
            self.segments_types.append([segment_nbr, style])
            self.channels_of[str(segment_nbr)] = channel_nodes(group)
        self.nbr_segments = len(self.segments_types)
        self.segments = {"extend": None, "retract": None, "pause": None,
                         "modulation": []}
        for segment_nbr, style in self.segments_types:
            if style in ("extend", "retract", "pause") and \
                    self.segments[style] is None:
                self.segments[style] = str(segment_nbr)
            elif style == "modulation":
                self.segments["modulation"].append(str(segment_nbr))
                self.modulation_seg_count += 1
        if self.segments["extend"] is None or \
                self.segments["retract"] is None:
            raise ValueError("No approach or retraction in the file")

        # Number of points of each curve, in each segment
        nbr_points = {}
        for segment_nbr, _ in self.segments_types:
            group = self.groups[segment_nbr]
            nbr_points[str(segment_nbr)] = numpy.asarray(
                group._f_get_child("meta-data")._f_get_child("num-points")
                [:, 0], int)
        approach = nbr_points[self.segments["extend"]]
        retraction = nbr_points[self.segments["retract"]]
        pause = nbr_points.get(self.segments["pause"])

        # Curves with data (corrupted_pixel is 1 for a valid curve)
        valid_curves = (approach > 0) & (retraction > 0)
        for segment in self.segments["modulation"] + [self.segments["pause"]]:
            if segment is not None:
                valid_curves &= nbr_points[segment] > 0
        self.corrupted_pixel = [int(flag) for flag in valid_curves]
        if not valid_curves.any():
            raise ValueError("No curve in the file")
        self.nbr_points_array_approach = \
            [int(n) if ok else 1 for n, ok in zip(approach, valid_curves)]
        self.nbr_points_array_retraction = \
            [int(n) if ok else 1 for n, ok in zip(retraction, valid_curves)]
        if pause is not None:
            self.nbr_points_array_pause = \
                [int(n) if ok else 1 for n, ok in zip(pause, valid_curves)]
        # Number of points of each modulation segment (from the first valid
        # curve, as for the other JPK files)
        first = int(numpy.argmax(valid_curves))
        self.nbr_points_array_modulation = [
            int(nbr_points[segment][first])
            for segment in self.segments["modulation"]]

        self.data.nbr_points_per_curve_approach = \
            int(max(self.nbr_points_array_approach))
        self.data.nbr_points_per_curve_retraction = \
            int(max(self.nbr_points_array_retraction))
        if self.nbr_points_array_pause:
            self.max_nbr_points_pause = \
                int(max(self.nbr_points_array_pause))
        self.data.nbr_points_per_pause_curve = self.max_nbr_points_pause
        if self.nbr_points_array_modulation:
            self.max_nbr_points_modulation = \
                int(max(self.nbr_points_array_modulation))
            self.data.modulation_nbr_points = \
                list(self.nbr_points_array_modulation)
        self.data.nbr_points_per_modulation_curve = \
            self.max_nbr_points_modulation

        # Duration of each segment (from the first valid curve)
        self.segment_durations = []
        for segment_nbr, _ in self.segments_types:
            durations = self.groups[segment_nbr]._f_get_child("meta-data") \
                ._f_get_child("duration")
            self.segment_durations.append(float(durations[first, 0]))
        self.data.segment_durations = self.segment_durations

        # Channels: the combined height if it is there (see HEIGHT_CHANNELS)
        channels = self.channels_of[self.segments["extend"]]
        for name in HEIGHT_CHANNELS:
            if all(name in self.channels_of[str(segment_nbr)]
                   for segment_nbr, _ in self.segments_types):
                self.height_channel_name = name
                break
        else:
            raise ValueError("No height channel in the file")
        if "vDeflection" not in channels:
            raise ValueError("No deflection channel in the file")
        self.data.zsensor_used = self.height_channel_name in SENSOR_CHANNELS

        # Conversion of the height (to m, net encoder of the channel)
        attributes = channels[self.height_channel_name][1]
        self.height_multiplier = float(
            attributes.get("net-encoder.scaling.multiplier", 1.0))
        self.height_offset = float(
            attributes.get("net-encoder.scaling.offset", 0.0))

        # Conversion of the deflection (to V, then the sensitivity and the
        # spring constant), as in the other JPK files
        index = None
        for key, name in header.items():
            if key.startswith("shared-data.lcd-info.") and \
                    key.endswith(".channel.name") and name == "vDeflection":
                index = key.split(".")[2]
        if index is None:
            raise ValueError("No conversion of the deflection in the file")
        pre = "shared-data.lcd-info." + index + "."
        conv = pre + "conversion-set.conversion."
        self.encoder_vDeflection_multiplier = number(
            pre + "encoder.scaling.multiplier")
        self.encoder_vDeflection_offset = number(
            pre + "encoder.scaling.offset")
        if value(conv + "distance.defined") == "true":
            self.vDeflection_distance_multiplier = \
                number(conv + "distance.scaling.multiplier") * 1e9
            self.vDeflection_distance_offset = \
                number(conv + "distance.scaling.offset", 0.0) * 1e9
        if value(conv + "force.defined") == "true":
            self.vDeflection_force_multiplier = \
                number(conv + "force.scaling.multiplier")
            self.vDeflection_force_offset = \
                number(conv + "force.scaling.offset", 0.0)

        self.data._original_deflection_sensitivity = \
            self.vDeflection_distance_multiplier
        self.data._original_spring_constant = self.vDeflection_force_multiplier
        # Temperatures are not saved for JPK, set -1 = None
        self.data._original_temperature = -1

        # Velocities
        if time_approach and time_retraction and self.data.ramp_size:
            self.data.approach_velocity = self.data.ramp_size / time_approach
            self.data.retraction_velocity = \
                self.data.ramp_size / time_retraction
            self.data.scan_rate = 1.0 / (time_approach + time_retraction)
        else:
            self.data.approach_velocity = 0
            self.data.retraction_velocity = 0

        # General data
        self.data.x_size = self.data.scan_size_x / self.data.nbr_pixels_x
        self.data.y_size = self.data.scan_size_y / self.data.nbr_pixels_y

        # Positions of the stage and of the scanners, grid on the sample
        frame = plant_tools.frame_from_header(header)
        self.data.frame = frame
        self.data.sample_grid = plant_tools.map_grid(
            {key: self.grid[key] for key in
             ("xcenter", "ycenter", "ulength", "vlength", "theta")}, frame)

        # Pixel of each curve on the map, and the region of the map (the
        # pixels without curve)
        self.place_curves()
        outside = numpy.ones((self.data.nbr_pixels_x,
                              self.data.nbr_pixels_y), int)
        for pixel in self.pixels:
            outside[pixel] = 0
        self.data.outside_region = outside if outside.any() else None

    def fetch_curves(self):
        """Loads the curves (see LoaderJPK), and the parts of the map
        between the moves of the z motor."""
        super().fetch_curves()
        self.data.motor_parts = self.motor_parts()

    def motor_parts(self):
        """Parts of the map between the moves of the z motor (number of the
        part of each pixel, -1 without curve; see math_tools.level_steps),
        None if the motor did not move or without combined height.

        The position of the motor is the combined height minus the height
        of the head (MOTOR_CHANNELS), at the end of the extend segment."""
        partner = MOTOR_CHANNELS.get(self.height_channel_name)
        channels = self.channels_of[self.segments["extend"]]
        if partner not in channels:
            return None
        attributes = channels[partner][1]
        multiplier = float(
            attributes.get("net-encoder.scaling.multiplier", 1.0))
        offset = float(attributes.get("net-encoder.scaling.offset", 0.0))

        parts = -numpy.ones((self.data.nbr_pixels_x, self.data.nbr_pixels_y),
                            int)
        part = -1
        last = None
        for index_nbr in range(self.data.real_nbr_pixels):
            if not self.corrupted_pixel[index_nbr]:
                continue
            last_point = self.nbr_points_array_approach[index_nbr] - 1
            combined = channels[self.height_channel_name][0][
                index_nbr, last_point]
            head = channels[partner][0][index_nbr, last_point]
            # nm
            motor = self.convert_height(combined) - \
                (head * multiplier + offset) * 1e9
            if last is None or abs(motor - last) > MOTOR_TOLERANCE:
                part += 1
            last = motor
            parts[self.pixels[index_nbr]] = part
        if part < 1:
            return None
        return parts

    def read_segment(self, index_nbr, segment, nbr_points):
        """Raw values of the height and of the deflection of a segment of
        a curve (see convert_height and convert_deflection)."""
        channels = self.channels_of[segment]
        return [numpy.asarray(channels[name][0][index_nbr, :nbr_points],
                              float)
                for name in (self.height_channel_name, "vDeflection")]

    def convert_height(self, values):
        """Height [nm] from the raw values of the height channel."""
        return (values * self.height_multiplier + self.height_offset) * 1e9
