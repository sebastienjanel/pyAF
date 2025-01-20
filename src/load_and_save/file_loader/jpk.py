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

"""Loads JPK files (unique curves, force maps and QI files)

Only well formated JPK files can be loaded. If for example you used the
experiment planner to make 5 forces curves per pixel with pauses and different
settings, PYAF will not be able to open it. PYAF can only open files with
one approach and one retraction curve per pixel (pauses are allowed).
"""
import struct

import numpy
from dateutil.parser import parse
from struct import unpack
import zipfile
from PyQt5 import QtWidgets
import pandas as pd
import h5py

from ... import widgets_list
from ... import shared


class LoaderJPK:
    """Class used to load JPK files."""

    def __init__(self, file_info, unique_id):
        self.unique_id = unique_id
        self.path = file_info["path"]
        self.file_type = file_info["file_type"]
        self.file_name = file_info["filename"]
        self.afm_file = zipfile.ZipFile(self.path, "r")
        # print(self.afm_file.namelist())

        # Some empty values
        self.vDeflection_base_multiplier = 1.0
        self.vDeflection_base_offset = 0.0
        self.vDeflection_distance_multiplier = 1.0
        self.vDeflection_distance_offset = 0.0
        self.vDeflection_force_multiplier = 1.0
        self.vDeflection_force_offset = 0.0
        self.capSensHeight_abs_mult = 1.0
        self.capSensHeight_abs_offset = 0.0
        self.capSensHeight_nom_mult = 1.0
        self.capSensHeight_nom_offset = 0.0
        self.segments_types = []
        self.segment_durations = []
        self.channels = []
        self.max_nbr_points = None
        self.max_nbr_points_pause = None
        self.max_nbr_points_modulation = None
        self.segments = None
        self.enc_capSensHeight_offset = None
        self.encoder_vDeflection_offset = None
        self.nbr_points_array_approach = None
        self.nbr_points_array_retraction = None
        self.nbr_points_array_pause = None
        self.encoder_vDeflection_multiplier = None
        self.enc_capSensHeight_mult = None
        self.nbr_segments = None
        self.height_channel_name = ""
        self.corrupted_pixel = []

        # Variables for modulation segments 
        self.segment_freqs = []
        self.segment_amps = []
        self.modulation_seg_count = 0

        # Pixels positions and way of parsing the data
        self.pos_x = 0
        self.pos_y = 0
        self.go_right = True

        # The data
        self.data = shared.exp.list[self.unique_id]

        # Fetch the data
        self.fetch_header()
        shared.exp.temp_file.create_data_tables(self.unique_id)
        self.fetch_curves()

    def fetch_header(self):
        """Saves all the informations in the header."""

        time_approach = None
        time_retraction = None

        if self.file_type == "JPK (Force Map)":
            prefix = "force-scan-map"
            pre_header = ".settings"
        elif self.file_type == "JPK (QI)":
            prefix = "quantitative-imaging-map"
            pre_header = ".settings"
        elif self.file_type == "JPK (Single File)":
            prefix = "force-scan-series"
            pre_header = ".header"
            self.data.nbr_pixels_x = 1
            self.data.nbr_pixels_y = 1
            self.data.scan_size_x = 0
            self.data.scan_size_y = 0

        ft = self.file_type

        # How many segments have been saved
        path = "shared-data/header.properties"
        for line in self.afm_file.read(path).splitlines():
            line = line.decode()
            if line.find("force-segment-header-infos.count") != -1:
                self.nbr_segments = extract_int(line)
                break

        # It no nbr_segments have been found in shared-data, look in the
        # main header (for single JPK files for example)
        path = "header.properties"
        for line in self.afm_file.read(path).splitlines():
            line = line.decode()
            if line.find("force-segments.count") != -1:
                self.nbr_segments = extract_int(line)
                break

        for line in self.afm_file.read("header.properties").splitlines():
            line = line.decode()
            if line.find(prefix + ".description.instrument") != -1:
                self.data.microscope_name = extract_str(line)
            if line.find(prefix + ".description.source-software") != -1:
                self.data.version = extract_str(line)
            path = prefix + ".settings.force-settings.retracted-pause-time"
            if line.find(path) != -1:
                self.data.retracted_delay = extract_float(line)
            path = prefix + ".settings.force-settings.extended-pause-time"
            if line.find(path) != -1:
                self.data.extended_delay = extract_float(line)
            if line.find(prefix + ".start-time") != -1:
                self.data.date = extract_jpk_date(line)
            if line.find(prefix + ".position-pattern.grid.theta") != -1:
                self.data.scan_angle = extract_float(line)
            if line.find(prefix + ".position-pattern.grid.ilength") != -1:
                self.data.nbr_pixels_x = extract_int(line)
            if line.find(prefix + ".position-pattern.grid.jlength") != -1:
                self.data.nbr_pixels_y = extract_int(line)
            if line.find(prefix + ".position-pattern.grid.ulength") != -1:
                self.data.scan_size_x = extract_float(line) * pow(10, 9)
            if line.find(prefix + ".position-pattern.grid.vlength") != -1:
                self.data.scan_size_y = extract_float(line) * pow(10, 9)
            path = prefix + ".settings.force-settings.closed-loop"
            if line.find(path) != -1:
                if extract_str(line) == "true":
                    self.data.z_closed_loop = "On"
                if extract_str(line) == "false":
                    self.data.z_closed_loop = "Off"
            if line.find(prefix + ".indexes.max") != -1:
                # indexes.min is (for all the files we had) set to 0, so the
                # number of pixels is max + 1, because the counting of the
                # pixels starts at 0.
                self.data.real_nbr_pixels = extract_int(line) + 1
            path = prefix + ".settings.force-settings.extended-pause-time"
            if line.find(path) != -1:
                self.data.pause_duration = extract_float(line)

            # Line added for loading extended pause when the file is created
            # using the Ramp ForceDesigner.
            # In this case the segments are identified using increasing
            # numbers: approach = 0, extended pause (if any) = 1 etc...
            # In this case we look to the extended pause as the segment # 1:
            # other cases, like multiple extended pauses or retracted pause,
            # are currently not supported.

            path = prefix + ".settings.force-settings.segment.1.duration"
            if line.find(path) != -1:
                self.data.pause_duration = extract_float(line)

            # Ramp size
            if ft == "JPK (Force Map)" or ft == "JPK (Single File)":
                if line.find("relative-z-start") != -1:
                    relative_z_start = extract_float(line) * pow(10, 9)
                if line.find("relative-z-end") != -1:
                    relative_z_end = extract_float(line) * pow(10, 9)
                    self.data.ramp_size = relative_z_start - relative_z_end

            elif ft == "JPK (QI)":
                if line.find("settings.force-settings.extend.z-start") != -1:
                    relative_z_start = extract_float(line) * pow(10, 9)
                if line.find("settings.force-settings.extend.z-end") != -1:
                    relative_z_end = extract_float(line) * pow(10, 9)
                    self.data.ramp_size = relative_z_start - relative_z_end

            # Duration of the ramp movement
            pr = prefix + pre_header + ".force-settings."
            if ft == "JPK (Force Map)" or ft == "JPK (Single File)":
                if line.find(pr + "extend-scan-time") != -1:
                    time_approach = extract_float(line)
                if line.find(pr + "retract-scan-time") != -1:
                    time_retraction = extract_float(line)

            elif ft == "JPK (QI)":
                if line.find(pr + "extend.duration") != -1:
                    time_approach = extract_float(line)
                if line.find(pr + "retract.duration") != -1:
                    time_retraction = extract_float(line)

            # Setpoint (Force)
            if line.find(pr + "relative-setpoint") != -1:
                self.data.trig_threshold = extract_float(line)

        if self.data.ramp_size is None:
            # Sometimes there is no relative-z-start, for example if you
            # used the experiment planner ...
            # Loop through the segments and find the first one which has
            # a ramp size
            found = False
            for line in self.afm_file.read("header.properties").splitlines():
                line = line.decode()
                for seg in range(self.nbr_segments):
                    if line.find(str(seg) + ".z-start") != -1:
                        if not found:
                            relative_z_start = extract_float(line) * pow(10, 9)
                    if line.find(str(seg) + ".z-end") != -1:
                        if not found:
                            relative_z_end = extract_float(line) * pow(10, 9)
                            self.data.ramp_size = \
                                relative_z_start - relative_z_end
                            found = True
                        else:
                            print("The ramp size was set as the ramp size of" +
                                  " the first segment. You have multiple " +
                                  "segments in your file and all the " +
                                  "features are not completely supported " +
                                  "by pyAF.")

        if ft == "JPK (Single File)":
            # In case of a single file there is only one pixel
            self.data.real_nbr_pixels = 1

            # The date is stored in the first segment for single force curves
            path = "segments/0/segment-header.properties"
            afile = self.afm_file.read(path)
            for line in afile.splitlines():
                line = line.decode()
                val = "force-segment-header.time-stamp"
                if line.find(val) != -1:
                    self.data.date = extract_jpk_date(line)
                    break

        # Find the segments
        if self.file_type != "JPK (Single File)":
            # For maps and QI, find the informations in
            # shared-data/header.properties
            for segment_nbr in range(self.nbr_segments):
                values = self.afm_file.read(
                    "shared-data/header.properties").splitlines()
                for line in values:
                    line = line.decode()
                    header_info_path = "force-segment-header-info." + str(segment_nbr) + ".settings.segment-settings."
                    val = "force-segment-header-info." + str(segment_nbr) + \
                          ".settings.segment-settings.style"

                    if line.find(val) != -1:
                        # Check if the files exists, at least for the first
                        # pixel. Sometimes you can have more segments in the
                        # header (for example 4), but only 2 segments in the
                        # folders. In this case only add the valid segments

                        try:
                            path = "index/0/segments/" + str(segment_nbr)
                            valinfo = path + "/channels/vDeflection.dat"

                            # Check if the file is present, this may fail
                            # (See explanation above)
                            self.afm_file.getinfo(valinfo)

                            # Do the actual segment type extraction
                            seg = extract_str(line)
                            self.segments_types.append([segment_nbr, seg])
                                
                        except KeyError:
                            print("You changed the number of points per "
                                  "curve during the scan. The file was "
                                  "loaded but data may be missing.")

                    val_2 = header_info_path + "duration"
                    if line.find(val_2) != -1:
                        duration = extract_float(line)
                        self.segment_durations.append(duration)

                    val_3 = header_info_path + "frequency"
                    if line.find(val_3) != -1: 
                        frequency = extract_float(line)
                        self.segment_freqs.append(frequency)
                    
                    val_4 = header_info_path + "amplitude"
                    if line.find(val_4) != -1:
                        amplitude = extract_float(line)
                        self.segment_amps.append(amplitude)

        else:
            # For single files, this data is stored in
            # segments/0/segment-header.properties
            for segment_nbr in range(self.nbr_segments):
                # Sometimes we have more segments than real segments folders
                # for example 4 segments but only 2 segments saved
                path = ("segments/" + str(segment_nbr) +
                        "/segment-header.properties")

                try:
                    # Try to open the segment to see if its present.
                    afile = self.afm_file.read(path)

                    for line in afile.splitlines():
                        line = line.decode()
                        header_info_path = "force-segment-header.settings.segment-settings."
                        
                        val = header_info_path + "style"
                        if line.find(val) != -1:
                            seg = extract_str(line)
                            self.segments_types.append([segment_nbr, seg])

                        val_2 = header_info_path + "duration"
                        if line.find(val_2) != -1:
                            duration = extract_float(line)
                            self.segment_durations.append(duration)
                        
                        val_3 = header_info_path + "frequency"
                        if line.find(val_3) != -1:
                                frequency = extract_float(line)
                                self.segment_freqs.append(frequency)
                        
                        val_4 = header_info_path + "amplitude"
                        if line.find(val_4) != -1:
                                amplitude = extract_float(line)
                                self.segment_amps.append(amplitude)

                except KeyError:
                    pass

        self.data.segment_durations = self.segment_durations

        self.segments = {"extend": None, "retract": None, "pause": None, "modulation": []}

        if len(self.segments_types) != self.nbr_segments:
            shared.exp.segment_handling.append([self.file_name, self.nbr_segments, len(self.segments_types)])

            # Reassign number of segments to the number of segments found.
            self.nbr_segments = len(self.segments_types)

        # print(self.segments_types)
        # print(self.segment_durations)
        # print(self.segment_freqs)
        # print(self.segment_amps)
        for item in self.segments_types:
            if item[1] == "extend":
                self.segments[item[1]] = str(item[0])
            elif item[1] == "retract":
                self.segments[item[1]] = str(item[0])
            elif item[1] == "pause":
                self.segments[item[1]] = str(item[0])
            elif item[1] == "modulation":
                self.segments[item[1]].append(str(item[0]))
                self.modulation_seg_count += 1
        # print(self.segments)
        self.nbr_points_array_approach = []
        self.nbr_points_array_retraction = []
        self.nbr_points_array_pause = []
        self.nbr_points_array_modulation = []

        # Dividing by 4 because data is written as 4 bytes values
        # Make a list with the file sizes, we cant rely on the headers because
        # there are sometimes errors with the number of points per curve

        for index_nbr in range(self.data.real_nbr_pixels):
            segments_per_pixel = self.nbr_segments
            if self.file_type != "JPK (Single File)":
                path = "index/" + str(index_nbr) + "/header.properties"
                afile = self.afm_file.read(path)
                for line in afile.splitlines():
                    line = line.decode()
                    if line.find("force-segments.count") != -1:
                        segments_per_pixel = extract_int(line)

            if segments_per_pixel == self.nbr_segments:
                self.corrupted_pixel.append(1)

                if self.file_type != "JPK (Single File)":
                    segment_path = "index/" + str(index_nbr) + "/segments/" + \
                                   self.segments["extend"] + "/"
                else:
                    segment_path = "segments/" + self.segments["extend"] + "/"

                self.nbr_points_array_approach.append(self.afm_file.getinfo(
                    segment_path + "channels/vDeflection.dat").file_size // 4)

                if self.file_type != "JPK (Single File)":
                    segment_path = "index/" + str(index_nbr) + "/segments/" + \
                                   self.segments["retract"] + "/"
                else:
                    segment_path = "segments/" + self.segments["retract"] + "/"

                self.nbr_points_array_retraction.append(self.afm_file.getinfo(
                    segment_path + "channels/vDeflection.dat").file_size // 4)

                if self.segments["pause"] is not None:
                    if self.file_type != "JPK (Single File)":
                        segment_path = "index/" + str(index_nbr) + "/segments/" + \
                                       self.segments["pause"] + "/"
                    else:
                        segment_path = "segments/" + self.segments["pause"] + "/"

                    self.nbr_points_array_pause.append(self.afm_file.getinfo(
                        segment_path + "channels/vDeflection.dat").file_size // 4)

                if self.segments["modulation"] is not []:
                    for i in range(len(self.segments["modulation"])):
                        if self.file_type != "JPK (Single File)":
                            segment_path = "index/" + str(index_nbr) + "/segments/" + \
                                           self.segments["modulation"][i] + "/"
                        else:
                            segment_path = "segments/" + self.segments["modulation"][i] + "/"

                        self.nbr_points_array_modulation.append(self.afm_file.getinfo(
                            segment_path + "channels/vDeflection.dat").file_size // 4)

                    # print(self.nbr_points_array_modulation)

            else:
                self.corrupted_pixel.append(0)

                self.nbr_points_array_approach.append(1)

                self.nbr_points_array_retraction.append(1)

                if self.segments["pause"] is not None:
                    self.nbr_points_array_pause.append(1)

                if self.segments["modulation"] != []:
                    self.nbr_points_array_modulation.append(1)

        max_approach = numpy.amax(self.nbr_points_array_approach)
        max_retraction = numpy.amax(self.nbr_points_array_retraction)

        # Fill the rest of the number of points
        max_curves = self.data.nbr_pixels_x * self.data.nbr_pixels_y
        for index_nbr in range(self.data.real_nbr_pixels, max_curves):
            self.nbr_points_array_approach.append(max_approach)
            self.nbr_points_array_retraction.append(max_retraction)

        self.data.nbr_points_per_curve_approach = int(
            numpy.amax(self.nbr_points_array_approach))
        self.data.nbr_points_per_curve_retraction = int(
            numpy.amax(self.nbr_points_array_retraction))

        # Number of points in the pause segment
        if self.nbr_points_array_pause != []:
            self.max_nbr_points_pause = int(
                numpy.amax(self.nbr_points_array_pause))
        self.data.nbr_points_per_pause_curve = self.max_nbr_points_pause

        # Number of points in the modulation segment
        if self.nbr_points_array_modulation != []:
            self.max_nbr_points_modulation = int(
                numpy.amax(self.nbr_points_array_modulation))
        self.data.nbr_points_per_modulation_curve = self.max_nbr_points_modulation

        header = self.afm_file.read(
            "shared-data/header.properties").splitlines()

        # Get the channels / infos
        # How many channels have been saved
        path = "shared-data/header.properties"
        for line in self.afm_file.read(path).splitlines():
            line = line.decode()
            if line.find("lcd-infos.count") != -1:
                nbr_channels = extract_int(line)
                break

        # In the older versions of the software, the Z sensor axis was called
        # capacitiveSensorHeight.
        # In an updated version of the software, following a changing in the
        # name, 2 identical channels of the Z sensor were saved:
        # capacitiveSensorHeight and measuredHeight.
        # In the last software (update date 22.04.2015), only the
        # measuredHeight is left. In the following I then check for the Z
        # sensor names and, in the case the 2 are there, I just load the
        # measuredHeight one, so I eliminate the line corresponding to
        # capacitiveSensorHeight from the list self.channels.

        capacitive_sensor_height_channel_found = False
        measured_height_channel_found = False
        height_channel_found = False

        for channel_id in range(nbr_channels):
            shared_header = "shared-data/header.properties"
            for line in self.afm_file.read(shared_header).splitlines():
                line = line.decode()
                channel = "lcd-info." + str(channel_id) + ".channel.name"
                if line.find(channel) != -1:
                    channel_name = extract_str(line)
                    if channel_name == "capacitiveSensorHeight":
                        capacitive_sensor_height_channel_found = True
                    if channel_name == "measuredHeight":
                        measured_height_channel_found = True
                    if channel_name == "height":
                        height_channel_found = True  # For working with felix data, no height sensor.

                    # Check if the channel file is present in the first
                    # segment of the first pixel, because this is not always
                    # the case ...
                    segment_nbr = "0"  # The first segment
                    if self.file_type != "JPK (Single File)":
                        # Look in self.segments["extend"] to find the path
                        # for the extend segment, which is always present.
                        path = (
                                "index/" + segment_nbr +
                                "/segments/" + self.segments["extend"] +
                                "/channels/" + channel_name + ".dat")
                    else:
                        path = ("segments/" + segment_nbr + "/channels/" +
                                channel_name + ".dat")

                    try:
                        self.afm_file.read(path)
                        self.channels.append([channel_id, channel_name])
                    except KeyError:
                        pass

        if capacitive_sensor_height_channel_found:
            self.height_channel_name = "capacitiveSensorHeight"
        # If there is a measuredHeight channel; use it instead of the old
        # capacitiveSensorHeight channel.
        if measured_height_channel_found:
            self.height_channel_name = "measuredHeight"

        if height_channel_found and not (capacitive_sensor_height_channel_found or measured_height_channel_found):
            self.height_channel_name = "height"  # For working with felix data, no height sensor.

        # Two height channels have been found; keep only the measured height
        if capacitive_sensor_height_channel_found and \
                measured_height_channel_found:
            for channel_id, channel_name in enumerate(self.channels):
                if channel_name == "capacitiveSensorHeight":
                    del self.channels[channel_id]

        # The encoders let you transform the raw data to floats
        for channel in self.channels:
            for line in header:
                line = line.decode()
                pre = "lcd-info." + str(channel[0])
                if channel[1] == "vDeflection":
                    if line.find(pre + ".encoder.scaling.offset") != -1:
                        self.encoder_vDeflection_offset = extract_float(line)
                    if line.find(pre + ".encoder.scaling.multiplier") != -1:
                        val = extract_float(line)
                        self.encoder_vDeflection_multiplier = val
                elif channel[1] == "capacitiveSensorHeight":
                    if line.find(pre + ".encoder.scaling.offset") != -1:
                        val = extract_float(line)
                        self.enc_capSensHeight_offset = val
                    if line.find(pre + ".encoder.scaling.multiplier") != -1:
                        val = extract_float(line)
                        self.enc_capSensHeight_mult = val
                elif channel[1] == "measuredHeight":
                    if line.find(pre + ".encoder.scaling.offset") != -1:
                        val = extract_float(line)
                        self.enc_capSensHeight_offset = val
                    if line.find(pre + ".encoder.scaling.multiplier") != -1:
                        val = extract_float(line)
                        self.enc_capSensHeight_mult = val
                elif channel[1] == "height":
                    if line.find(pre + ".encoder.scaling.offset") != -1:
                        val = extract_float(line)
                        self.enc_capSensHeight_offset = val
                    if line.find(pre + ".encoder.scaling.multiplier") != -1:
                        val = extract_float(line)
                        self.enc_capSensHeight_mult = val

        # Are there some conversions to do ?
        # print(self.channels)
        for channel in self.channels:
            if channel[1] == "vDeflection":
                v_deflection_index = str(channel[0])
            if channel[1] == "capacitiveSensorHeight":
                height_index = str(channel[0])
            if channel[1] == "measuredHeight":
                height_index = str(channel[0])
            if channel[1] == 'height':  # For working with felix data, no height sensor.
                height_index = str(channel[0])

        # vDeflection
        base_defined = False
        pre1 = "lcd-info." + v_deflection_index + ".conversion-set"
        pre2 = "lcd-info." + height_index + ".conversion-set"
        convd = ".conversion.distance"
        convf = ".conversion.force"
        conva = ".conversion.absolute"
        convn = ".conversion.nominal"
        for line in header:
            line = line.decode()
            if line.find(pre1 + ".conversions.base") != -1:
                base = extract_str(line)
                if line.find(pre1 + ".conversion." + base + ".defined") != -1:
                    if extract_str(line) == "true":
                        base_defined = True

            if base_defined:
                print("Loading error : The base conversion is set to true, "
                      "you should check the values loaded.")

        distance_defined = False
        for line in header:
            line = line.decode()
            if line.find(pre1 + convd + ".defined") != -1:
                if extract_str(line) == "true":
                    distance_defined = True
        if distance_defined:
            for line in header:
                line = line.decode()
                if line.find(pre1 + convd + ".scaling.offset") != -1:
                    val = extract_float(line) * 1e9  # in nm/V
                    self.vDeflection_distance_offset = val
                if line.find(pre1 + convd + ".scaling.multiplier") != -1:
                    val = extract_float(line) * 1e9  # in nm/V
                    self.vDeflection_distance_multiplier = val

        force_defined = False
        for line in header:
            line = line.decode()
            if line.find(pre1 + convf + ".defined") != -1:
                if extract_str(line) == "true":
                    force_defined = True
        if force_defined:
            for line in header:
                line = line.decode()
                if line.find(pre1 + convf + ".scaling.offset") != -1:
                    self.vDeflection_force_offset = extract_float(line)
                if line.find(pre1 + convf + ".scaling.multiplier") != -1:
                    self.vDeflection_force_multiplier = extract_float(line)

        # capacitiveSensorHeight_index or measuredHeight_index
        base_defined = False
        for line in header:
            line = line.decode()
            if line.find(pre2 + ".conversions.base") != -1:
                base = extract_str(line)
                if line.find(pre2 + ".conversion." + base + ".defined") != -1:
                    if extract_str(line) == "true":
                        base_defined = True

            if base_defined:
                print("Loading error : The base conversion is set to true, "
                      "you should check the values loaded.")

        absolute_defined = False
        for line in header:
            line = line.decode()
            if line.find(pre2 + conva + ".defined") != -1:
                if extract_str(line) == "true":
                    absolute_defined = True
            if absolute_defined:
                for line in header:
                    line = line.decode()
                    if line.find(pre2 + conva + ".scaling.offset") != -1:
                        self.capSensHeight_abs_offset = extract_float(line)
                    if line.find(pre2 + conva + ".scaling.multiplier") != -1:
                        self.capSensHeight_abs_mult = extract_float(line)

        nominal_defined = False
        for line in header:
            line = line.decode()
            if line.find(pre2 + convn + ".defined") != -1:
                if extract_str(line) == "true":
                    nominal_defined = True
        if nominal_defined:
            for line in header:
                line = line.decode()
                if line.find(pre2 + convn + ".scaling.offset") != -1:
                    self.capSensHeight_nom_offset = extract_float(line)
                if line.find(pre2 + convn + ".scaling.multiplier") != -1:
                    self.capSensHeight_nom_mult = extract_float(line)

        self.data._original_deflection_sensitivity = \
            self.vDeflection_distance_multiplier
        self.data._original_spring_constant = self.vDeflection_force_multiplier
        # Temperatures are not saved for JPK, set -1 = None
        self.data._original_temperature = -1

        # Check for missing times
        if time_approach is None:
            # If you used the experiment planner
            for line in self.afm_file.read("header.properties").splitlines():
                line = line.decode()
                for seg in self.segments_types:
                    if seg[1] == "extend":
                        if line.find(str(seg[0]) + ".duration") != -1:
                            time_approach = extract_float(line)
                    if seg[1] == "retract":
                        if line.find(str(seg[0]) + ".duration") != -1:
                            time_retraction = extract_float(line)

        # Check for missing threshold
        if self.data.trig_threshold is None:
            for line in self.afm_file.read("header.properties").splitlines():
                line = line.decode()
                for seg in self.segments_types:
                    if seg[1] == "extend":
                        if line.find(str(seg[0]) + ".setpoint") != -1:
                            self.data.trig_threshold = extract_float(line)
            print(("The loaded setpoint is the setpoint of the extend " +
                   "segment. You saved multiple segments and this is not " +
                   "fully supported"))

        # Velocities
        if time_approach is not None and time_retraction is not None:
            self.data.approach_velocity = self.data.ramp_size / time_approach
            val = self.data.ramp_size / time_retraction
            self.data.retraction_velocity = val
        else:
            self.data.approach_velocity = 0
            self.data.retraction_velocity = 0

        # General data
        if self.data.nbr_pixels_x != 1:
            self.data.x_size = self.data.scan_size_x / self.data.nbr_pixels_x
        if self.data.nbr_pixels_y != 1:
            self.data.y_size = self.data.scan_size_y / self.data.nbr_pixels_y
        # Scan rate
        self.data.scan_rate = 1.0 / (time_approach + time_retraction)

    def fetch_curves(self):
        """Loads the curves from the zip file to the tmp hdf5 file.

        Curves indexes are stored like this in jpk files::
            6 7 8
            5 4 3
            0 1 2

        Curves indexes are stored like this in jpk QI files::
            6 7 8
            3 4 5
            0 1 2

        I store them as (i, j) values in the curves_approach and
        curves_retraction nodes of the hdf5 file::
            (2, 0)     (2, 1)     (2, 2)
            (1, 0)     (1, 1)     (1, 2)
            (0, 0)     (0, 1)     (0, 2)

        The piezo_image is the last position of the piezo (fully extended).
        It is like the piezo_image in the nanoscope files, but here the z
        position of the curves is relative to one starting z position.
        In nanoscope files each curve has to be adjusted with the piezo_image.

        The positions array is only an array to store the lenght of every curve
        to check if the curve is corrupted at load in ../load.py.
        For each pixel you have a list: [start, end], which gives you the
        position of the first and last valid points of the curve (because
        points may be missing, especially on the retraction curve).
        """
        tf = shared.exp.temp_file.file
        dt = "/data/_" + str(self.unique_id)
        curves_approach = tf.get_node(dt + "/curves", "approach")
        curves_retraction = tf.get_node(dt + "/curves", "retraction")
        piezo_image = tf.get_node(dt + "/piezo_image", "piezo_image")
        app_positions = tf.get_node(dt + "/positions", "approach_positions")
        ret_positions = tf.get_node(dt + "/positions", "retraction_positions")
        if self.max_nbr_points_pause is not None:
            curves_pause = tf.get_node(dt + "/curves", "pause")

        if self.max_nbr_points_modulation is not None:
            curves_modulation = tf.get_node(dt + "/curves", "modulation")

        shape = [self.data.nbr_pixels_x, self.data.nbr_pixels_y, 1]
        temp_piezo_image = numpy.zeros(shape, numpy.float64)

        ch_1 = "channels/" + self.height_channel_name + ".dat"
        ch_2 = "channels/vDeflection.dat"

        time_segments = self.get_time_segments()

        """ 
        # Decide if we want to implement data export to HDF5
        # hdf5 file name
        hdf_name = f"{self.file_name}.hdf5"
        hdf_file = h5py.File(hdf_name, "w")

        grp_metadata = hdf_file.create_group("metadata")

        # Append metadata to dataframe
        metadata = {
            "deflection_sens": self.data._original_deflection_sensitivity,
            "spring_cte": self.data._original_spring_constant,
            "exp_temp": self.data._original_temperature,
            "nbr_points_per_curve_approach": self.data.nbr_points_per_curve_approach,
            "nbr_points_pause": self.nbr_points_array_pause,
            "nbr_points_per_curve_modulation": self.max_nbr_points_modulation,
            "nbr_points_retract": self.data.nbr_points_per_curve_retraction,
            "app_positions": app_positions,
            "ret_positions": ret_positions
        }

        for key, value in metadata.items():
            if value is not None:
                grp_metadata.attrs[key] = value

        grp_data = hdf_file.create_group("data")
        """

        for index_nbr in range(self.data.real_nbr_pixels):
            if self.file_type != "JPK (Single File)":
                seg_path = "index/" + str(index_nbr) + "/segments/"
            else:
                seg_path = "segments/"

            # Approach ----------------------------
            segment_path = seg_path + self.segments["extend"] + "/"

            max_nbr_points = self.data.nbr_points_per_curve_approach

            nbr_points = self.nbr_points_array_approach[index_nbr]
            startpos = max_nbr_points - nbr_points
            sc_approach = numpy.zeros([3, max_nbr_points])

            approach_startpos = startpos

            if self.corrupted_pixel[index_nbr]:
                sc_approach[0][startpos:] = unpack(
                    ">" + str(nbr_points) + "i",
                    self.afm_file.read(segment_path + ch_1))

                sc_approach[1][startpos:] = unpack(
                    ">" + str(nbr_points) + "i",
                    self.afm_file.read(segment_path + ch_2))

                # Conversion of the values
                base_val = (sc_approach[0] * self.enc_capSensHeight_mult +
                            self.enc_capSensHeight_offset)
                absolute_val = (base_val * self.capSensHeight_abs_mult +
                                self.capSensHeight_abs_offset)

                sc_approach[0] = (absolute_val * self.capSensHeight_nom_mult +
                                  self.capSensHeight_nom_offset) * 1e9

                # No multiplication through
                # self.vDeflection_distance_multiplier
                # (defl sensitivity).
                # This is done during the displaying of the
                # curves, and during the calculations
                sc_approach[1] = (
                        sc_approach[1] * self.encoder_vDeflection_multiplier +
                        self.encoder_vDeflection_offset +
                        self.vDeflection_distance_offset)

                approach_startpos = startpos

                sc_approach[2] = time_segments["approach"]

            # Retraction ----------------------------
            segment_path = seg_path + self.segments["retract"] + "/"

            max_nbr_points = self.data.nbr_points_per_curve_retraction

            nbr_points = self.nbr_points_array_retraction[index_nbr]
            startpos = self.data.nbr_points_per_curve_retraction - nbr_points
            sc_retraction = numpy.zeros([3, max_nbr_points])

            retraction_startpos = startpos

            if self.corrupted_pixel[index_nbr]:
                sc_retraction[0][:nbr_points] = unpack(
                    ">" + str(nbr_points) + "i",
                    self.afm_file.read(segment_path + ch_1))
                sc_retraction[1][:nbr_points] = unpack(
                    ">" + str(nbr_points) + "i",
                    self.afm_file.read(segment_path + ch_2))

                # Conversion of the values
                base_val2 = (sc_retraction[0] * self.enc_capSensHeight_mult +
                             self.enc_capSensHeight_offset)
                absolute_val2 = (base_val2 * self.capSensHeight_abs_mult +
                                 self.capSensHeight_abs_offset)

                sc_retraction[0] = (absolute_val2 * self.capSensHeight_nom_mult +
                                    self.capSensHeight_nom_offset) * 1e9

                # No multiplication through
                # self.vDeflection_distance_multiplier
                # (defl sensitivity).
                # This is done during the displaying of the
                # curves, and during the calculations
                sc_retraction[1] = (
                        sc_retraction[1] * self.encoder_vDeflection_multiplier +
                        self.encoder_vDeflection_offset +
                        self.vDeflection_distance_offset)

                # Turn the retraction curve around to store it the right way
                sc_retraction[0] = sc_retraction[0][::-1]
                sc_retraction[1] = sc_retraction[1][::-1]

                sc_retraction[2] = time_segments["retraction"]

                # print(sc_approach[0])

            # Pause segment ----------------------------
            if self.max_nbr_points_pause is not None:
                segment_path = seg_path + self.segments["pause"] + "/"

                sc_pause = numpy.zeros([3, self.max_nbr_points_pause])
                nbr_points = self.nbr_points_array_pause[index_nbr]

                if self.corrupted_pixel[index_nbr]:
                    # Capacitive Sensor Height
                    sc_pause[0][:nbr_points] = unpack(
                        ">" + str(nbr_points) + "i",
                        self.afm_file.read(segment_path + ch_1))
                    # V deflection
                    sc_pause[1][:nbr_points] = unpack(
                        ">" + str(nbr_points) + "i",
                        self.afm_file.read(segment_path + ch_2))


                    # Conversion of the values
                    base_val2 = (sc_pause[0] * self.enc_capSensHeight_mult +
                                 self.enc_capSensHeight_offset)
                    absolute_val2 = (base_val2 * self.capSensHeight_abs_mult +
                                     self.capSensHeight_abs_offset)
                    sc_pause[0] = (absolute_val2 * self.capSensHeight_nom_mult +
                                   self.capSensHeight_nom_offset) * 1e9

                    val_mult = self.encoder_vDeflection_multiplier
                    sc_pause[1] = (sc_pause[1] * val_mult +
                                   self.encoder_vDeflection_offset +
                                   self.vDeflection_distance_offset)

                    sc_pause[2] = time_segments["pause"]

            # Modulation segment ----------------------------
            if self.max_nbr_points_modulation is not None:

                modulation_segments = []

                for i in range(len(self.segments["modulation"])):
                    segment_path = seg_path + self.segments["modulation"][i] + "/"

                    nbr_points = self.nbr_points_array_modulation[i]
                    sc_modulation = numpy.zeros([3, nbr_points])

                    if self.corrupted_pixel[index_nbr]:

                        # Capacitive Sensor Height
                        sc_modulation[0][:nbr_points] = unpack(
                            ">" + str(nbr_points) + "i",
                            self.afm_file.read(segment_path + ch_1))

                        # V deflection
                        sc_modulation[1][:nbr_points] = unpack(
                            ">" + str(nbr_points) + "i",
                            self.afm_file.read(segment_path + ch_2))

                        # Conversion of the values
                        base_val2 = (sc_modulation[0] * self.enc_capSensHeight_mult +
                                     self.enc_capSensHeight_offset)
                        absolute_val2 = (base_val2 * self.capSensHeight_abs_mult +
                                         self.capSensHeight_abs_offset)
                        sc_modulation[0] = (absolute_val2 * self.capSensHeight_nom_mult +
                                            self.capSensHeight_nom_offset) * 1e9

                        val_mult = self.encoder_vDeflection_multiplier
                        sc_modulation[1] = (sc_modulation[1] * val_mult +
                                            self.encoder_vDeflection_offset +
                                            self.vDeflection_distance_offset)

                        sc_modulation[2] = time_segments["modulation"][i]

                        modulation_segments.append(sc_modulation)

            # Store data ----------------------------
            # Check redundant x positions at start of the curve
            # Some points which are also not at the rigth place (x(i) < x(i+1))
            # print(sc_approach[0])
            i = approach_startpos
            if i < len(sc_approach[0]) - 1:
                while sc_approach[0][i] <= sc_approach[0][i + 1]:
                    approach_startpos = i + 1
                    i += 1
                    if i == len(sc_approach[0]) - 1:
                        break
            # print(approach_startpos)
            # Check redundant x positions at start of the retraction curve
            i = retraction_startpos
            if i < len(sc_retraction[0]) - 1:
                while sc_retraction[0][i] <= sc_retraction[0][i + 1]:
                    retraction_startpos = i + 1
                    i += 1
                    if i == len(sc_retraction[0]) - 1:
                        break

            app_positions[self.pos_x, self.pos_y, :] = \
                [approach_startpos, len(sc_approach[0])]
            ret_positions[self.pos_x, self.pos_y, :] = \
                [retraction_startpos, len(sc_retraction[0])]

            # Store piezo_image (last extend position)
            temp_piezo_image[self.pos_x, self.pos_y] = \
                sc_approach[0][len(sc_approach[0]) - 1]

            # Xzero is substracted from the curves, so that they are correctly
            # displayed (with positive values)
            # The xzero value is the first value on the retraction curve which
            # is not corrupted.
            xzero = sc_retraction[0][retraction_startpos]
            curves_approach[self.pos_x, self.pos_y, :, :] = \
                [xzero - sc_approach[0], sc_approach[1], sc_approach[2]]
            curves_retraction[self.pos_x, self.pos_y, :, :] = \
                [xzero - sc_retraction[0], sc_retraction[1], sc_retraction[2]]
            if self.max_nbr_points_pause is not None:
                # For pause segments we consider that there are no corrupted
                # parts (but this may be wrong, it was not tested with a lot
                # of files for the moment)
                xzero = sc_pause[0][0]
                curves_pause[self.pos_x, self.pos_y, :, :] = \
                    [xzero - sc_pause[0], sc_pause[1], sc_pause[2]]

            if self.max_nbr_points_modulation is not None:
                # For modulation segments we consider that there are no corrupted
                # parts (but this may be wrong, it was not tested with a lot
                # of files for the moment)
                xzero = sc_modulation[0][0]
                # print(xzero)
                # print(xzero.shape)
                curves_modulation[self.pos_x, self.pos_y] = \
                    [xzero - sc_modulation[0], sc_modulation[1], sc_modulation[2]]

            self.update_positions()

            # Update progressbar every curve
            widgets_list.widget_progressbar.update()

            """
            grp_curve = grp_data.create_group(f"{index_nbr}")
                
            xzero_1 = sc_retraction[0][retraction_startpos]
            approach = {
                "Height": xzero_1 -sc_approach[0],
                "vDeflection": sc_approach[1],
                "time": sc_approach[2]
            }

            subgrp_app = grp_curve.create_group("extend")
            subgrp_app.create_dataset("height", data=approach["Height"])
            subgrp_app.create_dataset("vDeflection", data=approach["vDeflection"])
            subgrp_app.create_dataset("time", data=approach["time"])

            if self.data.nbr_points_per_pause_curve is not None:
                pause = {
                    "Height": xzero_1 -sc_pause[0],
                    "vDeflection": sc_pause[1],
                    "time": sc_pause[2]
                }

                subgrp_pause = grp_curve.create_group("pause")
                subgrp_pause.create_dataset("height", data=pause["Height"])
                subgrp_pause.create_dataset("vDeflection", data=pause["vDeflection"])
                subgrp_pause.create_dataset("time", data=pause["time"])

            if self.max_nbr_points_modulation is not None:
                subgrp_modulation = grp_curve.create_group("modulation")
                for i in range(len(self.segments["modulation"])):
                    modulation = {
                        "Height": xzero_1 -modulation_segments[i][0],
                        "vDeflection": modulation_segments[i][1],
                        "time": modulation_segments[i][2]
                    }

                    subgrp = subgrp_modulation.create_group(self.segments["modulation"][i])
                    # Save frequency and amplitude for every segment as an attribute
                    subgrp.attrs["frequency"] = self.segment_freqs[i]
                    subgrp.attrs["amplitude"] = self.segment_amps[i]
                    # Create datasets for every modulation segment
                    subgrp.create_dataset("height", data=modulation["Height"])
                    subgrp.create_dataset("vDeflection", data=modulation["vDeflection"])
                    subgrp.create_dataset("time", data=modulation["time"])

            retract = {
                "segment_type": "retract",
                "segment_index": self.segments["retract"],
                "Height": xzero_1 -sc_retraction[0],
                "vDeflection": sc_retraction[1],
                "time": sc_retraction[2]
            }

            subgrp_retract = grp_curve.create_group("retract")
            subgrp_retract.create_dataset("height", data=retract["Height"])
            subgrp_retract.create_dataset("vDeflection", data=retract["vDeflection"])
            subgrp_retract.create_dataset("time", data=retract["time"])

        hdf_file.close()
        """

        # Fill empty curves with zeros
        max_curves = self.data.nbr_pixels_x * self.data.nbr_pixels_y

        for index_nbr in range(self.data.real_nbr_pixels, max_curves):
            # The number of points is updated for apporach and
            # retraction in order to account for the case in which the
            # sampling frequency and or the velocity is different in
            # approach and retraction.

            # Approach
            nbr_points = self.nbr_points_array_approach[index_nbr]
            defl = numpy.zeros(nbr_points).tolist()
            ext = list(range(nbr_points))
            curves_approach[self.pos_x, self.pos_y, :, :] = [ext, defl]

            # Retraciton
            nbr_points = self.nbr_points_array_retraciton[index_nbr]
            defl = numpy.zeros(nbr_points).tolist()
            ext = list(range(nbr_points))
            curves_retraction[self.pos_x, self.pos_y, :, :] = [ext, defl]

            self.update_positions()

        # Go to 0 in z for piezo_image
        piezo_min = numpy.amin(temp_piezo_image)
        for pos_x in range(self.data.nbr_pixels_x):
            for pos_y in range(self.data.nbr_pixels_y):
                piezo_image[pos_x, pos_y] = \
                    temp_piezo_image[pos_x][pos_y] - piezo_min

        # Force write to disk
        shared.exp.temp_file.flush_file()

    def update_positions(self):
        """Update the x and y pixel positions."""
        # Go to next curve
        if self.file_type == "JPK (Force Map)":
            # Move on to next curve (x axis)
            if self.go_right:
                self.pos_x = self.pos_x + 1
            else:
                self.pos_x = self.pos_x - 1
            # Move on to next line
            if self.pos_x == self.data.nbr_pixels_x and self.go_right:
                shared.exp.temp_file.flush_file()
                self.pos_x = self.pos_x - 1
                self.pos_y = self.pos_y + 1
                self.go_right = False
            elif self.pos_x < 0 and self.go_right is False:
                self.pos_x = 0
                self.pos_y = self.pos_y + 1
                self.go_right = True
        elif self.file_type == "JPK (QI)":
            if self.pos_x <= self.data.nbr_pixels_x - 1:
                # Move on to next curve (x axis)
                self.pos_x = self.pos_x + 1
            # Move on to next line
            if self.pos_x > self.data.nbr_pixels_x - 1:
                shared.exp.temp_file.flush_file()
                self.pos_x = 0
                self.pos_y = self.pos_y + 1

    def get_time_segments(self):
        time_segments = {"approach": None, "pause": None, "modulation": [], "retraction": None}

        # Extract durations for each segment and calculate start and end time for each segment.
        sd = self.data.segment_durations
        # print(sd)

        # print(sd)

        times = [sum(sd[:i]) for i in range(len(sd) + 1)]

        # print(times)

        # Get number of points for each segment.
        pa = self.data.nbr_points_per_curve_approach
        pr = self.data.nbr_points_per_curve_retraction
        pp = self.data.nbr_points_per_pause_curve
        pm = self.nbr_points_array_modulation

        # Distribute time evently for every segment.
        index = 0
        approach_t = numpy.linspace(times[index], times[index + 1], pa, endpoint=False)
        time_segments["approach"] = approach_t
        index = index + 1

        if pp is not None:
            pause_t = numpy.linspace(times[index], times[index + 1], pp, endpoint=False)
            time_segments["pause"] = pause_t
            index = index + 1

        if pm is not None or pm != []:
            for i in range(len(self.segments["modulation"])):
                modulation_t = numpy.linspace(times[index], times[index + 1], pm[i], endpoint=False)
                time_segments["modulation"].append(modulation_t)
                index = index + 1

        retraction_t = numpy.linspace(times[index], times[index + 1], pr, endpoint=False)
        time_segments["retraction"] = retraction_t

        return time_segments


def extract_str(string):
    """Extract str from line in jpk file."""
    return str(string.split("=")[1])


def extract_int(string):
    """Extract int from line in jpk file."""
    return int(string.split("=")[1])


def extract_float(string):
    """Extract float from line in jpk file."""
    return float(string.split("=")[1])


def extract_jpk_date(line):
    """Format the JPK timestamp.

    Format the JPK timestamp correctly so that it can be nicely
    displayed in PYAF.
    """
    date = parse(extract_str(line).replace("\\", ""), ignoretz=True)
    return str(date).split(".")[0]
