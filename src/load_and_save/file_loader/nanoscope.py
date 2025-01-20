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
Loads Nanoscope files (unique curves, force Volume and Peak Force)

Originally written by Charles Roduit (See OpenFovea's source code), heavily
modified by Michka Popoff (added support for Peak Force, Nanoscope 8, and
integration in PYAF).

"""

import re
import numpy
import struct
from PyQt5 import QtWidgets
from dateutil.parser import parse
from struct import unpack
from ... import widgets_list
from ... import shared


class LoaderNanoscope:
    """Class used to load Nanoscope files."""

    def __init__(self, file_info, unique_id):
        self.unique_id = unique_id
        self.path = file_info["path"]
        self.file_type = file_info["file_type"]
        self.file_name = file_info["filename"]
        self.afm_file = open(self.path, "rb")

        # Some empty values
        self.matrix_length = 0
        self.bpp_used_force = None
        self.bpp_used_sensor = None
        self.force_data_length = 0
        self.sensor_data_length = 0
        self.force_data_offset = 0
        self.one_image = 0
        self.hard_z_scale = 0
        self.image_scale = 0
        self.image_offset = 0
        self.line = None
        self.peakforce = None

        # Fetch the data
        self.fetch_header()
        shared.exp.temp_file.create_data_tables(self.unique_id)
        self.fetch_image()
        self.fetch_curves()

    def fetch_header(self):
        """Saves all the informations in the header."""
        data = shared.exp.list[self.unique_id]
        position = ""
        # Sometimes when you change the spring constant manually in
        # nanoscope, there a two data offsets. The first is the
        # right one. (Another Nanoscope Bug). This variable lets you select the
        # first one. (For the force curves AND the piezo_image)
        piezo_image_data_offset_found = False
        force_data_offset_found = False
        # Initialization of 2 regular expressions used for finding first
        # instance of sync distances
        reqnm = sync_distance_re("qnm")
        renew = sync_distance_re("new")

        self.line = self.afm_file.readline().decode("latin_1")

        while "\\*File list end" not in self.line:
            # Position in file
            if self.line.find("\\Version:") + 1:
                data.version = extract_num(self.line, -1)
            elif "\\*Force file list" in self.line:
                position = "File Description"
            elif "\\*Equipment list" in self.line:
                position = "Equipment list"
            elif "\\*Ciao scan list" in self.line:
                position = "Scan list"
            elif "\\*Ciao force list" in self.line:
                position = "Force list"
            elif "\\*Ciao image list" in self.line:
                position = "Image list"
            elif "\\*Ciao force image list" in self.line:
                position = "Force image"
            elif position == "File Description":
                if self.line.find("\\Date:") + 1:
                    data.date = str(parse(self.line[7:-2]))
            # Description of the equipment
            elif position == "Equipment list":
                if self.line.find("\\@Sens. Zscan") + 1:
                    data.sens_z_scan = float(extract_num(self.line))
                elif self.line.find("\\@Sens. Zsens: V") + 1:
                    data.sens_z_scan = float(extract_num(self.line))
                elif self.line.find("\\Microscope") + 1:
                    data.microscope_name = extract_text(self.line)
                elif self.line.find("\\Scanner file") + 1:
                    data.scanner_file = extract_text(self.line)

            # Parameters of the scan
            elif position == "Scan list":
                if self.line.find("\\Scan Size") + 1:
                    # Sometimes there is no scan size
                    # (scanning on the same pixel multiple times).
                    # We set a default value for this (scan size = 1 um)
                    scan_size = extract_num(self.line)
                    if scan_size == 0.0:
                        data.scan_size_x = 0
                        data.scan_size_y = 0
                    else:
                        data.scan_size_x = scan_size
                        data.scan_size_y = scan_size
                elif self.line.find("\\Rotate Ang.") + 1:
                    data.scan_angle = extract_num(self.line)
                elif self.line.find("\\@Sens. Deflection: V") + 1:
                    # Not sure this value still exists in recent
                    # nanoscope software ?
                    val = extract_num(self.line)
                    data._original_deflection_sensitivity = val
                elif self.line.find("\\@Sens. DeflSens: V") + 1:
                    val = extract_num(self.line)
                    data._original_deflection_sensitivity = val
                elif self.line.find("\\XY Closed Loop") + 1:
                    data.xy_closed_loop = extract_text(self.line)
                elif self.line.find("\\Z Closed Loop") + 1:
                    data.z_closed_loop = extract_text(self.line)
                elif self.line.find("\\PeakForce Capture") + 1:
                    val = extract_text(self.line)
                    if val == "Allow":
                        self.peakforce = True
                elif self.line.find("\\Peak Force Amplitude") + 1:
                    # data.ramp_size = float(extract_num(self.line)) * 2.0
                    ramp_size_g = float(extract_num(self.line)) * 2.0
                elif self.line.find("\\PFT Freq") + 1:
                    # Get the Peak Force Frequency in Hz
                    self.pft_freq = 1000 * float(extract_num(self.line))
                elif self.line.find("\\Sample Points") + 1:
                    self.sample_curve_pft = int(extract_num(self.line))
                elif reqnm.match(self.line):
                    self.sync_dist_qnm = int(extract_num(self.line))
                elif renew.match(self.line):
                    self.sync_dist_new = int(extract_num(self.line))
                elif self.line.find("\\Samps/line") + 1:
                    data.piezo_image_samps_per_line = extract_num(self.line)
                # Get sensibility of Zsensor
                elif self.line.find("\\@Sens. ZsensSens") + 1:
                    self.sens_Z_sensor = extract_num(self.line)

            # Parameters of the acquisition of the force curves
            elif position == "Force list":
                if self.line.find("\\force/line") + 1:
                    self.matrix_length = extract_num(self.line)
                    if self.file_type == "Nanoscope (Single File)":
                        data.nbr_pixels_x = 1
                        data.nbr_pixels_y = 1
                    else:
                        data.nbr_pixels_x = extract_num(self.line)
                        data.nbr_pixels_y = extract_num(self.line)
                elif self.line.find("\\Scan rate") + 1:
                    data.scan_rate = extract_num(self.line)
                elif self.line.find("\\Forward vel") + 1:
                    val = extract_num(self.line) * data.sens_z_scan
                    data.approach_velocity = val
                elif self.line.find("\\Reverse vel") + 1:
                    val = extract_num(self.line) * data.sens_z_scan
                    data.retraction_velocity = val
                elif self.line.find("\\Hold Time") + 1:
                    # Starting from software version 9.20b43 the Ramp Delay
                    # parameter is called Hold Time
                    data.extended_delay = extract_num(self.line)
                elif self.line.find("\\Ramp delay") + 1:
                    data.extended_delay = extract_num(self.line)
                elif self.line.find("\\Reverse delay") + 1:
                    data.retracted_delay = extract_num(self.line)
                elif self.line.find("\\Temperature (Celsius)") + 1:
                    val = extract_num(self.line) + 273.15
                    data._original_temperature = val
                elif self.line.find("\\@4:Trig threshold Deflection: V") + 1 or \
                        self.line.find("\\@4:Trig Threshold Deflection: V") + 1:
                    val = extract_num(self.line, special="trig_threshold")
                    data.trig_threshold = val

            # Parameters of the image
            elif position == "Image list":
                # Go inside the first "Image list"
                if not self.one_image:
                    if self.line.find("\\Data length") + 1:
                        # pass
                        self.image_length = extract_num(self.line)
                        # Not needed ?
                    elif self.line.find("\\Plane fit") + 1:
                        data.plane_fit = extract_plane_fit(self.line,
                                                           "plane_fit")
                    elif self.line.find("\\Frame direction") + 1:
                        data.frame_direction = extract_text(self.line)
                    elif self.line.find("\\Realtime Planefit") + 1:
                        data.realtime_plane_fit = extract_plane_fit(self.line)
                    elif self.line.find("\\@2:Image Data: S") + 1:
                        self.map_channel = extract_bracket(self.line, 3)
                        # Move 2 lines further
                        self.line = self.afm_file.readline().decode("latin_1")
                        self.line = self.afm_file.readline().decode("latin_1")
                        if self.line.find("\\@2:Z scale:") + 1:
                            self.image_scale = extract_num(self.line, -1)
                    elif self.line.find("\\Samps/line") + 1:
                        data.image_samps_per_line = extract_num(self.line)
                    elif self.line.find("\\Number of lines") + 1:
                        data.image_number_lines = extract_num(self.line)
                    elif self.line.find("\\Data offset") + 1:
                        # Sometimes when you change the spring constant
                        # manually in nanoscope, there are a two data offsets
                        # for the piezo_image.
                        # The first is the right one
                        if not piezo_image_data_offset_found:
                            if not self.one_image:
                                self.image_offset = extract_num(self.line)
                            else:
                                self.image_offset = 0
                            piezo_image_data_offset_found = True
                    elif self.line.find("\\Bytes/pixel") + 1:
                        # This value is the actual number of bytes in which
                        # masured numbers (piezo voltage, zsensor readings,
                        # deflection etc...) are encoded.
                        # For the data storage anyway, the allocated memory
                        # for each pixel can be equal (e.g. 2 or 4) or greater
                        # if for example data are encoded as 2 bytes but stored
                        # as 4.
                        # The bytes allocated for each image are calculated
                        # from the data lenght for each channel
                        # (in the following).
                        self.bpp_used_image = extract_num(self.line, -1)
                    # Find the last line of \*Ciao force image list and stop
                    # reading the other one
                    elif self.line.find("\\@2:Z offset: V") + 1:
                        self.one_image = 1

            # Parameters of the force curves
            elif position == "Force image":
                if self.line.find("\\Offline Planefit") + 1:
                    offline_plane_fit_g = extract_plane_fit(self.line)
                elif self.line.find("\\Data length:") + 1:
                    force_data_length_g = extract_num(self.line, -1)
                elif self.line.find("\\Data offset:") + 1:
                    # Sometimes when you change the spring constant manually in
                    # nanoscope, there a two data offsets. The first is the
                    # right one
                    if not force_data_offset_found:
                        force_data_offset_g = extract_num(self.line, -1)
                        force_data_offset_found_g = True

                elif self.line.find("\\Bytes/pixel:") + 1:
                    # This number is related to data encoding and NOT to data
                    # storage, that can differ
                    force_bytes_per_pixel_g = extract_num(self.line, -1)

                elif self.line.find("\\Samps/line") + 1:
                    # For Nanoscope file, the number of points per curve
                    # are the same for approach and retraction.

                    val = extract_num(self.line, 0)
                    nbr_points_per_curve_approach_g = val
                    nbr_points_per_curve_retraction_g = val

                    # Real number of points per curve can be less than nbr
                    # of points per curve.

                    if data.version >= 20:
                        # Starting from software version 9.40b20 there was a change on this parameter.
                        index = 0
                    else:
                        index = 1

                    val = extract_num(self.line, index)
                    nbr_points_per_curve_approach_real_g = val
                    nbr_points_per_curve_retraction_real_g = val

                elif self.line.find("\\Data type: FORCE") + 1:
                    pass
                    # In the original version there was a NFJflag=1 ?
                    # This one is from Charles Roduit, don't know if it's
                    # important. Michka
                elif self.line.find("\\@4:Z scale: V [Sens. DeflSens]") + 1:
                    hard_z_scale_g = extract_num(self.line, pos=1)
                elif self.line.find("\\@4:FV scale: V [Sens. ZsensSens]") + 1:
                    hard_z_scale_g = extract_num(self.line, -1)
                elif (self.line.find("\\@4:Ramp size: V") + 1 or
                      self.line.find("\\@4:Ramp Size: V") + 1) and \
                        self.file_type != "Nanoscope (Peak Force)":
                    val = extract_num(self.line, -1) * data.sens_z_scan
                    ramp_size_g = val
                elif self.line.find("\\Spring Constant") + 1:
                    _original_spring_constant_g = extract_num(self.line, -1)

                elif self.line.find("\\@4:Image Data") + 1:
                    channel = extract_bracket(self.line, 3)

                # Find the last line of \*Ciao force image list
                elif self.line.find("\\@4:Z Display") + 1 or \
                        self.line.find("\\@4:Z display:") + 1 or \
                        self.line.find("\\@4:Ramp End:") + 1:
                    # Stock value of the different channel
                    if channel == "ZSensor":
                        data.zsensor_used = True
                        self.force_data_offset_Zs = force_data_offset_g
                        self.hard_z_scale_Zs = hard_z_scale_g
                        self.bpp_used_sensor = force_bytes_per_pixel_g
                        self.sensor_data_length = force_data_length_g

                    if channel == "DeflectionError":
                        data.offline_plane_fit = offline_plane_fit_g
                        self.force_data_length = force_data_length_g
                        self.force_data_offset = force_data_offset_g
                        self.force_data_offset_found = \
                            force_data_offset_found_g
                        data.nbr_points_per_curve_approach = \
                            nbr_points_per_curve_approach_g
                        data.nbr_points_per_curve_retraction = \
                            nbr_points_per_curve_retraction_g
                        data.nbr_points_per_curve_approach_real = \
                            nbr_points_per_curve_approach_real_g
                        data.nbr_points_per_curve_retraction_real = \
                            nbr_points_per_curve_retraction_real_g
                        data.hard_z_scale = hard_z_scale_g
                        data.ramp_size = ramp_size_g
                        data._original_spring_constant = \
                            _original_spring_constant_g
                        self.bpp_used_force = \
                            force_bytes_per_pixel_g

            # Go to next line
            self.line = self.afm_file.readline().decode("latin_1")

        # General data
        if data.nbr_pixels_x != 1:
            data.x_size = data.scan_size_x / data.nbr_pixels_x
        if data.nbr_pixels_y != 1:
            data.y_size = data.scan_size_y / data.nbr_pixels_y

        # There is a bug in nanoscope 7.3.
        # Sometimes the Bytes per pixel is not stored in the file
        if data.nbr_points_per_curve_approach is None:
            self.bpp_used_force = 2

        data.matrix_length = self.matrix_length
        data.number_curves = self.matrix_length ** 2

    def fetch_image(self):
        """Method to return the piezo height image.

        The piezo height is the last position of the piezo on the approach
        curve. It is used to get the relative positions of the curves.

        if headerimage_samps_per_line is different than "nbr_pixels_x",
        there are too many points in the piezo image. We skip these points with
        the seek function in the while loop below. From the nanoscope doc :
        <<Thus, the first image pixel matches the first force curve.>>
        """
        data = shared.exp.list[self.unique_id]

        error_positions = []

        # Get the piezo image array from the temp file
        tf = shared.exp.temp_file.file
        dt = "/data/_" + str(self.unique_id)
        piezo_image = tf.get_node(dt + "/piezo_image", "piezo_image")
        shape = [data.matrix_length, data.matrix_length, 1]
        temp_piezo_image = numpy.zeros(shape, numpy.float64)

        # seek(pos, mode) with mode = 0 = absolute pos
        self.afm_file.seek(self.image_offset, 0)

        # Check if we have the same number of pixels and force curves
        # Skip(bytes) = (nbr of pixels to skip/2)
        skip = int(numpy.round(((data.piezo_image_samps_per_line / data.nbr_pixels_x) - 1) * 2))

        if self.image_offset != 0:
            pos_y = 0
            pos_x = 0

            # Use the correct sensibility : Piezo or Sensor sensibility (depend
            # of the map channel)
            if self.map_channel == "ZSensor":
                z_sens = self.sens_Z_sensor
            else:
                z_sens = data.sens_z_scan

            # Use the good parameters to unpack the binary file and rescale data
            # bpp_allocated is the bytes/pixel allocated for data storage and is
            # used for a correct data reading.
            # bpp_used is the bytes/pixel in which data are actually encoded and
            # is used for rescaling the data in nm.
            bpp_allocated_image = self.image_length // \
                                  (data.image_number_lines * data.image_samps_per_line)

            bpp_used = self.bpp_used_image

            if bpp_allocated_image == 2:
                fmt = '<h'
            elif bpp_allocated_image == 4:
                fmt = '<i'

            scaling_factor = self.image_scale * z_sens / (2. ** (bpp_used * 8))

            while pos_y < data.nbr_pixels_y:
                while pos_x < data.nbr_pixels_x:

                    try:
                        hdata = unpack(fmt, self.afm_file.read(bpp_allocated_image))
                        temp_piezo_image[pos_x, pos_y] = hdata[0] * scaling_factor

                    except struct.error:
                        # For handling when the software encounters an error decoding a byte
                        # of the height map the while loop finishes and only the data successfully
                        # unpacked till that point is loaded. The software currently does not handle
                        # missing values.
                        shared.exp.missing_z_positions.append((pos_x, pos_y))
                        break


                    if skip != 0:
                        self.afm_file.seek(skip, 1)
                    pos_x += 1
                pos_x = 0
                pos_y += 1

                # Update progressbar every line
                widgets_list.widget_progressbar.update()

            if shared.exp.missing_z_positions:
                # Create a message box
                str_error_pos = ""
                for a, b in shared.exp.missing_z_positions:
                    str_error_pos = str_error_pos + "<li> " + f"{a}, {b}" + "</li>"

                text = f"While unpacking the piezo height image from file {self.file_name} there was an error unpacking data on positions: <ul>" \
                       + str_error_pos + "</ul>These positions will be skiped."

                msg_corr = QtWidgets.QMessageBox()
                msg_corr.setText("Error (Failed to unpack data)")
                msg_corr.setInformativeText(text)
                msg_corr.setIcon(QtWidgets.QMessageBox.Critical)
                msg_corr.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
                msg_corr.exec_()

            # Go to 0
            piezo_min = numpy.amin(temp_piezo_image)
            pos_y = 0
            pos_x = 0
            while pos_y < data.nbr_pixels_y:
                while pos_x < data.nbr_pixels_x:
                    piezo_image[pos_x, pos_y] = \
                        temp_piezo_image[pos_x][pos_y] - piezo_min
                    pos_x += 1
                pos_x = 0
                pos_y += 1

                # Update progressbar every line
                widgets_list.widget_progressbar.update()

        else:
            # Save empty piezo image
            piezo_image = temp_piezo_image

    def fetch_curves(self):
        """Method to get the approach and retraction curves.

        When a peak force file is read, the two curves have to be cut in 2
        (to get in total 4 segments), which have to be stitched together to get
        a classical force curve (See below).
        """
        data = shared.exp.list[self.unique_id]

        # Get the curves from the temp file
        tf = shared.exp.temp_file.file
        dt = "/data/_" + str(self.unique_id)
        curves_approach = tf.get_node(dt + "/curves", "approach")
        curves_retraction = tf.get_node(dt + "/curves", "retraction")
        app_positions = tf.get_node(dt + "/positions", "approach_positions")
        ret_positions = tf.get_node(dt + "/positions", "retraction_positions")
        scaling_factor = data.hard_z_scale
        f_samples = data.nbr_points_per_curve_approach

        # Get sensor values if we have acces to Z-Sensor data
        if data.zsensor_used:
            curves_x_app, curves_x_ret = self.fetch_approach()
        else:
            # Create a new calculated curve_x
            curve_x = numpy.arange(f_samples) * (data.ramp_size / f_samples)

        self.afm_file.seek(self.force_data_offset, 0)

        pos_y = 0
        pos_x = 0

        curve_av = numpy.zeros([f_samples])
        curve_re = numpy.zeros([f_samples])

        # Define provisory arrays to store the raw curve
        if self.peakforce:
            prov_app = numpy.zeros([f_samples])
            prov_re = numpy.zeros([f_samples])
            # I am forced to impose this condition, since the sync_distance value
            # that is measured in pixel, is referred to the theoretical samples
            # per curve. Anyway the samples per curve are limited to a max
            # of 1024 for memory reasons, so for frequencies lower than 0.5 kHz,
            # given a sampling frequency of 500 kHz (on Resolve), there is a
            # discrepancy.
            if f_samples != self.sample_curve_pft:
                pft_factor = self.sample_curve_pft / (2 * f_samples)
                self.sync_dist_new = self.sync_dist_new / pft_factor
                self.sync_dist_qnm = self.sync_dist_qnm / pft_factor

        # Check the real start value for single force curves
        if data.nbr_points_per_curve_approach_real != f_samples:
            # Single force curves can start at a random place
            end = data.nbr_points_per_curve_approach_real
        else:
            # Force Volume
            end = f_samples - 1

        # Use the good parameters to unpack the binary file and rescale data.
        # bpp_allocated is the bytes/pixel allocated for data storage and is
        # used for a correct data reading.
        # bpp_used is the bytes/pixel in which data are actually encoded and
        # is used for rescaling the data in nm. In the case of Deflection signal
        # the scaling factor is already stored into the header file, so bpp_used
        # is not necessary.
        bpp_allocated_force = self.force_data_length // \
                              (2 * f_samples * self.matrix_length ** 2)

        if self.file_type == "Nanoscope (Single File)":
            bpp_allocated_force = self.force_data_length // \
                                  (2 * f_samples)

        if bpp_allocated_force == 2:
            fmt = 'h'
        elif bpp_allocated_force == 4:
            fmt = 'i'

        while pos_y < data.nbr_pixels_y:
            while pos_x < data.nbr_pixels_x:
                # Load the x, y curve
                # Unpack (< means little endian, h means short)
                # Endianess is important if you are on a different platform
                # (Veeco/Bruker files are recorded on Win XP)

                if self.peakforce:
                    # Peak-force curves have to be cut in 4 and stitched
                    # together to produce a normal approach/retration curve

                    prov_re[:] = unpack("<" + str(f_samples) + fmt,
                                        self.afm_file.read(bpp_allocated_force
                                                           * f_samples))
                    prov_app[:] = unpack("<" + str(f_samples) + fmt,
                                         self.afm_file.read(bpp_allocated_force
                                                            * f_samples))

                    # For peakforce files I first read the binary data and then
                    # I build 2 arrays: the sinusoidal z movement (curve_x) that
                    # I calculate and the deflection (curve_pft) using the data
                    # I read. These curves are approach + retraction.
                    # I then locate the maximum on the curve_pft and cut the 2
                    # arrays in 2 parts: approach and retract.
                    # The peak is usually not centered so 1 of the 2 curves is
                    # normally longer than the other. Since pyAF does not manage
                    # cases in which the lenght of the 2 curves is not the same,
                    # I add constant values to the sortest curve, that are then
                    # removed by the filter checking for corrupted parts.
                    # Since, after calibration, the Z movement is syncronized
                    # using the sync_distance_qnm, this parameter is used to
                    # reconstruct a curve_x with the proper phase shift deltat.

                    curve_pft = numpy.zeros([2 * f_samples])
                    curve_pft = numpy.concatenate((prov_re[::-1], prov_app))
                    max_force_index = numpy.argmax(curve_pft)

                    # sd = self.sync_dist_new / (self.pft_freq * 2 * f_samples)
                    sd = self.sync_dist_qnm / (self.pft_freq * 2 * f_samples)
                    deltat = sd - 1 / (self.pft_freq * 4)
                    curve_t = numpy.arange(2 * f_samples) * \
                              ((0.5 / self.pft_freq) / f_samples)

                    curve_x = (data.ramp_size / 2.) * \
                              numpy.sin(2 * numpy.pi * self.pft_freq *
                                        (curve_t - deltat))

                    # This array is used to add dummy points to the shorter of
                    # the 2 curves, since they need to have an equal number of
                    # points because of the structure of the curves array.
                    defl_temp = numpy.empty(abs(f_samples - max_force_index))

                    # I made these 2 conditions (and not just the first one),
                    # just in case the peak force may shift.
                    if max_force_index < f_samples:
                        curve_x_ret_pft = curve_x[(max_force_index):
                                                  (max_force_index + f_samples)]
                        curve_re = curve_pft[(max_force_index):
                                             (max_force_index + f_samples)]
                        cut_x = curve_x[0:(max_force_index)]
                        x_temp = range(abs(f_samples - max_force_index))
                        x_temp = x_temp + (cut_x[0] - x_temp[-1])
                        curve_x_app_pft = numpy.concatenate((x_temp, cut_x))
                        cut_pft = curve_pft[0:(max_force_index)]

                        # This condition is because of 1 case in which a curve
                        # showed no interaction peak and a maximum at the second
                        # point
                        if len(cut_pft) == 1:
                            cut_pft[0] = curve_re[-1]
                        else:
                            cut_pft[0] = cut_pft[1]
                        curve_re = curve_re - cut_pft[0]
                        cut_pft = cut_pft - cut_pft[0]
                        defl_temp.fill(0)
                        curve_av = numpy.concatenate((defl_temp, cut_pft))
                    elif max_force_index == f_samples:
                        curve_x_ret_pft = curve_x[f_samples: 2 * f_samples]
                        curve_re = curve_pft[f_samples: 2 * f_samples]
                        curve_x_app_pft = curve_x[0:f_samples]
                        curve_av = curve_pft[0:f_samples]
                        curve_av[0] = curve_av[1]
                        curve_re = curve_re - curve_av[0]
                        curve_av = curve_av - curve_av[0]
                    else:
                        curve_x_app_pft = curve_x[(max_force_index - f_samples):
                                                  max_force_index]
                        curve_av = curve_pft[(max_force_index - f_samples):
                                             (max_force_index)]
                        cut_x = curve_x[(max_force_index)::]
                        x_temp = range(abs(f_samples - max_force_index))
                        x_temp = x_temp + (cut_x[-1] - x_temp[0])
                        curve_x_ret_pft = numpy.concatenate((cut_x, x_temp))
                        cut_pft = curve_pft[(max_force_index)::]
                        cut_pft = cut_pft - curve_av[0]
                        curve_av = curve_av - curve_av[0]
                        defl_temp.fill(0)
                        curve_re = numpy.concatenate((cut_pft, defl_temp))

                    curve_av = curve_av[::-1]
                    curve_re = curve_re[::-1]

                    # Check for flat parts in the curves
                    # I use this routine to do not consider the added dummy
                    # points into the calculations.
                    end_pos = end
                    for i in range(end, 1, -1):
                        if curve_av[i] == curve_av[i - 1]:
                            end_pos = i
                        else:
                            break
                    start_pos = 0
                    for i in range(end):
                        if curve_av[i] == curve_av[i + 1]:
                            start_pos = i
                        else:
                            break

                else:
                    curve_av[:] = unpack("<" + str(f_samples) + fmt,
                                         self.afm_file.read(bpp_allocated_force
                                                            * f_samples))

                    curve_re[:] = unpack("<" + str(f_samples) + fmt,
                                         self.afm_file.read(bpp_allocated_force
                                                            * f_samples))

                    # Check for flat parts in the curves
                    end_pos = end
                    for i in range(end, 1, -1):
                        if curve_av[i] == curve_av[i - 1]:
                            end_pos = i
                        else:
                            break
                    start_pos = 0
                    for i in range(end):
                        if curve_av[i] == curve_av[i + 1]:
                            start_pos = i
                        else:
                            break

                # Store the start and end positions for the flat parts of the
                # curves.
                app_positions[pos_x, pos_y, :] = \
                    [f_samples - end_pos + 1, f_samples - start_pos]
                ret_positions[pos_x, pos_y, :] = [0, f_samples]

                if data.zsensor_used:
                    curve_x_app = curves_x_app[pos_x, pos_y]
                    curve_x_ret = curves_x_ret[pos_x, pos_y]
                elif self.peakforce:
                    curve_x_app = curve_x_app_pft
                    curve_x_ret = curve_x_ret_pft
                else:
                    curve_x_app = curve_x
                    curve_x_ret = curve_x

                # Fill arrays
                curves_approach[pos_x, pos_y, :, :] = \
                    [curve_x_app, curve_av[::-1] * scaling_factor]
                curves_retraction[pos_x, pos_y, :, :] = \
                    [curve_x_ret, curve_re[::-1] * scaling_factor]

                pos_x += 1

                # Update progressbar every line
                widgets_list.widget_progressbar.update()

            pos_x = 0
            pos_y += 1

        # Force write to disk
        shared.exp.temp_file.flush_file()

    def fetch_approach(self):
        """Method to get the x coordinate for approach and retraction curves when
        we have acces to the Z-sensor data
        """
        data = shared.exp.list[self.unique_id]

        f_samples = data.nbr_points_per_curve_approach

        # Make empty vectors
        curve_x_appraoch = numpy.zeros(
            [data.matrix_length, data.matrix_length, f_samples])
        curve_x_retract = numpy.zeros(
            [data.matrix_length, data.matrix_length, f_samples])

        # Use the good parameters to unpack the binary file and rescale data
        # bpp_allocated is the bytes/pixel allocated for data storage and is
        # used for a correct data reading.
        # bpp_used is the bytes/pixel in which data are actually encoded and
        # is used for rescaling the data in nm.
        bpp_allocated_sensor = self.sensor_data_length // \
                               (2 * f_samples * self.matrix_length ** 2)
        bpp_used = self.bpp_used_sensor

        if self.file_type == "Nanoscope (Single File)":
            bpp_allocated_sensor = self.force_data_length // \
                                   (2 * f_samples)

        if bpp_allocated_sensor == 2:
            fmt = 'h'
        elif bpp_allocated_sensor == 4:
            fmt = 'i'

        # Zsensor scaling_factor

        zsens_scaling_factor = self.hard_z_scale_Zs * \
                               self.sens_Z_sensor / (2. ** (bpp_used * 8))

        # Goto to the good offset in the file
        self.afm_file.seek(self.force_data_offset_Zs, 0)

        prov_app = numpy.zeros([f_samples])
        prov_ret = numpy.zeros([f_samples])

        pos_y = 0
        pos_x = 0

        while pos_y < data.nbr_pixels_y:
            while pos_x < data.nbr_pixels_x:
                prov_app[:] = unpack("<" + str(f_samples) + fmt,
                                     self.afm_file.read(bpp_allocated_sensor
                                                        * f_samples))

                prov_ret[:] = unpack("<" + str(f_samples) + fmt,
                                     self.afm_file.read(bpp_allocated_sensor
                                                        * f_samples))

                min_ret = numpy.amin(prov_ret)
                min_app = numpy.amin(prov_app)

                if min_app < min_ret:
                    prov_app = prov_app - min_ret
                    prov_ret = prov_ret - min_ret

                elif min_ret < min_app:
                    prov_app = prov_app - min_app
                    prov_ret = prov_ret - min_app

                # Fill arrays
                curve_x_retract[pos_x, pos_y] = prov_ret[::-1] * \
                                                zsens_scaling_factor
                curve_x_appraoch[pos_x, pos_y] = prov_app[::-1] * \
                                                 zsens_scaling_factor

                pos_x += 1

            pos_x = 0
            pos_y += 1

        return curve_x_appraoch, curve_x_retract


def extract_num(line, pos=0, special=None):
    """Extract number from a string."""
    if special is None:
        number = re.findall("([0-9]+[.0-9]*)", line)[pos]
    elif special == "trig_threshold":
        number = line.rstrip("\r\n").split()[-2]
    try:
        number = int(number)
    except ValueError:
        number = float(number)

    return number


def extract_text(line):
    """Extract text from line."""
    line = line.split(":")  # Get the whole line
    text = line[1]  # After the split get the text
    text = text[1:].rstrip("\r\n")  # Strip the left space and \n\r

    return text


def extract_plane_fit(line, name=None):
    """Extract plane fit data from line."""
    if name == "plane_fit":
        line = line.split()
        del line[0:2]
        newline = []
        for i in range(len(line)):
            if line[i][0] == "-":
                newline.append(-float(line[i][1:]))
            else:
                newline.append(float(line[i]))
        line = newline

    return str(line)


def extract_bracket(line, position):
    """Extract text into extract_brackets"""
    text = line.split(' ')
    text = text[position][1:-1]

    return text


def sync_distance_re(sync_name):
    """Uses regular expressions to find first istance of sync distances (peakforce)
    into header file
    """
    if sync_name == "new":
        resync = re.compile(r"\\Sync Distance New: [0-9]*")
    if sync_name == "qnm":
        resync = re.compile(r"\\Sync Distance QNM: [0-9]*")

    return resync
