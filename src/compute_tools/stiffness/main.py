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

"""
Computes the stiffness.

This class will call the methods to get the point of contact on the approach
curve, and then segment the indentation part of the curve to extract the
Young modulus (depending on the model).

"""

import multiprocessing
import numpy
from ...compute_tools import multiproc

if multiprocessing.current_process().name == "MainProcess":
    # Do not reimport this in the child processes.
    import os
    from ...tools import misc_tools
    from ...tools import math_tools
    from ...widgets.progressbar import Progressbar
    from ... import widgets_list
    from ... import shared
    import math


class Stiffness:
    """Computes the stiffness for all the curves.

    The code takes advantage of the multiprocessing module from python to be
    able to use multiple cores. Great care must be taken while accessing the
    data from the hdf5 file. The multiprocessing code was adapted from the
    multiprocess_access_queues.py example from pytables.

    The data can be read and written in blocks to reduce the memory footprint
    of the computation.

    In each block, queues are set up with make_queues. Then, DataProcessors are
    created (one per core), and the computation is run.

    The progressbar needs to be updated in this module (it can not be called
    from one of the child processes.

    To be compatible with windows, all the parameters for the fits are passed
    to the processors as arguments.

    Please read the documentation to get more details about the multiprocessing
    part.
    """

    def __init__(self):
        self.temp_file = shared.exp.temp_file
        self.max_indentation_index = 0
        self.list_of_smoothing_errors = []
        self.error = False

    def get_stiffness(self, calc_id):
        """Method used to launch the computation."""
        data = shared.exp.list[calc_id]
        # Be sure that the data is updated (when computing "all", the hdf5
        # file is closed en reopened many times, so that we need to call an
        # update here)
        data.update()

        misc_tools.save_used_params(data, "stiffness")

        # Create a progressbar which is displayed during the computations
        Progressbar(os.path.basename(data.filename))
        if data.stiffness_model_selected == 3:
            widgets_list.widget_progressbar.set_label("Getting stiffness")
        else:
            widgets_list.widget_progressbar.set_label("Getting elasticity")
        nbr_pixels = data.nbr_pixels_x * data.nbr_pixels_y
        widgets_list.widget_progressbar.set_range(0, nbr_pixels)

        # If we recalculate the stiffness, it is better to remove the data
        # and recreate the tables.
        if data.stiffness_calculated:
            self.temp_file.delete_tables_for_results(calc_id, "stiffness")
            self.temp_file.create_tables_for_results(calc_id, "stiffness")

        # Calculate the coefficients only once for all curves
        if data.stiffness_model_selected == 0:
            # Hertz (sphere)
            div = numpy.sqrt(data.tip_radius * pow(10, -9))
            coeff = (3.0 / 4.0) * ((1.0 - data.poisson_ratio ** 2) / div)
        elif data.stiffness_model_selected == 1:
            # Sneddon (cone)
            div = 2.0 * math.tan(math.radians(data.tip_angle))
            coeff = math.pi / div * (1.0 - data.poisson_ratio ** 2)
        elif data.stiffness_model_selected == 2:
            # Bilodeau (Pyramid)
            div = math.tan(math.radians(data.tip_angle))
            coeff = math.sqrt(2) * ((1.0 - data.poisson_ratio ** 2) / div)
        elif data.stiffness_model_selected == 3:
            # Slope
            coeff = 1.0
        elif data.stiffness_model_selected == 4:
            # Flat punch
            div = 2 * data.tip_radius * pow(10, -9)
            coeff = (1.0 - data.poisson_ratio ** 2) / div

        # Piezo image (last position of the extended piezo)
        piezo_image = \
            numpy.zeros([data.nbr_pixels_x, data.nbr_pixels_y], float)
        if data.flatten_applied:
            piezo_image = math_tools.flatten(data.piezo_image,
                                             data.applied_flatten_order)
        else:
            piezo_image[:, :] = data.piezo_image

        # Close the file, it will be re-opened only in the FileAcces method
        # from the multiprocessing procedure.
        self.temp_file.close_file()

        # Dictionnary containing parameters for the computation
        dt = {
            "nbr_pixels_x": data.nbr_pixels_x,
            "nbr_pixels_y": data.nbr_pixels_y,
            "piezo_image": piezo_image,
            "discarded_curves": data.discarded_curves,
            "tilt_limit_1": data.used_tilt_limit_1,
            "tilt_limit_2": data.used_tilt_limit_2,
            "deflection_sensitivity": data.deflection_sensitivity,
            "spring_constant": data.spring_constant,
            "tilt_applied": data.used_tilt_applied,
            "sg_smoothing_enabled": data.used_sg_smoothing_enabled,
            "sg_smoothing_order": data.used_sg_smoothing_order,
            "sg_smoothing_width": data.used_sg_smoothing_width,
            "sg_smoothing_uniform": data.used_sg_smoothing_uniform,
            "stretch_applied_app": data.used_stretch_applied_app,
            "stretch_applied_ret": data.used_stretch_applied_ret,
            "stretch_app_lim1": data.used_stretch_app_lim1,
            "stretch_app_lim2": data.used_stretch_app_lim2,
            "stretch_len_app": data.used_stretch_len_app,
            "stretch_ret_lim1": data.used_stretch_ret_lim1,
            "stretch_ret_lim2": data.used_stretch_ret_lim2,
            "stretch_len_ret": data.used_stretch_len_ret,
            "poc_skip_start": data.fitparam_poc_skip_start,
            "poc_fit_length": data.fitparam_poc_fit_length,
            "poc_refit_option": data.fitparam_poc_refit_option,
            "poc_noise_multiplicator": data.fitparam_poc_noise_multiplicator,
            "poc_refit_times": data.fitparam_poc_refit_times,
            "indentation_start": data.indentation_start,
            "indentation_stop": data.indentation_stop,
            "indentation_step": data.indentation_step,
            "fit_range_type": data.fit_range_type,
            "force_start": data.force_start,
            "force_stop": data.force_stop,
            "strict_stop": data.used_strict_stop,
            "model_selected": data.stiffness_model_selected,
            "tomography": data.tomography,
            "trig_threshold": data.trig_threshold,
            "tip_radius": data.tip_radius,
            "tip_angle": data.tip_angle,
            "poisson_ratio": data.poisson_ratio,
            "perform_fit": data.perform_fit,
            "coeff": coeff
            }

        # Compute and get the list of errors if there are some
        filepath = self.temp_file.filepath
        self.list_of_smoothing_errors = multiproc.compute(
            "stiffness", calc_id, dt, filepath)

        # Reopen the hdf5 file
        self.temp_file.open_file()

        st = "/data/_" + str(calc_id)
        tf = self.temp_file.file
        stiffness = tf.get_node(st + "/results", "stiffness_array")
        topography = tf.get_node(st + "/results", "topography")

        # Check for strict stop (NOT USED ANYMORE)
        strict_stop_error = False
        if data.used_strict_stop and data.indentation_stop != 0:
            dist = data.used_indentation_stop - data.used_indentation_start
            if data.used_indentation_step == 0:
                segments_per_curve = 1
            else:
                segments_per_curve = dist // data.used_indentation_step
            for i in range(len(stiffness)):
                if len(stiffness[i]) != segments_per_curve:
                    strict_stop_error = True
                    # No need to check further, we found an error
                    break

            # Delete the stiffness
            if strict_stop_error:
                widgets_list.widget_progressbar.close()
                self.error = "Error: strict stop"

                # Delete result table
                name = "stiffness"
                self.temp_file.delete_tables_for_results(str(calc_id), name)
                self.temp_file.create_tables_for_results(str(calc_id), name)

                return True

        # Calculate the maximum indentation index
        z_array = []
        max_ind = 0
        nosegment = True

        row = 0
        col = 0
        # Stiffness_array is a Vlarray !!! (2D and not 3D)
        # Stiffness array is stored as : [row, col]
        # Pocs topo array is stored as : [col, row]
        # So we have to shift the indices so that the correct length is
        # substracted from topography
        # Store the base to get a minimal value to correct the pocs array
        base = []
        step = data.indentation_step
        for i in range(len(stiffness)):
            k = len(stiffness[i])
            z_array.append(k * data.indentation_step)
            # At least one value was found
            if k != 0 and nosegment:
                nosegment = False
            if k > max_ind:
                max_ind = k
            val = topography[row][col] - (len(stiffness[i])) * step
            base.append(val)
            col = col + 1
            if col == data.nbr_pixels_y:
                col = 0
                row += 1

        self.max_indentation_index = max_ind
        if nosegment:
            # The step was too big or we wanted to calculate too far
            # z_array is only populated wit 0, 0 segments
            widgets_list.widget_progressbar.close()
            self.error = "Error: no segments"

            # Delete result table
            self.temp_file.delete_tables_for_results(str(calc_id), "stiffness")
            self.temp_file.create_tables_for_results(str(calc_id), "stiffness")

            return True

        # Correct POC with piezo height (piezo_image) for topography and
        # stiffness display. Get a corrected POC array (when 3D is
        # reconstructed, the last point in z is = 0)
        minimum = numpy.amin(base)
        for i in range(data.nbr_pixels_x):
            for j in range(data.nbr_pixels_y):
                topography[i, j] = topography[i][j] - minimum

        # Flush the file and update the values in the experiment class
        self.temp_file.flush_file()
        data.stiffness_calculated = True

        # Close the progressbar and return
        widgets_list.widget_progressbar.close()

        return True
