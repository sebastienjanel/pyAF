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

"""Module to compute the work and rupture fore on the retraction curves."""

import multiprocessing
from .. import multiproc

if multiprocessing.current_process().name == "MainProcess":
    # Do not reimport this in the child processes.
    import os
    from ...tools import misc_tools
    from ...widgets.progressbar import Progressbar
    from ... import widgets_list
    from ... import shared


class WorkAndRuptureForce:
    """Will fit the retraction curve and get the work and maximum rupture force.

    A linear fit will be used, defined by parameters inputed by the user in the
    GUI. A progressbar is displayed during the computation. The results are
    saved to the temporary file.
    """

    def __init__(self):
        self.temp_file = shared.exp.temp_file

    def get_work_and_rupture_force(self, calc_id):
        """Method to calculate the work and the rupture force."""
        data = shared.exp.list[calc_id]
        # Be sure that the data is updated (when computing "all", the hdf5
        # file is closed en reopened many times, so that we need to call an
        # update here)
        data.update()

        misc_tools.save_used_params(data, "work_and_rupture_force")

        # Create a progressbar which is displayed during the computations
        Progressbar(os.path.basename(data.filename))
        wl = widgets_list
        wl.widget_progressbar.set_label("Getting work and rupture force")
        nbr_pixels = data.nbr_pixels_x * data.nbr_pixels_y
        wl.widget_progressbar.set_range(0, nbr_pixels)

        # If we recalculate the work, it is better to remove the data
        # and recreate the tables.
        if data.work_and_rupture_force1_calculated:
            self.temp_file.delete_tables_for_results(calc_id, "work")
            self.temp_file.create_tables_for_results(calc_id, "work")

        # Close the file, it will be re-opened only in the FileAcces method
        # from the multiprocessing procedure.
        self.temp_file.close_file()

        # Dictionnary containing parameters for the computation
        dt = {
            "nbr_pixels_x": data.nbr_pixels_x,
            "nbr_pixels_y": data.nbr_pixels_y,
            "discarded_curves": data.discarded_curves,
            "tilt_limit_1": data.used_tilt_limit_1,
            "tilt_limit_2": data.used_tilt_limit_2,
            "deflection_sensitivity": data.deflection_sensitivity,
            "spring_constant": data.spring_constant,
            "tilt_applied": data.used_tilt_applied,
            "trig_threshold": data.trig_threshold,
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
            "joc_skip_start": data.fitparam_joc_skip_start,
            "joc_fit_length": data.fitparam_joc_fit_length,
            "joc_refit_option": data.fitparam_joc_refit_option,
            "joc_noise_multiplicator": data.fitparam_joc_noise_multiplicator,
            "joc_refit_times": data.fitparam_joc_refit_times}

        multiproc.compute("work", calc_id, dt, self.temp_file.filepath)

        # Reopen the hdf5 file
        self.temp_file.open_file()

        data.work_and_rupture_force1_calculated = True

        # Close progress bar
        wl.widget_progressbar.close()
