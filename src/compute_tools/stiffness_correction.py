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

"""Tool used to correct the stiffness by different models."""

import os
import numpy
import math
import logging
from .. import shared
from .. import widgets_list
from ..tools import curve_tools
from ..widgets.progressbar import Progressbar


class StiffnessCorrection:
    """Corrects the stiffness.

    The stiffness is multiplied by a coefficient depending on the distance
    between the topography and the substrate. One has to define a ROI in the
    roi manager, and set it as a glass substrate.

    The coefficient chosen depends on the model chosen (Dimitriadis, Chadwick),
    depending if the stiffness was computed with the Hertz model or the Sneddon
    model.
    2024 -> BEC models have been updated to Garcia Garcia 2018.

    Note on the units:
    indentation, tip_radius, height (h) are all kept in nm as they are divided
    to compute the X value. No need to convert them to meters, X has no
    dimension.
    """

    def __init__(self):
        self.temp_file = shared.exp.temp_file
        self.logger = logging.getLogger()

    def correct_stiffness(self, calc_id):
        """Method called by PYAF to do the actual stiffness correction."""
        data = shared.exp.list[calc_id]

        self.logger.debug("Correctiong stiffness, file %s", calc_id)
        self.logger.debug("Using ROI %s", data.roi_glass_id)

        # Create a progressbar which is displayed during the computations
        Progressbar(os.path.basename(data.filename))
        widgets_list.widget_progressbar.set_label("Correcting stiffness")
        nbr_pts = data.nbr_pixels_x * data.nbr_pixels_y
        widgets_list.widget_progressbar.set_range(0, nbr_pts)

        # If we recalculate the stiffness, it is better to remove the data
        # and recreate the tables.
        if data.stiffness_corrected:
            self.temp_file.delete_tables_for_results(
                calc_id, "stiffness_corrected")
            self.temp_file.create_tables_for_results(
                calc_id, "stiffness_corrected")

        # Get the table
        stiffness_array = self.temp_file.file.get_node(
            "/data/_" + str(calc_id) + "/results", "stiffness_corrected")

        # Create empty arrays to store intermediate results
        resstiffness_array = numpy.zeros([data.nbr_pixels_y], list)

        tan_angle = math.tan(math.radians(data.used_tip_angle))
        tip_radius = data.used_tip_radius

        # Get the coefficients of the glass surface if there is no user defined sample height
        if data.user_h is None or data.user_h == 0:
            d, a, b, _ = data.roi_list[data.roi_glass_id - 1].glass_coeffs

        x_size = data.scan_size_x / data.nbr_pixels_x
        y_size = data.scan_size_x / data.nbr_pixels_x

        for i in range(data.nbr_pixels_x):
            for j in range(data.nbr_pixels_y):

                if data.user_h is None or data.user_h == 0:
                    z_glass = a * (j * y_size) + b * (i * x_size) + d

                E_new = []
                for z in range(len(data.stiffness_array[i][j])):
                    E_original = data.stiffness_array[i][j][z]

                    if data.user_h is None or data.user_h == 0:
                        h = data.topography[i][j] - z_glass

                    else:
                        h = data.user_h

                    # First indentation must be at 1, not at 0
                    z += 1

                    if not data.used_tomography:
                        if len(data.stiffness_array[i][j]) == 1:
                            poc_pos = data.pocs_real[i][j]

                            # Get the deflection - extension curve
                            app, _, _, _, _, _ = curve_tools.get_curve(
                                data, [i, j], mode="curve_results")
                            # Get the force curve
                            force_curve = curve_tools.get_force_curves(
                                app[0], app[1], poc_pos, data.spring_constant)

                            # The last value of the force-distance curve is
                            # the maximal indentation
                            # (The point of contact is at the position 0, 0)
                            indentation = force_curve[0][-1]
                        else:
                            indentation = z * data.used_indentation_step
                    else:
                        # Kasas model, not sure what to do here ...
                        indentation = data.used_indentation_step
                        h = h - z * data.used_indentation_step

                    # For very rare cases on glass, where the indentation
                    # is negative
                    if indentation < 0:
                        indentation = 0

                    model = data.used_stiffness_model_selected

                    if model == 0:
                        # Hertz (paraboloid)

                        if h >= 0:
                            # Old Chadwick/Dimitriadis correction
                            # X = math.sqrt(indentation * tip_radius) / h
                            # o2 = 1.133 * X + 1.283 * X ** 2
                            # o4 = 0.769 * X ** 3 + 0.0975 * X ** 4
                            # Updated correction according to Garcia 2018
                            o1 = (1.133 * math.sqrt(indentation * tip_radius)) / h
                            o2 = (1.497 * indentation * tip_radius) / h ** 2
                            o3 = (1.469 * indentation * tip_radius * math.sqrt(indentation * tip_radius)) / h ** 3
                            o4 = (0.755 * indentation ** 2 * tip_radius ** 2) / h ** 4
                            coeff = 1.0 / (1 + o1 + o2 + o3 + o4)

                        else:
                            # Don't correct if we are "under" the glass
                            coeff = 1.0

                    elif model == 1 or model == 2:
                        # Sneddon (Cone and Pyramid)

                        if h >= 0:
                            # Old Chadwick/Dimitriadis correction
                            # X = (indentation) / h
                            # o1 = 1.7795 * (2 * tan_angle / math.pi ** 2) * X
                            # o2 = 16.0 * (1.7795) ** 2 * tan_angle ** 2 * X ** 2
                            # Updated correction according to Garcia 2018
                            o1 = (0.721 * indentation * tan_angle) / h
                            o2 = (0.650 * indentation ** 2 * tan_angle ** 2) / h ** 2
                            o3 = (0.491 * indentation ** 3 * tan_angle ** 3) / h ** 3
                            o4 = (0.225 * indentation ** 4 * tan_angle ** 4) / h ** 4
                            coeff = 1.0 / (1 + o1 + o2 + o3 + o4)
                        else:
                            # Don't correct if we are "under" the glass
                            coeff = 1.0

                    # print(f"Pos: ({i},{j}), Ind:{indentation}, tanangle{tan_angle}, h{h}")

                    # Correct the stiffness
                    E_new.append(E_original * coeff)

                # Update the progressbar
                widgets_list.widget_progressbar.update()

                resstiffness_array[j] = E_new

            # Stiffness_array is a VLarray with only one dimension
            for k in range(len(resstiffness_array)):
                stiffness_array.append(resstiffness_array[k])

        # Flush the file to be sure everything has been written on the disk
        self.temp_file.flush_file()

        # Flush the file and update the values in the experiment class
        data.stiffness_corrected = True
        data.update()

        # Close the progressbar and return
        widgets_list.widget_progressbar.close()

        self.logger.debug("Height correction finished")
