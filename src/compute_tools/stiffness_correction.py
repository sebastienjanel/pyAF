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
import logging
from .. import shared
from .. import widgets_list
from ..widgets.progressbar import Progressbar
from .stiffness import bec
from .stiffness import contact_models


class StiffnessCorrection:
    """Corrects the stiffness.

    The stiffness is multiplied by a coefficient depending on the distance
    between the topography and the substrate. One has to define a ROI in the
    roi manager, and set it as a glass substrate.

    Garcia & Garcia 2018 corrections (paraboloid, cone and pyramid, flat
    punch), applied inside the fit model (see stiffness.bec). The height
    of the sample is the user height, or the topography above the plane of
    the substrate ROI.
    """

    def __init__(self):
        self.temp_file = shared.exp.temp_file
        self.logger = logging.getLogger()

    def correct_stiffness(self, calc_id):
        """Method called by PYAF to do the actual stiffness correction."""
        data = shared.exp.list[calc_id]

        self.logger.debug("Correcting stiffness, file %s", calc_id)
        self.logger.debug("Using ROI %s", data.roi_glass_id)

        # Create a progressbar which is displayed during the computations
        Progressbar(os.path.basename(data.filename))
        widgets_list.widget_progressbar.set_label("Correcting elasticity")
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

        model = data.used_stiffness_model_selected
        tip_radius = data.used_tip_radius
        # Height of the sample [nm] (nan below the noise floor: no
        # correction)
        heights = bec.get_height_map(data)

        for i in range(data.nbr_pixels_x):
            for j in range(data.nbr_pixels_y):
                h = heights[i, j]
                values = data.stiffness_array[i][j]
                E_new = []
                if model not in contact_models.HAS_BEC:
                    # No bottom effect correction for this model (the
                    # computation is not available, see the compute tab)
                    E_new = [numpy.nan] * len(values)
                elif data.used_perform_fit and len(values) == 1:
                    # Fitting: the B.E.C. factor is inside the fit model
                    # (elasticity fitted again with the point of contact of
                    # the elasticity fit)
                    E_new = [values[0] *
                             bec.correction_factor(data, i, j, h)]
                else:
                    # Previous Calculus mode: factor at the maximal
                    # indentation (one value) or at the end of each segment
                    for z in range(len(values)):
                        indentation = (z + 1) * data.used_indentation_step
                        if len(values) == 1:
                            points = bec.approach_points(data, i, j)
                            indentation = 0 if points is None else \
                                max(float(numpy.max(points[0])) * 1e9, 0)
                        factor = contact_models.bec_factor(
                            indentation, model, tip_radius,
                            data.used_tip_angle, h)
                        E_new.append(values[z] / float(factor))

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
