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

from pyAF.src import shared
from pyAF.src import widgets_list
from pyAF.src.single import SingleFile

from base_test import pyAFTestCase
from base_test_files import HELA


class ComputeTabTest(pyAFTestCase):
    """Test the elements related to the compute tab."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([HELA])
        cls.base_setUpClass()

    def test_stiffness_parameters(self):
        """See if the stiffness parameters were not modified."""
        data = shared.exp.current_data
        default = SingleFile(None, None, None, None)
        self.assertEqual(data.fitparam_poc_skip_start,
                         default.fitparam_poc_skip_start)
        self.assertEqual(data.fitparam_poc_fit_length,
                         default.fitparam_poc_fit_length)
        self.assertEqual(data.fitparam_poc_refit_option,
                         default.fitparam_poc_refit_option)
        self.assertEqual(data.fitparam_poc_noise_multiplicator,
                         default.fitparam_poc_noise_multiplicator)
        self.assertEqual(data.fitparam_poc_refit_times,
                         default.fitparam_poc_refit_times)
        self.assertEqual(data.stiffness_model_selected,
                         default.stiffness_model_selected)
        self.assertEqual(data.tip_radius,
                         default.tip_radius)
        self.assertEqual(data.tip_angle,
                         default.tip_angle)
        self.assertEqual(data.poisson_ratio,
                         default.poisson_ratio)
        self.assertEqual(data.indentation_start,
                         default.indentation_start)
        self.assertEqual(data.indentation_stop,
                         default.indentation_stop)
        self.assertEqual(data.indentation_step,
                         default.indentation_step)
        self.assertEqual(data.strict_stop,
                         default.strict_stop)
        self.assertEqual(data.tomography,
                         default.tomography)

    def test_fit_pocs_inputs(self):
        """Test the Poc inputs in the compute widget.

        Only done for the first file, which is enough.
        """
        data = shared.exp.current_data
        W_pocs = widgets_list.widget_fit_param_pocs
        self.assertIsNotNone(W_pocs)
        self.assertEqual(W_pocs.IN_skip_start.get_int_value(),
                         data.fitparam_poc_skip_start)
        self.assertEqual(W_pocs.IN_fit_length.get_int_value(),
                         data.fitparam_poc_fit_length)
        self.assertEqual(W_pocs.IN_noise_multi.get_float_value(),
                         data.fitparam_poc_noise_multiplicator)
        self.assertEqual(W_pocs.IN_refit.get_int_value(),
                         data.fitparam_poc_refit_option)
        self.assertEqual(W_pocs.IN_refit_times.get_int_value(),
                         data.fitparam_poc_refit_times)

    def test_fit_jocs_inputs(self):
        """Test the Joc inputs in the compute widget.

        Only done for the first file, which is enough.
        """
        data = shared.exp.current_data
        W_jocs = widgets_list.widget_fit_param_jocs
        self.assertIsNotNone(W_jocs)
        self.assertEqual(W_jocs.IN_skip_start.get_int_value(),
                         data.fitparam_joc_skip_start)
        self.assertEqual(W_jocs.IN_fit_length.get_int_value(),
                         data.fitparam_joc_fit_length)
        self.assertEqual(W_jocs.IN_noise_multi.get_float_value(),
                         data.fitparam_joc_noise_multiplicator)
        self.assertEqual(W_jocs.IN_refit.get_int_value(),
                         data.fitparam_joc_refit_option)
        self.assertEqual(W_jocs.IN_refit_times.get_int_value(),
                         data.fitparam_joc_refit_times)

    def test_fit_joc_events_inputs(self):
        """Test the Joc events inputs in the compute widget.

        Only done for the first file, which is enough.
        """
        data = shared.exp.current_data
        w_jocs_events = widgets_list.widget_fit_param_jocs_events
        self.assertIsNotNone(w_jocs_events)
        self.assertEqual(w_jocs_events.IN_skip_start.get_int_value(),
                         data.fitparam_events_joc_skip_start)
        self.assertEqual(w_jocs_events.IN_fit_length.get_int_value(),
                         data.fitparam_events_joc_fit_length)
        self.assertEqual(w_jocs_events.IN_noise_multi.get_float_value(),
                         data.fitparam_events_joc_noise_multiplicator)
        self.assertEqual(w_jocs_events.IN_refit.get_int_value(),
                         data.fitparam_events_joc_refit_option)
        self.assertEqual(w_jocs_events.IN_refit_times.get_int_value(),
                         data.fitparam_events_joc_refit_times)
