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

from base_test import pyAFTestCase
from base_test_files import HELA, NS_SINGLE, MAP16, JPK_SINGLE


class FileImportTest(pyAFTestCase):
    """Test the elements related to file import."""

    @classmethod
    def setUpClass(cls):
        """Load files and set-up the system."""
        cls.use_files([HELA, NS_SINGLE, MAP16, JPK_SINGLE])
        cls.base_setUpClass()

    def test_approach_velocity_Nanoscope(self):
        """Check the NS approach velocity extraction."""
        val = shared.exp.list[0].approach_velocity
        self.assertEqual(val, 20559.2727368)

    def test_retraction_velocity_Nanoscope(self):
        """Check the NS retraction velocity extraction."""
        val = shared.exp.list[0].retraction_velocity
        self.assertEqual(val, 20559.2727368)

    def test_scan_rate_Nanoscope(self):
        """Check the NS scan rate extraction."""
        self.assertEqual(shared.exp.list[0].scan_rate, 2.05592)

    def test_approach_velocity_JPK_force_map(self):
        """Check the JPK Map approach velocity extraction."""
        val = shared.exp.list[2].approach_velocity
        self.assertEqual(val, 16000)

    def test_retraction_velocity_JPK_force_map(self):
        """Check the JPK Map retraction velocity extraction."""
        val = shared.exp.list[2].retraction_velocity
        self.assertEqual(val, 16000)

    def test_scan_rate_JPK_force_map(self):
        """Check the JPK Map scan rate extraction."""
        self.assertEqual(shared.exp.list[2].scan_rate, 2)

    def test_approach_velocity_JPK_single_curve(self):
        """Check the JPK Single File approach velocity extraction."""
        val = shared.exp.list[3].approach_velocity
        self.assertEqual(val, 5000.0)

    def test_retraction_velocity_JPK_single_curve(self):
        """Check the JPK Single File approach velocity extraction."""
        val = shared.exp.list[3].retraction_velocity
        self.assertEqual(val, 5000.0)

    def test_scan_rate_JPK_single_curve(self):
        """Check the JPK Single File approach velocity extraction."""
        self.assertEqual(shared.exp.list[3].scan_rate, 1.25)

    def test_NS_velocity_consistency(self):
        """Check if the imported velocities are correct.

        Just compare the two methods used to compute the velocities :
        2*ramp_size*scan_rate or velocity*sens_z_scan
        """
        for i in range(2):
            data = shared.exp.list[i]

            velocity1 = data.approach_velocity
            velocity2 = 2 * data.ramp_size * data.scan_rate
            self.assertAlmostEqual(velocity1, velocity2, places=0)

            velocity1 = data.retraction_velocity
            velocity2 = 2 * data.ramp_size * data.scan_rate
            self.assertAlmostEqual(velocity1, velocity2, places=0)
