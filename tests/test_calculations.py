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


"""Checks of the calculations of pyAF on synthetic data with known results
(no file, no window): flatten and levelling of the z motor steps,
zero-force topography after a flatten, contact models and bottom effect
correction (docs/models.md, sections 3 and 5), plant cells (section 13).

Usage (from the pyAF folder):
    .venv/bin/python -m unittest discover tests
"""

import math
import os
import sys
import unittest

import numpy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from src.compute_tools.stiffness import contact_models  # noqa: E402
from src.single import SingleFile, UNKNOWN_FLATTEN  # noqa: E402
from src.tools import math_tools  # noqa: E402
from src.tools import plant_tools  # noqa: E402

NX, NY = 40, 30


def grids(nx=NX, ny=NY):
    """Pixel coordinates x, y of a map [x][y]."""
    return numpy.meshgrid(numpy.arange(nx, dtype=float),
                          numpy.arange(ny, dtype=float), indexing="ij")


def surface(nx=NX, ny=NY):
    """Smooth sample: a bump on a gentle slope [nm]."""
    x, y = grids(nx, ny)
    return 800.0 * numpy.exp(-((x - nx / 2) ** 2 + (y - ny / 2) ** 2) /
                             (2 * 8.0 ** 2)) + 3.0 * x + 2.0 * y


def scan_parts(nx, ny, lengths):
    """Parts of a map between the moves of the z motor: runs of curves of
    the given lengths, in the order of a back-and-forth scan (lines along
    x)."""
    parts = numpy.zeros((nx, ny), int)
    order = []
    for j in range(ny):
        line = range(nx) if j % 2 == 0 else range(nx - 1, -1, -1)
        order += [(i, j) for i in line]
    position = 0
    for number, length in enumerate(lengths):
        for pixel in order[position:position + length]:
            parts[pixel] = number
        position += length
    for pixel in order[position:]:
        parts[pixel] = len(lengths)
    return parts


class Flatten(unittest.TestCase):
    """math_tools.flatten: along the scan lines (x, for each y)."""

    def test_order_0_removes_the_offset_of_each_line(self):
        data = surface() + numpy.linspace(0, 5000, NY)[None, :]
        flat = math_tools.flatten(data, 0)
        medians = numpy.median(flat, axis=0)
        self.assertLess(numpy.ptp(medians), 1e-9)

    def test_lines_are_along_x(self):
        # Offsets between the lines (y) are removed, not between the
        # columns (x)
        offsets = numpy.random.default_rng(1).normal(0, 1000, NY)
        flat = math_tools.flatten(numpy.zeros((NX, NY)) + offsets, 0)
        self.assertLess(numpy.ptp(flat), 1e-9)

    def test_order_1_removes_the_tilt_of_each_line(self):
        x, y = grids()
        data = (2.0 + 0.1 * y) * x + 50.0 * y + 7.0
        self.assertLess(numpy.ptp(math_tools.flatten(data, 1)), 1e-6)

    def test_order_2_removes_the_curvature_of_each_line(self):
        x, y = grids()
        data = 0.3 * x ** 2 - 4.0 * x + 10.0 * y
        self.assertLess(numpy.ptp(math_tools.flatten(data, 2)), 1e-6)

    def test_pixels_outside_the_region(self):
        x, y = grids()
        data = 2.0 * x + 10.0 * y
        inside = numpy.ones((NX, NY), bool)
        inside[:5, :] = False
        data[~inside] = 1e9
        flat = math_tools.flatten(data, 1, inside)
        # The fit ignores them, and they are set to the minimum (0)
        self.assertLess(numpy.ptp(flat[inside]), 1e-6)
        self.assertTrue(numpy.all(flat[~inside] == 0))
        self.assertEqual(flat.min(), 0)


class MotorSteps(unittest.TestCase):
    """math_tools.level_steps: offsets of the parts between the moves of
    the z motor, also in the middle of the lines."""

    def setUp(self):
        rng = numpy.random.default_rng(2)
        self.surface = surface()
        lengths = rng.integers(15, 90, 25)
        self.parts = scan_parts(NX, NY, lengths)
        self.offsets = rng.normal(0, 3000, self.parts.max() + 1)
        self.data = self.surface + self.offsets[self.parts]

    def test_mid_line_steps(self):
        # At least one part starts in the middle of a line
        starts = numpy.diff(self.parts, axis=0) != 0
        self.assertTrue(starts.any())

    def test_offsets_removed(self):
        levelled = math_tools.level_steps(self.data, self.parts)
        error = levelled - self.surface
        self.assertLess(numpy.ptp(error), 5.0)

    def test_reference_is_the_largest_part(self):
        levelled = math_tools.level_steps(self.data, self.parts)
        largest = numpy.argmax(numpy.bincount(self.parts.ravel()))
        mask = self.parts == largest
        numpy.testing.assert_allclose(levelled[mask], self.data[mask])

    def test_pixels_without_curve_are_kept(self):
        parts = self.parts.copy()
        parts[:3, :3] = -1
        levelled = math_tools.level_steps(self.data, parts)
        numpy.testing.assert_array_equal(levelled[:3, :3],
                                         self.data[:3, :3])

    def test_flatten_after_levelling(self):
        flat = math_tools.flatten(self.data, 0, None, self.parts)
        expected = math_tools.flatten(self.surface, 0)
        self.assertLess(numpy.abs(flat - expected).max(), 5.0)


class FakeFile:
    """The attributes of a file (SingleFile) used by the flatten of the
    piezo height and by the zero-force topography."""

    flatten_state = SingleFile.flatten_state
    piezo_for = SingleFile.piezo_for
    topography_state = SingleFile.topography_state
    sync_topography = SingleFile.sync_topography

    def __init__(self, piezo, parts=None):
        self.piezo_image = piezo
        self.motor_parts = parts
        self.flatten_applied = False
        self.applied_flatten_order = None
        self.applied_flatten_motor = False
        self.stiffness_calculated = False
        self.topography = None
        self.topography_flatten = UNKNOWN_FLATTEN
        self.max_topo = None
        self.outside_region = None
        self.scan_size_x = self.scan_size_y = 0.0
        self.nbr_pixels_x, self.nbr_pixels_y = piezo.shape

    @staticmethod
    def inside_region():
        return None

    def compute_topography(self, contact):
        """Topography as computed by the elasticity: piezo height (with the
        flatten applied) plus the height of the point of contact of each
        curve, minimum at 0."""
        topo = self.piezo_for(self.flatten_state()) + contact
        self.topography = topo - topo.min()
        self.topography_flatten = self.flatten_state()
        self.stiffness_calculated = True

    def set_flatten(self, order, motor=False):
        self.flatten_applied = order is not None
        self.applied_flatten_order = order
        self.applied_flatten_motor = motor


class TopographyFollowsTheFlatten(unittest.TestCase):
    """SingleFile.sync_topography: the same topography as a new computation
    with the flatten applied."""

    def setUp(self):
        rng = numpy.random.default_rng(3)
        parts = scan_parts(NX, NY, rng.integers(20, 80, 15))
        piezo = surface() + rng.normal(0, 2000, parts.max() + 1)[parts]
        self.contact = rng.normal(0, 50, (NX, NY))
        self.file = FakeFile(piezo, parts)
        self.recomputed = FakeFile(piezo, parts)

    def check(self, order, motor):
        self.file.set_flatten(order, motor)
        self.assertTrue(self.file.sync_topography())
        self.recomputed.set_flatten(order, motor)
        self.recomputed.compute_topography(self.contact)
        numpy.testing.assert_allclose(self.file.topography,
                                      self.recomputed.topography,
                                      atol=1e-6)
        self.assertEqual(self.file.topography_state(),
                         self.file.flatten_state())

    def test_flatten_after_the_elasticity(self):
        self.file.compute_topography(self.contact)
        for order, motor in ((0, False), (1, True), (2, True), (None, False)):
            with self.subTest(order=order, motor=motor):
                self.check(order, motor)

    def test_unchanged_flatten(self):
        self.file.set_flatten(1)
        self.file.compute_topography(self.contact)
        before = self.file.topography.copy()
        self.assertFalse(self.file.sync_topography())
        numpy.testing.assert_array_equal(self.file.topography, before)

    def test_unknown_flatten_is_the_current_one(self):
        # Topography of an older save: the flatten in place is assumed
        self.file.set_flatten(0)
        self.file.compute_topography(self.contact)
        self.file.topography_flatten = UNKNOWN_FLATTEN
        self.assertFalse(self.file.sync_topography())
        self.assertEqual(self.file.topography_flatten, (0, False))


class ContactModels(unittest.TestCase):
    """contact_models.unit_force and bec_factor (docs/models.md, sections 3
    and 5)."""

    R = 50e-9
    NU = 0.5

    def test_paraboloid(self):
        delta = numpy.array([10e-9, 40e-9])
        expected = 4.0 / 3.0 * math.sqrt(self.R) * delta ** 1.5 / \
            (1 - self.NU ** 2)
        numpy.testing.assert_allclose(contact_models.unit_force(
            contact_models.PARABOLOID, self.R, 0, self.NU, delta), expected)

    def test_four_sided_pyramid(self):
        delta = numpy.array([100e-9])
        expected = math.tan(math.radians(35)) / math.sqrt(2) * delta ** 2 / \
            (1 - self.NU ** 2)
        numpy.testing.assert_allclose(contact_models.unit_force(
            contact_models.PYRAMID, self.R, 35, self.NU, delta), expected)

    def test_no_force_before_the_contact(self):
        force = contact_models.unit_force(contact_models.CONE, self.R, 20,
                                          self.NU, numpy.array([-1e-9, 0.0]))
        numpy.testing.assert_array_equal(force, 0)

    def test_sphere_tends_to_hertz(self):
        delta = numpy.array([0.5e-9, 1e-9])
        sphere = contact_models.unit_force(contact_models.SPHERE, self.R, 0,
                                           self.NU, delta)
        hertz = contact_models.unit_force(contact_models.PARABOLOID, self.R,
                                          0, self.NU, delta)
        numpy.testing.assert_allclose(sphere, hertz, rtol=0.02)

    def test_needle_and_nanowire_are_continuous(self):
        tan = math.tan(math.radians(20))
        for model, critical in ((contact_models.NEEDLE,
                                 math.pi * self.R / (2 * tan)),
                                (contact_models.NANOWIRE, self.R)):
            with self.subTest(model=model):
                force = contact_models.unit_force(
                    model, self.R, 20, self.NU,
                    critical * numpy.array([1 - 1e-6, 1 + 1e-6]))
                self.assertAlmostEqual(force[1] / force[0], 1.0, places=4)

    def test_bottom_effect_polynomial(self):
        # Paraboloid, chi = sqrt(R delta) / h = 0.5
        height = 1e-6
        delta = (0.5 * height) ** 2 / self.R
        x = 0.5
        expected = 1 + 1.133 * x + 1.497 * x ** 2 + 1.469 * x ** 3 + \
            0.755 * x ** 4
        factor = contact_models.bec_factor(numpy.array([delta]),
                                           contact_models.PARABOLOID,
                                           self.R, 0, height)
        self.assertAlmostEqual(float(factor[0]), expected, places=9)

    def test_bottom_effect_without_height(self):
        for height in (None, 0.0, float("nan")):
            factor = contact_models.bec_factor(
                numpy.array([1e-8]), contact_models.CONE, self.R, 20, height)
            self.assertEqual(float(factor[0]), 1.0)

    def test_nanowire_correction_of_its_cap(self):
        # Below the critical indentation: the correction of the paraboloid
        delta = numpy.array([0.3 * self.R])
        height = 2e-6
        self.assertAlmostEqual(
            float(contact_models.bec_factor(
                delta, contact_models.NANOWIRE, self.R, 0, height)[0]),
            float(contact_models.bec_factor(
                delta, contact_models.PARABOLOID, self.R, 0, height)[0]),
            places=12)


class PlantCells(unittest.TestCase):
    """plant_tools: apparent stiffness, geometry and shell models
    (docs/models.md, section 13)."""

    GEOMETRY = {"L_major": 30e-6, "L_minor": 12e-6, "k_major": 1 / 60e-6,
                "k_minor": 1 / 15e-6}
    THICKNESS = 0.5e-6
    NU = 0.5

    def geometry(self):
        geometry = dict(self.GEOMETRY)
        geometry["kM"] = (geometry["k_major"] + geometry["k_minor"]) / 2
        kappa_g = geometry["k_major"] * geometry["k_minor"]
        geometry["f"] = plant_tools.structural_factor(geometry["kM"],
                                                      kappa_g)
        return geometry

    def test_structural_factor(self):
        # 1 for a sphere, 3/4 for a cylinder
        self.assertAlmostEqual(plant_tools.structural_factor(2.0, 4.0), 1.0)
        self.assertAlmostEqual(plant_tools.structural_factor(2.0, 0.0), 0.75)

    def test_wall_strain(self):
        product = 0.4
        self.assertAlmostEqual(plant_tools.wall_strain(1.0, product),
                               2 * (math.asin(product) / product - 1))
        # Small L kappa: (L kappa)^2 / 3
        self.assertAlmostEqual(plant_tools.wall_strain(1.0, 1e-3) / 1e-6,
                               1 / 3, places=4)
        for product in (0.0, -0.1, 1.2):
            self.assertTrue(math.isnan(plant_tools.wall_strain(1.0,
                                                               product)))

    def test_yamasaki(self):
        geometry = self.geometry()
        k_as = 2.0
        modulus, pressure, tau = plant_tools.pressure_yamasaki(
            k_as, geometry, self.THICKNESS, self.NU)
        strain = plant_tools.wall_strain(geometry["L_major"],
                                         geometry["k_major"]) + \
            plant_tools.wall_strain(geometry["L_minor"], geometry["k_minor"])
        c = math.sqrt(3 * (1 - self.NU ** 2))
        argument = strain * c / (self.THICKNESS * geometry["kM"])
        self.assertAlmostEqual(tau, argument / 2)
        self.assertAlmostEqual(
            modulus, k_as * math.log(argument) /
            (math.pi * self.THICKNESS * strain))
        # P / E = t sum(epsilon) kM / f (eqs. 10 and 11)
        self.assertAlmostEqual(
            pressure / modulus,
            self.THICKNESS * strain * geometry["kM"] / geometry["f"])

    def test_yamasaki_without_solution(self):
        geometry = self.geometry()
        geometry["kM"] = -1e4
        self.assertTrue(all(math.isnan(value) for value in
                            plant_tools.pressure_yamasaki(
                                2.0, geometry, self.THICKNESS, self.NU)))

    def test_beauzamy_with_the_modulus_of_yamasaki(self):
        # Eq. 1 of Beauzamy with the modulus of Yamasaki gives the pressure
        # of Yamasaki (same shell equations, Lambert W on its branch -1)
        geometry = self.geometry()
        modulus, pressure, tau = plant_tools.pressure_yamasaki(
            2.0, geometry, self.THICKNESS, self.NU)
        pressure_b, tau_b = plant_tools.pressure_beauzamy(
            2.0, modulus, geometry, self.THICKNESS, self.NU)
        self.assertAlmostEqual(pressure_b / pressure, 1.0, places=9)
        self.assertAlmostEqual(tau_b / tau, 1.0, places=9)

    def test_stiffness_from_a_force_range(self):
        delta = numpy.linspace(0, 1e-6, 500)
        force = 0.8 * delta + 1e-9
        result = plant_tools.stiffness_force_range(delta, force)
        self.assertAlmostEqual(result["k_as"], 0.8, places=9)

    def test_stiffness_from_the_derivative(self):
        # Force whose derivative is k_as (1 - exp(-alpha (delta - beta)))
        k_as, alpha, beta = 0.5, 5e6, 0.2e-6
        delta = numpy.linspace(0, 2e-6, 2000)
        d = numpy.maximum(delta - beta, 0)
        force = k_as * (d - (1 - numpy.exp(-alpha * d)) / alpha)
        result = plant_tools.stiffness_derivative(delta, force)
        self.assertAlmostEqual(result["k_as"] / k_as, 1.0, delta=0.03)

    def test_geometry_of_a_spherical_cap(self):
        # Cap of a sphere of radius 50 um, ROI of radius 15 um, pixels of
        # 1 um: curvatures 1 / R, half-lengths 15 um, f = 1
        radius, roi_radius, n = 50.0, 15.0, 60
        x, y = grids(n, n)
        r2 = (x + 0.5 - n / 2) ** 2 + (y + 0.5 - n / 2) ** 2
        heights = numpy.sqrt(numpy.maximum(radius ** 2 - r2, 0)) * 1000.0
        cell = FakeFile(heights)
        cell.scan_size_x = cell.scan_size_y = n * 1000.0
        pixels = numpy.argwhere(r2 <= roi_radius ** 2)
        geometry = plant_tools.cell_geometry(cell, pixels.tolist())
        centre = geometry["points"][0]
        for name in ("k_major", "k_minor"):
            self.assertAlmostEqual(centre[name] * radius, 1.0, delta=0.05)
        self.assertAlmostEqual(centre["f"], 1.0, delta=0.01)
        for name in ("L_major", "L_minor"):
            self.assertAlmostEqual(geometry[name], roi_radius, delta=0.75)


if __name__ == "__main__":
    unittest.main()
