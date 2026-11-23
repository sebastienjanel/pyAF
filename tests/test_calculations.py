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
correction (docs/models.md, sections 3 and 5), adhesion (section 9),
dynamic force spectroscopy (section 11), plant cells (section 13), the
order of the files of the loading window, the pixels sorted into the ROIs
by ROI assistance, the curves sorted from a few marked curves, and the
automatic tilt correction.

Usage (from the pyAF folder):
    .venv/bin/python -m unittest discover tests
"""

import math
import os
import sys
import unittest
from unittest import mock

import numpy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from src.compute_tools.stiffness import bec  # noqa: E402
from src.compute_tools.stiffness import contact_models  # noqa: E402
from src.compute_tools.work_and_rupture_force import \
    work_and_rupture_force_compute as work  # noqa: E402
from src.load_and_save.pre_check import sort_infos  # noqa: E402
from src.single import SingleFile, UNKNOWN_FLATTEN  # noqa: E402
from src.tools import auto_roi  # noqa: E402
from src.tools import curve_sorting  # noqa: E402
from src.tools import dfs_tools  # noqa: E402
from src.tools import math_tools  # noqa: E402
from src.tools import search_joc  # noqa: E402
from src.tools import tilt_tools  # noqa: E402
from src.tools.results_sorting import summary  # noqa: E402
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

    def test_sphere_limit(self):
        # Contact radius equal to the radius: about 14 R
        self.assertAlmostEqual(
            contact_models.sphere_max_indentation(1.0), 14.16, delta=0.01)

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


class Adhesion(unittest.TestCase):
    """Jumps of contact, work and rupture force of a retraction curve
    (deflection [nm] vs distance [nm])."""

    PARAMS = {"joc_skip_start": 20, "joc_fit_length": 200,
              "joc_refit_option": 0, "joc_noise_multiplicator": 5,
              "joc_refit_times": 0}
    K = 0.05e-9   # N/nm

    def curve(self, dip):
        """Baseline, adhesion of 20 nm from 600 to 800 nm (dip), contact
        after 800 nm."""
        rng = numpy.random.default_rng(0)
        x = numpy.linspace(0, 1000, 1001)
        y = rng.normal(0, 0.05, x.size)
        if dip:
            part = (x > 600) & (x <= 800)
            y[part] -= 20 * (x[part] - 600) / 200
        contact = x > 800
        y[contact] += x[contact] - 800
        return [x, y]

    def compute(self, curve):
        joc = search_joc.get_JOC(curve, [0, len(curve[0])], self.PARAMS)
        if not search_joc.joc_found(joc[1], curve):
            return None
        return work.get_surface_and_force(curve[0], curve[1], joc[0],
                                          joc[1], joc[2], joc[3], self.K)

    def test_rupture_force(self):
        found = self.compute(self.curve(True))
        self.assertIsNotNone(found)
        work_value, force = found
        self.assertGreater(work_value, 0)
        # 20 nm of deflection x 0.05 N/m = 1 nN
        self.assertAlmostEqual(force / 1e-9, 1.0, delta=0.05)

    def test_no_jump_found_is_not_a_work_of_0(self):
        self.assertIsNone(self.compute(self.curve(False)))


class ForceSpectroscopy(unittest.TestCase):
    """dfs_tools.get_points: the events with a loading rate <= 0 (not
    fitted) are not used, with both methods."""

    def items(self):
        forces = numpy.array([100.0, 110.0, 120.0, 130.0, 10.0, 10.0, 10.0])
        rates = numpy.array([1e4, 2e4, 3e4, 4e4, 0.0, 0.0, 0.0])
        return [("file", forces, rates, "k", [300.0])]

    def test_medians_per_file(self):
        with mock.patch.object(dfs_tools, "get_items",
                               return_value=self.items()):
            points_r, points_f, used, negative, *_ = dfs_tools.get_points(
                "single", "medians")
        self.assertEqual(float(points_r[0]), 2.5e4)
        self.assertEqual(float(points_f[0]), 115.0)
        self.assertEqual((used, negative), (4, 3))

    def test_median_point_of_the_plot(self):
        _, forces, rates, _, _ = self.items()[0]
        r, f, _, _ = dfs_tools.median_point(forces, rates)
        self.assertEqual((float(r), float(f)), (2.5e4, 115.0))
        self.assertIsNone(dfs_tools.median_point(forces[4:], rates[4:]))


class CantileverOfTheBEC(unittest.TestCase):
    """bec.cantilever_changed: the B.E.C. corrects the elasticity with its
    spring constant and deflection sensitivity."""

    def file(self, k, used_k, sensitivity=10.0, used_sensitivity=10.0):
        return mock.Mock(spring_constant=k, used_spring_constant=used_k,
                         deflection_sensitivity=sensitivity,
                         used_deflection_sensitivity=used_sensitivity)

    def test_same_values(self):
        self.assertFalse(bec.cantilever_changed(self.file(5e-11, 5e-11)))

    def test_changed_values(self):
        self.assertTrue(bec.cantilever_changed(self.file(6e-11, 5e-11)))
        self.assertTrue(bec.cantilever_changed(
            self.file(5e-11, 5e-11, 12.0, 10.0)))

    def test_not_computed(self):
        self.assertFalse(bec.cantilever_changed(
            self.file(5e-11, None, 10.0, None)))


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


class FileOrder(unittest.TestCase):
    """pre_check.sort_infos: order of the files of the loading window."""

    INFOS = [{"filename": "map10.jpk-force-map", "date": "2026-01-03 10:00:00"},
             {"filename": "project.pyaf", "date": ""},
             {"filename": "Map2.jpk-force-map", "date": "2026-01-01 10:00:00"},
             {"filename": "map1.jpk-force-map", "date": "2026-01-02 10:00:00"}]

    def names(self, infos):
        return [info["filename"].split(".")[0] for info in infos]

    def test_natural_order_of_the_names(self):
        self.assertEqual(self.names(sort_infos(self.INFOS, "name")),
                         ["map1", "Map2", "map10", "project"])
        self.assertEqual(self.names(sort_infos(self.INFOS, "name", True)),
                         ["project", "map10", "Map2", "map1"])

    def test_dates_and_files_without_date(self):
        self.assertEqual(self.names(sort_infos(self.INFOS, "date")),
                         ["Map2", "map1", "map10", "project"])
        # Without date: at the end in both orders
        self.assertEqual(self.names(sort_infos(self.INFOS, "date", True)),
                         ["map10", "map1", "Map2", "project"])


def domains_map(nx=64, ny=64, seed=0):
    """Supported lipid bilayer: round domains 1 nm higher and stiffer than
    the rest, noise on both maps; the topography, the elasticity [Pa] and
    the domains (boolean [x][y])."""
    rng = numpy.random.default_rng(seed)
    x, y = grids(nx, ny)
    domains = numpy.zeros((nx, ny), dtype=bool)
    for cx, cy, radius in ((16, 16, 9), (45, 20, 7), (30, 47, 11)):
        domains |= (x - cx) ** 2 + (y - cy) ** 2 < radius ** 2
    topo = numpy.where(domains, 1.0, 0.0) + rng.normal(0, 0.4, (nx, ny))
    stiffness = numpy.where(domains, 40e3, 20e3) * \
        rng.lognormal(0, 0.3, (nx, ny))
    return topo, stiffness, domains


class AutoRoi(unittest.TestCase):
    """Pixels sorted into the ROIs from a few marked pixels."""

    def test_classify_threshold_and_size(self):
        proba = numpy.zeros((6, 6, 2))
        proba[..., 0] = 0.9
        proba[..., 1] = 0.1
        # One pixel of the second ROI (too small), one uncertain
        proba[1, 1] = (0.2, 0.8)
        proba[4, 4] = (0.55, 0.45)
        proba[5, 5] = numpy.nan
        marks = numpy.full((6, 6), -1)
        marks[0, 0] = 1
        result = auto_roi.classify(proba, marks, 0.7, 2)
        self.assertEqual(result[1, 1], -1)
        self.assertEqual(result[4, 4], -1)
        self.assertEqual(result[5, 5], -2)
        # The marks stay in their ROI
        self.assertEqual(result[0, 0], 1)
        self.assertEqual(result[3, 3], 0)
        result = auto_roi.classify(proba, marks, 0.5, 1)
        self.assertEqual(result[1, 1], 1)
        self.assertEqual(result[4, 4], 0)

    def test_smoothed_ignores_nan(self):
        values = numpy.ones((20, 20))
        values[5:8, 5:8] = numpy.nan
        result = auto_roi.smoothed(values, 2.0)
        self.assertTrue(numpy.allclose(result, 1.0))

    def test_normalized(self):
        self.assertIsNone(auto_roi.normalized(numpy.full((4, 4), 3.0)))
        values = auto_roi.normalized([[1.0, 10.0], [100.0, numpy.nan]],
                                     log=True)
        self.assertAlmostEqual(values[0, 1], 0.0)
        self.assertTrue(numpy.isnan(values[1, 1]))

    @unittest.skipUnless(auto_roi.available(), "scikit-learn not installed")
    def test_domains_from_few_marks(self):
        topo, stiffness, domains = domains_map()
        # Curves without values (nan) in both maps
        topo[10:13, 40:43] = numpy.nan
        stiffness[10:13, 40:43] = numpy.nan
        feats = numpy.stack(
            auto_roi.channel_features(auto_roi.normalized(topo)) +
            auto_roi.channel_features(auto_roi.normalized(stiffness, True)),
            axis=-1)
        # A few strokes in a domain and outside
        marks = [{(i, 16) for i in range(12, 21)},
                 {(i, 60) for i in range(5, 60)} |
                 {(60, j) for j in range(5, 40)}]
        labels = auto_roi.marks_labels(domains.shape, marks)
        proba = auto_roi.probabilities(feats, labels)
        result = auto_roi.classify(proba, labels, 0.7, 5)
        # Away from the borders of the domains, nearly all right
        inner = ~(numpy.logical_xor(
            domains, auto_roi.ndimage.binary_erosion(domains, iterations=2)) |
            numpy.logical_xor(domains, auto_roi.ndimage.binary_dilation(
                domains, iterations=2)))
        sure = inner & (result >= 0)
        right = (result[sure] == 0) == domains[sure]
        self.assertGreater(right.mean(), 0.97)
        self.assertGreater(sure.sum(), 0.8 * inner.sum())
        # Region of the map: nan outside
        inside = numpy.ones(domains.shape, dtype=bool)
        inside[:, :4] = False
        proba = auto_roi.probabilities(feats, labels, inside)
        self.assertTrue(numpy.isnan(proba[:, :4]).all())
        self.assertEqual(auto_roi.classify(proba, labels, 0.7, 5)[0, 0], -2)

    def test_cells_of_a_tissue(self):
        """Plant tissue: 3 x 3 domed cells of different heights, grooves
        between them, noise; one stroke in a cell gives the whole cell."""
        rng = numpy.random.default_rng(1)
        x, y = grids(60, 60)
        centers = [(10 + 20 * a, 10 + 20 * b) for a in range(3)
                   for b in range(3)]
        nearest = numpy.argmin([(x - cx) ** 2 + (y - cy) ** 2
                                for cx, cy in centers], axis=0)
        height = numpy.zeros((60, 60))
        for index, (cx, cy) in enumerate(centers):
            dome = (3000 + 300 * index) * numpy.clip(
                1 - ((x - cx) ** 2 + (y - cy) ** 2) / 16.0 ** 2, 0, None)
            height = numpy.where(nearest == index, dome, height)
        height += rng.normal(0, 50, height.shape) + 2.0 * x
        height[5, 5] = numpy.nan
        dome_labels, number = auto_roi.domes(height, 0.05)
        self.assertEqual(number, 9)
        # A stroke in the central cell (index 4), and a click in a corner
        labels = auto_roi.marks_labels(
            height.shape, [{(i, 30) for i in range(27, 33)}, {(8, 50)}])
        result, shared = auto_roi.dome_classes(dome_labels, labels)
        self.assertEqual(shared, 0)
        central = nearest == 4
        self.assertGreater(numpy.mean(result[central] == 0), 0.9)
        # At most a pixel more along its edge (80 px): between domes of
        # different heights the groove is a step, within a pixel
        self.assertLess(numpy.count_nonzero(result[~central] == 0), 80)
        self.assertGreater(numpy.mean(result[nearest == 2] == 1), 0.9)
        # A cell marked in two ROIs: uncertain, the marks kept
        labels = auto_roi.marks_labels(height.shape, [{(30, 30)},
                                                      {(31, 31)}])
        result, shared = auto_roi.dome_classes(dome_labels, labels)
        self.assertEqual(shared, 1)
        self.assertEqual(result[30, 30], 0)
        self.assertEqual(result[28, 33], -1)

    def test_erased_substrate(self):
        """Cell on a flat substrate, a shallow groove between them: the
        substrate is no dome (too low above the groove), it is flooded from
        the cell; erased from the ROI with one stroke, it is flooded from
        the stroke up to the groove."""
        x, y = grids(50, 40)
        height = 300.0 + 2000.0 * numpy.clip(
            1 - ((x - 15) ** 2 + (y - 20) ** 2) / 12.0 ** 2, 0, None)
        height = numpy.where(numpy.abs(x - 30) < 1.5, 250.0, height)
        cell = (x - 15) ** 2 + (y - 20) ** 2 < 11.0 ** 2
        labels = auto_roi.marks_labels(height.shape,
                                       [{(i, 20) for i in range(12, 18)}])
        dome_labels, number = auto_roi.domes(height, 0.05)
        result = auto_roi.dome_classes(dome_labels, labels)[0]
        # The substrate is taken with the cell
        self.assertTrue((result[x > 32] == 0).all())
        stroke = {(45, j) for j in range(10, 30)}
        dome_labels, number = auto_roi.domes(height, 0.05, [stroke])
        result = auto_roi.dome_classes(dome_labels, labels,
                                banned={0: {number + 1}})[0]
        self.assertTrue((result[x > 32] == -1).all())
        self.assertTrue((result[cell] == 0).all())

    @unittest.skipUnless(auto_roi.available(), "scikit-learn not installed")
    def test_one_roi_marked(self):
        feats = numpy.zeros((8, 8, 3))
        labels = auto_roi.marks_labels((8, 8), [{(1, 1), (2, 2)}])
        with self.assertRaises(auto_roi.AutoRoiError):
            auto_roi.probabilities(feats, labels)


def synthetic_curves(kinds, seed=0, points=300):
    """Curves [curve][z, deflection][point] (approach, retraction) of each
    kind: "indentation" (contact, then indentation), "no contact" (flat),
    "spikes" (indentation with jumps), "unbinding" (indentation, then an
    adhesion peak on the retraction), "no unbinding" (no peak)."""
    rng = numpy.random.default_rng(seed)
    z = numpy.linspace(0, 1000, points)
    approach, retraction = [], []
    for kind in kinds:
        contact = rng.uniform(550, 750)
        push = numpy.clip(z - contact, 0, None)
        app = 0.002 * push ** 1.5 * rng.uniform(0.8, 1.2)
        if kind == "no contact":
            app = numpy.zeros(points)
        if kind == "spikes":
            for _ in range(3):
                start = rng.integers(20, points - 20)
                app[start:start + 3] += rng.uniform(5, 10)
        ret = app[::-1].copy()
        if kind == "unbinding":
            peak = rng.integers(points // 2, points - 30)
            ret[peak - 15:peak] -= numpy.linspace(0, rng.uniform(4, 8), 15)
        app = app + rng.normal(0, 0.2, points)
        ret = ret + rng.normal(0, 0.2, points)
        approach.append([z, app])
        retraction.append([z[::-1], ret])
    return numpy.array(approach), numpy.array(retraction)


class CurveSorting(unittest.TestCase):
    """Bad curves learned from a few marked curves, whatever they are."""

    def sort(self, kinds, bad_kinds, marked):
        """Probability to be bad of the curves of kinds, the curves of
        marked (indexes) marked as bad if of bad_kinds."""
        approach, retraction = synthetic_curves(kinds)
        shape, numbers = curve_sorting.measures(approach, retraction)
        # As a map of one row
        all_features = curve_sorting.features(
            [(shape[None], numbers[None])])
        marks = [(0, 0, index, kinds[index] in bad_kinds)
                 for index in marked]
        return curve_sorting.bad_probability(all_features, marks)[0][0]

    @unittest.skipUnless(curve_sorting.available(),
                         "scikit-learn not installed")
    def test_indentation(self):
        kinds = ["indentation"] * 150 + ["no contact"] * 30 + \
            ["spikes"] * 30
        rng = numpy.random.default_rng(1)
        kinds = list(rng.permutation(kinds))
        # 10 marks of each kind
        marked = [index for kind in ("indentation", "no contact", "spikes")
                  for index in [i for i, k in enumerate(kinds)
                                if k == kind][:10]]
        proba = self.sort(kinds, ("no contact", "spikes"), marked)
        bad = numpy.array([kind != "indentation" for kind in kinds])
        self.assertGreater(numpy.mean((proba >= 0.5) == bad), 0.95)

    @unittest.skipUnless(curve_sorting.available(),
                         "scikit-learn not installed")
    def test_unbinding(self):
        """SMFS: here the good curves have an unbinding peak."""
        kinds = list(numpy.random.default_rng(2).permutation(
            ["unbinding"] * 80 + ["no unbinding"] * 120))
        # 5 marks of each kind
        marked = [index for index, kind in enumerate(kinds)
                  if kind == "unbinding"][:5] + \
            [index for index, kind in enumerate(kinds)
             if kind == "no unbinding"][:5]
        proba = self.sort(kinds, ("no unbinding",), marked)
        bad = numpy.array([kind == "no unbinding" for kind in kinds])
        self.assertGreater(numpy.mean((proba >= 0.5) == bad), 0.95)

    def test_missing_curve(self):
        approach, retraction = synthetic_curves(["indentation"] * 2)
        approach[1] = 0.0
        retraction[1] = 0.0
        shape, numbers = curve_sorting.measures(approach, retraction)
        self.assertTrue(numpy.isfinite(shape[0]).all())
        self.assertTrue(numpy.isnan(shape[1]).all())
        self.assertTrue(numpy.isnan(numbers[1]).all())

    def test_curves_without_points(self):
        """Files with curves of one point, or without retraction: missing
        curves, no error."""
        approach, retraction = synthetic_curves(["indentation"] * 3)
        for app, ret in ((approach, retraction[:, :, :1]),
                         (approach, retraction[:, :, :0]),
                         (approach[:, :, :1], retraction)):
            shape, numbers = curve_sorting.measures(app, ret)
            self.assertEqual(shape.shape, (3, 2 * curve_sorting.POINTS))
            self.assertEqual(numbers.shape, (3, 6))
            self.assertTrue(numpy.isnan(shape).all())

    def test_resampled(self):
        values = numpy.arange(64.0).reshape(1, 64)
        self.assertTrue(numpy.allclose(
            curve_sorting.resampled(values, 32)[0],
            numpy.arange(0.5, 64, 2)))
        # Fewer points than asked: interpolated
        self.assertEqual(curve_sorting.resampled(values[:, :10], 32).shape,
                         (1, 32))

    @unittest.skipUnless(curve_sorting.available(),
                         "scikit-learn not installed")
    def test_too_few_marks(self):
        approach, retraction = synthetic_curves(["indentation"] * 20)
        shape, numbers = curve_sorting.measures(approach, retraction)
        all_features = curve_sorting.features([(shape[None],
                                                numbers[None])])
        with self.assertRaises(curve_sorting.CurveSortingError):
            curve_sorting.bad_probability(all_features,
                                          [(0, 0, 0, True), (0, 0, 1, False)])


def tilted_curve(contact, rng, slope=0.004, offset=3.0, noise=0.08,
                 points=1000, events=False):
    """Approach (far from the surface first) of 2 um with a tilt, contact at
    the height contact [nm] (Hertz-like indentation after it); with
    events: two unbinding-like dips before the contact (retraction)."""
    z = numpy.linspace(0, 2000.0, points)
    d = slope * z + offset + rng.normal(0, noise, points)
    d += 2e-4 * numpy.clip(z - contact, 0, None) ** 1.5
    if events:
        for position in (0.6 * contact, 0.85 * contact):
            dip = (z > position - 60) & (z < position)
            d[dip] -= numpy.linspace(0, 1.5, dip.sum())
    return z, d


class AutoTilt(unittest.TestCase):
    """Baseline of each curve found without limits, slope of the file from
    all the baselines (tools/tilt_tools)."""

    def file_slope(self, contacts, **options):
        rng = numpy.random.default_rng(2)
        baselines, ends = [], []
        for contact in contacts:
            z, d = tilted_curve(contact, rng, **options)
            slope, _, first, end, _ = tilt_tools.baseline(z, d)
            baselines.append((slope, end - first, z[end - 1] - z[first]))
            ends.append(z[end - 1] / contact)
        return tilt_tools.file_slope(baselines), numpy.array(ends)

    def test_varying_baselines(self):
        # Contact from 3 % to 90 % of the curve (cells of different
        # heights): the baseline ends at the contact, the slope is found
        contacts = numpy.random.default_rng(3).uniform(60, 1800, 150)
        (slope, _, _), ends = self.file_slope(contacts)
        self.assertAlmostEqual(slope, 0.004, delta=0.0001)
        self.assertGreater(numpy.median(ends), 0.9)
        self.assertLess(numpy.median(ends), 1.1)

    def test_events(self):
        # Retraction with unbinding events: the baseline stops before them
        contacts = numpy.random.default_rng(4).uniform(800, 1800, 50)
        (slope, _, _), ends = self.file_slope(contacts, events=True)
        self.assertAlmostEqual(slope, 0.004, delta=0.0001)
        self.assertLess(ends.max(), 0.6)

    def test_turnaround(self):
        # Noise of the turnaround of the piezo at the start of the curves
        # (old instruments), with a small tilt (0.2 nm over the curve, as in
        # hela.003.spm): the baseline is not cut by it
        rng = numpy.random.default_rng(6)
        baselines = []
        for contact in rng.uniform(1200, 1700, 30):
            z, d = tilted_curve(contact, rng, slope=0.0001)
            d[:60] += 0.6 * numpy.exp(-numpy.arange(60) / 15.0)
            slope, _, first, end, _ = tilt_tools.baseline(z, d)
            baselines.append((slope, end - first, z[end - 1] - z[first]))
            self.assertGreater(z[end - 1], 0.8 * contact)
        slope, short, _ = tilt_tools.file_slope(baselines)
        self.assertAlmostEqual(slope, 0.0001, delta=0.00005)
        self.assertEqual(short, [])

    def test_turnaround_strong_tilt(self):
        # With a strong tilt, a baseline found after the turnaround is not
        # flat: the first one is kept (not a part of the contact)
        rng = numpy.random.default_rng(7)
        z, d = tilted_curve(300.0, rng)
        _, _, first, end, _ = tilt_tools.baseline(z, d)
        self.assertEqual(first, 0)
        self.assertLess(z[end - 1], 1.2 * 300.0)

    def test_correction(self):
        # Slope of the file, offset of each curve: the baseline at 0
        rng = numpy.random.default_rng(5)
        z, d = tilted_curve(700.0, rng)
        retraction = [z, d - 1.0]
        approach, corrected_retraction, info = math_tools.correct_tilt(
            [z, d], retraction, 0, 0, "trace", [0], [0], slope=0.004)
        self.assertTrue(info["is_corrected"])
        base = z < 650
        self.assertLess(abs(numpy.median(approach[1][base])), 0.02)
        self.assertLess(numpy.abs(approach[1][base]).max(), 0.5)
        # Same line removed from the retraction
        self.assertAlmostEqual(
            numpy.median(corrected_retraction[1][base]), -1.0, delta=0.02)
        # Fitted on the retraction
        corrected_approach, retraction, _ = math_tools.correct_tilt(
            [z, d + 1.0], [z, d], 0, 0, "retrace", [0], [0], slope=0.004)
        self.assertLess(abs(numpy.median(retraction[1][base])), 0.02)
        self.assertAlmostEqual(
            numpy.median(corrected_approach[1][base]), 1.0, delta=0.02)

    def test_outliers(self):
        # A few curves with another slope are reported, not followed
        baselines = [(0.004 + 0.0001 * (k % 3), 500, 1000.0)
                     for k in range(40)]
        baselines += [(0.05, 500, 1000.0), (0.004, 10, 20.0)]
        slope, short, far = tilt_tools.file_slope(baselines)
        self.assertAlmostEqual(slope, 0.004, delta=0.0003)
        self.assertEqual(far, [40])
        self.assertEqual(short, [41])

    def test_no_baseline(self):
        self.assertIsNone(tilt_tools.baseline([0, 1, 2], [0, 1, 2]))
        self.assertEqual(tilt_tools.file_slope([]), (None, [], []))


class Statistics(unittest.TestCase):
    """Mean, median and standard deviation of the tables of the plots."""

    def test_curves_without_value(self):
        # Curves without value (nan, e.g. a fit which failed) are not
        # counted, as in the violins
        values = numpy.r_[numpy.arange(1.0, 6.0), [numpy.nan] * 3]
        mean, median, std = summary(values)
        self.assertEqual((mean, median), (3.0, 3.0))
        self.assertAlmostEqual(std, numpy.std(numpy.arange(1.0, 6.0)))
        self.assertTrue(all(numpy.isnan(summary([numpy.nan]))))


if __name__ == "__main__":
    unittest.main()
