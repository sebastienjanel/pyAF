"""Contact models of the elasticity fit: force of an elastic half-space
indented by a tip, F = E * g(delta).

Models (index of the model in pyAF):

0. Paraboloid (Hertz): F = 4/3 E/(1 - nu^2) sqrt(R) delta^1.5, valid for
   delta << R.
1. Cone (Sneddon): F = 2/pi E/(1 - nu^2) tan(theta) delta^2.
2. 4-sided pyramid (Barber & Billings): n-sided pyramid
   F = 2 tan(theta) / (n sin(pi/n)) E/(1 - nu^2) delta^2 with n = 4,
   i.e. F = 1/sqrt(2) E/(1 - nu^2) tan(theta) delta^2 (Barber & Billings
   1990, Rico et al. 2005; within 6 % of the numerical solution of
   Bilodeau 1992, see Lekka et al. 2022, table 2.2.1).
3. Slope (Tomography, no contact model).
4. Flat punch (Sneddon): F = 2 R E/(1 - nu^2) delta.
5. Sphere (Sneddon): exact solution for a sphere of radius R (Sneddon,
   Int. J. Eng. Sci. 1965), parametrized by the contact radius a:
   delta = a/2 ln((R + a)/(R - a)),
   F = E/(2 (1 - nu^2)) [(a^2 + R^2) ln((R + a)/(R - a)) - 2 a R].
6. Blunted cone (Briscoe): cone of half-angle theta with a paraboloid
   apex of radius R tangent to the faces at the radius b = R / tan(theta)
   (Briscoe et al., J. Phys. D 1994; b as recommended in Lekka et al.
   2022, table 2.2.1, for a smooth transition). Model 7 with c = 1.
7. Blunted 4-sided pyramid (Rico): n-sided pyramid (n = 4) with a
   paraboloid apex (Rico et al., Phys. Rev. E 2005, as written in Lekka
   et al. 2022, table 2.2.1), with c = n / pi sin(pi / n) and b = c R / tan(theta)
   (smooth transition; with the b = R / tan(theta) of the book, the
   indentation decreases after the transition for c < 1):
   a <= b (apex): F = 4/3 E/(1 - nu^2) sqrt(R) delta^1.5,
   a > b, the contact radius a given by
   delta - c a / tan(theta) (pi/2 - arcsin(b/a))
       + a / R (sqrt(a^2 - b^2) - a) = 0,
   F = 2 E/(1 - nu^2) [delta a
       - c a^2 / (2 tan(theta)) (pi/2 - arcsin(b/a)) - a^3 / (3 R)
       + sqrt(a^2 - b^2) (c b / (2 tan(theta)) + (a^2 - b^2) / (3 R))].
   Sharp limit (b = 0) equal to model 2.
8. Needle (Garcia & Garcia): flat punch of radius R capped with a cone
   of half-angle theta (Garcia & Garcia, Biophys. J. 2018, eq. 18; the
   names needle and nanowire are swapped in Lekka et al. 2022, table
   2.2.1): cone (model 1) up to delta_c where the contact radius
   2 tan(theta) delta / pi reaches R, then the flat punch (model 4):
   F = 2/pi E/(1 - nu^2) tan(theta) delta^2 for delta <= delta_c,
   F = 2/pi E/(1 - nu^2) tan(theta) delta_c^2
       + 2 R E/(1 - nu^2) (delta - delta_c) beyond,
   delta_c = pi R / (2 tan(theta)).
9. Nanowire (Garcia & Garcia): flat punch of radius R capped with a
   hemisphere of radius R (Garcia & Garcia 2018, eq. 17): paraboloid
   (model 0) up to delta_c where the contact radius sqrt(R delta) reaches
   R, then the flat punch:
   F = 4/3 E/(1 - nu^2) sqrt(R) delta^1.5 for delta <= delta_c = R,
   F = 4/3 E/(1 - nu^2) R^2 + 2 R E/(1 - nu^2) (delta - R) beyond.
   In both models dF/d delta = 2 a E/(1 - nu^2) (a: contact radius) is
   continuous at delta_c.
10. Thin shell (Reissner): point load on a thin spherical shell of
   radius R (the vesicle, not the tip) and thickness h, linear for an
   indentation small compared to R (Reissner 1946; Calò et al.,
   Nanoscale 2014; Vorselen et al., Front. Mol. Biosci. 2020):
   F = 4 E h^2 / (R sqrt(3 (1 - nu^2))) delta.

The models 5 to 7 are tabulated (force per unit modulus vs indentation)
and interpolated. Units: SI (m, N, Pa); angles in degrees (half-angle
between the axis and the face of the tip).
"""

import math

import numpy

PARABOLOID = 0
CONE = 1
PYRAMID = 2
SLOPE = 3
FLAT_PUNCH = 4
SPHERE = 5
BLUNTED_CONE = 6
BLUNTED_PYRAMID = 7
NEEDLE = 8
NANOWIRE = 9
THIN_SHELL = 10

# Number of sides of the pyramids (models 2 and 7)
PYRAMID_SIDES = 4

MODEL_NAMES = ["Paraboloid (Hertz)", "Cone (Sneddon)",
               "4-sided pyramid (Barber & Billings)", "Slope",
               "Flat punch (Sneddon)", "Sphere (Sneddon)",
               "Blunted cone (Briscoe)",
               "Blunted 4-sided pyramid (Rico)",
               "Needle (Garcia & Garcia)", "Nanowire (Garcia & Garcia)",
               "Thin shell (Reissner)"]

# Parameters of the tips: radius and angle
USES_RADIUS = (PARABOLOID, FLAT_PUNCH, SPHERE, BLUNTED_CONE,
               BLUNTED_PYRAMID, NEEDLE, NANOWIRE, THIN_SHELL)
USES_ANGLE = (CONE, PYRAMID, BLUNTED_CONE, BLUNTED_PYRAMID, NEEDLE)

# Meaning of the radius for each model (label of the input)
RADIUS_LABELS = {FLAT_PUNCH: "Punch radius [nm]",
                 NEEDLE: "Needle radius [nm]",
                 NANOWIRE: "Wire radius [nm]",
                 THIN_SHELL: "Vesicle radius [nm]"}
RADIUS_LABEL = "Tip radius [nm]"

# Default radius of the tip of each model [nm] (MLCT-Bio probes for the
# blunted tips, smallest nanowire of Garcia & Garcia 2018; radius of the
# vesicle for the thin shell), used until the user changes it for the model
DEFAULT_RADIUS = {PARABOLOID: 75, FLAT_PUNCH: 100, SPHERE: 75,
                  BLUNTED_CONE: 10, BLUNTED_PYRAMID: 20, NEEDLE: 5,
                  NANOWIRE: 50, THIN_SHELL: 50}

# Models with a shell thickness (lipid bilayer of the vesicles: 4 nm by
# default, as in Calò et al. 2014 and Vorselen et al. 2020)
USES_THICKNESS = (THIN_SHELL,)
DEFAULT_THICKNESS = 4.0

# Default half-opening angle of each model [degrees]: 15 for the cones
# (also the end of the needle), 35 for the pyramids (MLCT-Bio probes)
DEFAULT_ANGLE = {CONE: 15, PYRAMID: 35, BLUNTED_CONE: 15,
                 BLUNTED_PYRAMID: 35, NEEDLE: 15}

# Models with a bottom effect correction (Garcia & Garcia, Biophys. J.
# 2018: paraboloid, cone, cylinder; the cone factors are used for the
# pyramid)
HAS_BEC = (PARABOLOID, CONE, PYRAMID, FLAT_PUNCH)

# Models whose validity is limited by the radius of the tip (the paraboloid
# needs delta << R)
RADIUS_LIMITED = (PARABOLOID,)

# Power laws: F = coeff * E * delta ** exponent
EXPONENTS = {PARABOLOID: 1.5, CONE: 2.0, PYRAMID: 2.0, FLAT_PUNCH: 1.0,
             THIN_SHELL: 1.0}

_TABLES = {}


def power_law(model, tip_radius, tip_angle, poisson_ratio,
              thickness=DEFAULT_THICKNESS * 1e-9):
    """(coeff, exponent) of the power-law models (None for the others).
    thickness: thickness of the shell [m] (thin shell)."""
    nu = 1.0 - poisson_ratio ** 2
    if model == PARABOLOID:
        return 4.0 / 3.0 * math.sqrt(tip_radius) / nu, 1.5
    if model == CONE:
        return 2.0 / math.pi * math.tan(math.radians(tip_angle)) / nu, 2.0
    if model == PYRAMID:
        # n-sided pyramid: 2 tan(theta) / (n sin(pi/n)), 1/sqrt(2) tan(theta)
        # for n = 4
        n = PYRAMID_SIDES
        return (2.0 * math.tan(math.radians(tip_angle)) /
                (n * math.sin(math.pi / n)) / nu, 2.0)
    if model == FLAT_PUNCH:
        return 2.0 * tip_radius / nu, 1.0
    if model == THIN_SHELL:
        # Radius of the vesicle
        return 4.0 * thickness ** 2 / (tip_radius * math.sqrt(3.0 * nu)), 1.0
    return None


def sphere_exact(a, radius):
    """Indentation and force per unit reduced modulus E/(1 - nu^2) of a
    sphere, for the contact radii a (< radius)."""
    log = numpy.log((radius + a) / (radius - a))
    delta = a / 2.0 * log
    force = 0.5 * ((a ** 2 + radius ** 2) * log - 2.0 * a * radius)
    return delta, force


def pyramid_factor(n):
    """Factor c = n / pi sin(pi / n) of the n-sided pyramid in the Rico
    model (c = 1 for the cone, n -> infinity)."""
    return n / math.pi * math.sin(math.pi / n)


def transition_radius(tip_radius, tip_angle, c):
    """Radius b = c R / tan(theta) where the paraboloid apex is tangent to
    the faces (of the equivalent cone for the pyramid). It is the only b
    with a smooth transition: for a larger b, the indentation decreases
    after the transition (two contact radii for one indentation)."""
    return c * tip_radius / math.tan(math.radians(tip_angle))


def blunted_exact(a, radius, angle, c):
    """Indentation and force per unit reduced modulus of a blunted cone
    (c = 1) or pyramid (c = pyramid_factor(n)) with a paraboloid apex
    (Briscoe 1994, Rico 2005), for the contact radii a > b."""
    a = numpy.asarray(a, dtype=float)
    tan = math.tan(math.radians(angle))
    b = transition_radius(radius, angle, c)
    arc = numpy.pi / 2 - numpy.arcsin(b / a)
    root = numpy.sqrt(a ** 2 - b ** 2)
    delta = c * a / tan * arc - a / radius * (root - a)
    force = 2.0 * (delta * a - c * a ** 2 / (2 * tan) * arc -
                   a ** 3 / (3 * radius) +
                   root * (c * b / (2 * tan) + root ** 2 / (3 * radius)))
    return delta, force


def _table(model, tip_radius, tip_angle):
    """Indentation and force per unit reduced modulus (log-spaced
    tables)."""
    key = (model, tip_radius, tip_angle)
    if key not in _TABLES:
        if model == SPHERE:
            # a/R from 1e-5 to 1 - 1e-12 (indentation up to about 14 R),
            # increasing
            u = numpy.r_[numpy.logspace(-5, -1, 200, endpoint=False),
                         1.0 - numpy.logspace(math.log10(0.9), -12, 600)]
            delta, force = sphere_exact(u * tip_radius, tip_radius)
        else:
            c = 1.0 if model == BLUNTED_CONE else \
                pyramid_factor(PYRAMID_SIDES)
            tan = math.tan(math.radians(tip_angle))
            b = transition_radius(tip_radius, tip_angle, c)
            # Contact radius for an indentation of 100 um (sharp tip)
            a_max = b + tan * 100e-6 / c
            a = b * (1.0 + numpy.logspace(-6, math.log10(a_max / b - 1.0),
                                          800))
            delta, force = blunted_exact(a, tip_radius, tip_angle, c)
        _TABLES[key] = (numpy.log(delta), numpy.log(force))
    return _TABLES[key]


def needle_force(delta, radius, angle):
    """Force per unit reduced modulus of the needle (flat punch of radius
    capped with a cone of half-angle)."""
    tan = math.tan(math.radians(angle))
    delta_c = math.pi * radius / (2.0 * tan)
    out = 2.0 / math.pi * tan * delta ** 2
    beyond = delta > delta_c
    out[beyond] = 2.0 / math.pi * tan * delta_c ** 2 + \
        2.0 * radius * (delta[beyond] - delta_c)
    return out


def nanowire_force(delta, radius):
    """Force per unit reduced modulus of the nanowire (flat punch of radius
    capped with a hemisphere of the same radius)."""
    out = 4.0 / 3.0 * math.sqrt(radius) * delta ** 1.5
    beyond = delta > radius
    out[beyond] = 4.0 / 3.0 * radius ** 2 + \
        2.0 * radius * (delta[beyond] - radius)
    return out


def _table_force(model, tip_radius, tip_angle, d):
    """Force per unit reduced modulus of the tabulated models, for the
    indentations d > 0."""
    log_delta, log_force = _table(model, tip_radius, tip_angle)
    values = numpy.exp(numpy.interp(numpy.log(d), log_delta, log_force))
    # Small indentations: Hertz (sphere, or contact on the paraboloid apex
    # of the blunted tips)
    small = d < math.exp(log_delta[0])
    values[small] = 4.0 / 3.0 * math.sqrt(tip_radius) * d[small] ** 1.5
    return values


def bec_factor(delta, model, tip_radius, tip_angle, height):
    """Bottom effect correction of Garcia & Garcia (Biophys. J. 2018, eqs.
    12, 14 and 16, cell bonded to the substrate): factor (>= 1) of the
    force of a layer of thickness height, for the models of HAS_BEC (the
    cone factors are used for the pyramid).

    delta, tip_radius and height in the same unit. Returns 1 if the height
    is not positive (or nan).
    """
    delta = numpy.clip(numpy.asarray(delta, dtype=float), 0, None)
    if height is None or not height > 0:
        return numpy.ones_like(delta)
    if model == PARABOLOID:
        x = numpy.sqrt(delta * tip_radius) / height
        coeffs = (1.133, 1.497, 1.469, 0.755)
    elif model in (CONE, PYRAMID):
        x = delta * math.tan(math.radians(tip_angle)) / height
        coeffs = (0.721, 0.650, 0.491, 0.225)
    elif model == FLAT_PUNCH:
        # -0.291: sign of eq. 12a (+0.291 in Lekka et al. 2022, table 2.2.2)
        x = numpy.full_like(delta, tip_radius / height)
        coeffs = (1.133, 1.283, 0.598, -0.291)
    else:
        raise ValueError("No bottom effect correction for the model " +
                         str(model))
    return 1.0 + sum(c * x ** (k + 1) for k, c in enumerate(coeffs))


def unit_force(model, tip_radius, tip_angle, poisson_ratio, delta,
               thickness=DEFAULT_THICKNESS * 1e-9):
    """Force per unit modulus g(delta) (F = E g(delta)), 0 for delta <= 0.

    tip_radius in m, tip_angle in degrees, delta in m (array), thickness of
    the shell in m (thin shell).
    """
    delta = numpy.asarray(delta, dtype=float)
    out = numpy.zeros(delta.shape)
    pos = delta > 0
    law = power_law(model, tip_radius, tip_angle, poisson_ratio, thickness)
    if law is not None:
        coeff, exponent = law
        out[pos] = coeff * delta[pos] ** exponent
        return out
    d = delta[pos]
    if model == NEEDLE:
        values = needle_force(d, tip_radius, tip_angle)
    elif model == NANOWIRE:
        values = nanowire_force(d, tip_radius)
    elif model in (SPHERE, BLUNTED_CONE, BLUNTED_PYRAMID):
        values = _table_force(model, tip_radius, tip_angle, d)
    else:
        raise ValueError("No contact model for the model " + str(model))
    out[pos] = values / (1.0 - poisson_ratio ** 2)
    return out
