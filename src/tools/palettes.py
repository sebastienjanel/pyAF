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


"""Color palettes for the plots of the results (the ones used with ggplot2).

The colors are lists of [r, g, b, 1.0] values between 0 and 1, as the
colors of the results (see experiment.Result).
"""

import math

from matplotlib import colormaps
from matplotlib.colors import to_rgba


def _hex(*colors):
    return [list(to_rgba("#" + color)) for color in colors]


# Fixed palettes (repeated when there are more files than colors)
FIXED = {
    # Colorblind safe (Okabe & Ito, 2008)
    "Okabe-Ito": _hex("E69F00", "56B4E9", "009E73", "0072B2", "D55E00",
                      "CC79A7", "F0E442", "999999"),
    # ggsci palettes of the journals
    "Nature (NPG)": _hex("E64B35", "4DBBD5", "00A087", "3C5488", "F39B7F",
                         "8491B4", "91D1C2", "DC0000", "7E6148", "B09C85"),
    "Science (AAAS)": _hex("3B4992", "EE0000", "008B45", "631879", "008280",
                           "BB0021", "5F559B", "A20056", "808180", "1B1919"),
    "Lancet": _hex("00468B", "ED0000", "42B540", "0099B4", "925E9F",
                   "FDAF91", "AD002A", "ADB6B6", "1B1919"),
    "NEJM": _hex("BC3C29", "0072B5", "E18727", "20854E", "7876B1",
                 "6F99AD", "FFDC91", "EE4C97"),
    "JCO": _hex("0073C2", "EFC000", "868686", "CD534C", "7AA6DC",
                "003C67", "8F7700", "3B3B3B", "A73030", "4A6990"),
}

# ColorBrewer palettes (scale_color_brewer)
BREWER = ("Set1", "Set2", "Dark2", "Paired")

# The name of the palette when the colors were chosen file by file
CUSTOM = "Custom"

NAMES = ["ggplot2", "Viridis"] + list(FIXED) + list(BREWER)


def _hcl_to_rgb(h, c, l):
    """Converts a HCL (polar CIE Luv) color to sRGB, as R does."""
    # D65 white point
    xn, yn, zn = 95.047, 100.000, 108.883
    un = 4 * xn / (xn + 15 * yn + 3 * zn)
    vn = 9 * yn / (xn + 15 * yn + 3 * zn)

    u = c * math.cos(math.radians(h))
    v = c * math.sin(math.radians(h))
    y = yn * (((l + 16) / 116) ** 3 if l > 8 else l / 903.3)
    up = u / (13 * l) + un
    vp = v / (13 * l) + vn
    x = 9 * y * up / (4 * vp)
    z = -x / 3 - 5 * y + 3 * y / vp

    def gamma(t):
        t = t / 100
        t = 12.92 * t if t <= 0.0031308 else 1.055 * t ** (1 / 2.4) - 0.055
        return min(max(t, 0.0), 1.0)

    r = gamma(3.240479 * x - 1.537150 * y - 0.498535 * z)
    g = gamma(-0.969256 * x + 1.875992 * y + 0.041556 * z)
    b = gamma(0.055648 * x - 0.204043 * y + 1.057311 * z)
    return [r, g, b, 1.0]


def get_colors(name, number):
    """Returns number colors of the palette name."""
    if number <= 0:
        return []
    if name == "ggplot2":
        # scales::hue_pal(): hues evenly spaced from 15, l = 65, c = 100
        step = 360 / number
        return [_hcl_to_rgb((15 + i * step) % 360, 100, 65)
                for i in range(number)]
    if name == "Viridis":
        # scale_color_viridis_d(): from purple to yellow
        cmap = colormaps["viridis"]
        if number == 1:
            return [list(cmap(0.0))]
        return [list(cmap(i / (number - 1))) for i in range(number)]
    if name in BREWER:
        colors = [list(c) for c in colormaps[name].colors]
    else:
        colors = FIXED[name]
    return [list(colors[i % len(colors)]) for i in range(number)]
