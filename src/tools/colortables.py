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

"""Definition of the color scales."""

import matplotlib.colors as mpl_colors
import matplotlib.pyplot as mpl_plt
import matplotlib.cm as mpl_cm
import numpy


class ColorTables:
    """Creates a new colortable."""

    def __init__(self, colortableid=None, color_saturation=None,
                 color_negative=None, saturation_limit=None,
                 middle_value=None, mode=None, color_nan=None):
        self.colortables_list = ["Copper",
                                 "Rainbow 1",
                                 "Viridis",
                                 "Inferno",
                                 "Magma",
                                 "Jet",
                                 "Cool-Warm",
                                 "Blue to white",
                                 "Orange to white",
                                 "Green to white",
                                 "Red to white",
                                 "Red to blue",
                                 "Gray",
                                 "Binary"]
        self.colortables_mpl_list = ["copper",
                                     "Rainbow 1",
                                     "viridis",
                                     "inferno",
                                     "magma",
                                     "jet",
                                     "Cool-Warm",
                                     "Blues",
                                     "Oranges",
                                     "Greens",
                                     "Reds",
                                     "RdBu",
                                     "gray",
                                     "binary"]

        self.colortableid = colortableid
        self.color_saturation = color_saturation
        self.color_negative = color_negative
        self.color_nan = color_nan
        if saturation_limit is not None:
            self.saturation_limit = float(saturation_limit)
        if middle_value is not None:
            self.middle_value = float(middle_value)
        if mode == "scalarMap" and colortableid != 1 and colortableid != 6:
            name = self.colortables_mpl_list[self.colortableid]
            cNorm = mpl_colors.Normalize(vmin=0,
                                         vmax=self.saturation_limit)
            self.scalarMap = mpl_cm.ScalarMappable(norm=cNorm,
                                                   cmap=mpl_plt.get_cmap(name))

    def get_color_as_list(self, value):
        """Returns the color as a RGBA list."""
        colortable = self.colortables_list[self.colortableid]
        if colortable == "Rainbow 1":
            if numpy.isnan(value):
                return self.color_nan
            elif value < 0:
                return self.color_negative
            elif value > self.saturation_limit:
                return self.color_saturation
            elif value <= self.saturation_limit and value > self.middle_value:
                # Hard, from yellow to red
                Erel = value - self.middle_value
                Ea = float(Erel) / (self.saturation_limit - self.middle_value)
                return [1.0, -1.0 * Ea + 1.0, 0.0, 1.0]
            elif value <= self.middle_value:
                # Soft, from dark blue to light blue
                Erel = self.middle_value - value
                Ea = float(Erel) / self.middle_value
                return [0.0, -1.0 * Ea + 1.0, 1.0, 1.0]
            else:
                print("error")
        elif colortable == "Cool-Warm":
            if numpy.isnan(value):
                yr = yg = yb = ya = self.color_nan[0], self.color_nan[1], \
                    self.color_nan[2], self.color_nan[3]
            else:
                value = value / self.saturation_limit
                ya = 1.0
                if value < 0:
                    yr = self.color_negative[0]
                    yg = self.color_negative[1]
                    yb = self.color_negative[2]
                elif value >= 0.0 and value < 0.03125:
                    yr = get_val(value, 0.0, 0.03125, 59.0, 68.0)
                    yg = get_val(value, 0.0, 0.03125, 76.0, 90.0)
                    yb = get_val(value, 0.0, 0.03125, 192.0, 204.0)
                elif value >= 0.03125 and value < 0.0625:
                    yr = get_val(value, 0.03125, 0.0625, 68.0, 77.0)
                    yg = get_val(value, 0.03125, 0.0625, 90.0, 104.0)
                    yb = get_val(value, 0.03125, 0.0625, 204.0, 215.0)
                elif value >= 0.0625 and value < 0.09375:
                    yr = get_val(value, 0.0625, 0.09375, 77.0, 87.0)
                    yg = get_val(value, 0.0625, 0.09375, 104.0, 117.0)
                    yb = get_val(value, 0.0625, 0.09375, 215.0, 225.0)
                elif value >= 0.09375 and value < 0.125:
                    yr = get_val(value, 0.09375, 0.125, 87.0, 98.0)
                    yg = get_val(value, 0.09375, 0.125, 117.0, 130.0)
                    yb = get_val(value, 0.09375, 0.125, 225.0, 234.0)
                elif value >= 0.125 and value < 0.15625:
                    yr = get_val(value, 0.125, 0.15625, 98.0, 108.0)
                    yg = get_val(value, 0.125, 0.15625, 130.0, 142.0)
                    yb = get_val(value, 0.125, 0.15625, 234.0, 241.0)
                elif value >= 0.15625 and value < 0.1875:
                    yr = get_val(value, 0.15625, 0.1875, 108.0, 119.0)
                    yg = get_val(value, 0.15625, 0.1875, 142.0, 154.0)
                    yb = get_val(value, 0.15625, 0.1875, 241.0, 247.0)
                elif value >= 0.1875 and value < 0.21875:
                    yr = get_val(value, 0.1875, 0.21875, 119.0, 130.0)
                    yg = get_val(value, 0.1875, 0.21875, 154.0, 165.0)
                    yb = get_val(value, 0.1875, 0.21875, 247.0, 251.0)
                elif value >= 0.21875 and value < 0.25:
                    yr = get_val(value, 0.21875, 0.25, 130.0, 141.0)
                    yg = get_val(value, 0.21875, 0.25, 165.0, 176.0)
                    yb = get_val(value, 0.21875, 0.25, 251.0, 254.0)
                elif value >= 0.25 and value < 0.28125:
                    yr = get_val(value, 0.25, 0.28125, 141.0, 152.0)
                    yg = get_val(value, 0.25, 0.28125, 176.0, 185.0)
                    yb = get_val(value, 0.25, 0.28125, 254.0, 255.0)
                elif value >= 0.28125 and value < 0.3125:
                    yr = get_val(value, 0.28125, 0.3125, 152.0, 163.0)
                    yg = get_val(value, 0.28125, 0.3125, 185.0, 194.0)
                    yb = get_val(value, 0.28125, 0.3125, 255.0, 255.0)
                elif value >= 0.3125 and value < 0.34375:
                    yr = get_val(value, 0.3125, 0.34375, 163.0, 174.0)
                    yg = get_val(value, 0.3125, 0.34375, 194.0, 201.0)
                    yb = get_val(value, 0.3125, 0.34375, 255.0, 253.0)
                elif value >= 0.34375 and value < 0.375:
                    yr = get_val(value, 0.34375, 0.375, 174.0, 184.0)
                    yg = get_val(value, 0.34375, 0.375, 201.0, 208.0)
                    yb = get_val(value, 0.34375, 0.375, 253.0, 249.0)
                elif value >= 0.375 and value < 0.40625:
                    yr = get_val(value, 0.375, 0.40625, 184.0, 194.0)
                    yg = get_val(value, 0.375, 0.40625, 208.0, 213.0)
                    yb = get_val(value, 0.375, 0.40625, 249.0, 244.0)
                elif value >= 0.40625 and value < 0.4375:
                    yr = get_val(value, 0.40625, 0.4375, 194.0, 204.0)
                    yg = get_val(value, 0.40625, 0.4375, 213.0, 217.0)
                    yb = get_val(value, 0.40625, 0.4375, 244.0, 238.0)
                elif value >= 0.4375 and value < 0.46875:
                    yr = get_val(value, 0.4375, 0.46875, 204.0, 213.0)
                    yg = get_val(value, 0.4375, 0.46875, 217.0, 219.0)
                    yb = get_val(value, 0.4375, 0.46875, 238.0, 230.0)
                elif value >= 0.46875 and value < 0.5:
                    yr = get_val(value, 0.46875, 0.5, 213.0, 221.0)
                    yg = get_val(value, 0.46875, 0.5, 219.0, 221.0)
                    yb = get_val(value, 0.46875, 0.5, 230.0, 221.0)
                elif value >= 0.5 and value < 0.53125:
                    yr = get_val(value, 0.5, 0.53125, 221.0, 229.0)
                    yg = get_val(value, 0.5, 0.53125, 221.0, 216.0)
                    yb = get_val(value, 0.5, 0.53125, 221.0, 209.0)
                elif value >= 0.53125 and value < 0.5625:
                    yr = get_val(value, 0.53125, 0.5625, 229.0, 236.0)
                    yg = get_val(value, 0.53125, 0.5625, 216.0, 211.0)
                    yb = get_val(value, 0.53125, 0.5625, 209.0, 197.0)
                elif value >= 0.5625 and value < 0.59375:
                    yr = get_val(value, 0.5625, 0.59375, 236.0, 241.0)
                    yg = get_val(value, 0.5625, 0.59375, 211.0, 204.0)
                    yb = get_val(value, 0.5625, 0.59375, 197.0, 185.0)
                elif value >= 0.59375 and value < 0.625:
                    yr = get_val(value, 0.59375, 0.625, 241.0, 245.0)
                    yg = get_val(value, 0.59375, 0.625, 204.0, 196.0)
                    yb = get_val(value, 0.59375, 0.625, 185.0, 173.0)
                elif value >= 0.625 and value < 0.65625:
                    yr = get_val(value, 0.625, 0.65625, 245.0, 247.0)
                    yg = get_val(value, 0.625, 0.65625, 196.0, 187.0)
                    yb = get_val(value, 0.625, 0.65625, 173.0, 160.0)
                elif value >= 0.65625 and value < 0.6875:
                    yr = get_val(value, 0.65625, 0.6875, 247.0, 247.0)
                    yg = get_val(value, 0.65625, 0.6875, 187.0, 177.0)
                    yb = get_val(value, 0.65625, 0.6875, 160.0, 148.0)
                elif value >= 0.6875 and value < 0.71875:
                    yr = get_val(value, 0.6875, 0.71875, 247.0, 247.0)
                    yg = get_val(value, 0.6875, 0.71875, 177.0, 166.0)
                    yb = get_val(value, 0.6875, 0.71875, 148.0, 135.0)
                elif value >= 0.71875 and value < 0.75:
                    yr = get_val(value, 0.71875, 0.75, 247.0, 244.0)
                    yg = get_val(value, 0.71875, 0.75, 166.0, 154.0)
                    yb = get_val(value, 0.71875, 0.75, 135.0, 123.0)
                elif value >= 0.75 and value < 0.78125:
                    yr = get_val(value, 0.75, 0.78125, 244.0, 241.0)
                    yg = get_val(value, 0.75, 0.78125, 154.0, 141.0)
                    yb = get_val(value, 0.75, 0.78125, 123.0, 111.0)
                elif value >= 0.78125 and value < 0.8125:
                    yr = get_val(value, 0.78125, 0.8125, 241.0, 236.0)
                    yg = get_val(value, 0.78125, 0.8125, 141.0, 127.0)
                    yb = get_val(value, 0.78125, 0.8125, 111.0, 99.0)
                elif value >= 0.8125 and value < 0.84375:
                    yr = get_val(value, 0.8125, 0.84375, 236.0, 229.0)
                    yg = get_val(value, 0.8125, 0.84375, 127.0, 112.0)
                    yb = get_val(value, 0.8125, 0.84375, 99.0, 88.0)
                elif value >= 0.84375 and value < 0.875:
                    yr = get_val(value, 0.84375, 0.875, 229.0, 222.0)
                    yg = get_val(value, 0.84375, 0.875, 112.0, 96.0)
                    yb = get_val(value, 0.84375, 0.875, 88.0, 77.0)
                elif value >= 0.875 and value < 0.90625:
                    yr = get_val(value, 0.875, 0.90625, 222.0, 213.0)
                    yg = get_val(value, 0.875, 0.90625, 96.0, 80.0)
                    yb = get_val(value, 0.875, 0.90625, 77.0, 66.0)
                elif value >= 0.90625 and value < 0.9375:
                    yr = get_val(value, 0.90625, 0.9375, 213.0, 203.0)
                    yg = get_val(value, 0.90625, 0.9375, 80.0, 62.0)
                    yb = get_val(value, 0.90625, 0.9375, 66.0, 56.0)
                elif value >= 0.9375 and value < 0.96875:
                    yr = get_val(value, 0.9375, 0.96875, 203.0, 192.0)
                    yg = get_val(value, 0.9375, 0.96875, 62.0, 40.0)
                    yb = get_val(value, 0.9375, 0.96875, 56.0, 47.0)
                elif value >= 0.96875 and value <= 1.0:
                    yr = get_val(value, 0.96875, 1.0, 192.0, 180.0)
                    yg = get_val(value, 0.96875, 1.0, 40.0, 4.0)
                    yb = get_val(value, 0.96875, 1.0, 47.0, 38.0)
                elif value > 1.0:
                    yr = self.color_saturation[0]
                    yg = self.color_saturation[1]
                    yb = self.color_saturation[2]

            return [yr, yg, yb, ya]

        else:
            if numpy.isnan(value):
                return [0, 0, 0, 0]
            else:
                return self.scalarMap.to_rgba(value)

    def get_color_as_cdict(self, freq=256):
        """Returns the colorscale as a cdict."""
        colortable = self.colortables_list[self.colortableid]
        if colortable == "Rainbow 1":
            # Normalize middle value if needed
            if self.middle_value > 1.0:
                middle_value = self.middle_value / self.saturation_limit
            else:
                middle_value = self.middle_value
            cdict = {"red": ((0.0, 0.0, 0.0),
                             (middle_value, 0.0, 1.0),
                             (1.0, 1.0, 0.0)),
                     "green": ((0.0, 0.0, 0.0),
                               (middle_value, 1.0, 1.0),
                               (1.0, 0.0, 0.0)),
                     "blue": ((0.0, 0.0, 1.0),
                              (middle_value, 1.0, 0.0),
                              (1.0, 0.0, 0.0))}
            colormap = mpl_colors.LinearSegmentedColormap("mycolormap",
                                                          cdict, freq)
            colormap.set_over(self.color_saturation)
            colormap.set_under(self.color_negative)
            colormap.set_bad(self.color_nan)
            norm = mpl_colors.Normalize(vmin=0, vmax=self.saturation_limit)
        elif colortable == "Cool-Warm":
            cdict = {"red": ((0.0, 0.0, 59 / 256.0),
                             (0.03125, 59 / 256.0, 68 / 256.0),
                             (0.0625, 68 / 256.0, 77 / 256.0),
                             (0.09375, 77 / 256.0, 87 / 256.0),
                             (0.125, 87 / 256.0, 98 / 256.0),
                             (0.15625, 98 / 256.0, 108 / 256.0),
                             (0.1875, 108 / 256.0, 119 / 256.0),
                             (0.21875, 119 / 256.0, 130 / 256.0),
                             (0.25, 130 / 256.0, 141 / 256.0),
                             (0.28125, 141 / 256.0, 152 / 256.0),
                             (0.3125, 152 / 256.0, 163 / 256.0),
                             (0.34375, 163 / 256.0, 174 / 256.0),
                             (0.375, 174 / 256.0, 184 / 256.0),
                             (0.40625, 184 / 256.0, 194 / 256.0),
                             (0.4375, 194 / 256.0, 204 / 256.0),
                             (0.46875, 204 / 256.0, 213 / 256.0),
                             (0.5, 213 / 256.0, 221 / 256.0),
                             (0.53125, 221 / 256.0, 229 / 256.0),
                             (0.5625, 229 / 256.0, 236 / 256.0),
                             (0.59375, 236 / 256.0, 241 / 256.0),
                             (0.625, 241 / 256.0, 245 / 256.0),
                             (0.65625, 245 / 256.0, 247 / 256.0),
                             (0.6875, 247 / 256.0, 247 / 256.0),
                             (0.71875, 247 / 256.0, 247 / 256.0),
                             (0.75, 247 / 256.0, 244 / 256.0),
                             (0.78125, 244 / 256.0, 241 / 256.0),
                             (0.8125, 241 / 256.0, 236 / 256.0),
                             (0.84375, 236 / 256.0, 229 / 256.0),
                             (0.875, 229 / 256.0, 222 / 256.0),
                             (0.90625, 222 / 256.0, 213 / 256.0),
                             (0.9375, 213 / 256.0, 203 / 256.0),
                             (0.96875, 203 / 256.0, 192 / 256.0),
                             (1.0, 192 / 256.0, 180 / 256.0)),
                     "green": ((0.0, 0.0, 76 / 256.0),
                               (0.03125, 76 / 256.0, 90 / 256.0),
                               (0.0625, 90 / 256.0, 104 / 256.0),
                               (0.09375, 104 / 256.0, 117 / 256.0),
                               (0.125, 117 / 256.0, 130 / 256.0),
                               (0.15625, 130 / 256.0, 142 / 256.0),
                               (0.1875, 142 / 256.0, 154 / 256.0),
                               (0.21875, 154 / 256.0, 165 / 256.0),
                               (0.25, 165 / 256.0, 176 / 256.0),
                               (0.28125, 176 / 256.0, 185 / 256.0),
                               (0.3125, 185 / 256.0, 194 / 256.0),
                               (0.34375, 194 / 256.0, 201 / 256.0),
                               (0.375, 201 / 256.0, 208 / 256.0),
                               (0.40625, 208 / 256.0, 213 / 256.0),
                               (0.4375, 213 / 256.0, 217 / 256.0),
                               (0.46875, 217 / 256.0, 219 / 256.0),
                               (0.5, 219 / 256.0, 221 / 256.0),
                               (0.53125, 221 / 256.0, 216 / 256.0),
                               (0.5625, 216 / 256.0, 211 / 256.0),
                               (0.59375, 211 / 256.0, 204 / 256.0),
                               (0.625, 204 / 256.0, 196 / 256.0),
                               (0.65625, 196 / 256.0, 187 / 256.0),
                               (0.6875, 187 / 256.0, 177 / 256.0),
                               (0.71875, 177 / 256.0, 166 / 256.0),
                               (0.75, 166 / 256.0, 154 / 256.0),
                               (0.78125, 154 / 256.0, 141 / 256.0),
                               (0.8125, 141 / 256.0, 127 / 256.0),
                               (0.84375, 127 / 256.0, 112 / 256.0),
                               (0.875, 112 / 256.0, 96 / 256.0),
                               (0.90625, 96 / 256.0, 80 / 256.0),
                               (0.9375, 80 / 256.0, 62 / 256.0),
                               (0.96875, 40 / 256.0, 40 / 256.0),
                               (1.0, 4 / 256.0, 4 / 256.0)),
                     "blue": ((0.0, 0.0, 192 / 256.0),
                              (0.03125, 192 / 256.0, 204 / 256.0),
                              (0.0625, 204 / 256.0, 215 / 256.0),
                              (0.09375, 215 / 256.0, 225 / 256.0),
                              (0.125, 225 / 256.0, 234 / 256.0),
                              (0.15625, 234 / 256.0, 241 / 256.0),
                              (0.1875, 241 / 256.0, 247 / 256.0),
                              (0.21875, 247 / 256.0, 251 / 256.0),
                              (0.25, 251 / 256.0, 254 / 256.0),
                              (0.28125, 254 / 256.0, 255 / 256.0),
                              (0.3125, 255 / 256.0, 255 / 256.0),
                              (0.34375, 255 / 256.0, 253 / 256.0),
                              (0.375, 253 / 256.0, 249 / 256.0),
                              (0.40625, 249 / 256.0, 244 / 256.0),
                              (0.4375, 244 / 256.0, 238 / 256.0),
                              (0.46875, 238 / 256.0, 230 / 256.0),
                              (0.5, 230 / 256.0, 221 / 256.0),
                              (0.53125, 221 / 256.0, 209 / 256.0),
                              (0.5625, 209 / 256.0, 197 / 256.0),
                              (0.59375, 197 / 256.0, 185 / 256.0),
                              (0.625, 185 / 256.0, 173 / 256.0),
                              (0.65625, 173 / 256.0, 160 / 256.0),
                              (0.6875, 160 / 256.0, 148 / 256.0),
                              (0.71875, 148 / 256.0, 135 / 256.0),
                              (0.75, 135 / 256.0, 123 / 256.0),
                              (0.78125, 123 / 256.0, 111 / 256.0),
                              (0.8125, 111 / 256.0, 99 / 256.0),
                              (0.84375, 99 / 256.0, 88 / 256.0),
                              (0.875, 88 / 256.0, 77 / 256.0),
                              (0.90625, 77 / 256.0, 66 / 256.0),
                              (0.9375, 66 / 256.0, 56 / 256.0),
                              (0.96875, 56 / 256.0, 47 / 256.0),
                              (1.0, 47 / 256.0, 38 / 256.0))}

            colormap = mpl_colors.LinearSegmentedColormap("mycolormap",
                                                          cdict, 256)
            colormap.set_over(self.color_saturation)
            colormap.set_under(self.color_negative)
            colormap.set_bad(self.color_nan)
            norm = mpl_colors.Normalize(vmin=0, vmax=self.saturation_limit)
        return colormap, norm


def get_val(value, xa, xb, ya, yb):
    """Utilitary function to compute the values for the Cool-Warm colorscale."""
    ya = ya / 256.0
    yb = yb / 256.0

    a = (yb - ya) / (xb - xa)
    b = ya - a * xa

    return a * value + b
