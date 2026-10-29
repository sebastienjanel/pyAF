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

"""Plot the curve on the data tab."""
"""Force vs distance curve of the data tab, to check the deflection
sensitivity: on a hard substrate (glass), the contact part is vertical."""

import numpy

from ..plot_curve import PlotCurve


class PlotDataCurveForce(PlotCurve):
    """Force (spring constant x deflection, from the baseline) vs distance
    (Z position - deflection, the tip-sample separation). With the right
    deflection sensitivity, the deflection follows the piezo on a hard
    substrate: the contact is vertical."""

    def __init__(self, parent):
        super().__init__(parent)

        self.xlabel = "Distance (Z - deflection) [nm]"
        self.ylabel = "Force [nN]"

        # Current spring constant (N/m) and deflection (sensitivity applied
        # in get_curve)
        k = self.data.spring_constant * 1e9
        base = None
        if self.approach is not None:
            first = int(self.app_pos[0])
            base = self.approach[1][first]

        self.curve_approach = self.to_force(self.approach, k, base)
        self.curve_retraction = self.to_force(self.retraction, k, base)

        self.plot_fig()

    def set_lims(self):
        """Zoom on the contact of the approach curve, if asked (Zoom on
        contact in the data tab): a slope of the contact due to a wrong
        deflection sensitivity is visible at this scale (a few nm), not on
        the whole curve (microns)."""
        zoom = self.settings.value("ForceDistanceZoom", False)
        if zoom in ("false", False, None) or self.curve_approach is None:
            super().set_lims()
            return
        seg = self.valid_seg_app
        x = numpy.asarray(self.curve_approach[0][seg], dtype=float)
        f = numpy.asarray(self.curve_approach[1][seg], dtype=float)
        if x.size < 2 or not numpy.isfinite(f).any():
            return
        fmax = numpy.nanmax(f)
        # Contact: above the middle between the baseline (median of the
        # curve, mostly out of contact) and the maximum force (the baseline
        # can drift)
        base = numpy.nanmedian(f)
        if fmax <= base:
            return
        contact = f > base + 0.5 * (fmax - base)
        if not contact.any():
            return
        xc = x[contact]
        width = max(30.0, 5 * float(numpy.ptp(xc)))
        center = float(numpy.median(xc))
        self.axes.set_xlim(center - width, center + width)
        shown = (x >= center - width) & (x <= center + width)
        low = float(numpy.nanmin(f[shown])) if shown.any() else 0.0
        margin = 0.05 * (fmax - low)
        self.axes.set_ylim(low - margin, fmax + margin)

    @staticmethod
    def to_force(curve, k, base):
        """Force vs distance curve, from a deflection vs Z position curve."""
        if curve is None:
            return None
        z = numpy.asarray(curve[0], dtype=float)
        defl = numpy.asarray(curve[1], dtype=float)
        if base is None:
            base = defl[0]
        return [z - defl, k * (defl - base)]
