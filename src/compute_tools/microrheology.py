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

"""Microrheology: complex modulus G*(f) = G'(f) + i G''(f) from sine
modulations of the piezo in contact (JPK force maps with modulation
segments), as in PyFMRheo (routines MicrorheologyFFT and
MicrorheologySine, github.com/DyNaMo-INSERM/PyFMlab_DyNaMo; Alcaraz et
al., Biophys. J. 2003).

For each curve and each modulation segment:
1. Slow drift of the height and of the deflection removed (rolling
   average over one period, forward and backward, as in PyFMRheo).
2. Transfer function H = F / delta at the frequency of the segment (N/m):
   ratio of the Fourier transforms (FFT method) or of the sines fitted at
   the frequency (sine method), F = k d, delta = z - d.
   The lag and the gain of the piezo, measured on a rigid sample (piezo
   characterization), correct the height.
3. Viscous drag of the liquid on the cantilever removed:
   H - i 2 pi f b(0).
4. Linearization of the contact model of the elasticity fit around the
   working indentation delta0 (F = E g(delta)): E* = H / g'(delta0), and
   G* = E* / (2 (1 + nu)). For the paraboloid,
   G* = (1 - nu) / (4 sqrt(R delta0)) H, as in PyFMRheo.

The working indentation is the maximal indentation of the approach of the
elasticity fit (point of contact of the elasticity fit), as in PyFMRheo.
Units: SI (m, N, s, Pa).
"""

import math
import os

import numpy

from .. import shared
from .. import widgets_list
from ..tools import curve_tools
from ..widgets.progressbar import Progressbar
from .stiffness import bec
from .stiffness import contact_models

METHODS = ("fft", "sine")
METHOD_NAMES = ("FFT (transfer function)", "Sine fit")

# Results for each frequency (index in data.microrheo_results[i, j, f])
VALUES = ("storage", "loss", "r2", "indentation")

# Models of the elasticity fit without microrheology (no contact model, or
# not the indentation of the sample)
UNSUPPORTED_MODELS = (contact_models.SLOPE, contact_models.THIN_SHELL)


def in_contact(data):
    """For each modulation segment, True if it is in contact (not above the
    surface: viscous drag measurement, see viscous_drag)."""
    freqs = data.modulation_frequencies or []
    heights = data.modulation_heights
    if heights is None or len(heights) != len(freqs):
        return [True] * len(freqs)
    return [height <= 0 for height in heights]


def available(data):
    """True if the file has modulation segments in contact."""
    return bool(data.modulation_frequencies) and \
        data.modulation_nbr_points is not None and any(in_contact(data))


def check_file(data):
    """Returns None if the microrheology can be computed for the file, else
    the reason."""
    if not available(data):
        return "no sine modulation of the piezo in contact in the file"
    if not data.stiffness_calculated:
        return "the elasticity has not been computed"
    if not data.used_perform_fit:
        return "the elasticity was computed with Tomography (Fitting needed)"
    if data.used_stiffness_model_selected in UNSUPPORTED_MODELS:
        return "no microrheology with the model " + contact_models.MODEL_NAMES[
            data.used_stiffness_model_selected]
    return None


def get_settings(data, used=False):
    """Parameters of the microrheology (the parameters used for the
    computation if used is True): method, viscous drag coefficient
    b(0) [uN s/m] and piezo characterization (None: no correction)."""
    prefix = "used_microrheo_" if used else "microrheo_"
    return {name: getattr(data, prefix + name)
            for name in ("method", "drag", "piezo")}


def detrend(frequency, z, d, t):
    """Height, deflection and time with the slow drift removed: rolling
    average over one period subtracted forward then backward (PyFMRheo,
    detrend_rolling_average); the half periods at the ends are dropped.
    Linear detrending if the segment is shorter than two periods."""
    n = len(z)
    points = int(round(1.0 / frequency / ((t[-1] - t[0]) / (n - 1))))
    if points < 2 or 2 * points >= n:
        z = z - numpy.polyval(numpy.polyfit(t, z, 1), t)
        d = d - numpy.polyval(numpy.polyfit(t, d, 1), t)
        return z, d, t
    kernel = numpy.ones(points) / points
    half = points // 2

    def rolling(values):
        # Centered rolling mean (pandas rolling(center=True)), nan where the
        # window is not complete
        out = numpy.full(n, numpy.nan)
        means = numpy.convolve(values, kernel, mode="valid")
        out[half:half + len(means)] = means
        return out

    z1 = z - rolling(z)
    d1 = d - rolling(d)
    z2 = z1 - rolling(z1[::-1])[::-1]
    d2 = d1 - rolling(d1[::-1])[::-1]
    keep = numpy.isfinite(z2) & numpy.isfinite(d2)
    return z2[keep], d2[keep], t[keep]


def sine_fit(t, values, omega):
    """Sine of angular frequency omega fitted on the values (linear least
    squares of a sin + b cos + c): complex amplitude (amplitude e^(i
    phase), values ~ |A| sin(omega t + phase)) and R2."""
    design = numpy.column_stack((numpy.sin(omega * t),
                                 numpy.cos(omega * t),
                                 numpy.ones(len(t))))
    coeffs, _, _, _ = numpy.linalg.lstsq(design, values, rcond=None)
    model = design @ coeffs
    ss_tot = numpy.sum((values - numpy.mean(values)) ** 2)
    r2 = 1.0 - numpy.sum((values - model) ** 2) / ss_tot \
        if ss_tot > 0 else numpy.nan
    return complex(coeffs[0], coeffs[1]), r2


def amplitudes(frequency, z, d, t, method):
    """Complex amplitudes of the height and of the deflection at the
    frequency (Fourier transforms at the bin of the frequency, PyFMRheo
    TransferFunction, or sines fitted at the frequency), and R2 of the sine
    fitted on the deflection (part of the force at the frequency). nan if
    the bin of the frequency is outside the transform."""
    omega = 2.0 * math.pi * frequency
    t = t - t[0]
    d_amp, r2 = sine_fit(t, d, omega)
    if method == "fft":
        n = len(t)
        step = (t[-1] - t[0]) / (n - 1)
        idx = int(round(frequency * step * n))
        if idx <= 0 or idx >= n // 2:
            nan = complex(numpy.nan, numpy.nan)
            return nan, nan, r2
        return complex(numpy.fft.fft(z)[idx]), \
            complex(numpy.fft.fft(d)[idx]), r2
    z_amp, _ = sine_fit(t, z, omega)
    return z_amp, d_amp, r2


def transfer_function(frequency, z, d, t, method, k, piezo=(0.0, 1.0)):
    """Transfer function H = F / delta at the frequency [N/m] (complex),
    F = k d, delta = z - d, and R2 of the sine fitted on the deflection.

    piezo: phase [degrees] and amplitude quotient of the deflection over the
    height measured on a rigid sample (piezo characterization): the motion
    of the piezo is q e^(i phi) z (PyFMRheo corrects the height amplitude
    by q and H by e^(-i phi), the same for a small deflection)."""
    z_amp, d_amp, r2 = amplitudes(frequency, z, d, t, method)
    phase, quotient = piezo
    delta = quotient * numpy.exp(1j * math.radians(phase)) * z_amp - d_amp
    if not abs(delta) > 0:
        return complex(numpy.nan, numpy.nan), r2
    return complex(k * d_amp / delta), r2


def piezo_correction(data, frequency, used=True):
    """Phase [degrees] and amplitude quotient of the piezo at the frequency
    (piezo characterization of the settings, interpolated linearly between
    its frequencies as in PyFMRheo), (0, 1) without correction."""
    piezo = get_settings(data, used)["piezo"]
    if not piezo:
        return 0.0, 1.0
    freqs = numpy.asarray(piezo["frequencies"], float)
    return (float(numpy.interp(frequency, freqs, piezo["phase"])),
            float(numpy.interp(frequency, freqs, piezo["quotient"])))


def piezo_characterization(data, i, j):
    """Phase [degrees] and amplitude quotient of the deflection over the
    height at each frequency in contact, for a curve on a rigid sample
    (substrate): the deflection follows the piezo, the phase and the
    quotient are the lag and the gain of the piezo (PyFMRheo,
    PiezoCharacterization: Fourier transforms at the frequency). Returns a
    dict (frequencies, phase, quotient, file, curve), None without
    modulation."""
    segments = curve_tools.get_modulation_segments(data, i, j)
    if segments is None:
        return None
    result = {"frequencies": [], "phase": [], "quotient": [],
              "file": os.path.basename(str(data.filename)),
              "curve": [i + 1, j + 1]}
    for (frequency, z, d, t), contact in sorted(
            zip(segments, in_contact(data)), key=lambda s: s[0][0]):
        if not contact or len(z) <= 8:
            continue
        zd, dd, td = detrend(frequency, z, d, t)
        z_amp, d_amp, _ = amplitudes(frequency, zd, dd, td, "fft")
        if not (abs(z_amp) > 0 and numpy.isfinite(abs(d_amp))):
            continue
        ratio = d_amp / z_amp
        result["frequencies"].append(float(frequency))
        result["phase"].append(float(numpy.degrees(numpy.angle(ratio))))
        result["quotient"].append(float(abs(ratio)))
    return result if result["frequencies"] else None


def unit_slope(data, delta0):
    """Derivative g'(delta0) of the force per unit modulus of the contact
    model of the elasticity fit (F = E g(delta)) [N/(Pa m)]."""
    step = 1e-3 * delta0
    g = bec.unit_force(data, numpy.array([delta0 - step, delta0 + step]))
    return (g[1] - g[0]) / (2.0 * step)


def complex_modulus(data, h, frequency, delta0, drag):
    """Shear complex modulus G* [Pa] from the transfer function h [N/m] at
    the working indentation delta0 [m], drag b(0) in N s/m."""
    slope = unit_slope(data, delta0)
    if not slope > 0:
        return complex(numpy.nan, numpy.nan)
    h = h - 1j * 2.0 * math.pi * frequency * drag
    nu = data.used_poisson_ratio
    return h / slope / (2.0 * (1.0 + nu))


def working_indentation(data, i, j):
    """Maximal indentation of the approach [m] (point of contact and
    baseline of the elasticity fit), None without point of contact."""
    points = bec.approach_points(data, i, j)
    if points is None:
        return None
    delta = points[0]
    delta = delta[numpy.isfinite(delta)]
    if not delta.size or not numpy.max(delta) > 0:
        return None
    return float(numpy.max(delta))


def curve(data, i, j, used=True, full=False):
    """Microrheology of a curve: list of results for the frequencies in
    increasing order, each a dict with frequency, storage, loss, r2,
    indentation [m] (and, with full, the detrended time, indentation and
    force). None if the curve can not be analysed."""
    segments = curve_tools.get_modulation_segments(data, i, j)
    delta0 = working_indentation(data, i, j)
    if segments is None or delta0 is None:
        return None
    segments = [segment for segment, contact in
                zip(segments, in_contact(data)) if contact]
    settings = get_settings(data, used)
    drag = (settings["drag"] or 0.0) * 1e-6
    k = data.used_spring_constant
    if k is None:
        k = data.spring_constant
    # Spring constant in N/m, heights and deflections in nm
    k = k * 1e9
    results = []
    for frequency, z, d, t in sorted(segments, key=lambda s: s[0]):
        out = {"frequency": frequency, "indentation": delta0,
               "storage": numpy.nan, "loss": numpy.nan, "r2": numpy.nan}
        if len(z) > 8 and numpy.ptp(t) > 0:
            zd, dd, td = detrend(frequency, z, d, t)
            if len(zd) > 8:
                h, r2 = transfer_function(
                    frequency, zd, dd, td, settings["method"], k,
                    piezo_correction(data, frequency, used))
                g = complex_modulus(data, h, frequency, delta0, drag)
                out.update(storage=g.real, loss=g.imag, r2=r2)
                if full:
                    out.update(time=td - td[0], delta=(zd - dd) * 1e-9,
                               force=dd * 1e-9 * k)
        results.append(out)
    return results


def frequencies(data):
    """Frequencies of the modulation segments in contact in increasing
    order [Hz]."""
    return sorted(f for f, contact in
                  zip(data.modulation_frequencies or [], in_contact(data))
                  if contact)


def compute(calc_id):
    """Microrheology of all the curves of a file, stored in
    data.microrheo_results [i, j, frequency, VALUES]. Returns the number of
    curves analysed."""
    data = shared.exp.list[calc_id]
    nx, ny = data.nbr_pixels_x, data.nbr_pixels_y
    freqs = frequencies(data)
    for name, value in get_settings(data).items():
        setattr(data, "used_microrheo_" + name, value)

    Progressbar(data.filename.split("/")[-1])
    widgets_list.widget_progressbar.set_label("Getting microrheology")
    widgets_list.widget_progressbar.set_range(0, nx * ny)
    results = numpy.full((nx, ny, len(freqs), len(VALUES)), numpy.nan)
    nbr = 0
    try:
        for i in range(nx):
            for j in range(ny):
                if not data.discarded_curves[i][j]:
                    values = curve(data, i, j)
                    if values is not None:
                        nbr += 1
                        for f, value in enumerate(values):
                            results[i, j, f] = [value[name]
                                                for name in VALUES]
                widgets_list.widget_progressbar.update()
    finally:
        widgets_list.widget_progressbar.close()

    data.microrheo_results = results.astype(numpy.float32)
    data.microrheo_calculated = True
    if data.microrheo_frequency not in freqs and freqs:
        data.microrheo_frequency = freqs[0]
    data.update()
    return nbr
