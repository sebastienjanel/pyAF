"""Viscoelasticity: Ting fit of the approach and the retraction curves.

The power-law rheology of the sample (modulus E0 at 1 s and fluidity) is
fitted on the approach and the retraction (Efremov et al.; Tejedor &
Garcia, Adv. Intell. Syst. 2025, with the bottom effect correction), see
ting.py.

The elasticity fit (Fitting mode) is needed: its point of contact and its
baseline give the indentation and the force, its modulus the initial value
of the fit. The curves are prepared in the main process (they are read in
the HDF5 file) and fitted in the processes of the pool.
"""

import concurrent.futures
import logging
import math
import os
import sys

import numpy

from ... import shared
from ... import widgets_list
from ...single import VISCO_RESULTS
from ...tools import curve_tools
from ...widgets.progressbar import Progressbar
from .. import multiproc
from . import ting
from ..stiffness import contact_models
from ..stiffness.bec import has_height, get_height_map

# Number of curves fitted in a task of the pool
CURVES_PER_TASK = 32


def get_settings(data, used=False):
    """Parameters of the viscoelasticity fit of a file (the parameters used
    for the computation if used is True)."""
    prefix = "used_visco_" if used else "visco_"
    names = ("solver", "bec", "contact_offset", "nbr_points", "beta_auto",
             "beta_init", "drag_correction")
    return {name: getattr(data, prefix + name) for name in names}


def check_file(data):
    """Returns None if the viscoelasticity can be computed for the file,
    else the reason."""
    if not data.stiffness_calculated:
        return "the elasticity has not been computed"
    if not data.used_perform_fit:
        return "the elasticity was computed with Tomography (Fitting needed)"
    if data.used_stiffness_model_selected not in ting.SUPPORTED_MODELS:
        return "the model of the elasticity fit is not supported"
    if data.visco_bec and data.used_stiffness_model_selected not in \
            contact_models.HAS_BEC:
        return "no bottom effect correction for this model"
    if data.visco_bec and not has_height(data):
        return "no substrate ROI or sample height for the B.E.C."
    if data.file_type == "Nanoscope (Peak Force)" and not data.pft_freq:
        return "PeakForce frequency unknown (file loaded with a previous " \
            "version of pyAF: load the .spm file again)"
    has_time = data.file_type in ("JPK (Single File)", "JPK (Force Map)",
                                  "JPK (QI)")
    if not has_time and curve_tools.get_nanoscope_times(data, 2, 2) is None:
        return "no time and no ramp velocities in the file"
    return None


def get_time(data, i, j, nbr_app, nbr_ret):
    """Time of the approach and of the retraction [s], in the convention of
    the JPK files: t_app for the stored approach points, t_ret for the
    stored retraction points in reverse order.

    JPK files: time of the file. Nanoscope files: computed as in PyFMReader
    (see curve_tools.get_nanoscope_times). None if the time is not known.
    """
    if data.file_type in ("JPK (Single File)", "JPK (Force Map)",
                          "JPK (QI)"):
        t_app = numpy.asarray(data.curves_approach[i][j][2], float)
        t_ret = numpy.asarray(data.curves_retraction[i][j][2], float)
        if len(t_app) < 2 or len(t_ret) < 2:
            return None
        dt = t_app[1] - t_app[0]
        pause = data.pause_duration or 0
        # Pause between the segments, not included in the times
        if pause > 0 and t_ret[0] - t_app[-1] < pause / 2:
            t_ret = t_ret + pause
        if t_ret[0] < t_app[-1] + dt / 2:
            # Times of the segments from 0
            t_ret = t_ret - t_ret[0] + t_app[-1] + dt + pause
        return t_app, t_ret

    return curve_tools.get_nanoscope_times(data, nbr_app, nbr_ret)


def downsample(values, nbr):
    """Mean of blocks of the values (nbr blocks, at most the number of
    values)."""
    values = numpy.asarray(values, float)
    nbr = max(1, min(nbr, len(values)))
    edges = numpy.linspace(0, len(values), nbr + 1).astype(int)
    return numpy.add.reduceat(values, edges[:-1]) / numpy.diff(edges)


def prepare_curve(data, i, j, settings):
    """Force, indentation and time of a curve for the Ting fit (SI units).

    Returns (time, delta, force, idx_tm, tc_init, f0_init, has_pause), or
    None if the curve can not be fitted (no point of contact).
    """
    poc = data.pocs_real[i][j]
    z_poc, d_poc = float(poc[0]), float(poc[1])
    if not (numpy.isfinite(z_poc) and numpy.isfinite(d_poc)):
        return None

    app, ret, _, _, _, _ = curve_tools.get_curve(data, [i, j],
                                                mode="curve_results")
    z_app = numpy.asarray(app[0], float)
    d_app = numpy.asarray(app[1], float)
    z_ret = numpy.asarray(ret[0], float)
    d_ret = numpy.asarray(ret[1], float)

    times = get_time(data, i, j, len(z_app), len(z_ret))
    if times is None:
        return None
    t_app = numpy.asarray(times[0], float)
    # Time of the stored retraction points (stored from the end of the
    # curve for the ramps, in the time order for PeakForce)
    t_ret = numpy.asarray(times[1], float)[::-1]

    # Valid parts of the curves (corrupted data, dummy points)
    pos = data.approach_positions[i][j]
    valid = slice(int(pos[0]), int(pos[1]))
    z_app, d_app, t_app = z_app[valid], d_app[valid], t_app[valid]
    pos = data.retraction_positions[i][j]
    valid = slice(int(pos[0]), int(pos[1]))
    z_ret, d_ret, t_ret = z_ret[valid], d_ret[valid], t_ret[valid]
    if len(z_app) < 4 or len(z_ret) < 4:
        return None
    # Time order
    order = numpy.argsort(t_app, kind="stable")
    z_app, d_app, t_app = z_app[order], d_app[order], t_app[order]
    order = numpy.argsort(t_ret, kind="stable")
    z_ret, d_ret, t_ret = z_ret[order], d_ret[order], t_ret[order]

    # Baseline of the search of the point of contact (as in the elasticity
    # fit): its slope is removed, the force is 0 at the contact
    slope = float(data.fits_poc[i][j][0])
    if not numpy.isfinite(slope):
        slope = 0.0
    d_app = d_app - slope * (z_app - z_poc)
    d_ret = d_ret - slope * (z_ret - z_poc)

    # Offset of the retraction baseline (viscous drag): difference of the
    # retraction and the approach far from the contact
    if settings["drag_correction"]:
        z_far = z_app[0] + (z_poc - z_app[0]) / 2
        far_app = z_app < z_far
        far_ret = z_ret < z_far
        if far_app.sum() >= 10 and far_ret.sum() >= 10:
            d_ret = d_ret - (numpy.median(d_ret[far_ret]) -
                             numpy.median(d_app[far_app]))

    k = data.used_spring_constant
    if k is None:
        k = data.spring_constant
    # Time of the point of contact (approach), pause between the segments
    tc_init = float(numpy.interp(z_poc, numpy.maximum.accumulate(z_app),
                                 t_app))
    dt = numpy.median(numpy.diff(t_app))
    has_pause = bool(t_ret[0] - t_app[-1] > 2 * dt)

    # Whole curve. The maximum indentation (end of the loading, tm) is at
    # the maximum of the piezo extension: it can be after the end of the
    # approach segment (the piezo continues for some points)
    z = numpy.r_[z_app, z_ret]
    t = numpy.r_[t_app, t_ret]
    # Indentation [m] and force [N] (spring constant in N/nm)
    defl = numpy.r_[d_app, d_ret] - d_poc
    delta = (z - z_poc - defl) * 1e-9
    force = defl * k
    k_tm = int(numpy.argmax(z))

    # Fitted part: indentation > -contact offset, from the approach to
    # the retraction
    offset = -settings["contact_offset"] * 1e-9
    inside = numpy.nonzero(delta > offset)[0]
    first, last = inside[0] if len(inside) else 0, \
        inside[-1] if len(inside) else 0
    if k_tm - first < 3 or last - k_tm < 4:
        return None
    sa = slice(first, k_tm + 1)
    sr = slice(k_tm + 1, last + 1)

    # Downsampled (mean of blocks), points shared between the loading and
    # the unloading
    nbr_app = sa.stop - sa.start
    nbr_ret = sr.stop - sr.start
    nbr = settings["nbr_points"]
    n_app = max(4, int(round(nbr * nbr_app / (nbr_app + nbr_ret))))
    n_ret = max(4, nbr - n_app)
    t_load = downsample(t[sa], n_app)
    time = numpy.r_[t_load, downsample(t[sr], n_ret)]
    delta = numpy.r_[downsample(delta[sa], n_app),
                     downsample(delta[sr], n_ret)]
    force = numpy.r_[downsample(force[sa], n_app),
                     downsample(force[sr], n_ret)]
    idx_tm = len(t_load) - 1

    if numpy.any(numpy.diff(time) <= 0):
        return None

    f0_init = float(numpy.interp(tc_init, time, force))

    return time, delta, force, idx_tm, tc_init, f0_init, has_pause


def fit_curves(task):
    """Fits the curves of a task (in a process of the pool).

    task: dictionary with "common" (model and settings) and "curves" (list
    of (i, j, time, delta, force, idx_tm, tc_init, f0_init, has_pause,
    e_init, height [m] or nan)).

    Returns the list of (i, j, result) (result: dictionary of fit_ting).
    """
    common = task["common"]
    results = []
    for (i, j, time, delta, force, idx_tm, tc_init, f0_init, has_pause,
         e_init, height) in task["curves"]:
        solver = "numerical" if common["solver"] == 1 or has_pause \
            else "analytical"
        model = ting.TingModel(
            common["model"], common["coeff"], common["exponent"], solver,
            tip_radius=common["tip_radius"], tip_angle=common["tip_angle"],
            height=height if numpy.isfinite(height) else None,
            poisson_ratio=common["poisson"])
        curve = ting.TingCurve(time, delta, force, idx_tm)
        beta_init = None if common["beta_auto"] else common["beta_init"]
        result = ting.fit_ting(curve, model, e_init, tc_init, f0_init,
                               beta_init)
        results.append((i, j, result))
    return results


def get_model(data, settings):
    """Common parameters of the fits of a file."""
    model = data.used_stiffness_model_selected
    tip_radius = data.used_tip_radius * 1e-9
    coeff, exponent = ting.get_geometry(model, tip_radius,
                                        data.used_tip_angle,
                                        data.used_poisson_ratio)
    # Numerical solver with the B.E.C., for PeakForce (sinusoidal
    # modulation: no constant velocities) and the tabulated tips
    numerical = settings["bec"] or exponent is None or \
        data.file_type == "Nanoscope (Peak Force)"
    return {"model": model, "coeff": coeff, "exponent": exponent,
            "tip_radius": tip_radius, "tip_angle": data.used_tip_angle,
            "solver": 1 if numerical else settings["solver"],
            "beta_auto": settings["beta_auto"],
            "beta_init": settings["beta_init"],
            "poisson": data.used_poisson_ratio}


def set_frequency_window(data, contact_times, steps):
    """Frequencies of E' and E'' supported by the curves [Hz]: from
    1 / contact time to 1 / (2 x time between two fitted points) (Nyquist),
    medians of the curves. The frequency of the maps is set to
    1 / contact time (time scale of the measurement), unless the user
    entered a frequency supported by the curves."""
    contact = numpy.asarray(contact_times, dtype=float)
    step = numpy.asarray(steps, dtype=float)
    contact = contact[numpy.isfinite(contact) & (contact > 0)]
    step = step[numpy.isfinite(step) & (step > 0)]
    if contact.size == 0 or step.size == 0:
        data.visco_frequency_window = None
        return
    low = 1.0 / float(numpy.median(contact))
    high = 1.0 / (2.0 * float(numpy.median(step)))
    if high <= low:
        data.visco_frequency_window = None
        return
    data.visco_frequency_window = [low, high]
    if getattr(data, "visco_frequency_user", False) and \
            low <= data.visco_frequency <= high:
        return
    # 1 / contact time rounded up to 2 significant digits (in the window)
    unit = 10.0 ** (math.floor(math.log10(low)) - 1)
    data.visco_frequency = min(math.ceil(low / unit - 1e-9) * unit, high)
    data.visco_frequency_user = False


def clamp_frequency(data, frequency):
    """Frequency of E' and E'' within the window of the file (no limit
    without window: files computed before)."""
    window = getattr(data, "visco_frequency_window", None)
    if window is None:
        return frequency
    return min(max(frequency, window[0]), window[1])


def get_fit_inputs(data, i, j, settings, heights=None):
    """Inputs of the fit of a curve (see fit_curves), or None if the curve
    is not fitted."""
    if data.discarded_curves[i][j]:
        return None
    stiffness = data.stiffness_array[i][j]
    e_init = float(stiffness[0]) if len(stiffness) else numpy.nan
    if not numpy.isfinite(e_init) or e_init <= 0:
        return None
    prepared = prepare_curve(data, i, j, settings)
    if prepared is None:
        return None
    height = numpy.nan
    if heights is not None:
        height = heights[i, j] * 1e-9
    return (i, j) + tuple(prepared) + (e_init, height)


def fit_curve(data, i, j, used=True):
    """Fits a single curve with the parameters used for the map (for the
    display of the fit). Returns (inputs, result) or None."""
    settings = get_settings(data, used)
    heights = None
    if settings["bec"]:
        heights = get_height_map(data)
    inputs = get_fit_inputs(data, i, j, settings, heights)
    if inputs is None:
        return None
    task = {"common": get_model(data, settings), "curves": [inputs]}
    return inputs, fit_curves(task)[0][2]


class ViscoCompute:
    """Computes the viscoelasticity of the curves of a file."""

    def __init__(self):
        self.logger = logging.getLogger()
        # Number of curves which could not be fitted (last computation)
        self.nbr_failed = 0
        self.nbr_fitted = 0

    def compute(self, calc_id):
        """Fits all the curves of a file. The results are stored in
        data.visco_results."""
        data = shared.exp.list[calc_id]
        settings = get_settings(data)
        nx, ny = data.nbr_pixels_x, data.nbr_pixels_y

        Progressbar(os.path.basename(data.filename))
        widgets_list.widget_progressbar.set_label("Fitting viscoelasticity")
        widgets_list.widget_progressbar.set_range(0, nx * ny)

        common = get_model(data, settings)
        heights = get_height_map(data) if settings["bec"] else None

        results = numpy.full((nx, ny, len(VISCO_RESULTS)), numpy.nan)
        if heights is not None:
            results[:, :, VISCO_RESULTS.index("height")] = heights

        nbr_workers = max(1, min(shared.exp.cores, nx * ny))
        if sys.platform == "win32":
            nbr_workers = min(nbr_workers, 61)
        executor = multiproc.get_executor(nbr_workers) \
            if nbr_workers > 1 else None

        self.nbr_failed = 0
        self.nbr_fitted = 0
        futures = []

        def store(task_results):
            for i, j, result in task_results:
                if result["success"]:
                    self.nbr_fitted += 1
                    for name, key in (("e0", "e0"), ("fluidity", "beta"),
                                      ("tc", "tc"), ("f0", "f0"),
                                      ("r2", "r2")):
                        results[i, j, VISCO_RESULTS.index(name)] = \
                            result[key]
                else:
                    self.nbr_failed += 1
            widgets_list.widget_progressbar.update(len(task_results))

        def submit(curves):
            task = {"common": common, "curves": curves}
            if executor is None:
                store(fit_curves(task))
            else:
                futures.append(executor.submit(fit_curves, task))

        # Time scales probed by the curves: contact time (from the contact
        # to the end of the fitted part) and time between two fitted points
        contact_times = []
        steps = []
        try:
            curves = []
            for i in range(nx):
                for j in range(ny):
                    inputs = get_fit_inputs(data, i, j, settings, heights)
                    if inputs is None:
                        widgets_list.widget_progressbar.update()
                        continue
                    time, tc_init = inputs[2], inputs[6]
                    contact_times.append(time[-1] - tc_init)
                    steps.append(numpy.median(numpy.diff(time)))
                    curves.append(inputs)
                    if len(curves) == CURVES_PER_TASK:
                        submit(curves)
                        curves = []
                    # Results of the finished tasks
                    done = [f for f in futures if f.done()]
                    for future in done:
                        futures.remove(future)
                        store(future.result())
            if curves:
                submit(curves)
            for future in concurrent.futures.as_completed(futures):
                store(future.result())
        except BaseException:
            for future in futures:
                future.cancel()
            multiproc.shutdown_pool(wait=False)
            widgets_list.widget_progressbar.close()
            raise

        data.visco_results = results.astype(numpy.float32)
        for name, value in settings.items():
            setattr(data, "used_visco_" + name, value)
        set_frequency_window(data, contact_times, steps)
        data.visco_calculated = True
        data.update()

        widgets_list.widget_progressbar.close()
        self.logger.debug("Viscoelasticity: %d curves fitted, %d failed",
                          self.nbr_fitted, self.nbr_failed)
