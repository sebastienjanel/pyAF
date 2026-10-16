"""Ting model of a power-law viscoelastic material.

The relaxation modulus of the material is a power law (Efremov et al.,
Sci. Rep. 2017, Soft Matter 2019):

    E(t) = E0 * (t / t0) ** -beta        (t0 = 1 s)

where beta is the fluidity (0: elastic solid, 1: viscous liquid). The
force during the approach and the retraction is given by Ting's solution
(Ting, J. Appl. Mech. 1966):

    approach:   F(t) = C * int_0^t E(t - t') dU/dt' dt'
    retraction: F(t) = C * int_0^t1(t) E(t - t') dU/dt' dt'

with U = delta ** n (C, n: Hertz coefficient and exponent of the tip) and
t1(t) defined by int_t1^t E(t - t') d(delta)/dt' dt' = 0. With a bottom
effect correction (Garcia & Garcia, Biophys. J. 2018), U is multiplied by
the polynomial correction of the indentation over the height of the
sample, inside the integral (Tejedor & Garcia, Adv. Intell. Syst. 2025).

Two solvers:
- analytical: constant approach and retraction velocities, no bottom
  effect correction (closed forms with the hypergeometric function).
- numerical: the measured indentation is used (linear between the
  points); the integrals are exact for this piecewise linear history.

The fit (E0, tc, beta, F0 free; tc is the time of contact) follows the
TingFit routine of PyFMRheo (DyNaMo, INSERM, MIT license).

Units: SI (m, N, s, Pa).
"""

import math
import warnings

import numpy
from scipy.optimize import least_squares
from scipy.special import beta as beta_function
from scipy.special import hyp2f1

T0 = 1.0

# Models of the elasticity fit (index of the model in pyAF)
PARABOLOID = 0
CONE = 1
PYRAMID = 2
FLAT_PUNCH = 4
SUPPORTED_MODELS = (PARABOLOID, CONE, PYRAMID, FLAT_PUNCH)

BETA_MIN = 0.01
BETA_MAX = 0.99


def get_geometry(model, tip_radius, tip_angle, poisson_ratio):
    """Returns the coefficient C and the exponent n of the tip
    (F = C * E * delta ** n for an elastic material).

    tip_radius in m, tip_angle (half angle) in degrees.
    """
    nu = 1.0 - poisson_ratio ** 2
    if model == PARABOLOID:
        return 4.0 / 3.0 * math.sqrt(tip_radius) / nu, 1.5
    if model == CONE:
        return 2.0 / math.pi * math.tan(math.radians(tip_angle)) / nu, 2.0
    if model == PYRAMID:
        return 1.0 / math.sqrt(2) * math.tan(math.radians(tip_angle)) / nu, 2.0
    if model == FLAT_PUNCH:
        return 2.0 * tip_radius / nu, 1.0
    raise ValueError("Model not supported by the Ting fit: " + str(model))


def bec_factor(delta, model, tip_radius, tip_angle, height):
    """Bottom effect correction of Garcia & Garcia (Biophys. J. 2018):
    factor of the force of an elastic material (>= 1).

    delta, tip_radius and height in the same unit. Returns 1 if the height
    is not positive.
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
        x = numpy.full_like(delta, tip_radius / height)
        coeffs = (1.133, 1.283, 0.598, -0.291)
    else:
        raise ValueError("Model not supported by the Ting fit: " + str(model))
    return 1.0 + sum(c * x ** (k + 1) for k, c in enumerate(coeffs))


class TingCurve:
    """Data of a curve for the Ting fit.

    time, delta (indentation) and force: approach followed by the
    retraction, in the order of the time. idx_tm is the index of the end
    of the approach (maximum indentation).
    """

    def __init__(self, time, delta, force, idx_tm):
        self.time = numpy.asarray(time, dtype=float)
        self.delta = numpy.asarray(delta, dtype=float)
        self.force = numpy.asarray(force, dtype=float)
        self.idx_tm = int(idx_tm)


class TingModel:
    """Ting model for a tip and a solver.

    model: index of the elasticity model in pyAF, coeff and exponent: from
    get_geometry. With a height (bottom effect correction), the numerical
    solver is used. tip_radius in m, height in m.
    """

    def __init__(self, model, coeff, exponent, solver="analytical",
                 tip_radius=None, tip_angle=None, height=None):
        self.model = model
        self.coeff = coeff
        self.exponent = exponent
        self.tip_radius = tip_radius
        self.tip_angle = tip_angle
        self.height = height
        if height is not None and height > 0:
            solver = "numerical"
        self.solver = solver

    def u_function(self, delta):
        """U(delta) = delta ** n, with the bottom effect correction."""
        delta = numpy.clip(delta, 0, None)
        u = delta ** self.exponent
        if self.height is not None and self.height > 0:
            u = u * bec_factor(delta, self.model, self.tip_radius,
                               self.tip_angle, self.height)
        return u

    def evaluate(self, curve, e0, tc, beta, f0, velocities=None):
        """Force of the model at the times of the curve."""
        if self.solver == "analytical":
            vt, vr = velocities if velocities is not None else \
                get_velocities(curve, tc)
            force = self.solve_analytical(curve, e0, tc, beta, vt, vr)
        else:
            force = self.solve_numerical(curve, e0, tc, beta)
        return force + f0

    def solve_analytical(self, curve, e0, tc, beta, vt, vr):
        """Closed forms for constant velocities vt (approach) and vr
        (retraction): delta = vt * t, then vt * tm - vr * (t - tm)."""
        n = self.exponent
        time = curve.time - tc
        tm = time[curve.idx_tm]
        force = numpy.zeros(len(time))
        if tm <= 0 or vt <= 0:
            return force
        pref = self.coeff * e0 * T0 ** beta * vt ** n

        # Approach: F = C E0 t0^b vt^n n B(n, 1 - b) t^(n - b)
        app = (time > 0) & (numpy.arange(len(time)) <= curve.idx_tm)
        force[app] = pref * n * beta_function(n, 1 - beta) * \
            time[app] ** (n - beta)

        # Retraction: t1 = t - (1 + vr / vt) ^ (1 / (1 - b)) (t - tm),
        # F = C E0 t0^b vt^n t1^n t^-b 2F1(b, n; n + 1; t1 / t)
        ret = numpy.arange(len(time)) > curve.idx_tm
        t = time[ret]
        t1 = t - (1 + vr / vt) ** (1 / (1 - beta)) * (t - tm)
        values = numpy.zeros(len(t))
        ok = t1 > 0
        z = t1[ok] / t[ok]
        values[ok] = pref * t1[ok] ** n * t[ok] ** -beta * \
            hyp2f1(beta, n, n + 1, z)
        force[ret] = values
        return force

    def solve_numerical(self, curve, e0, tc, beta):
        """Ting's integrals for the measured indentation, linear between
        the points (exact integration of the power law)."""
        time = curve.time
        nbr = len(time)
        force = numpy.zeros(nbr)
        idx_tm = curve.idx_tm
        if not time[0] <= tc < time[idx_tm]:
            return force

        # Contact at tc: point inserted, indentation from delta(tc)
        k0 = numpy.searchsorted(time, tc, side="right")
        delta_tc = numpy.interp(tc, time, curve.delta)
        tau = numpy.r_[0.0, time[k0:] - tc]
        delta = numpy.r_[0.0, curve.delta[k0:] - delta_tc]
        m = idx_tm - k0 + 1   # Index of tm in tau
        u = self.u_function(delta)

        dtau = numpy.diff(tau)
        dtau[dtau <= 0] = numpy.nan
        s_u = numpy.nan_to_num(numpy.diff(u) / dtau)       # dU/dt
        s_d = numpy.nan_to_num(numpy.diff(delta) / dtau)   # d(delta)/dt

        # w[i, k] = int_tau_k^tau_k+1 (t0 / (tau_i - t')) ^ b dt'
        # for k < i
        a = 1 - beta
        diff = tau[:, None] - tau[None, :]
        pw = numpy.where(diff > 0, diff, 0.0) ** a
        w = (pw[:, :-1] - pw[:, 1:]) * T0 ** beta / a
        lower = numpy.tril(numpy.ones((len(tau), len(tau) - 1), bool), -1)
        w = numpy.where(lower, w, 0.0)

        result = numpy.zeros(len(tau))
        # Approach
        result[:m + 1] = w[:m + 1] @ s_u

        # Retraction: t1 such as int_t1^t E(t - t') d(delta) = 0. The
        # cumulated integral from the end, g[i, j] = int_tau_j^tau_i, is
        # negative during the retraction and increases with the approach
        # segments: t1 is in the last approach segment where g >= 0.
        if m + 1 < len(tau) and m > 0:
            wr = w[m + 1:]
            g = numpy.cumsum((wr * s_d[None, :])[:, ::-1], axis=1)[:, ::-1]
            g = numpy.c_[g, numpy.zeros(len(g))]
            ga = g[:, :m + 1] >= 0
            found = ga.any(axis=1)
            rows = numpy.nonzero(found)[0]
            # Segment of t1 (jj = m: t1 = tm, beginning of the retraction)
            jj = m - numpy.argmax(ga[rows, ::-1], axis=1)
            jj = numpy.minimum(jj, m - 1)
            at_tm = g[rows, m] >= 0
            ti = tau[m + 1:][rows]
            c = T0 ** beta / a
            # t1 in [tau_jj, tau_jj+1]: int_t1^tau_jj+1 = -g[jj + 1]
            rest = numpy.where(at_tm, 0.0, -g[rows, jj + 1])
            v = s_d[jj]
            with numpy.errstate(divide="ignore", invalid="ignore"):
                x = (ti - tau[jj + 1]) ** a + rest / (v * c)
                t1 = ti - numpy.clip(x, 0, None) ** (1 / a)
            t1 = numpy.where(v > 0, t1, tau[jj + 1])
            t1 = numpy.clip(t1, tau[jj], tau[jj + 1])
            t1 = numpy.where(at_tm, tau[m], t1)
            # Integral of the complete segments before jj, and of the part
            # of the segment jj before t1
            before = numpy.cumsum(wr * s_u[None, :], axis=1)
            before = numpy.c_[numpy.zeros(len(wr)), before]
            partial = s_u[jj] * c * ((ti - tau[jj]) ** a - (ti - t1) ** a)
            values = numpy.zeros(len(wr))
            values[rows] = before[rows, jj] + partial
            result[m + 1:] = values

        force[k0:] = self.coeff * e0 * result[1:]
        return force


def get_velocities(curve, tc):
    """Indentation velocities of the approach and of the retraction
    (linear fits of the indentation in contact)."""
    time = curve.time
    idx_tm = curve.idx_tm
    app = numpy.arange(len(time)) <= idx_tm
    app &= time >= tc
    vt = numpy.nan
    if app.sum() >= 2:
        vt = numpy.polyfit(time[app], curve.delta[app], 1)[0]
    # Retraction while the indentation is positive
    ret = numpy.arange(len(time)) > idx_tm
    delta_tc = numpy.interp(tc, time, curve.delta)
    ret &= curve.delta - delta_tc > 0
    vr = numpy.nan
    if ret.sum() >= 2:
        vr = -numpy.polyfit(time[ret], curve.delta[ret], 1)[0]
    elif numpy.isfinite(vt):
        vr = vt
    return vt, vr


def solve_linear(g, force, f0_min, f0_max, e0_min, e0_max):
    """E0 and F0 minimizing |E0 g + F0 - force| (F0 and E0 in their
    bounds)."""
    n = len(g)
    sg, sf = g.sum(), force.sum()
    sgg, sgf = g @ g, g @ force
    det = n * sgg - sg * sg
    if det > 0:
        e0 = (n * sgf - sg * sf) / det
        f0 = (sf - e0 * sg) / n
    else:
        e0, f0 = numpy.nan, numpy.nan
    if not (numpy.isfinite(f0) and f0_min <= f0 <= f0_max and
            e0_min <= e0 <= e0_max):
        # Bounded solution: best F0 for each bound of E0, best E0 for each
        # bound of F0
        candidates = []
        for f0_b in (f0_min, f0_max):
            e = (sgf - f0_b * sg) / sgg if sgg > 0 else e0_min
            candidates.append((min(max(e, e0_min), e0_max), f0_b))
        for e_b in (e0_min, e0_max):
            f = min(max((sf - e_b * sg) / n, f0_min), f0_max)
            candidates.append((e_b, f))
        e0, f0 = min(candidates, key=lambda c: numpy.sum(
            (c[0] * g + c[1] - force) ** 2))
    return e0, f0


def fit_ting(curve, model, e_init, tc_init, f0_init, beta_init=None,
             f0_range=100e-12):
    """Fits the Ting model to a curve.

    e_init: initial modulus (from the Hertz fit), tc_init: time of the
    point of contact, f0_init: force at the contact. beta_init: None for
    the automatic initial value (0.05 if E > 10 kPa, else 0.25).

    The force is linear in E0 and F0: they are computed by linear least
    squares for each (tc, beta), which are fitted (variable projection,
    same minimum as the fit of the four parameters).

    Returns a dictionary (e0, tc, beta, f0, r2, tan_delta, success); the
    values are nan if the fit failed.
    """
    nan = {"e0": numpy.nan, "tc": numpy.nan, "beta": numpy.nan,
           "f0": numpy.nan, "r2": numpy.nan, "tan_delta": numpy.nan,
           "success": False}

    time = curve.time
    force = curve.force
    tm = time[curve.idx_tm]
    if not (numpy.isfinite(e_init) and e_init > 0 and
            time[0] <= tc_init < tm and len(time) > 8):
        return nan

    if beta_init is None:
        beta_init = 0.05 if e_init > 10e3 else 0.25
    beta_init = min(max(beta_init, BETA_MIN), BETA_MAX)

    # Bounds of E0 (at t0 = 1 s: E(t) = E0 (t / t0)^-b, the Hertz modulus
    # corresponds to the duration of the indentation)
    e0_ref = e_init * ((tm - tc_init) / T0) ** beta_init
    e0_min, e0_max = e0_ref / 1000, e0_ref * 1000
    f0_min, f0_max = f0_init - f0_range, f0_init + f0_range

    tc_min = time[0]
    tc_max = time[curve.idx_tm - 1]
    tc_init = min(max(tc_init, tc_min), tc_max)

    scale_f = numpy.max(numpy.abs(force - f0_init))
    if not scale_f > 0:
        return nan

    # Velocities: from the initial point of contact (analytical solver)
    velocities = None
    if model.solver == "analytical":
        velocities = get_velocities(curve, tc_init)
        if not (numpy.all(numpy.isfinite(velocities)) and
                velocities[0] > 0 and velocities[1] > 0):
            return nan

    def linear_part(x):
        g = model.evaluate(curve, 1.0, x[0], x[1], 0.0, velocities)
        return g, solve_linear(g, force, f0_min, f0_max, e0_min, e0_max)

    def residuals(x):
        g, (e0, f0) = linear_part(x)
        return (e0 * g + f0 - force) / scale_f

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = least_squares(residuals, [tc_init, beta_init],
                                bounds=([tc_min, BETA_MIN],
                                        [tc_max, BETA_MAX]),
                                x_scale=[(tm - time[0]) / 10, 0.1],
                                method="trf")
    except (ValueError, FloatingPointError, numpy.linalg.LinAlgError):
        return nan

    if not res.success or not numpy.all(numpy.isfinite(res.x)):
        return nan

    tc, beta = res.x
    g, (e0, f0) = linear_part(res.x)
    if not (numpy.isfinite(e0) and e0 > 0):
        return nan
    fitted = e0 * g + f0
    ss_res = numpy.sum((force - fitted) ** 2)
    ss_tot = numpy.sum((force - numpy.mean(force)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else numpy.nan

    return {"e0": float(e0), "tc": float(tc), "beta": float(beta),
            "f0": float(f0), "r2": float(r2),
            "tan_delta": float(math.tan(math.pi * beta / 2)),
            "success": True}
