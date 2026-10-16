"""Closed-form parametric cycle analysis with losses - Mattingly,
*Elements of Propulsion*, ch. 7 ("PARA" equivalents).

These are the textbook equations, written out in the textbook's own
non-dimensional variables, with the "modified" two-gamma gas model
(gamma_c, cp_c upstream of the burner; gamma_t, cp_t downstream).
They serve three purposes:
  1. reproduce the hand-calculations of design problems 7.D1 / 7.D2;
  2. generate the carpet plots (TSFC vs specific thrust);
  3. act as an independent check of the station-by-station solver in
     engine.py, which must give the same answer when run with
     TwoGammaGas.

Notation:  tau = total-temperature ratio across a component,
           pi  = total-pressure ratio across a component,
           tau_r, pi_r = free-stream stagnation/static ratios ("ram"),
           tau_lambda  = cp_t Tt4 / (cp_c T0)  (burner exit enthalpy / ambient enthalpy).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict

from scipy.optimize import brentq, minimize_scalar

from .components import ram_recovery


@dataclass
class GasConstants:
    gamma_c: float = 1.4
    cp_c: float = 1004.0      # J/(kg K)
    gamma_t: float = 1.3
    cp_t: float = 1235.0
    h_PR: float = 42.8e6      # J/kg

    @property
    def R_c(self):
        # R = cp (gamma - 1)/gamma   from cp - cv = R and gamma = cp/cv
        return (self.gamma_c - 1) / self.gamma_c * self.cp_c

    @property
    def R_t(self):
        return (self.gamma_t - 1) / self.gamma_t * self.cp_t


def _ram(M0, T0, gc: GasConstants):
    """Free-stream ratios.  Energy eq. with constant cp:
        Tt0/T0 = 1 + (g-1)/2 M0^2 = tau_r;   isentropic: pi_r = tau_r^(g/(g-1))."""
    g = gc.gamma_c
    a0 = math.sqrt(g * gc.R_c * T0)
    tau_r = 1 + 0.5 * (g - 1) * M0 ** 2
    pi_r = tau_r ** (g / (g - 1))
    return a0, tau_r, pi_r


def _poly_tau(pi, gamma, e):
    """Polytropic compression: tau = pi^((g-1)/(g e)).  For a turbine,
    call with e -> 1/e_t: tau_t = pi_t^((g-1) e_t / g)."""
    return pi ** ((gamma - 1) / (gamma * e))


def _isentropic_eff_comp(pi, tau, gamma):
    """eta = (pi^((g-1)/g) - 1)/(tau - 1): ideal over actual temperature rise."""
    return (pi ** ((gamma - 1) / gamma) - 1) / (tau - 1)


def _exit_mach(Pt_P, gamma):
    """Isentropic Pt/P = (1 + (g-1)/2 M^2)^(g/(g-1)) solved for M."""
    if Pt_P < 1.0:
        raise ValueError(f"nozzle pressure ratio {Pt_P:.4f} < 1: no jet (turbine over-loaded)")
    return math.sqrt(2 / (gamma - 1) * (Pt_P ** ((gamma - 1) / gamma) - 1))


@dataclass
class SeparateFlowResult:
    tau_r: float; pi_r: float; eta_r: float; pi_d: float; tau_lambda: float
    tau_c: float; eta_c: float; tau_f: float; eta_f: float; f: float
    tau_t: float; pi_t: float; eta_t: float
    Pt9_P9: float; M9: float; T9_T0: float; V9_a0: float
    Pt19_P19: float; M19: float; T19_T0: float; V19_a0: float
    specific_thrust: float   # N/(kg/s)
    tsfc: float              # kg/(N s)
    thrust_ratio: float      # core/fan-stream thrust ratio  FR
    eta_P: float; eta_T: float; eta_O: float

    def as_dict(self):
        return asdict(self)


def separate_flow_turbofan(*, M0, T0, Tt4, pi_c, pi_f, alpha, gc: GasConstants = GasConstants(),
                           pi_d_max=0.995, pi_b=0.96, pi_n=0.995, pi_fn=0.995,
                           e_c=0.90, e_f=0.89, e_t=0.89, eta_b=0.999, eta_m=0.99,
                           P0_P9=1.0, P0_P19=1.0) -> SeparateFlowResult:
    """Mattingly ch. 7, turbofan with separate exhaust streams and losses."""
    gcm, gt = gc.gamma_c, gc.gamma_t
    a0, tau_r, pi_r = _ram(M0, T0, gc)
    eta_r = ram_recovery(M0)
    pi_d = pi_d_max * eta_r
    # Burner-exit enthalpy relative to ambient enthalpy (cp_t Tt4)/(cp_c T0)
    tau_l = gc.cp_t * Tt4 / (gc.cp_c * T0)
    tau_c = _poly_tau(pi_c, gcm, e_c)
    eta_c = _isentropic_eff_comp(pi_c, tau_c, gcm)
    tau_f = _poly_tau(pi_f, gcm, e_f)
    eta_f = _isentropic_eff_comp(pi_f, tau_f, gcm)
    # Burner energy balance, per kg air:  cp_c Tt3 + f eta_b h_PR = (1+f) cp_t Tt4
    f = (tau_l - tau_r * tau_c) / (eta_b * gc.h_PR / (gc.cp_c * T0) - tau_l)
    # Power balance: turbine drives HPC (core air) + fan (core+bypass air)
    #   eta_m (1+f) cp_t Tt4 (1 - tau_t) = cp_c T0 tau_r [(tau_c - 1) + alpha (tau_f - 1)]
    tau_t = 1 - 1 / (eta_m * (1 + f)) * tau_r / tau_l * (tau_c - 1 + alpha * (tau_f - 1))
    # Turbine polytropic: pi_t = tau_t^(g_t/((g_t-1) e_t))
    pi_t = tau_t ** (gt / ((gt - 1) * e_t))
    eta_t = (1 - tau_t) / (1 - tau_t ** (1 / e_t))
    # Core nozzle: chain the total-pressure ratios from P0 to Pt9
    Pt9_P9 = P0_P9 * pi_r * pi_d * pi_c * pi_b * pi_t * pi_n
    M9 = _exit_mach(Pt9_P9, gt)
    # T9/T0 = (Tt9/T0)/(Tt9/T9), with Tt9/T0 = tau_lambda tau_t cp_c/cp_t
    T9_T0 = tau_l * tau_t / Pt9_P9 ** ((gt - 1) / gt) * gc.cp_c / gc.cp_t
    V9_a0 = M9 * math.sqrt(gt * gc.R_t * T9_T0 / (gcm * gc.R_c))
    # Fan nozzle
    Pt19_P19 = P0_P19 * pi_r * pi_d * pi_f * pi_fn
    M19 = _exit_mach(Pt19_P19, gcm)
    T19_T0 = tau_r * tau_f / Pt19_P19 ** ((gcm - 1) / gcm)
    V19_a0 = M19 * math.sqrt(T19_T0)
    # Uninstalled thrust per unit total air flow (momentum + pressure terms):
    #  F/mdot0 = a0/(1+alpha) [ (1+f) V9/a0 - M0 + (1+f) (R_t/R_c) (T9/T0)/(V9/a0) (1 - P0/P9)/g_c ]
    #          + a0 alpha/(1+alpha) [ V19/a0 - M0 + (T19/T0)/(V19/a0) (1 - P0/P19)/g_c ]
    core = (1 + f) * V9_a0 - M0 + (1 + f) * gc.R_t / gc.R_c * T9_T0 / V9_a0 * (1 - P0_P9) / gcm
    fan = V19_a0 - M0 + T19_T0 / V19_a0 * (1 - P0_P19) / gcm
    F = a0 / (1 + alpha) * core + a0 * alpha / (1 + alpha) * fan
    S = f / ((1 + alpha) * F)
    FR = core / fan / alpha if alpha > 0 else float("inf")
    # Efficiencies (Mattingly eqs. for eta_P, eta_T with losses)
    dKE = (a0 ** 2 / 2) * ((1 + f) * V9_a0 ** 2 + alpha * V19_a0 ** 2 - (1 + alpha) * M0 ** 2)
    eta_T = dKE / (f * gc.h_PR)
    eta_P = (1 + alpha) * F * M0 * a0 / dKE if M0 > 0 else 0.0
    return SeparateFlowResult(tau_r, pi_r, eta_r, pi_d, tau_l, tau_c, eta_c, tau_f, eta_f, f,
                              tau_t, pi_t, eta_t, Pt9_P9, M9, T9_T0, V9_a0,
                              Pt19_P19, M19, T19_T0, V19_a0, F, S, FR,
                              eta_P, eta_T, eta_P * eta_T)


def optimum_bypass_separate(**kw):
    """Bypass ratio that minimises TSFC at fixed pi_c, pi_f (Mattingly's
    'optimum alpha').  As alpha grows, more turbine work goes into the
    bypass stream; S falls until the core jet is drained so far that its
    exit velocity drops below the fan jet - beyond that, extra work moved
    to the fan is wasted.  Bounded by the alpha at which the core nozzle
    pressure ratio Pt9/P9 falls to 1 (no core jet at all)."""
    kw = dict(kw)
    kw.pop("alpha", None)

    def ok(a):
        try:
            r = separate_flow_turbofan(alpha=a, **kw)
            return r.Pt9_P9 > 1.0001 and r.tau_t > 0
        except (ValueError, ZeroDivisionError):
            return False

    hi = 0.1
    while ok(hi * 1.25) and hi < 100:
        hi *= 1.25
    a_max = brentq(lambda a: (1.0 if ok(a) else -1.0), hi, hi * 1.25, xtol=1e-6) if ok(hi) else hi
    res = minimize_scalar(lambda a: separate_flow_turbofan(alpha=a, **kw).tsfc,
                          bounds=(0.01, a_max * 0.9999), method="bounded", options={"xatol": 1e-7})
    return float(res.x), separate_flow_turbofan(alpha=float(res.x), **kw)


# ---------------------------------------------------------------- mixed
def phi_mixer(M, g):
    """Mattingly's mixer function  Phi(M, g) = M^2 (1 + (g-1)/2 M^2) / (1 + g M^2)^2.

    It is (mdot sqrt(R Tt/g) / I)^2 with I = P A (1 + g M^2) the impulse
    function: mass flow x sqrt(total temperature) over impulse.  Mass,
    energy and momentum conservation across a constant-area mixer
    therefore fix Phi of the mixed stream directly.
    """
    return M * M * (1 + 0.5 * (g - 1) * M * M) / (1 + g * M * M) ** 2


def mach_from_phi(phi, g):
    """Subsonic inverse of Phi:  M^2 = 2 Phi / (1 - 2 g Phi + sqrt(1 - 2 (g+1) Phi))."""
    disc = 1 - 2 * (g + 1) * phi
    if disc < 0:
        raise ValueError("mixer would choke (Phi above sonic value)")
    return math.sqrt(2 * phi / (1 - 2 * g * phi + math.sqrt(disc)))


def mfp(M, g, R):
    """Mass-flow parameter  mdot sqrt(Tt)/(Pt A) = M sqrt(g/R) (1+(g-1)/2 M^2)^(-(g+1)/(2(g-1)))."""
    return M * math.sqrt(g / R) * (1 + 0.5 * (g - 1) * M * M) ** (-(g + 1) / (2 * (g - 1)))


@dataclass
class MixedFlowResult:
    tau_r: float; pi_r: float; eta_r: float; pi_d: float; tau_lambda: float
    tau_c: float; eta_c: float; tau_f: float; pi_f: float; eta_f: float; f: float
    tau_t: float; pi_t: float; eta_t: float
    M16: float; alpha_prime: float; cp_6A: float; R_6A: float; gamma_6A: float
    tau_M: float; A16_A6: float; M6A: float; pi_M_ideal: float; pi_M: float
    Pt9_P9: float; M9: float; T9_T0: float; V9_a0: float
    specific_thrust: float; tsfc: float; f0: float
    eta_P: float; eta_T: float; eta_O: float

    def as_dict(self):
        return asdict(self)


def mixed_flow_turbofan(*, M0, T0, Tt4, pi_c, alpha, gc: GasConstants = GasConstants(),
                        pi_d_max=0.96, pi_b=0.96, pi_n=0.97, e_c=0.90, e_f=0.89, e_t=0.89,
                        eta_b=0.999, eta_m=0.99, M6=0.4, pi_M_max=0.96, P0_P9=1.0,
                        pi_f=None) -> MixedFlowResult:
    """Mixed-flow turbofan with losses (Mattingly ch. 7 and Problem 7.27).

    With alpha and pi_c given, the fan pressure ratio is found from the
    mixer design condition Pt16 = Pt6, i.e.  pi_f = pi_c pi_b pi_t
    (bypass air at fan exit has the same total pressure as core gas at LPT
    exit).  pi_t depends on tau_f through the power balance, so the system
    is solved for tau_f with a bracketing root-finder.
    """
    gcm, gt = gc.gamma_c, gc.gamma_t
    a0, tau_r, pi_r = _ram(M0, T0, gc)
    eta_r = ram_recovery(M0)
    pi_d = pi_d_max * eta_r
    tau_l = gc.cp_t * Tt4 / (gc.cp_c * T0)
    tau_c = _poly_tau(pi_c, gcm, e_c)
    eta_c = _isentropic_eff_comp(pi_c, tau_c, gcm)
    f = (tau_l - tau_r * tau_c) / (eta_b * gc.h_PR / (gc.cp_c * T0) - tau_l)

    def turbine(tau_f):
        tau_t = 1 - 1 / (eta_m * (1 + f)) * tau_r / tau_l * (tau_c - 1 + alpha * (tau_f - 1))
        return tau_t, tau_t ** (gt / ((gt - 1) * e_t))

    if pi_f is None:
        # Problem 7.27: find tau_f such that  pi_f(tau_f) = pi_c pi_b pi_t(tau_f)
        def resid(tau_f):
            tau_t, pi_t = turbine(tau_f)
            return tau_f ** (gcm * e_f / (gcm - 1)) - pi_c * pi_b * pi_t
        tau_f_max = tau_c
        while turbine(tau_f_max)[0] <= 0:
            tau_f_max = 1 + 0.9 * (tau_f_max - 1)
        tau_f = brentq(resid, 1.0 + 1e-9, tau_f_max, xtol=1e-14)
        pi_f = tau_f ** (gcm * e_f / (gcm - 1))
    else:
        tau_f = _poly_tau(pi_f, gcm, e_f)
    eta_f = _isentropic_eff_comp(pi_f, tau_f, gcm)
    tau_t, pi_t = turbine(tau_f)
    eta_t = (1 - tau_t) / (1 - tau_t ** (1 / e_t))

    # ---- mixer
    # Total pressure ratio of the two streams entering the mixer
    Pt16_Pt6 = pi_f / (pi_c * pi_b * pi_t)
    # Equal static pressure P16 = P6 fixes M16 from M6:
    #   Pt16/P16 = (Pt16/Pt6)(Pt6/P6)
    Pt6_P6 = (1 + 0.5 * (gt - 1) * M6 ** 2) ** (gt / (gt - 1))
    M16 = _exit_mach(Pt16_Pt6 * Pt6_P6, gcm)
    alpha_p = alpha / (1 + f)                 # bypass mass per unit core gas mass
    # Mass-weighted properties of the mixed gas
    cp6A = (gc.cp_t + alpha_p * gc.cp_c) / (1 + alpha_p)
    R6A = (gc.R_t + alpha_p * gc.R_c) / (1 + alpha_p)
    g6A = cp6A / (cp6A - R6A)
    # Energy: cp6A Tt6A (1+a') = cp_t Tt6 + a' cp_c Tt16   ->  tau_M = Tt6A/Tt6
    Tt16_Tt6 = tau_r * tau_f / (tau_l * tau_t) * gc.cp_t / gc.cp_c
    tau_M = gc.cp_t / cp6A * (1 + alpha_p * gc.cp_c / gc.cp_t * Tt16_Tt6) / (1 + alpha_p)
    # Areas from continuity at equal static pressure (A = mdot sqrt(R T/g)/(P M))
    T16_T6 = Tt16_Tt6 * (1 + 0.5 * (gt - 1) * M6 ** 2) / (1 + 0.5 * (gcm - 1) * M16 ** 2)
    A16_A6 = alpha_p * math.sqrt(gc.R_c * T16_T6 * gt / (gc.R_t * gcm)) * M6 / M16
    # Momentum + mass + energy -> Phi of mixed stream:
    #   sqrt(Phi6A) = m6A sqrt(R6A Tt6A/g6A) / (I6 + I16),
    #   I = m sqrt(R Tt/g) / sqrt(Phi)  for each inlet stream
    I6 = math.sqrt(gc.R_t / gt) / math.sqrt(phi_mixer(M6, gt))
    I16 = alpha_p * math.sqrt(gc.R_c * Tt16_Tt6 / gcm) / math.sqrt(phi_mixer(M16, gcm))
    sqrt_phi6A = (1 + alpha_p) * math.sqrt(R6A * tau_M / g6A) / (I6 + I16)
    M6A = mach_from_phi(sqrt_phi6A ** 2, g6A)
    # Ideal mixer total-pressure ratio from continuity, mdot = Pt A MFP / sqrt(Tt):
    #   Pt6A/Pt6 = (1+a') sqrt(tau_M) (A6/A6A) MFP(M6)/MFP(M6A)
    pi_M_ideal = (1 + alpha_p) * math.sqrt(tau_M) / (1 + A16_A6) * mfp(M6, gt, gc.R_t) / mfp(M6A, g6A, R6A)
    pi_M = pi_M_max * pi_M_ideal

    # ---- nozzle
    Pt9_P9 = P0_P9 * pi_r * pi_d * pi_c * pi_b * pi_t * pi_M * pi_n
    M9 = _exit_mach(Pt9_P9, g6A)
    # Tt9/T0 = Tt6A/T0 = tau_M tau_t Tt4/T0 = tau_M tau_t tau_lambda cp_c/cp_t
    Tt9_T0 = tau_M * tau_t * tau_l * gc.cp_c / gc.cp_t
    T9_T0 = Tt9_T0 / Pt9_P9 ** ((g6A - 1) / g6A)
    V9_a0 = M9 * math.sqrt(g6A * R6A * T9_T0 / (gcm * gc.R_c))
    f0 = f / (1 + alpha)                      # fuel per unit TOTAL air
    F = a0 * ((1 + f0) * V9_a0 - M0 + (1 + f0) * R6A / gc.R_c * T9_T0 / V9_a0 * (1 - P0_P9) / gcm)
    S = f0 / F
    dKE = 0.5 * a0 ** 2 * ((1 + f0) * V9_a0 ** 2 - M0 ** 2)
    eta_T = dKE / (f0 * gc.h_PR)
    eta_P = F * M0 * a0 / dKE if M0 > 0 else 0.0
    return MixedFlowResult(tau_r, pi_r, eta_r, pi_d, tau_l, tau_c, eta_c, tau_f, pi_f, eta_f, f,
                           tau_t, pi_t, eta_t, M16, alpha_p, cp6A, R6A, g6A, tau_M, A16_A6, M6A,
                           pi_M_ideal, pi_M, Pt9_P9, M9, T9_T0, V9_a0, F, S, f0,
                           eta_P, eta_T, eta_P * eta_T)
