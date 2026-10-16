"""Component models used by both the design-point and off-design solvers.

Each function implements one physical statement about one component;
the cycle modules chain them station by station.  Station numbering
follows Mattingly / SAE ARP 755:

    0   free stream              13  fan exit, bypass side
    2   fan face (inlet exit)    16  bypass duct exit = mixer bypass inlet
    25  fan exit, core side      6   LPT exit = mixer core inlet
    3   HPC exit                 6A  mixed-out (mixer exit)
    4   burner exit / HPT inlet  8   nozzle throat
    45  HPT exit / LPT inlet     9   nozzle exit
    5   LPT exit
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from scipy.optimize import brentq

from . import flow


@dataclass
class Station:
    """Stagnation state + mass flow at one station."""

    name: str
    Tt: float
    Pt: float
    comp: object
    mdot: float = float("nan")
    static: flow.StaticState | None = None
    area: float | None = None


# ---------------------------------------------------------------- inlet
def ram_recovery(M0: float) -> float:
    """MIL-E-5008B inlet recovery eta_r (Mattingly eq. 6.6).

    Subsonic: no shocks, eta_r = 1.  Supersonic: the external/internal shock
    system loses total pressure; the spec curve is
        eta_r = 1 - 0.075 (M0 - 1)^1.35      for 1 < M0 < 5.
    Total recovery is pi_d = pi_d,max * eta_r, where pi_d,max captures wall
    friction (set by the diffuser 'type' in Mattingly Table 6.2).
    """
    if M0 <= 1.0:
        return 1.0
    return 1.0 - 0.075 * (M0 - 1.0) ** 1.35


def free_stream(gas, amb, M0):
    """Station 0: static ambient + flight velocity -> stagnation state.

    V0 = M0 a0 with a0 = sqrt(gamma(T0) R T0); then Tt0 from the energy
    equation ht0 = h(T0) + V0^2/2 and Pt0 along the isentrope.  For constant
    gamma this is exactly tau_r = 1 + (g-1)/2 M0^2,  pi_r = tau_r^(g/(g-1)).
    """
    air = gas.air()
    a0 = math.sqrt(gas.gamma(amb.T, air) * gas.R(air) * amb.T)
    V0 = M0 * a0
    Tt0, Pt0 = flow.total_from_static(gas, air, amb.T, amb.P, V0)
    return Station("0", Tt0, Pt0, air), V0, a0


def inlet(st0: Station, M0: float, pi_d_max: float) -> Station:
    """Adiabatic diffuser: Tt2 = Tt0, Pt2 = pi_d Pt0."""
    pi_d = pi_d_max * ram_recovery(M0)
    return Station("2", st0.Tt, st0.Pt * pi_d, st0.comp, st0.mdot)


# --------------------------------------------------- compressors / fan
def compress(gas, st: Station, pi: float, e: float, name: str) -> Station:
    """Adiabatic compression at polytropic efficiency e.

    Polytropic definition: every infinitesimal stage has isentropic
    efficiency e, i.e. dh = dh_s / e.  With dh_s = dp/rho = R T dp/p,
        cp dT/T = (R/e) dp/p   =>   phi(Tt_out) = phi(Tt_in) + (R/e) ln(pi).
    (Constant cp: tau = pi^((g-1)/(g e)), Mattingly eq. 6.23.)
    """
    R = gas.R(st.comp)
    phi_out = gas.phi(st.Tt, st.comp) + R / e * math.log(pi)
    Tt = gas.T_from_phi(phi_out, st.comp, T_guess=st.Tt * pi ** 0.3)
    return Station(name, Tt, st.Pt * pi, st.comp, st.mdot)


def compress_by_work(gas, st: Station, dh: float, e: float, name: str):
    """Compressor whose specific work dh is fixed by its turbine (spool
    power balance).  Exit enthalpy h_out = h_in + dh gives Tt_out; the
    pressure ratio then follows from the polytropic relation above solved
    for pi:  ln(pi) = e [phi(Tt_out) - phi(Tt_in)] / R."""
    Tt = gas.T_from_h(gas.h(st.Tt, st.comp) + dh, st.comp, T_guess=st.Tt * 1.5)
    pi = math.exp(e * (gas.phi(Tt, st.comp) - gas.phi(st.Tt, st.comp)) / gas.R(st.comp))
    return Station(name, Tt, st.Pt * pi, st.comp, st.mdot), pi


def isentropic_efficiency_compression(gas, st_in: Station, st_out: Station) -> float:
    """eta = (h_out,s - h_in)/(h_out - h_in): the 'real' efficiency implied by
    the polytropic efficiency at this pressure ratio (always < e)."""
    R = gas.R(st_in.comp)
    Ts = gas.T_from_phi(gas.phi(st_in.Tt, st_in.comp) + R * math.log(st_out.Pt / st_in.Pt), st_in.comp)
    h_in = gas.h(st_in.Tt, st_in.comp)
    return (gas.h(Ts, st_in.comp) - h_in) / (gas.h(st_out.Tt, st_in.comp) - h_in)


# ---------------------------------------------------------------- burner
def burner(gas, st3: Station, Tt4: float, pi_b: float, eta_b: float, h_PR: float):
    """Fuel/air ratio from the first law across the combustor.

    Per unit mass of air:  air enthalpy in + chemical energy released
    = products enthalpy out
        h(Tt3, air) + f eta_b h_PR = (1 + f) h(Tt4, products(f))
    =>  f = [h4 - h3] / [eta_b h_PR - h4]
    h4 itself depends on f through the product composition, so iterate
    (converges in ~4 passes because the dependence is weak).
    Constant-cp limit: Mattingly eq. 7.9, f = (tau_l - tau_r tau_c)/(eta_b h_PR/(cpc T0) - tau_l).
    Total-pressure loss pi_b represents friction + the 'Rayleigh' loss of
    heating a moving gas.
    """
    h3 = gas.h(st3.Tt, st3.comp)
    f = 0.02
    for _ in range(30):
        comp4 = gas.burn(st3.comp, f)
        h4 = gas.h(Tt4, comp4)
        f_new = (h4 - h3) / (eta_b * h_PR - h4)
        if abs(f_new - f) < 1e-12:
            break
        f = f_new
    comp4 = gas.burn(st3.comp, f_new)
    if f_new <= 0:
        raise ValueError("Tt4 below compressor exit temperature: no fuel can be burned")
    return Station("4", Tt4, st3.Pt * pi_b, comp4, st3.mdot * (1 + f_new)), f_new


# --------------------------------------------------------------- turbines
def expand_by_work(gas, st: Station, dh: float, e: float, name: str):
    """Turbine delivering specific work dh (J per kg of turbine gas).

    h_out = h_in - dh fixes Tt_out; for expansion the polytropic definition is
    dh = e dh_s, so  phi(Tt_out) - phi(Tt_in) = e R ln(pi_t)  =>
        ln(pi_t) = [phi(Tt_out) - phi(Tt_in)] / (e R).
    (Constant cp: pi_t = tau_t^(g/((g-1) e)), Mattingly eq. 7.12.)
    """
    h_out = gas.h(st.Tt, st.comp) - dh
    Tt = gas.T_from_h(h_out, st.comp, T_guess=st.Tt * 0.8)
    if Tt <= 0.2 * st.Tt:
        raise ValueError("turbine asked for more work than the gas holds")
    pi = math.exp((gas.phi(Tt, st.comp) - gas.phi(st.Tt, st.comp)) / (e * gas.R(st.comp)))
    return Station(name, Tt, st.Pt * pi, st.comp, st.mdot), pi


def expand_by_pressure_ratio(gas, st: Station, pi: float, e: float, name: str) -> Station:
    phi_out = gas.phi(st.Tt, st.comp) + e * gas.R(st.comp) * math.log(pi)
    Tt = gas.T_from_phi(phi_out, st.comp, T_guess=st.Tt * pi ** 0.25)
    return Station(name, Tt, st.Pt * pi, st.comp, st.mdot)


def choked_turbine_pressure_ratio(gas, st4: Station, A_in: float, A_out: float, e: float) -> float:
    """Pressure ratio of a turbine whose inlet AND exit throats are choked.

    Choked throat:  mdot = A * G*(Tt, Pt), and G* is proportional to Pt, so
    write G* = Pt * g*(Tt).  Continuity through the turbine (same mdot):
        A_in Pt_in g*(Tt_in) = A_out Pt_out g*(Tt_out)
     => pi_t = Pt_out/Pt_in = (A_in/A_out) g*(Tt_in)/g*(Tt_out)
    and Tt_out is tied to pi_t by the polytropic relation, so solve the pair.
    For constant gamma g* ~ 1/sqrt(Tt), giving the classic result that a
    turbine between two choked throats runs at fixed pi_t and tau_t.
    """
    g_in = flow.choked_mass_flux(gas, st4.comp, st4.Tt, 1.0)

    def resid(pi):
        out = expand_by_pressure_ratio(gas, st4, pi, e, "x")
        g_out = flow.choked_mass_flux(gas, st4.comp, out.Tt, 1.0)
        return pi - (A_in / A_out) * g_in / g_out

    return brentq(resid, 0.02, 0.9999, xtol=1e-13)


# ------------------------------------------------------------------ mixer
def mixed_out(gas, s6: Station, s16: Station, A6: float, A16: float, pi_M_max: float) -> Station:
    """Constant-area mixer, control volume from planes 6/16 to 6A.

    Conservation laws (no wall shear inside the ideal CV):
      mass      m6A = m6 + m16
      energy    m6A h(Tt6A) = m6 h(Tt6) + m16 h(Tt16)
      momentum  P6A A6A + m6A V6A = (P6 A6 + m6 V6) + (P16 A16 + m16 V16)
                (the impulse function I = P A + m V is conserved;  A6A = A6 + A16)
    With Tt6A known, pick static T6A: energy gives V6A, continuity gives
    rho6A = m6A/(A6A V6A), state gives P6A = rho R T; iterate T6A until the
    momentum equation balances (subsonic root).  Then Pt6A follows from the
    isentrope.  Wall friction is lumped into pi_M,max (Mattingly Table 6.2):
        Pt6A = pi_M,max * Pt6A,ideal.
    The ideal mixing loss (irreversible mixing of a hot and cold stream) is
    what the momentum balance produces automatically.
    """
    m6A = s6.mdot + s16.mdot
    A6A = A6 + A16
    comp = gas.mix([(s6.mdot, s6.comp), (s16.mdot, s16.comp)])
    h6A = (s6.mdot * gas.h(s6.Tt, s6.comp) + s16.mdot * gas.h(s16.Tt, s16.comp)) / m6A
    Tt6A = gas.T_from_h(h6A, comp, T_guess=s6.Tt)
    I_in = (s6.static.P * A6 + s6.mdot * s6.static.V
            + s16.static.P * A16 + s16.mdot * s16.static.V)
    R = gas.R(comp)

    def impulse(T):
        V = math.sqrt(max(2 * (h6A - gas.h(T, comp)), 1e-12))
        P = m6A / (A6A * V) * R * T
        return P * A6A + m6A * V, P, V

    Ts = flow.sonic_T(gas, comp, Tt6A)
    # impulse(T) -> infinity as V -> 0 and is minimal at M = 1, so the
    # subsonic root lies between the sonic temperature and Tt6A
    T6A = brentq(lambda T: impulse(T)[0] - I_in, Ts, Tt6A * (1 - 1e-10), xtol=1e-12 * Tt6A)
    _, P6A, V6A = impulse(T6A)
    Pt_ideal = P6A * math.exp((gas.phi(Tt6A, comp) - gas.phi(T6A, comp)) / R)
    st = flow.static_from_T(gas, comp, Tt6A, Pt_ideal, T6A)
    return Station("6A", Tt6A, Pt_ideal * pi_M_max, comp, m6A, static=st, area=A6A)
