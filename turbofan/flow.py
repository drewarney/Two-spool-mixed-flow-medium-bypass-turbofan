"""Compressible-flow relations for a thermally perfect gas.

Constant-gamma textbooks write closed forms such as
    Tt/T = 1 + (g-1)/2 M^2,   Pt/P = (Tt/T)^(g/(g-1)),
    MFP(M) = m sqrt(Tt) / (Pt A) = M sqrt(g/R) (1+(g-1)/2 M^2)^(-(g+1)/(2(g-1))).
With cp = cp(T) those closed forms do not exist, so every relation here is
built directly from the three statements they were derived from:

    energy     ht = h(T) + V^2/2            (adiabatic, no work)
    isentrope  phi(Tt) - phi(T) = R ln(Pt/P)
    state      rho = P / (R T),   a = sqrt(gamma(T) R T),   M = V/a

Pick a static temperature T between 0 and Tt and everything else follows,
so each function below is a 1-D root-find in T.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.optimize import brentq, minimize_scalar


@dataclass(frozen=True)
class StaticState:
    T: float      # static temperature  K
    P: float      # static pressure     Pa
    V: float      # velocity            m/s
    M: float      # Mach number
    rho: float    # density             kg/m^3

    @property
    def mass_flux(self) -> float:
        return self.rho * self.V


def static_from_T(gas, comp, Tt, Pt, T) -> StaticState:
    """All static properties at static temperature T on the isentrope (Tt, Pt)."""
    R = gas.R(comp)
    # Isentropic expansion from stagnation: ln(Pt/P) = [phi(Tt) - phi(T)]/R
    P = Pt * math.exp(-(gas.phi(Tt, comp) - gas.phi(T, comp)) / R)
    # Steady adiabatic energy equation: V^2/2 = ht - h
    dh = gas.h(Tt, comp) - gas.h(T, comp)
    V = math.sqrt(max(2.0 * dh, 0.0))
    a = math.sqrt(gas.gamma(T, comp) * R * T)
    return StaticState(T, P, V, V / a, P / (R * T))


def sonic_T(gas, comp, Tt) -> float:
    """Static temperature where M = 1, i.e. 2(ht - h(T)) = gamma(T) R T."""
    R = gas.R(comp)
    ht = gas.h(Tt, comp)
    f = lambda T: 2.0 * (ht - gas.h(T, comp)) - gas.gamma(T, comp) * R * T
    return brentq(f, 0.5 * Tt, Tt, xtol=1e-10 * Tt)


def sonic_state(gas, comp, Tt, Pt) -> StaticState:
    return static_from_T(gas, comp, Tt, Pt, sonic_T(gas, comp, Tt))


def state_from_mach(gas, comp, Tt, Pt, M) -> StaticState:
    if M <= 0:
        return StaticState(Tt, Pt, 0.0, 0.0, Pt / (gas.R(comp) * Tt))
    f = lambda T: static_from_T(gas, comp, Tt, Pt, T).M - M
    return static_from_T(gas, comp, Tt, Pt, brentq(f, 0.05 * Tt, Tt * (1 - 1e-12), xtol=1e-11 * Tt))


def state_from_pressure(gas, comp, Tt, Pt, P) -> StaticState:
    """Isentropic expansion to static pressure P (e.g. nozzle exit at P0)."""
    phi = gas.phi(Tt, comp) - gas.R(comp) * math.log(Pt / P)
    return static_from_T(gas, comp, Tt, Pt, gas.T_from_phi(phi, comp, T_guess=Tt * 0.8))


def choked_mass_flux(gas, comp, Tt, Pt) -> float:
    """(rho V)* - the largest mass flux an isentropic stream can pass.

    In a 1-D channel d(rho V)/dT = 0 exactly at M = 1 (the area-Mach
    relation dA/A = (M^2 - 1) dV/V), so the throat flux is the sonic flux.
    """
    return sonic_state(gas, comp, Tt, Pt).mass_flux


def state_from_mass_flux(gas, comp, Tt, Pt, G, supersonic=False) -> StaticState:
    """Static state carrying mass flux G = mdot/A on the chosen branch.

    Raises ValueError if G exceeds the choked flux (the duct would choke).
    """
    Ts = sonic_T(gas, comp, Tt)
    Gs = static_from_T(gas, comp, Tt, Pt, Ts).mass_flux
    if G > Gs * (1 + 1e-12):
        raise ValueError(f"mass flux {G:.4g} exceeds choked flux {Gs:.4g}")
    f = lambda T: static_from_T(gas, comp, Tt, Pt, T).mass_flux - G
    if supersonic:
        lo, hi = 0.05 * Tt, Ts
    else:
        lo, hi = Ts, Tt * (1 - 1e-12)
    if abs(G - Gs) <= 1e-12 * Gs:
        return static_from_T(gas, comp, Tt, Pt, Ts)
    return static_from_T(gas, comp, Tt, Pt, brentq(f, lo, hi, xtol=1e-12 * Tt))


def total_from_static(gas, comp, T, P, V):
    """Stagnation state of a moving static state (inverse of static_from_T)."""
    ht = gas.h(T, comp) + 0.5 * V * V
    Tt = gas.T_from_h(ht, comp, T_guess=T)
    Pt = P * math.exp((gas.phi(Tt, comp) - gas.phi(T, comp)) / gas.R(comp))
    return Tt, Pt
