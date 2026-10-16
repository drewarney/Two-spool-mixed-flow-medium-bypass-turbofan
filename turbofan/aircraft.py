"""Aircraft-side requirements: thrust needed at cruise and the minimum
uninstalled specific thrust an engine must deliver (design problems
7.D1 / 7.D2, 'Calculate minimum specific thrust at cruise').

The logic:
  1. Required installed thrust from the energy equation (Mattingly eq. 1.28):
        P_s = V (T - D) / W    =>    T = D + W P_s / V
     with drag from the parabolic polar
        D = q S C_D,   C_D = K1 C_L^2 + K2 C_L + C_D0,   C_L = W / (q S),
        q = 0.5 rho V^2 = 0.5 gamma P M^2.
  2. Installation losses phi_inlet + phi_nozzle are fractions of the
     uninstalled thrust F:  T = F (1 - phi)   =>   F = T / (1 - phi).
  3. The inlet can swallow at most the free-stream tube the size of its
     capture area:  mdot_max = rho0 A V0.
  4. Engine must give F/mdot0 >= F_required / mdot_max; a smaller inlet
     demands a higher specific thrust (and a hotter, less efficient jet).

AIRCRAFT DATA WARNING: the default weights / drag polars below are
placeholders of the right order of magnitude so the scripts run.  Replace
them with the HP-1 and HF-1 values from your Chapter 1 / Chapter 2 design
problems before quoting any result.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .atmosphere import Ambient


@dataclass
class Aircraft:
    name: str
    weight: float            # N, at the flight condition
    wing_area: float         # m^2
    CD0: float
    K1: float
    K2: float = 0.0
    n_engines: int = 2
    install_loss: float = 0.02   # phi_inlet + phi_nozzle
    tsfc_max: float | None = None  # kg/(N s), S_max from Chapter 1


# Placeholder data -- see module docstring.
HP1_PLACEHOLDER = Aircraft("HP-1 (placeholder data)", weight=0.95 * 1.6e6, wing_area=343.0,
                           CD0=0.016, K1=0.056, K2=-0.008, install_loss=0.02, tsfc_max=17.0e-6)
HF1_PLACEHOLDER = Aircraft("HF-1 (placeholder data)", weight=0.92 * 40000 * 4.4482216, wing_area=375 * 0.09290304,
                           CD0=0.028, K1=0.20, K2=0.0, install_loss=0.05, tsfc_max=31.0e-6)


def dynamic_pressure(amb: Ambient, M0: float, gamma=1.4) -> float:
    return 0.5 * gamma * amb.P * M0 ** 2


def drag(ac: Aircraft, amb: Ambient, M0: float, load_factor=1.0) -> float:
    q = dynamic_pressure(amb, M0)
    CL = load_factor * ac.weight / (q * ac.wing_area)
    CD = ac.K1 * CL ** 2 + ac.K2 * CL + ac.CD0
    return q * ac.wing_area * CD


def required_uninstalled_thrust_per_engine(ac: Aircraft, amb: Ambient, M0: float, Ps: float = 0.0) -> float:
    V = M0 * amb.a
    T_installed = drag(ac, amb, M0) + ac.weight * Ps / V
    return T_installed / (1.0 - ac.install_loss) / ac.n_engines


def capture_mass_flow(amb: Ambient, M0: float, area: float) -> float:
    return amb.rho * area * M0 * amb.a


def inlet_area_from_diameter(d: float) -> float:
    return math.pi * d * d / 4.0


def minimum_specific_thrust(F_required: float, mdot_max: float) -> float:
    return F_required / mdot_max
