"""Convergent-divergent exhaust nozzle with variable exit-to-throat area
ratio  eps = A9/A8, and the sweep that finds the thrust-maximising eps.

Physics summary (what to say in an interview)
--------------------------------------------
* Upstream of the throat nothing depends on A9 once the throat is choked:
  the throat passes  mdot = A8 * (rho V)*  regardless of what happens
  downstream, so the engine's gas generator is isolated from eps.
* Gross thrust from a control volume around the nozzle:
      Fg = mdot9 V9 + (P9 - P0) A9
  Differentiate along the isentropic expansion (mdot fixed):
      d(mdot V) = -A dP      (1-D momentum: rho V dV = -dP times A)
      dFg/dA9 = -A9 dP9/dA9 + (P9 - P0) + A9 dP9/dA9 = (P9 - P0)
  so Fg is maximised where P9 = P0 (ideal expansion).  Under-expanded
  (P9 > P0): opening A9 still adds thrust.  Over-expanded (P9 < P0): the
  divergent walls feel net suction and thrust drops.
* Real-nozzle effects that move the optimum off P9 = P0:
  - divergence (angularity) loss: the exit flow is conical, half-angle
    theta, so only lambda = (1 + cos theta)/2 of the momentum is axial
    (Malina's 1-D conical correction).  For a fixed-length flap, theta
    grows with eps, so the optimum moves slightly toward under-expansion.
  - total-pressure loss pi_n (friction) - applied upstream of the throat.
  - flow separation in strong over-expansion (Summerfield criterion): the
    boundary layer separates where wall pressure drops below ~0.37-0.4 P0,
    and the separated wall then sits near P0.  We model it by truncating
    the nozzle at the separation area ratio.
* If the throat is not choked (low nozzle pressure ratio) the nozzle is a
  subsonic venturi: exit pressure equals P0 and the exit area sets the flow.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from . import flow


@dataclass
class NozzleModel:
    pi_n: float = 0.97                  # total-pressure ratio Pt8/Pt6A
    flap_length_ratio: float | None = 3.0   # divergent flap length / throat radius; None = no angularity loss
    separation_pressure_ratio: float | None = 0.37  # P_wall/P0 at separation; None = never separates


@dataclass
class NozzleResult:
    area_ratio: float
    effective_area_ratio: float
    choked: bool
    separated: bool
    exit: flow.StaticState
    A8: float
    A9: float
    divergence_factor: float
    gross_thrust: float


def throat_capacity(gas, comp, Tt8, Pt8, P0, A8, area_ratio):
    """Mass flow the nozzle can swallow at these conditions.

    Subsonic everywhere when back pressure P0 exceeds the 'first critical'
    pressure for this eps (isentropic subsonic exit with sonic throat);
    then the exit, at P9 = P0, limits the flow.  Otherwise the throat is
    choked and mdot = A8 (rho V)*.
    """
    G_star = flow.choked_mass_flux(gas, comp, Tt8, Pt8)
    exit_sub = flow.state_from_mass_flux(gas, comp, Tt8, Pt8, G_star / area_ratio, supersonic=False)
    if P0 >= exit_sub.P:
        ex = flow.state_from_pressure(gas, comp, Tt8, Pt8, P0)
        return A8 * area_ratio * ex.mass_flux, False
    return A8 * G_star, True


def ideal_expansion_area_ratio(gas, comp, Tt8, Pt8, P0) -> float:
    """eps giving P9 = P0 exactly (the textbook optimum)."""
    G_star = flow.choked_mass_flux(gas, comp, Tt8, Pt8)
    s = flow.sonic_state(gas, comp, Tt8, Pt8)
    if s.P <= P0:
        return 1.0
    return G_star / flow.state_from_pressure(gas, comp, Tt8, Pt8, P0).mass_flux


def nozzle(gas, comp, Tt8, Pt8, mdot, A8, area_ratio, P0, model: NozzleModel) -> NozzleResult:
    """Exit state and gross thrust for throat area A8 and A9/A8 = area_ratio."""
    G_star = flow.choked_mass_flux(gas, comp, Tt8, Pt8)
    capacity, choked = throat_capacity(gas, comp, Tt8, Pt8, P0, A8, area_ratio)
    if not choked:
        ex = flow.state_from_pressure(gas, comp, Tt8, Pt8, P0)
        A9 = mdot / ex.mass_flux
        return NozzleResult(area_ratio, A9 / A8, False, False, ex, A8, A9, 1.0, mdot * ex.V)

    eps_eff, separated = area_ratio, False
    ex = flow.state_from_mass_flux(gas, comp, Tt8, Pt8, G_star / area_ratio, supersonic=area_ratio > 1.0)
    psep = model.separation_pressure_ratio
    if psep is not None and ex.P < psep * P0:
        # Separation: the jet leaves the wall where P = psep*P0; beyond that
        # the wall sees ~ambient pressure (no thrust, no drag).
        separated = True
        sonic = flow.sonic_state(gas, comp, Tt8, Pt8)
        if sonic.P <= psep * P0:
            eps_eff = 1.0
            ex = sonic
        else:
            ex = flow.state_from_pressure(gas, comp, Tt8, Pt8, psep * P0)
            eps_eff = G_star / ex.mass_flux
    A9_eff = A8 * eps_eff

    lam = 1.0
    if model.flap_length_ratio:
        # Exit radius / throat radius = sqrt(eps); straight flap of length L
        theta = math.atan((math.sqrt(eps_eff) - 1.0) / model.flap_length_ratio)
        lam = 0.5 * (1.0 + math.cos(theta))
    # Fg = lambda * mdot V9 + (P9 - P0) A9  (pressure term acts on the
    # projected exit area, so it is not reduced by angularity)
    Fg = lam * mdot * ex.V + (ex.P - P0) * A9_eff
    return NozzleResult(area_ratio, eps_eff, True, separated, ex, A8, A8 * area_ratio, lam, Fg)


@dataclass
class AreaSweep:
    area_ratios: np.ndarray
    gross_thrust: np.ndarray
    best_area_ratio: float
    best_gross_thrust: float
    ideal_area_ratio: float


def optimise_area_ratio(gas, comp, Tt8, Pt8, mdot, A8, P0, model: NozzleModel,
                        eps_min=1.0, eps_max=3.0, n=81) -> AreaSweep:
    """Brute-force sweep of eps followed by a bounded scalar refinement.

    The sweep is kept (and returned) because the thrust-vs-eps curve itself
    is worth plotting: it is flat near the optimum, which is why real
    nozzles tolerate schedule errors of +/-10 % in eps with < 0.5 % thrust loss.
    """
    eps = np.linspace(eps_min, eps_max, n)
    Fg = np.array([nozzle(gas, comp, Tt8, Pt8, mdot, A8, e, P0, model).gross_thrust for e in eps])
    i = int(np.argmax(Fg))
    lo, hi = eps[max(i - 1, 0)], eps[min(i + 1, n - 1)]
    if hi > lo:
        res = minimize_scalar(lambda e: -nozzle(gas, comp, Tt8, Pt8, mdot, A8, e, P0, model).gross_thrust,
                              bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
        best, Fbest = float(res.x), float(-res.fun)
        if Fbest < Fg[i]:
            best, Fbest = float(eps[i]), float(Fg[i])
    else:
        best, Fbest = float(eps[i]), float(Fg[i])
    return AreaSweep(eps, Fg, best, Fbest, ideal_expansion_area_ratio(gas, comp, Tt8, Pt8, P0))
