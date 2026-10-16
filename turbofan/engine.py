"""Two-spool, mixed-flow turbofan: design point and off-design operation.

Architecture
------------
    LP spool:  fan (both streams, pressure ratio pi_f)  <-  LP turbine
    HP spool:  HP compressor (core, pi_cH)               <-  HP turbine
    overall core pressure ratio  pi_c = pi_f * pi_cH
    bypass and core streams meet in a constant-area mixer, then a single
    variable-geometry convergent-divergent nozzle.

Design point ("parametric", Mattingly ch. 7)
--------------------------------------------
Given pi_c, bypass ratio alpha, Tt4 and component efficiencies, the fan
pressure ratio is NOT free: for a mixer to work the two streams must enter
at the same static pressure, and Mattingly's design rule is to also match
total pressure, Pt16 = Pt6.  The fan pressure ratio is found by root-finding
pi_f until that holds (this is the "system of equations of Problem 7.27").
The design then fixes the flow areas A4, A45 (turbine nozzle throats),
A6, A16 (mixer inlets, from the design mixer Mach number M6) and A8.

Off-design (Mattingly ch. 8 assumptions, but with cp(T))
--------------------------------------------------------
* component polytropic efficiencies stay at their design values;
* HPT and LPT nozzle-guide-vane throats are choked (A4, A45 fixed), so the
  HPT pressure ratio is fixed by A4/A45 (see components.choked_turbine_pressure_ratio);
* mixer inlet areas A6, A16 are fixed and the streams meet at equal static
  pressure, P6 = P16 (Kutta condition at the splitter) - this sets alpha;
* the nozzle throat A8 passes the mixed flow - this sets the fan operating point.
Unknowns (pi_f, alpha); residuals (P16 - P6, nozzle flow mismatch).
Everything else follows station by station.

Throttle schedule: Mattingly's throttle ratio TR = theta_0 at the design
point.  For theta_0 = Tt0/288.15 K <= TR the controller holds Tt4/Tt0
constant (so corrected speeds - and therefore all pressure ratios - stay at
design); above TR it holds Tt4 = Tt4,max and the compressors unload.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np
from scipy.optimize import brentq, root

from . import components as C
from . import flow
from .atmosphere import Ambient, isa, T_SL
from .gas import ThermallyPerfectGas
from .nozzle import NozzleModel, nozzle, optimise_area_ratio, ideal_expansion_area_ratio, throat_capacity


@dataclass
class DesignInputs:
    """Everything the cycle needs.  Defaults: Mattingly Table 6.2 level-4
    technology (cooled turbine, type C diffuser, type F nozzle) for a
    medium-bypass supercruise engine."""

    M0: float = 1.6
    altitude: float = 12192.0          # m (40 kft)
    mdot0: float = 100.0               # kg/s total airflow at design
    bypass_ratio: float = 0.8
    pi_c: float = 25.0                 # overall (fan x HPC) pressure ratio
    Tt4: float = 1805.0                # K  (3250 degR)
    h_PR: float = 42.8e6               # J/kg, Jet-A lower heating value
    pi_d_max: float = 0.96             # diffuser type C, level 4
    e_f: float = 0.89                  # fan polytropic efficiency
    e_cH: float = 0.90                 # HP compressor polytropic efficiency
    e_tH: float = 0.89                 # HP turbine (cooled) polytropic efficiency
    e_tL: float = 0.89                 # LP turbine polytropic efficiency
    pi_b: float = 0.96                 # burner total-pressure ratio
    eta_b: float = 0.999               # combustion efficiency
    eta_mH: float = 0.99               # HP spool mechanical efficiency
    eta_mL: float = 0.99               # LP spool mechanical efficiency
    pi_bypass_duct: float = 1.0        # Pt16/Pt13
    pi_M_max: float = 0.96             # mixer wall-friction loss
    M6: float = 0.4                    # mixer core-inlet Mach number at design
    nozzle: NozzleModel = field(default_factory=lambda: NozzleModel(pi_n=0.97))
    design_area_ratio: float | None = None   # None = ideally expanded at design


@dataclass
class CycleResult:
    """One converged operating point, station by station."""

    inputs: DesignInputs
    M0: float
    ambient: Ambient
    V0: float
    stations: dict
    f: float                 # burner fuel/air ratio (per kg core air)
    alpha: float
    pi_f: float
    pi_cH: float
    pi_tH: float
    pi_tL: float
    mdot0: float
    nozzle_result: object
    eta_r: float

    @property
    def fuel_flow(self) -> float:
        return self.f * self.mdot0 / (1.0 + self.alpha)

    @property
    def thrust(self) -> float:
        # Uninstalled net thrust: gross thrust minus ram drag mdot0 V0
        return self.nozzle_result.gross_thrust - self.mdot0 * self.V0

    @property
    def specific_thrust(self) -> float:
        return self.thrust / self.mdot0

    @property
    def tsfc(self) -> float:
        """kg/(N s)"""
        return self.fuel_flow / self.thrust

    @property
    def pi_c(self) -> float:
        return self.pi_f * self.pi_cH

    @property
    def Tt4(self) -> float:
        return self.stations["4"].Tt

    def thermal_efficiency(self) -> float:
        """Change in jet kinetic energy / fuel chemical power."""
        V9 = self.nozzle_result.exit.V
        m9 = self.mdot0 + self.fuel_flow
        dKE = 0.5 * (m9 * V9 ** 2 - self.mdot0 * self.V0 ** 2)
        return dKE / (self.fuel_flow * self.inputs.h_PR)

    def propulsive_efficiency(self) -> float:
        V9 = self.nozzle_result.exit.V
        m9 = self.mdot0 + self.fuel_flow
        dKE = 0.5 * (m9 * V9 ** 2 - self.mdot0 * self.V0 ** 2)
        return self.thrust * self.V0 / dKE if dKE > 0 else float("nan")

    def summary(self) -> str:
        lines = [f"M0={self.M0:.3f}  alt={self.ambient.altitude:.0f} m  T0={self.ambient.T:.2f} K  P0={self.ambient.P:.0f} Pa",
                 f"{'stn':>4} {'Tt [K]':>9} {'Pt [kPa]':>10} {'mdot [kg/s]':>12}"]
        for k, s in self.stations.items():
            lines.append(f"{k:>4} {s.Tt:9.2f} {s.Pt / 1e3:10.3f} {s.mdot:12.3f}")
        nr = self.nozzle_result
        lines += [f"pi_f={self.pi_f:.4f} pi_cH={self.pi_cH:.4f} pi_c={self.pi_c:.3f} alpha={self.alpha:.4f}",
                  f"pi_tH={self.pi_tH:.4f} pi_tL={self.pi_tL:.4f} f={self.f:.5f}",
                  f"A9/A8={nr.area_ratio:.4f}  P9/P0={nr.exit.P / self.ambient.P:.4f}  M9={nr.exit.M:.3f}  V9={nr.exit.V:.1f} m/s",
                  f"F={self.thrust / 1e3:.3f} kN  F/mdot0={self.specific_thrust:.2f} N/(kg/s)  TSFC={self.tsfc * 1e6:.3f} mg/(N s)"]
        return "\n".join(lines)


# ======================================================================
def _gas_path(gas, inp: DesignInputs, amb: Ambient, M0, Tt4, pi_f, alpha, mdot_core,
              pi_cH=None, A4=None, A45=None, pi_tH=None):
    """Stations 0 -> 5 and 13/16 for given fan ratio, bypass ratio, core flow.

    Design mode passes pi_cH; off-design mode passes the HPT throat areas
    and gets pi_cH from the HP-spool power balance instead.
    Returns (stations, f, pi_cH, pi_tH, pi_tL, V0, eta_r, mdot_core).
    """
    st0, V0, a0 = C.free_stream(gas, amb, M0)
    st2 = C.inlet(st0, M0, inp.pi_d_max)
    eta_r = C.ram_recovery(M0)

    # Fan: one pressure ratio for both streams (simplification; real fans
    # have a hub/tip split).  Fan work per kg of total air:
    st13 = C.compress(gas, st2, pi_f, inp.e_f, "13")
    dh_fan = gas.h(st13.Tt, st13.comp) - gas.h(st2.Tt, st2.comp)
    st25 = replace(st13, name="25")

    if pi_cH is not None:                       # ----- design mode
        st3 = C.compress(gas, st25, pi_cH, inp.e_cH, "3")
        st4, f = C.burner(gas, st3, Tt4, inp.pi_b, inp.eta_b, inp.h_PR)
        # HP spool power balance (per kg core air):
        #   h3 - h25 = eta_mH (1 + f) (h4 - h45)
        dh_c = gas.h(st3.Tt, st3.comp) - gas.h(st25.Tt, st25.comp)
        st45, pi_tH = C.expand_by_work(gas, st4, dh_c / (inp.eta_mH * (1 + f)), inp.e_tH, "45")
    else:                                       # ----- off-design mode
        # HPT between two choked throats -> pi_tH set by A4/A45 (and Tt4).
        # The HPC work then follows from HPT work; f depends on Tt3 so iterate.
        f = 0.025
        for _ in range(40):
            comp4 = gas.burn(gas.air(), f)
            s4 = C.Station("4", Tt4, 1.0, comp4)
            pi_tH = C.choked_turbine_pressure_ratio(gas, s4, A4, A45, inp.e_tH)
            s45 = C.expand_by_pressure_ratio(gas, s4, pi_tH, inp.e_tH, "45")
            dh_t = gas.h(Tt4, comp4) - gas.h(s45.Tt, comp4)
            st3, pi_cH = C.compress_by_work(gas, st25, inp.eta_mH * (1 + f) * dh_t, inp.e_cH, "3")
            st4, f_new = C.burner(gas, st3, Tt4, inp.pi_b, inp.eta_b, inp.h_PR)
            if abs(f_new - f) < 1e-11:
                break
            f = f_new
        f = f_new
        st45 = C.expand_by_pressure_ratio(gas, st4, pi_tH, inp.e_tH, "45")
        # Choked HPT nozzle: mdot4 = A4 (rho V)*(Tt4, Pt4)  -> core airflow
        mdot_core = A4 * flow.choked_mass_flux(gas, st4.comp, st4.Tt, st4.Pt) / (1 + f)

    # LP spool power balance (per kg core air):
    #   (1 + alpha) (h13 - h2) = eta_mL (1 + f) (h45 - h5)
    st5, pi_tL = C.expand_by_work(gas, st45, (1 + alpha) * dh_fan / (inp.eta_mL * (1 + f)), inp.e_tL, "5")
    st16 = C.Station("16", st13.Tt, st13.Pt * inp.pi_bypass_duct, st13.comp)

    # mass flows
    for s, m in ((st0, (1 + alpha) * mdot_core), (st2, (1 + alpha) * mdot_core),
                 (st13, alpha * mdot_core), (st16, alpha * mdot_core), (st25, mdot_core), (st3, mdot_core),
                 (st4, (1 + f) * mdot_core), (st45, (1 + f) * mdot_core), (st5, (1 + f) * mdot_core)):
        s.mdot = m
    stations = {"0": st0, "2": st2, "13": st13, "16": st16, "25": st25, "3": st3, "4": st4, "45": st45, "5": st5}
    return stations, f, pi_cH, pi_tH, pi_tL, V0, eta_r, mdot_core


def _exhaust(gas, inp, amb, stations, A6, A16, A8, area_ratio, model):
    st5, st16 = stations["5"], stations["16"]
    st6 = replace(st5, name="6", area=A6)
    st6.static = flow.state_from_mass_flux(gas, st6.comp, st6.Tt, st6.Pt, st6.mdot / A6)
    st16.static = flow.state_from_mass_flux(gas, st16.comp, st16.Tt, st16.Pt, st16.mdot / A16)
    st16.area = A16
    st6A = C.mixed_out(gas, st6, st16, A6, A16, inp.pi_M_max)
    st8 = C.Station("8", st6A.Tt, st6A.Pt * model.pi_n, st6A.comp, st6A.mdot, area=A8)
    nr = nozzle(gas, st8.comp, st8.Tt, st8.Pt, st8.mdot, A8, area_ratio, amb.P, model)
    st8.static = flow.sonic_state(gas, st8.comp, st8.Tt, st8.Pt) if nr.choked else None
    st9 = C.Station("9", st8.Tt, st8.Pt, st8.comp, st8.mdot, static=nr.exit, area=nr.A9)
    stations.update({"6": st6, "6A": st6A, "8": st8, "9": st9})
    return nr


@dataclass
class DesignPoint:
    """A sized engine: the design cycle plus its fixed flow areas."""

    inputs: DesignInputs
    gas: object
    cycle: CycleResult
    A4: float
    A45: float
    A6: float
    A16: float
    A8: float
    throttle_ratio: float

    def engine(self) -> "Engine":
        return Engine(self)


def design(inp: DesignInputs = None, gas=None, pi_f_bounds=(1.05, 8.0)) -> DesignPoint:
    """Design-point (parametric) cycle with Pt16 = Pt6 at the mixer."""
    inp = inp or DesignInputs()
    gas = gas or ThermallyPerfectGas()
    amb = isa(inp.altitude)
    alpha = inp.bypass_ratio
    mdot_core = inp.mdot0 / (1 + alpha)

    def path(pi_f):
        return _gas_path(gas, inp, amb, inp.M0, inp.Tt4, pi_f, alpha, mdot_core, pi_cH=inp.pi_c / pi_f)

    def mismatch(pi_f):
        st = path(pi_f)[0]
        return st["16"].Pt / st["5"].Pt - 1.0

    # Raising pi_f raises Pt16 but takes more LPT work, lowering Pt6: the
    # mismatch is monotone so a bracketing root-finder is safe.
    hi = min(pi_f_bounds[1], inp.pi_c * 0.999)
    pi_f = brentq(mismatch, pi_f_bounds[0], hi, xtol=1e-12)
    stations, f, pi_cH, pi_tH, pi_tL, V0, eta_r, _ = path(pi_f)

    # ---- size the flow areas
    st4, st45, st5, st16 = stations["4"], stations["45"], stations["5"], stations["16"]
    A4 = st4.mdot / flow.choked_mass_flux(gas, st4.comp, st4.Tt, st4.Pt)
    A45 = st45.mdot / flow.choked_mass_flux(gas, st45.comp, st45.Tt, st45.Pt)
    s6 = flow.state_from_mach(gas, st5.comp, st5.Tt, st5.Pt, inp.M6)
    A6 = st5.mdot / s6.mass_flux
    # equal static pressure at the splitter -> bypass Mach number M16
    s16 = flow.state_from_pressure(gas, st16.comp, st16.Tt, st16.Pt, s6.P)
    A16 = st16.mdot / s16.mass_flux

    # mixed-out state is needed to size the throat
    st6 = replace(st5, name="6", static=s6, area=A6)
    st16.static, st16.area = s16, A16
    st6A = C.mixed_out(gas, st6, st16, A6, A16, inp.pi_M_max)
    Pt8 = st6A.Pt * inp.nozzle.pi_n
    A8 = st6A.mdot / flow.choked_mass_flux(gas, st6A.comp, st6A.Tt, Pt8)
    eps = inp.design_area_ratio or ideal_expansion_area_ratio(gas, st6A.comp, st6A.Tt, Pt8, amb.P)
    nr = _exhaust(gas, inp, amb, stations, A6, A16, A8, eps, inp.nozzle)

    cyc = CycleResult(inp, inp.M0, amb, V0, stations, f, alpha, pi_f, pi_cH, pi_tH, pi_tL,
                      inp.mdot0, nr, eta_r)
    TR = stations["0"].Tt / T_SL
    return DesignPoint(inp, gas, cyc, A4, A45, A6, A16, A8, TR)


# ======================================================================
class Engine:
    """Fixed-geometry engine for off-design operation (variable A9 only,
    plus optionally a scaled A8)."""

    def __init__(self, dp: DesignPoint):
        self.dp = dp
        self.inp = dp.inputs
        self.gas = dp.gas
        self._last = (math.log(dp.cycle.pi_f), math.log(dp.cycle.alpha))

    def scheduled_Tt4(self, Tt0: float) -> float:
        """Throttle-ratio control law (Mattingly ch. 8, 'theta break')."""
        theta0 = Tt0 / T_SL
        return self.inp.Tt4 * min(1.0, theta0 / self.dp.throttle_ratio)

    def _residuals(self, x, amb, M0, Tt4, A8, eps):
        pi_f, alpha = math.exp(x[0]), math.exp(x[1])
        g = self.gas
        try:
            st, f, pi_cH, pi_tH, pi_tL, V0, eta_r, mdot_core = _gas_path(
                g, self.inp, amb, M0, Tt4, pi_f, alpha, None, A4=self.dp.A4, A45=self.dp.A45)
        except (ValueError, RuntimeError):
            return np.array([1.0, 1.0]), None
        st5, st16 = st["5"], st["16"]

        def p_static(s, A):
            # static pressure for this flow through area A; if it would choke,
            # extend continuously below P* so the solver is pushed back
            G = s.mdot / A
            Gs = flow.choked_mass_flux(g, s.comp, s.Tt, s.Pt)
            if G <= Gs:
                return flow.state_from_mass_flux(g, s.comp, s.Tt, s.Pt, G).P
            return flow.sonic_state(g, s.comp, s.Tt, s.Pt).P * (1 - (G / Gs - 1))

        r1 = (p_static(st16, self.dp.A16) - p_static(st5, self.dp.A6)) / st5.Pt
        if not (r1 == r1):
            return np.array([1.0, 1.0]), None
        try:
            st6 = replace(st5, name="6")
            st6.static = flow.state_from_mass_flux(g, st5.comp, st5.Tt, st5.Pt, st5.mdot / self.dp.A6)
            st16.static = flow.state_from_mass_flux(g, st16.comp, st16.Tt, st16.Pt, st16.mdot / self.dp.A16)
            st6A = C.mixed_out(g, st6, st16, self.dp.A6, self.dp.A16, self.inp.pi_M_max)
        except ValueError:
            return np.array([r1, 1.0]), None
        Pt8 = st6A.Pt * self.inp.nozzle.pi_n
        cap, _ = throat_capacity(g, st6A.comp, st6A.Tt, Pt8, amb.P, A8, eps)
        r2 = st6A.mdot / cap - 1.0
        return np.array([r1, r2]), (st, f, pi_cH, pi_tH, pi_tL, V0, eta_r, mdot_core)

    def operate(self, M0: float, altitude: float, Tt4: float | None = None,
                area_ratio: float | None = None, A8_scale: float = 1.0, guess=None) -> CycleResult:
        """Converge the engine at a flight condition.

        area_ratio None -> use the thrust-optimal A9/A8 (found by sweep).
        """
        amb = isa(altitude)
        st0, _, _ = C.free_stream(self.gas, amb, M0)
        Tt4 = Tt4 or self.scheduled_Tt4(st0.Tt)
        A8 = self.dp.A8 * A8_scale
        eps_for_match = area_ratio or 2.0  # flow capacity is independent of eps once choked
        x0 = np.array(guess if guess is not None else self._last)
        sol = root(lambda x: self._residuals(x, amb, M0, Tt4, A8, eps_for_match)[0], x0,
                   method="hybr", options={"xtol": 1e-12})
        r, data = self._residuals(sol.x, amb, M0, Tt4, A8, eps_for_match)
        if data is None or np.max(np.abs(r)) > 1e-7:
            raise RuntimeError(f"off-design did not converge at M0={M0}, alt={altitude}: {r}")
        self._last = tuple(sol.x)
        st, f, pi_cH, pi_tH, pi_tL, V0, eta_r, mdot_core = data
        pi_f, alpha = math.exp(sol.x[0]), math.exp(sol.x[1])

        model = self.inp.nozzle
        if area_ratio is None:
            # quick mixed-out state to optimise eps; _exhaust rebuilds it below
            tmp = dict(st)
            _exhaust(self.gas, self.inp, amb, tmp, self.dp.A6, self.dp.A16, A8, 1.5, model)
            s8 = tmp["8"]
            sweep = optimise_area_ratio(self.gas, s8.comp, s8.Tt, s8.Pt, s8.mdot, A8, amb.P, model)
            area_ratio = sweep.best_area_ratio
        nr = _exhaust(self.gas, self.inp, amb, st, self.dp.A6, self.dp.A16, A8, area_ratio, model)
        return CycleResult(self.inp, M0, amb, V0, st, f, alpha, pi_f, pi_cH, pi_tH, pi_tL,
                           (1 + alpha) * mdot_core, nr, eta_r)

    def area_sweep(self, M0, altitude, Tt4=None, **kw):
        """Converged gas generator + full A9/A8 sweep at one flight condition."""
        op = self.operate(M0, altitude, Tt4=Tt4, area_ratio=1.5)
        s8 = op.stations["8"]
        sweep = optimise_area_ratio(self.gas, s8.comp, s8.Tt, s8.Pt, s8.mdot, s8.area, op.ambient.P,
                                    self.inp.nozzle, **kw)
        return op, sweep
