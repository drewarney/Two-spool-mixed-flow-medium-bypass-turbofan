"""Working-fluid models.

Two interchangeable gas models share one small interface, so the *same*
station-by-station cycle code can run either with temperature-dependent
properties (the "real" model) or with Mattingly's two-constant-gamma
idealisation (used to cross-check against the textbook closed-form
equations).

Interface (all SI, per unit mass):

    comp = gas.air()                      composition of fresh air
    comp = gas.burn(comp_air, far)        composition after the burner
    comp = gas.mix([(m1, c1), (m2, c2)])  composition of a mixed stream
    gas.cp(T, comp), gas.R(comp), gas.gamma(T, comp)
    gas.h(T, comp)        static or total enthalpy  [J/kg]
    gas.T_from_h(h, comp)
    gas.phi(T, comp)      entropy function  phi(T) = integral cp(T)/T dT [J/kg/K]
    gas.T_from_phi(phi, comp)

Why the entropy function phi?
-----------------------------
For a thermally perfect gas (p = rho R T, but cp = cp(T)) the Gibbs
equation  T ds = dh - dp/rho  gives

    ds = cp(T) dT / T - R dp / p
    s2 - s1 = [phi(T2) - phi(T1)] - R ln(p2/p1)

so an *isentropic* process obeys  phi(T2) - phi(T1) = R ln(p2/p1).
That one line replaces the constant-cp relation  T2/T1 = (p2/p1)^((g-1)/g)
everywhere in the code, and reduces to it exactly when cp is constant
(phi = cp ln T).  Polytropic efficiency e is defined on the differential
process, so for compression  dh = (1/e) dh_s  =>  cp dT/T = (R/e) dp/p  =>

    phi(T2) - phi(T1) = (R / e_c) ln(pi)          (compressor / fan)
    phi(T2) - phi(T1) = (R * e_t) ln(pi)          (turbine, pi < 1)

which again collapses to tau = pi^((g-1)/(g e)) for constant cp.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

# Reference temperature for enthalpy.  The fuel heating value h_PR is
# defined with reactants and products at 298.15 K, so enthalpies are
# measured from there; the burner energy balance is then simply
#   (1 + f) h4 - h3 = f * eta_b * h_PR    (fuel enters at T_ref, h_fuel = 0)
T_REF = 298.15


@dataclass(frozen=True)
class AirFuel:
    """Composition for the thermally perfect model: fuel/air mass ratio."""

    far: float

    @property
    def fuel_mass_fraction(self) -> float:
        # mass of fuel-derived products per unit mass of mixture
        return self.far / (1.0 + self.far)


class ThermallyPerfectGas:
    """Air and kerosene combustion products with cp = cp(T, FAR).

    The polynomials are the widely used fits of Walsh & Fletcher,
    *Gas Turbine Performance* (2nd ed., ch. 3) for dry air plus a
    combustion-products correction that is linear in fuel mass fraction:

        cp(T, FAR) = A(Tz) + FAR/(1+FAR) * B(Tz),      Tz = T / 1000 K

    Valid roughly 200-2000 K.  Because the correction is linear in the
    fuel mass fraction, mixing two streams by mass-weighting their
    enthalpies is exactly consistent with re-evaluating h at the mixed FAR.

    The gas constant of kerosene (CnH2n) products differs from air by
    well under 1 % at any realistic FAR (CH2 + 1.5 O2 -> CO2 + H2O adds
    0.5 mol per 1.5 mol O2 consumed but also adds the fuel mass), so R is
    held at the air value.
    """

    name = "thermally-perfect (Walsh & Fletcher)"

    # cp of dry air, kJ/(kg K), polynomial in Tz = T/1000
    A = (0.992313, 0.236688, -1.852148, 6.083152, -8.893933,
         7.097112, -3.234725, 0.794571, -0.081873)
    # combustion-products increment, multiplied by FAR/(1+FAR)
    B = (-0.718874, 8.747481, -15.863157, 17.254096, -10.233795,
         3.081778, -0.361112, -0.003919)
    R_AIR = 287.05  # J/(kg K)

    def __init__(self):
        # Pre-build numpy polynomials (coefficients in kJ/kg/K -> J/kg/K).
        self._cpA = np.polynomial.Polynomial(np.array(self.A) * 1000.0)
        self._cpB = np.polynomial.Polynomial(np.array(self.B) * 1000.0)
        # h = integral cp dT = 1000 * integral cp(Tz) dTz  (dT = 1000 dTz)
        self._hA = self._cpA.integ() * 1000.0
        self._hB = self._cpB.integ() * 1000.0
        # phi = integral cp/T dT = a0 ln T + integral (cp - a0)/Tz dTz
        self._phiA = np.polynomial.Polynomial(self._cpA.coef[1:]).integ()
        self._phiB = np.polynomial.Polynomial(self._cpB.coef[1:]).integ()
        zr = T_REF / 1000.0
        self._hA0, self._hB0 = self._hA(zr), self._hB(zr)

    # ---- composition handling -------------------------------------------
    def air(self) -> AirFuel:
        return AirFuel(0.0)

    def burn(self, comp: AirFuel, far: float) -> AirFuel:
        # far is fuel added per unit mass of the *air* in the stream
        return AirFuel(comp.far + far)

    def mix(self, streams) -> AirFuel:
        # Total fuel / total air, both obtained from each stream's FAR.
        air = sum(m / (1.0 + c.far) for m, c in streams)
        fuel = sum(m * c.far / (1.0 + c.far) for m, c in streams)
        return AirFuel(fuel / air)

    # ---- properties --------------------------------------------------------
    def R(self, comp) -> float:
        return self.R_AIR

    def cp(self, T, comp) -> float:
        z = T / 1000.0
        return self._cpA(z) + comp.fuel_mass_fraction * self._cpB(z)

    def gamma(self, T, comp) -> float:
        cp = self.cp(T, comp)
        return cp / (cp - self.R(comp))

    def h(self, T, comp) -> float:
        z = T / 1000.0
        y = comp.fuel_mass_fraction
        return (self._hA(z) - self._hA0) + y * (self._hB(z) - self._hB0)

    def phi(self, T, comp) -> float:
        z = T / 1000.0
        y = comp.fuel_mass_fraction
        a0 = self._cpA.coef[0] + y * self._cpB.coef[0]
        return a0 * math.log(T) + self._phiA(z) + y * self._phiB(z)

    # ---- inversions (Newton; cp > 0 so both functions are monotone) --------
    def T_from_h(self, h, comp, T_guess=800.0) -> float:
        T = T_guess
        for _ in range(50):
            dT = (self.h(T, comp) - h) / self.cp(T, comp)
            dT = max(min(dT, 0.5 * T), -T)   # damp: never step to T <= 0
            T -= dT
            if abs(dT) < 1e-9 * T:
                return T
        raise RuntimeError("T_from_h did not converge")

    def T_from_phi(self, phi, comp, T_guess=800.0) -> float:
        lnT = math.log(T_guess)
        for _ in range(50):
            # d(phi)/d(ln T) = cp  ->  Newton step in ln T (keeps T > 0)
            d = (self.phi(math.exp(lnT), comp) - phi) / self.cp(math.exp(lnT), comp)
            lnT -= d
            if abs(d) < 1e-12:
                return math.exp(lnT)
        raise RuntimeError("T_from_phi did not converge")


@dataclass(frozen=True)
class HotFraction:
    """Composition for the two-gamma model: mass fraction that was burned."""

    x: float


class TwoGammaGas:
    """Mattingly's 'modified' constant-property model.

    Unburned air has (cp_c, gamma_c); everything downstream of the burner
    has (cp_t, gamma_t).  A mixed stream gets mass-weighted cp and R,
    exactly as in Mattingly's mixer equations
        cp6A = (cp_t + a' cp_c)/(1 + a'),   R6A = (R_t + a' R_c)/(1 + a').
    Enthalpy is h = cp T measured from 0 K, which makes the burner balance
    (1+f) cp_t Tt4 - cp_c Tt3 = f eta_b h_PR, i.e. Mattingly's
        f = (cp_t Tt4 - cp_c Tt3) / (eta_b h_PR - cp_t Tt4).
    """

    name = "two-gamma constant cp (Mattingly)"

    def __init__(self, cp_c=1004.0, gamma_c=1.4, cp_t=1235.0, gamma_t=1.3):
        self.cp_c, self.gamma_c = cp_c, gamma_c
        self.cp_t, self.gamma_t = cp_t, gamma_t
        self.R_c = cp_c * (gamma_c - 1.0) / gamma_c
        self.R_t = cp_t * (gamma_t - 1.0) / gamma_t

    def air(self):
        return HotFraction(0.0)

    def burn(self, comp, far):
        return HotFraction(1.0)

    def mix(self, streams):
        m = sum(mi for mi, _ in streams)
        return HotFraction(sum(mi * c.x for mi, c in streams) / m)

    def R(self, comp):
        return self.R_c + comp.x * (self.R_t - self.R_c)

    def cp(self, T, comp):
        return self.cp_c + comp.x * (self.cp_t - self.cp_c)

    def gamma(self, T, comp):
        cp = self.cp(T, comp)
        return cp / (cp - self.R(comp))

    def h(self, T, comp):
        return self.cp(T, comp) * T

    def T_from_h(self, h, comp, T_guess=None):
        return h / self.cp(0.0, comp)

    def phi(self, T, comp):
        return self.cp(T, comp) * math.log(T)

    def T_from_phi(self, phi, comp, T_guess=None):
        return math.exp(phi / self.cp(0.0, comp))
