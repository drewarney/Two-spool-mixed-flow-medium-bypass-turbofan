"""1976 US Standard / ISA atmosphere, 0-20 km (troposphere + lower stratosphere).

Troposphere (h < 11 km): temperature falls linearly with lapse rate L,
    T = T_sl - L h
and hydrostatic balance dp = -rho g dh with p = rho R T integrates to
    p/p_sl = (T/T_sl)^(g0 / (R L)).
Lower stratosphere (11-20 km): isothermal at 216.65 K, so
    p = p_11 exp(-g0 (h - 11 km) / (R T_11)).
Altitudes are geopotential metres (the difference from geometric altitude
is < 0.2 % at 11 km and irrelevant for cycle analysis).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

T_SL = 288.15      # K
P_SL = 101325.0    # Pa
RHO_SL = 1.225     # kg/m^3
LAPSE = 0.0065     # K/m
R_AIR = 287.05     # J/(kg K)
G0 = 9.80665       # m/s^2
H_TROP = 11000.0   # m
T_TROP = T_SL - LAPSE * H_TROP  # 216.65 K
P_TROP = P_SL * (T_TROP / T_SL) ** (G0 / (R_AIR * LAPSE))


@dataclass(frozen=True)
class Ambient:
    altitude: float
    T: float
    P: float

    @property
    def rho(self) -> float:
        return self.P / (R_AIR * self.T)

    @property
    def a(self) -> float:
        # speed of sound of cold air, gamma = 1.4
        return math.sqrt(1.4 * R_AIR * self.T)


def isa(altitude_m: float) -> Ambient:
    h = float(altitude_m)
    if h < 0 or h > 20000:
        raise ValueError("isa() valid for 0-20 km")
    if h <= H_TROP:
        T = T_SL - LAPSE * h
        P = P_SL * (T / T_SL) ** (G0 / (R_AIR * LAPSE))
    else:
        T = T_TROP
        P = P_TROP * math.exp(-G0 * (h - H_TROP) / (R_AIR * T_TROP))
    return Ambient(h, T, P)
