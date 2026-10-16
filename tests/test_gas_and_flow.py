"""Property model and compressible-flow utilities."""
import math

import pytest

from turbofan import flow
from turbofan.gas import ThermallyPerfectGas, TwoGammaGas, AirFuel


def test_air_cp_matches_tables():
    # Air tables (e.g. Cengel A-2 / NIST): 1.005 at 300 K, 1.141 at 1000 K
    g = ThermallyPerfectGas()
    assert g.cp(300, g.air()) == pytest.approx(1005, abs=2.5)
    assert g.cp(1000, g.air()) == pytest.approx(1141, abs=3)
    assert g.gamma(300, g.air()) == pytest.approx(1.400, abs=0.002)


def test_products_have_higher_cp_than_air():
    g = ThermallyPerfectGas()
    for T in (500, 1000, 1500):
        assert g.cp(T, AirFuel(0.03)) > g.cp(T, g.air())


def test_enthalpy_and_entropy_inversions_round_trip():
    g = ThermallyPerfectGas()
    c = AirFuel(0.025)
    for T in (250.0, 700.0, 1800.0):
        assert g.T_from_h(g.h(T, c), c) == pytest.approx(T, rel=1e-10)
        assert g.T_from_phi(g.phi(T, c), c) == pytest.approx(T, rel=1e-10)


def test_entropy_function_is_integral_of_cp_over_T():
    g = ThermallyPerfectGas()
    c = g.air()
    T1, T2, n = 400.0, 1200.0, 4000
    dT = (T2 - T1) / n
    integral = sum(g.cp(T1 + (i + 0.5) * dT, c) / (T1 + (i + 0.5) * dT) * dT for i in range(n))
    assert g.phi(T2, c) - g.phi(T1, c) == pytest.approx(integral, rel=1e-7)


def test_mixing_rule_conserves_fuel():
    g = ThermallyPerfectGas()
    mixed = g.mix([(1.03, AirFuel(0.03)), (0.8, g.air())])
    assert mixed.far == pytest.approx(0.03 / 1.8)


def test_constant_gamma_limits():
    """With constant gamma the general routines must reproduce the closed forms."""
    g = TwoGammaGas(1004.0, 1.4, 1235.0, 1.3)
    a = g.air()
    Tt, Pt = 400.0, 2e5
    s = flow.sonic_state(g, a, Tt, Pt)
    assert Tt / s.T == pytest.approx(1.2, rel=1e-10)                 # (g+1)/2
    assert Pt / s.P == pytest.approx(1.2 ** 3.5, rel=1e-9)            # ((g+1)/2)^(g/(g-1))
    mfp_star = math.sqrt(1.4 / g.R_c) * 1.2 ** -3                     # MFP(M=1)
    assert flow.choked_mass_flux(g, a, Tt, Pt) == pytest.approx(mfp_star * Pt / math.sqrt(Tt), rel=1e-9)
    # area ratio A/A* at M = 2, gamma = 1.4 is 1.6875 (isentropic tables)
    s2 = flow.state_from_mach(g, a, Tt, Pt, 2.0)
    assert s.mass_flux / s2.mass_flux == pytest.approx(1.6875, rel=1e-9)
    # supersonic branch of the mass-flux inversion recovers M = 2
    back = flow.state_from_mass_flux(g, a, Tt, Pt, s2.mass_flux, supersonic=True)
    assert back.M == pytest.approx(2.0, rel=1e-8)


def test_polytropic_relation_reduces_to_textbook_form():
    from turbofan.components import Station, compress
    g = TwoGammaGas(1004.0, 1.4, 1235.0, 1.3)
    st = Station("2", 288.15, 101325.0, g.air(), 1.0)
    out = compress(g, st, 20.0, 0.9, "3")
    assert out.Tt / st.Tt == pytest.approx(20.0 ** (0.4 / (1.4 * 0.9)), rel=1e-12)
