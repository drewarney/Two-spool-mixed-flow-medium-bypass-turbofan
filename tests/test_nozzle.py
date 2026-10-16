"""C-D nozzle area-ratio optimisation."""
import pytest

from turbofan import flow
from turbofan.gas import ThermallyPerfectGas, AirFuel, TwoGammaGas
from turbofan.nozzle import NozzleModel, nozzle, optimise_area_ratio, ideal_expansion_area_ratio

GAS = ThermallyPerfectGas()
COMP = AirFuel(0.015)
TT8, PT8, P0, A8 = 950.0, 260e3, 18.75e3, 0.30


def mdot():
    return A8 * flow.choked_mass_flux(GAS, COMP, TT8, PT8)


def test_lossless_optimum_is_perfect_expansion():
    """dFg/dA9 = P9 - P0, so the lossless optimum is exactly P9 = P0."""
    model = NozzleModel(pi_n=1.0, flap_length_ratio=None, separation_pressure_ratio=None)
    sw = optimise_area_ratio(GAS, COMP, TT8, PT8, mdot(), A8, P0, model, eps_max=4.0)
    assert sw.best_area_ratio == pytest.approx(sw.ideal_area_ratio, rel=1e-4)
    r = nozzle(GAS, COMP, TT8, PT8, mdot(), A8, sw.best_area_ratio, P0, model)
    assert r.exit.P / P0 == pytest.approx(1.0, abs=1e-3)


def test_thrust_derivative_equals_pressure_mismatch():
    model = NozzleModel(pi_n=1.0, flap_length_ratio=None, separation_pressure_ratio=None)
    m = mdot()
    for eps in (1.4, 2.0, 3.0):
        d = 1e-5
        F1 = nozzle(GAS, COMP, TT8, PT8, m, A8, eps - d, P0, model).gross_thrust
        F2 = nozzle(GAS, COMP, TT8, PT8, m, A8, eps + d, P0, model).gross_thrust
        dFdA9 = (F2 - F1) / (2 * d * A8)
        P9 = nozzle(GAS, COMP, TT8, PT8, m, A8, eps, P0, model).exit.P
        assert dFdA9 == pytest.approx(P9 - P0, rel=2e-4, abs=1.0)


def test_divergence_loss_moves_optimum_to_underexpansion():
    model = NozzleModel(pi_n=1.0, flap_length_ratio=3.0, separation_pressure_ratio=None)
    sw = optimise_area_ratio(GAS, COMP, TT8, PT8, mdot(), A8, P0, model, eps_max=4.0)
    assert sw.best_area_ratio < sw.ideal_area_ratio
    r = nozzle(GAS, COMP, TT8, PT8, mdot(), A8, sw.best_area_ratio, P0, model)
    assert r.exit.P > P0


def test_separation_limits_overexpansion_penalty():
    m = mdot()
    attached = NozzleModel(pi_n=1.0, flap_length_ratio=None, separation_pressure_ratio=None)
    sep = NozzleModel(pi_n=1.0, flap_length_ratio=None, separation_pressure_ratio=0.37)
    P0_sl = 101325.0  # NPR ~ 2.6: a big eps is badly over-expanded at sea level
    Fa = nozzle(GAS, COMP, TT8, PT8, m, A8, 2.5, P0_sl, attached)
    Fs = nozzle(GAS, COMP, TT8, PT8, m, A8, 2.5, P0_sl, sep)
    assert Fs.separated and not Fa.separated
    assert Fs.gross_thrust > Fa.gross_thrust


def test_ideal_area_ratio_matches_isentropic_tables():
    g = TwoGammaGas(1004.0, 1.4, 1004.0, 1.4)
    a = g.air()
    # Pt/P = 7.824 <-> M = 2 for gamma = 1.4, A/A* = 1.6875
    eps = ideal_expansion_area_ratio(g, a, 300.0, 7.824449 * 1e5, 1e5)
    assert eps == pytest.approx(1.6875, rel=1e-5)


def test_unchoked_nozzle_exits_at_ambient():
    # eps = 1.2: a subsonic exit with a sonic throat has P9/Pt ~ 0.79, so a
    # back pressure above that (here P0/Pt = 0.84) leaves the throat unchoked
    model = NozzleModel()
    m = 0.5 * mdot()
    r = nozzle(GAS, COMP, TT8, 120e3, m, A8, 1.2, 101325.0, model)
    assert not r.choked
    assert r.exit.P == pytest.approx(101325.0)
