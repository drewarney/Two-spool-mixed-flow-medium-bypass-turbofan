"""Station-by-station engine: conservation checks and off-design behaviour."""
import math

import pytest

from turbofan import flow
from turbofan.engine import design
from turbofan.gas import TwoGammaGas
from turbofan.atmosphere import T_SL

from conftest import hf1_inputs


def h(gas, st):
    return gas.h(st.Tt, st.comp)


def test_spool_power_balances(tp_design):
    dp, g = tp_design, tp_design.gas
    c = dp.cycle
    s = c.stations
    inp = dp.inputs
    hp_comp = s["3"].mdot * (h(g, s["3"]) - h(g, s["25"]))
    hp_turb = s["4"].mdot * (h(g, s["4"]) - h(g, s["45"]))
    assert hp_comp == pytest.approx(inp.eta_mH * hp_turb, rel=1e-9)
    fan = s["2"].mdot * (h(g, s["13"]) - h(g, s["2"]))
    lp_turb = s["45"].mdot * (h(g, s["45"]) - h(g, s["5"]))
    assert fan == pytest.approx(inp.eta_mL * lp_turb, rel=1e-9)


def test_burner_energy_balance(tp_design):
    g, c = tp_design.gas, tp_design.cycle
    s = c.stations
    lhs = s["3"].mdot * h(g, s["3"]) + c.f * s["3"].mdot * tp_design.inputs.eta_b * tp_design.inputs.h_PR
    assert lhs == pytest.approx(s["4"].mdot * h(g, s["4"]), rel=1e-10)


def test_mixer_conserves_mass_energy_momentum(tp_design):
    g, c, dp = tp_design.gas, tp_design.cycle, tp_design
    s6, s16, s6A = c.stations["6"], c.stations["16"], c.stations["6A"]
    assert s6A.mdot == pytest.approx(s6.mdot + s16.mdot, rel=1e-12)
    assert s6A.mdot * h(g, s6A) == pytest.approx(s6.mdot * h(g, s6) + s16.mdot * h(g, s16), rel=1e-10)
    I_in = s6.static.P * dp.A6 + s6.mdot * s6.static.V + s16.static.P * dp.A16 + s16.mdot * s16.static.V
    I_out = s6A.static.P * s6A.area + s6A.mdot * s6A.static.V
    assert I_out == pytest.approx(I_in, rel=1e-9)
    # design rules: equal static pressure and equal total pressure at the splitter
    assert s6.static.P == pytest.approx(s16.static.P, rel=1e-9)
    assert s6.Pt == pytest.approx(s16.Pt, rel=1e-9)
    # mixing is irreversible: Pt6A below the inlet total pressure even before friction
    assert s6A.Pt / tp_design.inputs.pi_M_max < s6.Pt


def test_design_has_ideal_expansion(tp_design):
    nr = tp_design.cycle.nozzle_result
    assert nr.exit.P == pytest.approx(tp_design.cycle.ambient.P, rel=1e-8)
    assert nr.choked


def test_offdesign_reproduces_design(tp_design):
    dp = tp_design
    eng = dp.engine()
    op = eng.operate(dp.inputs.M0, dp.inputs.altitude, area_ratio=dp.cycle.nozzle_result.area_ratio)
    assert op.pi_f == pytest.approx(dp.cycle.pi_f, rel=1e-7)
    assert op.alpha == pytest.approx(dp.cycle.alpha, rel=1e-7)
    assert op.mdot0 == pytest.approx(dp.inputs.mdot0, rel=1e-7)
    assert op.thrust == pytest.approx(dp.cycle.thrust, rel=1e-7)


def test_offdesign_similarity_at_equal_inlet_temperature(two_gamma):
    """Constant-gamma engine: two flight conditions with the same Tt0 get the
    same Tt4 from the schedule, so the engine behind the inlet differs only
    by the scale factor Pt2.  Every ratio must then be identical - a strong
    check that the off-design equations are dimensionally consistent."""
    dp = design(hf1_inputs(), two_gamma)
    eng = dp.engine()
    a = eng.operate(0.0, 0.0, area_ratio=1.3)                    # Tt0 = 288.15 K
    M = math.sqrt((T_SL / 216.65 - 1) / 0.2)                     # same Tt0 at 11 km
    b = eng.operate(M, 11000.0, area_ratio=1.3)
    assert b.stations["0"].Tt == pytest.approx(a.stations["0"].Tt, rel=1e-6)
    for attr in ("pi_f", "pi_cH", "pi_tH", "pi_tL", "alpha", "f"):
        assert getattr(b, attr) == pytest.approx(getattr(a, attr), rel=1e-6), attr
    ratio = b.stations["2"].Pt / a.stations["2"].Pt
    assert b.mdot0 / a.mdot0 == pytest.approx(ratio, rel=1e-6)   # corrected flow unchanged


def test_offdesign_near_similar_below_theta_break(two_gamma):
    """Below the theta break Tt4/Tt0 is held, so corrected operation is
    almost frozen; only the (1 + f) mass-addition terms (f does not scale
    with temperature because h_PR is fixed) move the pressure ratios."""
    dp = design(hf1_inputs(), two_gamma)
    eng = dp.engine()
    op = eng.operate(0.0, 0.0, area_ratio=1.3)
    assert op.Tt4 < dp.inputs.Tt4
    assert op.pi_f == pytest.approx(dp.cycle.pi_f, rel=0.02)
    assert op.alpha == pytest.approx(dp.cycle.alpha, rel=0.03)


def test_offdesign_unloads_above_theta_break(tp_design):
    """Above the theta break Tt4 is capped, Tt4/Tt2 falls, the compressors
    unload (pi_c down) and the bypass ratio rises - textbook behaviour."""
    eng = tp_design.engine()
    hot = eng.operate(1.4, 0.0)
    assert hot.Tt4 == pytest.approx(tp_design.inputs.Tt4)
    assert hot.pi_c < tp_design.cycle.pi_c - 3
    assert hot.alpha > tp_design.cycle.alpha
    assert hot.stations["0"].Tt / T_SL > tp_design.throttle_ratio
