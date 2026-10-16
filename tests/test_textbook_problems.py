"""Reproduce the hand-calculations of Mattingly design problems 7.D1 and 7.D2.

Two layers:
  * the individual hand-calc steps (tau_r, pi_r, eta_r, tau_lambda, tau_c ...)
    are checked against values evaluated directly from the problem data;
  * the full mixed-flow cycle is computed twice - by the closed-form
    textbook equations and by the station-by-station solver run with the
    same two-gamma gas - and the two must agree to round-off.
"""
import math

import pytest

from turbofan.atmosphere import isa
from turbofan.engine import design
from turbofan.parametric import (GasConstants, separate_flow_turbofan, mixed_flow_turbofan,
                                 optimum_bypass_separate)
from turbofan.tech_levels import hp1_losses, hf1_losses
from turbofan.units import specific_thrust_to_us, tsfc_to_us, LBM, FT
from turbofan import aircraft as AC

from conftest import HF1_GC, HF1_T0, HF1_TT4, hf1_inputs

HP1_GC = GasConstants(1.4, 1004.0, 1.3, 1235.0, 42.8e6)


# ------------------------------------------------------------- 7.D1 (HP-1)
@pytest.fixture(scope="module")
def hp1():
    return separate_flow_turbofan(M0=0.83, T0=isa(11000).T, Tt4=1560.0, pi_c=36, pi_f=1.8, alpha=10,
                                  gc=HP1_GC, **hp1_losses())


def test_hp1_hand_calc_steps(hp1):
    T0 = isa(11000).T
    assert T0 == pytest.approx(216.65, abs=0.01)
    assert hp1.tau_r == pytest.approx(1 + 0.2 * 0.83 ** 2, rel=1e-12)          # 1.13778
    assert hp1.pi_r == pytest.approx(hp1.tau_r ** 3.5, rel=1e-12)              # 1.5711
    assert hp1.eta_r == 1.0                                                    # subsonic
    assert hp1.tau_lambda == pytest.approx(1235 * 1560 / (1004 * T0), rel=1e-12)  # 8.857
    assert hp1.tau_c == pytest.approx(36 ** (0.4 / (1.4 * 0.90)), rel=1e-12)   # 3.1194
    assert hp1.tau_f == pytest.approx(1.8 ** (0.4 / (1.4 * 0.89)), rel=1e-12)  # 1.2077
    f = (hp1.tau_lambda - hp1.tau_r * hp1.tau_c) / (0.999 * 42800 / (1.004 * T0) - hp1.tau_lambda)
    assert hp1.f == pytest.approx(f, rel=1e-12)


def test_hp1_power_balance_closes(hp1):
    # turbine work = (compressor + fan work)/eta_m, in cp T0 units
    lhs = 0.99 * (1 + hp1.f) * hp1.tau_lambda * (1 - hp1.tau_t)
    rhs = hp1.tau_r * (hp1.tau_c - 1 + 10 * (hp1.tau_f - 1))
    assert lhs == pytest.approx(rhs, rel=1e-12)


def test_hp1_hand_calc_result(hp1):
    # Our hand-calculation (regression values; compare with PARA in your copy)
    assert hp1.specific_thrust == pytest.approx(143.25, abs=0.02)   # N/(kg/s)
    assert hp1.tsfc * 1e6 == pytest.approx(17.946, abs=0.002)        # mg/(N s)
    assert hp1.f == pytest.approx(0.028278, abs=2e-6)


def test_hp1_optimum_bypass_equalises_jets():
    """At the TSFC-optimal bypass ratio the core and fan jets leave at
    (almost) the same velocity - the classic result (Mattingly: V9 ~ V19)."""
    a_opt, r = optimum_bypass_separate(M0=0.83, T0=isa(11000).T, Tt4=1560.0, pi_c=36, pi_f=1.8,
                                       gc=HP1_GC, **hp1_losses())
    assert 9 < a_opt < 11
    # with losses the optimum is where d(KE)/d(alpha) balances, close to V9/V19 = eta_ratio
    assert r.V9_a0 / r.V19_a0 == pytest.approx(1.0, abs=0.12)


def test_hp1_capture_mass_flow():
    amb = isa(11000)
    m = AC.capture_mass_flow(amb, 0.83, AC.inlet_area_from_diameter(2.2))
    assert m == pytest.approx(amb.rho * math.pi * 1.1 ** 2 * 0.83 * amb.a)
    assert 330 < m < 345


# ------------------------------------------------------------- 7.D2 (HF-1)
@pytest.fixture(scope="module")
def hf1():
    return mixed_flow_turbofan(M0=1.6, T0=HF1_T0, Tt4=HF1_TT4, pi_c=25, alpha=0.5, gc=HF1_GC, **hf1_losses())


def test_hf1_hand_calc_steps(hf1):
    assert hf1.tau_r == pytest.approx(1.512, rel=1e-12)
    assert hf1.pi_r == pytest.approx(1.512 ** 3.5, rel=1e-12)                   # 4.2504
    assert hf1.eta_r == pytest.approx(1 - 0.075 * 0.6 ** 1.35, rel=1e-12)       # 0.96237 (MIL-E-5008B)
    assert hf1.pi_d == pytest.approx(0.96 * hf1.eta_r, rel=1e-12)
    assert hf1.tau_lambda == pytest.approx(0.296 * 3250 / (0.240 * 390), rel=1e-12)  # 10.2778
    assert hf1.tau_c == pytest.approx(25 ** (0.4 / 1.26), rel=1e-12)            # 2.7784
    f = (hf1.tau_lambda - hf1.tau_r * hf1.tau_c) / (0.999 * 18400 / (0.240 * 390) - hf1.tau_lambda)
    assert hf1.f == pytest.approx(f, rel=1e-12)                                 # 0.03265


def test_hf1_mixer_design_condition(hf1):
    """Problem 7.27 system: pi_f chosen so that Pt16 = Pt6."""
    assert hf1.pi_f == pytest.approx(25 * 0.96 * hf1.pi_t, rel=1e-10)
    assert hf1.M16 < hf1.M6A < 1 or hf1.M16 < 1
    # mixed-gas properties are the mass-weighted averages
    ap = 0.5 / (1 + hf1.f)
    assert hf1.alpha_prime == pytest.approx(ap)
    assert hf1.cp_6A == pytest.approx((HF1_GC.cp_t + ap * HF1_GC.cp_c) / (1 + ap))
    assert hf1.pi_M < 0.96                 # friction x mixing loss


def test_hf1_hand_calc_result(hf1):
    # Our hand-calculation in US units (regression; compare with PARA in your copy)
    assert specific_thrust_to_us(hf1.specific_thrust) == pytest.approx(66.27, abs=0.01)  # lbf/(lbm/s)
    assert tsfc_to_us(hf1.tsfc) == pytest.approx(1.1826, abs=0.0002)                      # (lbm/h)/lbf
    assert hf1.pi_f == pytest.approx(4.268, abs=0.001)


def test_hf1_given_inlet_mass_flow():
    # Problem statement: mdot = sigma rho_ref A M a = 146.3 lbm/s for a 5 ft^2 inlet
    amb = isa(40000 * FT)
    m = AC.capture_mass_flow(amb, 1.6, 5 * FT ** 2) / LBM
    assert m == pytest.approx(146.3, rel=0.006)


def test_station_solver_reproduces_closed_form(two_gamma):
    """THE key verification: run the general station-by-station code with
    Mattingly's two-gamma gas and compare with his closed-form equations.

    Mattingly compresses ALL core air at e_c (it lumps the fan's core side
    into tau_c); our two-spool model sends core air through the fan (e_f)
    and then the HPC (e_c).  Setting e_f = e_c makes the two identical.
    """
    inp = hf1_inputs(e_f=0.90)
    dp = design(inp, two_gamma)
    cyc = dp.cycle
    ref = mixed_flow_turbofan(M0=1.6, T0=cyc.ambient.T, Tt4=HF1_TT4, pi_c=25, alpha=0.5, gc=HF1_GC,
                              **{**hf1_losses(), "e_f": 0.90})
    assert cyc.pi_f == pytest.approx(ref.pi_f, rel=1e-9)
    assert cyc.f == pytest.approx(ref.f, rel=1e-9)
    assert cyc.pi_tH * cyc.pi_tL == pytest.approx(ref.pi_t, rel=1e-9)
    assert dp.A16 / dp.A6 == pytest.approx(ref.A16_A6, rel=1e-8)
    assert cyc.stations["6A"].static.M == pytest.approx(ref.M6A, rel=1e-8)
    assert cyc.stations["6A"].Pt / cyc.stations["6"].Pt == pytest.approx(ref.pi_M, rel=1e-9)
    assert cyc.nozzle_result.exit.M == pytest.approx(ref.M9, rel=1e-8)
    assert cyc.specific_thrust == pytest.approx(ref.specific_thrust, rel=1e-9)
    assert cyc.tsfc == pytest.approx(ref.tsfc, rel=1e-9)


def test_variable_cp_changes_answer_by_a_few_percent(two_gamma):
    """Sanity: the thermally perfect model should differ from the two-gamma
    idealisation, but only by a few percent (both describe the same engine)."""
    a = design(hf1_inputs(), two_gamma).cycle
    b = design(hf1_inputs()).cycle
    assert b.specific_thrust / a.specific_thrust == pytest.approx(1.0, abs=0.06)
    # TSFC differs more: with cp_t = 0.296 Btu/lbm/R applied from 0 K, the
    # two-gamma burner balance over-states the enthalpy rise to Tt4 and
    # so over-states the fuel flow (~13 % here).
    assert b.f < a.f
    assert b.tsfc / a.tsfc == pytest.approx(1.0, abs=0.15)
    assert b.specific_thrust != pytest.approx(a.specific_thrust, rel=1e-3)
