import pytest

from turbofan.engine import DesignInputs, design
from turbofan.gas import TwoGammaGas
from turbofan.nozzle import NozzleModel
from turbofan.parametric import GasConstants
from turbofan.units import BTU_PER_LBM, BTU_PER_LBM_R, rankine_to_kelvin, FT

# Problem 7.D2 data (HF-1 supercruise), converted to SI
HF1_GC = GasConstants(gamma_c=1.4, cp_c=0.240 * BTU_PER_LBM_R, gamma_t=1.3, cp_t=0.296 * BTU_PER_LBM_R,
                      h_PR=18400 * BTU_PER_LBM)
HF1_T0 = rankine_to_kelvin(390.0)
HF1_TT4 = rankine_to_kelvin(3250.0)


@pytest.fixture(scope="session")
def tp_design():
    """Default engine with the thermally perfect gas."""
    return design(DesignInputs())


@pytest.fixture(scope="session")
def two_gamma():
    return TwoGammaGas(HF1_GC.cp_c, 1.4, HF1_GC.cp_t, 1.3)


def hf1_inputs(**kw):
    base = dict(M0=1.6, altitude=40000 * FT, bypass_ratio=0.5, pi_c=25.0, Tt4=HF1_TT4, h_PR=HF1_GC.h_PR,
                pi_d_max=0.96, e_f=0.89, e_cH=0.90, e_tH=0.89, e_tL=0.89, pi_b=0.96, eta_b=0.999,
                eta_mH=0.99, eta_mL=0.99, pi_M_max=0.96, M6=0.4,
                nozzle=NozzleModel(pi_n=0.97, flap_length_ratio=None, separation_pressure_ratio=None))
    base.update(kw)
    return DesignInputs(**base)
