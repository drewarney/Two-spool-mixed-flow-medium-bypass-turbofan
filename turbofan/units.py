"""Unit conversions between SI (used internally) and US customary units
(used by Mattingly's HF-1 problems)."""

FT = 0.3048                    # m
LBM = 0.45359237               # kg
LBF = 4.4482216152605          # N
RANKINE = 5.0 / 9.0            # K per degR
BTU_PER_LBM = 2326.0           # J/kg  (1 Btu/lbm)
BTU_PER_LBM_R = 4186.8         # J/(kg K)  (1 Btu/(lbm degR))

# Specific thrust: 1 N/(kg/s)  ->  lbf/(lbm/s)
N_PER_KGS_TO_LBF_PER_LBMS = LBM / LBF          # 0.10197
# TSFC: 1 kg/(N s)  ->  (lbm/h)/lbf
KG_PER_NS_TO_LBM_PER_HR_LBF = 3600.0 * LBF / LBM  # 35303.9


def rankine_to_kelvin(T_R: float) -> float:
    return T_R * RANKINE


def specific_thrust_to_us(F_per_mdot_si: float) -> float:
    return F_per_mdot_si * N_PER_KGS_TO_LBF_PER_LBMS


def tsfc_to_us(S_si: float) -> float:
    """kg/(N s) -> (lbm/h)/lbf."""
    return S_si * KG_PER_NS_TO_LBM_PER_HR_LBF


def tsfc_to_mg_per_ns(S_si: float) -> float:
    return S_si * 1e6
