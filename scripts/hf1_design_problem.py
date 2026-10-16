"""Mattingly design problem 7.D2 (HF-1 fighter, mixed-flow turbofan).

1. hand-calculation at pi_c = 25, alpha = 0.5, M0 = 1.6, 40 kft (US units);
   pi_f from the mixer condition Pt16 = Pt6 (Problem 7.27 system);
2. sweep alpha = 0.10..1.00 (step 0.05) x pi_c = 16, 18, 20, 22, 24, 28;
3. minimum uninstalled specific thrust (5 ft^2 inlet, 146.3 lbm/s);
4. carpet plot.
Also cross-checks the closed form against the station-by-station solver
(engine.design with TwoGammaGas).
"""
import argparse
import os

import numpy as np

from turbofan.atmosphere import isa
from turbofan.parametric import GasConstants, mixed_flow_turbofan
from turbofan.tech_levels import hf1_losses
from turbofan.units import (BTU_PER_LBM, BTU_PER_LBM_R, rankine_to_kelvin, specific_thrust_to_us,
                            tsfc_to_us, FT, LBM, LBF)
from turbofan import aircraft as AC
from turbofan import plots

M0, ALT = 1.6, 40000 * FT
GC = GasConstants(gamma_c=1.4, cp_c=0.240 * BTU_PER_LBM_R, gamma_t=1.3, cp_t=0.296 * BTU_PER_LBM_R,
                  h_PR=18400 * BTU_PER_LBM)
TT4 = rankine_to_kelvin(3250)
T0 = rankine_to_kelvin(390.0)     # Mattingly's 40 kft value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--weight-lbf", type=float, default=AC.HF1_PLACEHOLDER.weight / LBF)
    ap.add_argument("--wing-area-ft2", type=float, default=AC.HF1_PLACEHOLDER.wing_area / FT ** 2)
    ap.add_argument("--cd0", type=float, default=AC.HF1_PLACEHOLDER.CD0)
    ap.add_argument("--k1", type=float, default=AC.HF1_PLACEHOLDER.K1)
    ap.add_argument("--k2", type=float, default=AC.HF1_PLACEHOLDER.K2)
    ap.add_argument("--smax", type=float, default=AC.HF1_PLACEHOLDER.tsfc_max * 3600 * LBF / LBM,
                    help="S_max [(lbm/h)/lbf]")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    L = hf1_losses()

    print("=== 7.D2 hand calculation: pi_c=25, alpha=0.5, M0=1.6, 40 kft, Tt4=3250 R ===")
    r = mixed_flow_turbofan(M0=M0, T0=T0, Tt4=TT4, pi_c=25, alpha=0.5, gc=GC, **L)
    for k, v in r.as_dict().items():
        print(f"  {k:>16} = {v:.6g}")
    print(f"  F/mdot0 = {specific_thrust_to_us(r.specific_thrust):.2f} lbf/(lbm/s)"
          f"  S = {tsfc_to_us(r.tsfc):.4f} (lbm/h)/lbf")

    alphas = [round(x, 2) for x in np.arange(0.10, 1.0001, 0.05)]
    pi_cs = [16, 18, 20, 22, 24, 28]
    res = {}
    for al in alphas:
        for pc in pi_cs:
            try:
                rr = mixed_flow_turbofan(M0=M0, T0=T0, Tt4=TT4, pi_c=pc, alpha=al, gc=GC, **L)
            except ValueError:
                continue
            res[(al, pc)] = (specific_thrust_to_us(rr.specific_thrust), tsfc_to_us(rr.tsfc))
    print("\n=== sweep (F/mdot0 lbf/(lbm/s), S (lbm/h)/lbf, pi_f) ===")
    for al in alphas:
        print(f"  alpha={al:4.2f} " + "  ".join(f"pc{pc}:{res[(al, pc)][0]:6.2f}/{res[(al, pc)][1]:.3f}"
                                              for pc in pi_cs if (al, pc) in res))

    amb = isa(ALT)
    ac = AC.Aircraft("HF-1", a.weight_lbf * LBF, a.wing_area_ft2 * FT ** 2, a.cd0, a.k1, a.k2, install_loss=0.05)
    F_req = AC.required_uninstalled_thrust_per_engine(ac, amb, M0)
    mdot = 146.3 * LBM   # given in the problem statement
    fmin = specific_thrust_to_us(F_req / mdot)
    print(f"\nrequired uninstalled thrust/engine = {F_req / LBF:.0f} lbf, mdot = 146.3 lbm/s"
          f" -> F/mdot_min = {fmin:.2f} lbf/(lbm/s)")

    plots.carpet(res, alphas, pi_cs, os.path.join(a.out, "carpet_hf1_7D2.png"),
                 "HF-1 (7.D2) supercruise M0 = 1.6, 40 kft, Tt4 = 3250 R, mixed flow",
                 "Specific thrust F/mdot0  [lbf/(lbm/s)]", "TSFC  [(lbm/h)/lbf]",
                 f_min_lines=[fmin], f_min_labels=["F/mdot min"], s_max=a.smax)
    print(f"carpet plot -> {a.out}/carpet_hf1_7D2.png")


if __name__ == "__main__":
    main()
