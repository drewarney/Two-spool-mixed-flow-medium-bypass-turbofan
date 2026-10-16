"""Mattingly design problem 7.D1 (HP-1 airliner, separate-exhaust turbofan).

1. hand-calculation at pi_c = 36, pi_f = 1.8, alpha = 10 (printed step by step);
2. PARA-style sweep pi_c = 20..40 (step 2), alpha = 4, 6, 8, 10, 12 and optimum;
3. minimum uninstalled specific thrust for inlet diameters 2.2 - 3.5 m;
4. carpet plot TSFC vs F/mdot0 with S_max and F/mdot_min lines.

Aircraft data are PLACEHOLDERS (see turbofan/aircraft.py) - pass your own
Chapter 1 numbers with --weight, --wing-area, --cd0, --k1, --k2, --smax.
"""
import argparse
import os

import numpy as np

from turbofan.atmosphere import isa
from turbofan.parametric import GasConstants, separate_flow_turbofan, optimum_bypass_separate
from turbofan.tech_levels import hp1_losses
from turbofan import aircraft as AC
from turbofan import plots

M0, ALT, TT4 = 0.83, 11000.0, 1560.0
GC = GasConstants(gamma_c=1.4, cp_c=1004.0, gamma_t=1.3, cp_t=1235.0, h_PR=42.8e6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pi-f", type=float, default=1.8, help="assigned fan pressure ratio for the sweep")
    ap.add_argument("--weight", type=float, default=AC.HP1_PLACEHOLDER.weight, help="cruise weight [N]")
    ap.add_argument("--wing-area", type=float, default=AC.HP1_PLACEHOLDER.wing_area)
    ap.add_argument("--cd0", type=float, default=AC.HP1_PLACEHOLDER.CD0)
    ap.add_argument("--k1", type=float, default=AC.HP1_PLACEHOLDER.K1)
    ap.add_argument("--k2", type=float, default=AC.HP1_PLACEHOLDER.K2)
    ap.add_argument("--smax", type=float, default=AC.HP1_PLACEHOLDER.tsfc_max * 1e6, help="S_max [mg/(N s)]")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    amb = isa(ALT)
    L = hp1_losses()

    print("=== 7.D1 hand calculation: pi_c=36, pi_f=1.8, alpha=10, M0=0.83, 11 km ===")
    r = separate_flow_turbofan(M0=M0, T0=amb.T, Tt4=TT4, pi_c=36, pi_f=1.8, alpha=10, gc=GC, **L)
    for k, v in r.as_dict().items():
        print(f"  {k:>16} = {v:.6g}")
    print(f"  F/mdot0 = {r.specific_thrust:.2f} N/(kg/s)   S = {r.tsfc * 1e6:.3f} mg/(N s)")

    pi_cs = list(range(20, 41, 2))
    alphas = [4, 6, 8, 10, 12]
    res, rows = {}, []
    for pc in pi_cs:
        for al in alphas:
            try:
                rr = separate_flow_turbofan(M0=M0, T0=amb.T, Tt4=TT4, pi_c=pc, pi_f=a.pi_f, alpha=al, gc=GC, **L)
            except ValueError:
                continue   # turbine cannot drive this much bypass flow
            res[(al, pc)] = (rr.specific_thrust, rr.tsfc * 1e6)
        a_opt, rr = optimum_bypass_separate(M0=M0, T0=amb.T, Tt4=TT4, pi_c=pc, pi_f=a.pi_f, gc=GC, **L)
        res[("alpha_opt", pc)] = (rr.specific_thrust, rr.tsfc * 1e6)
        rows.append((pc, a_opt, rr.specific_thrust, rr.tsfc * 1e6))
    print(f"\n=== optimum bypass ratio at pi_f = {a.pi_f} ===")
    for pc, ao, F, S in rows:
        print(f"  pi_c={pc:2d}  alpha*={ao:6.3f}  F/mdot0={F:7.2f}  S={S:6.3f}")

    ac = AC.Aircraft("HP-1", a.weight, a.wing_area, a.cd0, a.k1, a.k2, install_loss=0.02)
    F_req = AC.required_uninstalled_thrust_per_engine(ac, amb, M0, Ps=1.5)
    print(f"\n=== minimum specific thrust (required uninstalled thrust per engine = {F_req / 1e3:.2f} kN) ===")
    fmins, labels = [], []
    for d in (2.2, 2.5, 2.75, 3.0, 3.25, 3.5):
        m = AC.capture_mass_flow(amb, M0, AC.inlet_area_from_diameter(d))
        fm = AC.minimum_specific_thrust(F_req, m)
        fmins.append(fm)
        labels.append(f"D = {d:g} m")
        print(f"  D={d:4.2f} m  mdot_max={m:7.1f} kg/s  F/mdot_min={fm:7.2f} N/(kg/s)")

    keys = alphas + ["alpha_opt"]
    plots.carpet(res, keys, pi_cs, os.path.join(a.out, "carpet_hp1_7D1.png"),
                 f"HP-1 (7.D1) cruise M0 = 0.83, 11 km, pi_f = {a.pi_f:g}, Tt4 = 1560 K",
                 "Specific thrust F/mdot0  [N/(kg/s)]", "TSFC  [mg/(N s)]",
                 f_min_lines=fmins, f_min_labels=labels, s_max=a.smax)
    print(f"\ncarpet plot -> {a.out}/carpet_hp1_7D1.png  (aircraft data: {'PLACEHOLDER' if a.weight == AC.HP1_PLACEHOLDER.weight else 'user'})")


if __name__ == "__main__":
    main()
