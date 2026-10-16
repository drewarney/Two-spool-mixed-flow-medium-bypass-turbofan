"""Design the default two-spool mixed-flow turbofan, fly it from M0 = 0 to
1.7, optimise the C-D nozzle area ratio at every point and write figures.

    python scripts/run_engine_study.py [--out figures]
"""
import argparse
import csv
import os

import numpy as np

from turbofan.engine import DesignInputs, design
from turbofan.sweeps import flight_envelope
from turbofan import plots


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="figures")
    ap.add_argument("--dM", type=float, default=0.1)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    dp = design(DesignInputs())
    print("=== Design point ===")
    print(dp.cycle.summary())
    print(f"A4={dp.A4:.4f}  A45={dp.A45:.4f}  A6={dp.A6:.4f}  A16={dp.A16:.4f}  A8={dp.A8:.4f} m^2   TR={dp.throttle_ratio:.4f}")
    plots.ts_diagram(dp, os.path.join(args.out, "ts_diagram.png"))

    eng = dp.engine()
    machs = np.round(np.arange(0.0, 1.7 + 1e-9, args.dM), 3)
    # dynamic-pressure limits: no supersonic flight on the deck at sea level
    env = flight_envelope(eng, [0.0, 6000.0, 11000.0], machs, max_mach={0.0: 1.2, 6000.0: 1.5})

    with open(os.path.join(args.out, "optimal_area_schedule.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["altitude_m", "M0", "Tt4_K", "pi_f", "pi_c", "alpha", "mdot0_kg_s", "thrust_kN",
                    "tsfc_mg_Ns", "A9_A8_opt", "A9_A8_ideal_P9eqP0", "P9_P0_at_opt", "gain_vs_convergent_pct"])
        for alt, pts in env.items():
            for p in pts:
                op = p.op
                w.writerow([alt, p.M0, round(op.Tt4, 2), round(op.pi_f, 4), round(op.pi_c, 3), round(op.alpha, 4),
                            round(op.mdot0, 3), round(op.thrust / 1e3, 3), round(op.tsfc * 1e6, 3),
                            round(p.sweep.best_area_ratio, 4), round(p.sweep.ideal_area_ratio, 4),
                            round(op.nozzle_result.exit.P / op.ambient.P, 4), round(100 * p.gain_vs_convergent, 2)])

    print("\n=== Optimal A9/A8 schedule ===")
    print(f"{'alt':>6} {'M0':>4} {'F kN':>8} {'TSFC':>7} {'eps_opt':>8} {'eps_P9=P0':>9} {'gain %':>7}")
    for alt, pts in env.items():
        for p in pts:
            print(f"{alt:6.0f} {p.M0:4.1f} {p.op.thrust / 1e3:8.2f} {p.op.tsfc * 1e6:7.2f} "
                  f"{p.sweep.best_area_ratio:8.3f} {p.sweep.ideal_area_ratio:9.3f} {100 * p.gain_vs_convergent:7.2f}")

    plots.thrust_tsfc_vs_mach(env, os.path.join(args.out, "thrust_tsfc_vs_mach.png"))
    plots.optimal_area_schedule(env, os.path.join(args.out, "optimal_area_schedule.png"))
    plots.nozzle_sweep_curves(env, os.path.join(args.out, "nozzle_area_sweep_11km.png"), 11000.0,
                              [0.0, 0.8, 1.2, 1.7])
    print(f"\nfigures written to {args.out}/")


if __name__ == "__main__":
    main()
