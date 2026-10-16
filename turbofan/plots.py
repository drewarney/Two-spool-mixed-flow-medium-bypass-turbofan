"""Matplotlib figures.  Style: thin 2 px lines, recessive grid, fixed
categorical colour order (validated palette), direct labels, one y-axis
per panel (two measures -> two panels, never a twin axis)."""
from __future__ import annotations

import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .units import tsfc_to_mg_per_ns

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"


def _style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0,
        "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False,
        "axes.titlelocation": "left",
    })


def _alt_label(alt):
    return f"{alt / 1000:g} km"


def thrust_tsfc_vs_mach(envelope, path):
    _style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for i, (alt, pts) in enumerate(envelope.items()):
        M = [p.M0 for p in pts]
        c = SERIES[i]
        ax1.plot(M, [p.op.thrust / 1e3 for p in pts], color=c, label=_alt_label(alt))
        ax2.plot(M, [tsfc_to_mg_per_ns(p.op.tsfc) for p in pts], color=c, label=_alt_label(alt))
        ax1.annotate(_alt_label(alt), (M[-1], pts[-1].op.thrust / 1e3), xytext=(4, 0),
                     textcoords="offset points", va="center", fontsize=9, color=INK2)
    ax1.set(title="Uninstalled net thrust (optimal A9/A8)", xlabel="Flight Mach number M0", ylabel="Thrust  [kN]")
    ax2.set(title="Thrust-specific fuel consumption", xlabel="Flight Mach number M0", ylabel="TSFC  [mg/(N s)]")
    ax2.legend(title="Altitude", loc="upper left")
    ax1.set_ylim(bottom=0)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def optimal_area_schedule(envelope, path):
    _style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for i, (alt, pts) in enumerate(envelope.items()):
        M = [p.M0 for p in pts]
        c = SERIES[i]
        ax1.plot(M, [p.sweep.best_area_ratio for p in pts], color=c, label=f"{_alt_label(alt)}  max thrust")
        ax1.plot(M, [p.sweep.ideal_area_ratio for p in pts], color=c, lw=1.2, ls="--",
                 label=f"{_alt_label(alt)}  P9 = P0")
        ax2.plot(M, [100 * p.gain_vs_convergent for p in pts], color=c, label=_alt_label(alt))
        ax2.annotate(_alt_label(alt), (M[-1], 100 * pts[-1].gain_vs_convergent), xytext=(4, 0),
                     textcoords="offset points", va="center", fontsize=9, color=INK2)
    ax1.set(title="Thrust-maximising nozzle area ratio", xlabel="Flight Mach number M0",
            ylabel="A9 / A8")
    ax1.legend(fontsize=8, ncol=1, loc="upper left")
    ax2.set(title="Net-thrust gain over a convergent nozzle (A9/A8 = 1)",
            xlabel="Flight Mach number M0", ylabel="Gain  [%]")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def nozzle_sweep_curves(envelope, path, altitude, machs):
    _style()
    fig, ax = plt.subplots(figsize=(7.5, 4.4), constrained_layout=True)
    pts = {round(p.M0, 3): p for p in envelope[altitude]}
    for i, M in enumerate(machs):
        p = pts[round(M, 3)]
        s = p.sweep
        Fn = s.gross_thrust - p.op.mdot0 * p.op.V0
        Fn_best = s.best_gross_thrust - p.op.mdot0 * p.op.V0
        y = Fn / Fn_best
        ax.plot(s.area_ratios, y, color=SERIES[i], label=f"M0 = {M:g}")
        ax.plot([s.best_area_ratio], [1.0], "o", ms=8, color=SERIES[i], mec=SURFACE, mew=2)
    ax.set(title=f"Net thrust vs nozzle area ratio at {_alt_label(altitude)} (dot = optimum)",
           xlabel="A9 / A8", ylabel="F / F_max")
    lo = min(np.min((p.sweep.gross_thrust - p.op.mdot0 * p.op.V0) /
                    (p.sweep.best_gross_thrust - p.op.mdot0 * p.op.V0)) for p in (pts[round(M, 3)] for M in machs))
    ax.set_ylim(max(lo, 0.8), 1.01)
    ax.legend(loc="lower right")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def ts_diagram(dp, path):
    """T-s diagram of the design cycle, entropy relative to ambient air.

    s - s0 = [phi(T, comp) - phi(T0, air)] - R ln(P/P0)    (thermally perfect gas)
    Stagnation states are plotted for internal stations (they are the
    quantities the cycle tracks), static states for 0 and 9."""
    _style()
    gas, cyc = dp.gas, dp.cycle
    amb = cyc.ambient
    air = gas.air()
    phi0 = gas.phi(amb.T, air)

    def s_of(T, P, comp):
        return gas.phi(T, comp) - phi0 - gas.R(comp) * math.log(P / amb.P)

    st = cyc.stations
    pts = {"0": (amb.T, amb.P, air)}
    for k in ["2", "13", "16", "25", "3", "4", "45", "5", "6A", "8"]:
        pts[k] = (st[k].Tt, st[k].Pt, st[k].comp)
    ex = cyc.nozzle_result.exit
    pts["9"] = (ex.T, ex.P, st["9"].comp)
    xy = {k: (s_of(T, P, c), T) for k, (T, P, c) in pts.items()}

    fig, ax = plt.subplots(figsize=(8, 5.6), constrained_layout=True)

    def path_(keys, color, label):
        # A polytropic process with constant cp is a straight line in (s, ln T)
        # (ds = cp (1 - e*) dT/T); interpolate that way between computed states.
        xs, ys = [], []
        for a, b in zip(keys[:-1], keys[1:]):
            t = np.linspace(0, 1, 30)
            xs += list(xy[a][0] + t * (xy[b][0] - xy[a][0]))
            ys += list(np.exp(np.log(xy[a][1]) + t * (np.log(xy[b][1]) - np.log(xy[a][1]))))
        ax.plot(xs, ys, color=color, label=label)

    path_(["0", "2", "25", "3", "4", "45", "5"], SERIES[0], "core: 0-2-25-3-4-45-5")
    path_(["2", "13"], SERIES[1], "bypass: 2-13(=16)")
    path_(["6A", "8", "9"], SERIES[2], "mixed exhaust: 6A-8-9")
    # Mixing is irreversible, not a quasi-static path: show it as straight dashed joins
    for k, c in (("5", SERIES[0]), ("16", SERIES[1])):
        ax.plot([xy[k][0], xy["6A"][0]], [xy[k][1], xy["6A"][1]], color=c, lw=1.2, ls="--")
    # constant-pressure lines for orientation: ambient and burner
    for P, lab, Tr in ((amb.P, "P0", (200, 700)), (st["4"].Pt, "Pt4", (700, 1950))):
        Ts = np.linspace(*Tr, 100)
        comp = air if lab == "P0" else st["4"].comp
        ax.plot([s_of(T, P, comp) for T in Ts], Ts, color=MUTED, lw=0.8, ls=":")
        ax.annotate(lab, (s_of(Tr[1], P, comp), Tr[1]), fontsize=8, color=MUTED, xytext=(3, 0),
                    textcoords="offset points")
    for k, (x, y) in xy.items():
        if k in ("16", "25", "8"):
            continue
        ax.plot([x], [y], "o", ms=6, color=INK, mec=SURFACE, mew=1.5, zorder=5)
        ax.annotate(k, (x, y), xytext=(6, -3), textcoords="offset points", fontsize=9, color=INK)
    ax.set(title=f"T-s diagram, design point (M0 = {cyc.M0:g}, {amb.altitude / 1000:.1f} km)",
           xlabel="Specific entropy s - s0  [J/(kg K)]", ylabel="Temperature  [K]  (stagnation; static at 0, 9)")
    ax.set_ylim(150, max(2000, st["4"].Tt * 1.05))
    ax.plot([], [], color=MUTED, lw=1.2, ls="--", label="mixing (5, 16 -> 6A), irreversible")
    ax.legend(loc="upper left")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def carpet(results, family_a, family_b, path, title, xlabel, ylabel, f_min_lines=(), s_max=None,
           label_a="alpha", label_b="pi_c", f_min_labels=()):
    """Carpet plot of TSFC vs specific thrust.

    results[(a, b)] = (specific_thrust, tsfc) in plot units.  Lines of
    constant `a` in series colour 1, constant `b` in series colour 2,
    each family directly labelled at one end.
    """
    _style()
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for a in family_a:
        xy = [results[(a, b)] for b in family_b if (a, b) in results]
        if len(xy) < 2:
            continue
        x, y = zip(*xy)
        ax.plot(x, y, color=SERIES[0], lw=1.5)
        if isinstance(a, str):   # e.g. the optimum-alpha line: label its other end
            ax.annotate(a, (x[0], y[0]), xytext=(4, 2), textcoords="offset points", fontsize=8, color=INK2)
        else:
            ax.annotate(f"{label_a}={a:g}", (x[-1], y[-1]), xytext=(4, -8),
                        textcoords="offset points", fontsize=8, color=INK2)
    for b in family_b:
        xy = [results[(a, b)] for a in family_a if (a, b) in results]
        if len(xy) < 2:
            continue
        x, y = zip(*xy)
        ax.plot(x, y, color=SERIES[1], lw=1.0)
        ax.annotate(f"{label_b}={b:g}", (x[0], y[0]), xytext=(5, 0), textcoords="offset points",
                    fontsize=8, color=INK2, ha="left", va="center")
    for fm, lab in zip(f_min_lines, f_min_labels or [""] * len(f_min_lines)):
        ax.axvline(fm, color=INK2, lw=1.0, ls="--")
        ax.annotate(lab, (fm, 0), xycoords=("data", "axes fraction"), xytext=(3, 6),
                    textcoords="offset points", fontsize=8, color=INK2, rotation=90, va="bottom")
    if s_max is not None:
        ax.axhline(s_max, color=INK2, lw=1.0, ls="--")
        ax.annotate("S_max", (1, s_max), xycoords=("axes fraction", "data"), xytext=(-4, 3),
                    textcoords="offset points", fontsize=8, color=INK2, ha="right")
    ax.plot([], [], color=SERIES[0], label=f"constant {label_a}")
    ax.plot([], [], color=SERIES[1], lw=1.0, label=f"constant {label_b}")
    ax.legend(loc="upper center")
    ax.margins(x=0.08)
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    fig.savefig(path, dpi=150)
    plt.close(fig)
