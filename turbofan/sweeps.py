"""Flight-envelope and design-space sweeps that feed the plots."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .engine import Engine, CycleResult


@dataclass
class EnvelopePoint:
    altitude: float
    M0: float
    op: CycleResult
    sweep: object                 # nozzle.AreaSweep
    thrust_convergent: float      # net thrust with A9/A8 = 1 (for comparison)

    @property
    def gain_vs_convergent(self) -> float:
        return self.op.thrust / self.thrust_convergent - 1.0


def flight_envelope(engine: Engine, altitudes, mach_numbers, max_mach=None):
    """Converge the engine along each altitude line with continuation (each
    point starts from the previous solution), optimise A9/A8 at every point.

    max_mach: optional {altitude: M_max} to stop lines at a q-limit.
    """
    out = {}
    for alt in altitudes:
        engine._last = (math.log(engine.dp.cycle.pi_f), math.log(engine.dp.cycle.alpha))
        pts = []
        for M in mach_numbers:
            if max_mach and M > max_mach.get(alt, 99) + 1e-9:
                continue
            op, sweep = engine.area_sweep(M, alt)
            op = engine.operate(M, alt, area_ratio=sweep.best_area_ratio)
            conv = engine.operate(M, alt, area_ratio=1.0)
            pts.append(EnvelopePoint(alt, M, op, sweep, conv.thrust))
        out[alt] = pts
    return out
