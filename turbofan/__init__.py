"""Two-spool, mixed-flow, medium-bypass turbofan cycle analysis.

Modules
-------
gas         thermally perfect (cp(T, FAR)) and two-gamma gas models
atmosphere  ISA 0-20 km
flow        compressible-flow relations for a variable-cp gas
components  inlet, fan/compressor, burner, turbines, mixer
nozzle      variable C-D nozzle and A9/A8 optimisation
engine      design point (sizing) and off-design operation
parametric  Mattingly closed-form 'PARA' equations (7.D1 / 7.D2)
aircraft    required thrust and minimum specific thrust
plots       matplotlib figures
"""
from .engine import DesignInputs, design, Engine
from .gas import ThermallyPerfectGas, TwoGammaGas
from .nozzle import NozzleModel

__all__ = ["DesignInputs", "design", "Engine", "ThermallyPerfectGas", "TwoGammaGas", "NozzleModel"]
