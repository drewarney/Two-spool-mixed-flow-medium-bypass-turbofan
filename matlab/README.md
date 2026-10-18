# MATLAB versions

The design-problem calculations in MATLAB (also run unmodified in GNU Octave). Every script checks itself against the Python package and stops with an error if the two disagree.

| File | What it does |
|---|---|
| `hp1_hand_calc.m` | 7.D1 (HP-1, separate flow) hand calculation, one printed line per step. The steps match [`docs/hand_calcs.md`](../docs/hand_calcs.md). |
| `hf1_hand_calc.m` | 7.D2 (HF-1, mixed flow), including the P_t16 = P_t6 iteration for π_f and the mixer. |
| `hp1_carpet.m` | 7.D1 computer part: π_c sweep, optimum α (`fminbnd`), minimum specific thrust per inlet diameter, carpet plot. |
| `hf1_carpet.m` | 7.D2 computer part: α × π_c sweep, minimum specific thrust, carpet plot. |
| `separate_flow_turbofan.m` | Reusable function: Mattingly ch. 7 separate-flow turbofan with losses ("PARA"). |
| `mixed_flow_turbofan.m` | Reusable function: mixed-flow turbofan with losses (solves π_f with `fzero`). |
| `ram_recovery.m`, `tsfc_of.m` | Small helpers. |

```matlab
cd matlab
hp1_hand_calc      % prints Steps 1-10, ends with "Self-check passed"
hf1_hand_calc
hp1_carpet         % writes carpet_hp1_matlab.png
hf1_carpet         % writes carpet_hf1_matlab.png
```

**Units.** The functions work in any consistent unit system:
* SI: `gc = 1`, `J = 1`.
* US (lbm, lbf, ft, s, Btu): `gc = 32.174`, `J = 778.16`.

**Before using the carpet-plot lines:** the aircraft data in `hp1_carpet.m` and `hf1_carpet.m` (weight, wing area, drag polar, S_max) are placeholders. Edit the "Aircraft inputs" block at the top of each script.
