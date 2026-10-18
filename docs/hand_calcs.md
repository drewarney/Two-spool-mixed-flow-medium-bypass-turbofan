# Hand calculations: design problems 7.D1 and 7.D2

This is a step-by-step walkthrough of both hand calculations (Mattingly, *Elements of Propulsion*, ch. 7), written so you can follow it with a pencil and a calculator. Each step gives:

1. the equation,
2. one or two lines on **why** it holds (the part you'll be asked about),
3. the numbers plugged in and the result to check yourself against.

The same steps are coded line for line in [`matlab/hp1_hand_calc.m`](../matlab/hp1_hand_calc.m) and [`matlab/hf1_hand_calc.m`](../matlab/hf1_hand_calc.m), with the same step numbers. If your paper and the script disagree, find the first step where they diverge.

> **Check the component values against your book first.** The level-4 values used here (π_d,max, e_c, e_f, e_t, π_b, η_b, π_n) are the Table 6.2 values as transcribed in `turbofan/tech_levels.py`. If your edition differs, change them in Step 0. Every later number moves with them, but the method doesn't change.

---

## How to set up the page

* **Column 1: givens.** Copy every number from the problem statement and Table 6.2 into one block at the top of the page, with units. Most hand-calc mistakes are a wrong input.
* **Column 2: derived constants.** R_c, R_t, a₀ and h_PR/(c_pc T₀). You use them repeatedly.
* **Then one line per step**: symbol = formula = numbers = result. Box the result.
* **Carry 5 significant figures.** π_t for the turbofan is a small number raised to a large power (τ_t^9.6). Rounding τ_t to 3 figures moves F/ṁ₀ by about 1 %.
* **Work with ratios (τ, π) as long as possible.** Only convert to dimensional numbers (a₀, V₉) at the end. The ratios keep the arithmetic short and make the sanity checks obvious.
* **Run the sanity check after every step** (listed under each step below). Catching an error at Step 4 is much cheaper than at Step 9.

Notation:

| Symbol | Meaning |
|---|---|
| τ | total-temperature ratio across a component |
| π | total-pressure ratio across a component |
| τ_r, π_r | free-stream total/static ratios ("ram") |
| τ_λ | c_pt T_t4 / (c_pc T₀), burner-exit enthalpy over ambient enthalpy |
| α | bypass ratio |
| f | fuel/air ratio per unit **core** air |

---

# Part 1. Problem 7.D1, HP-1 airliner (separate exhausts)

## Step 0. Givens

| Quantity | Value | Source |
|---|---|---|
| M₀ | 0.83 | problem |
| altitude | 11 km → T₀ = 216.65 K | standard atmosphere |
| T_t4 | 1560 K | problem |
| π_c / π_f / α | 36 / 1.8 / 10 | problem |
| γ_c, c_pc | 1.4, 1004 J/(kg·K) | problem |
| γ_t, c_pt | 1.3, 1235 J/(kg·K) | problem |
| h_PR | 42 800 kJ/kg | problem |
| η_m | 0.99 | problem |
| π_d,max | 0.995 | Table 6.2, diffuser type A |
| π_b, η_b | 0.96, 0.999 | Table 6.2 |
| e_c, e_f | 0.90, 0.89 | Table 6.2 |
| e_t | 0.89 | Table 6.2, cooled turbine |
| π_n = π_fn | 0.995 | Table 6.2, nozzle type D |
| P₀/P₉ = P₀/P₁₉ | 1 | assume both nozzles fully expanded |

## Step 1. Gas constants and speed of sound

$$R_c = \frac{\gamma_c - 1}{\gamma_c}\,c_{pc}, \qquad R_t = \frac{\gamma_t - 1}{\gamma_t}\,c_{pt}, \qquad a_0 = \sqrt{\gamma_c R_c T_0}$$

**Why:** c_p − c_v = R and γ = c_p/c_v give R = c_p(γ−1)/γ. a₀ is the speed of sound of the ambient air.

R_c = (0.4/1.4)(1004) = **286.86 J/(kg·K)**. R_t = (0.3/1.3)(1235) = **285.00 J/(kg·K)**.
a₀ = √(1.4 × 286.86 × 216.65) = **294.97 m/s**.

*Check:* a₀ is about 295 m/s at 11 km. ✓

## Step 2. Free stream and inlet

$$\tau_r = 1 + \frac{\gamma_c - 1}{2}M_0^2, \qquad \pi_r = \tau_r^{\gamma_c/(\gamma_c - 1)}, \qquad \pi_d = \pi_{d,max}\,\eta_r$$

**Why:**
* τ_r is the energy equation for the free stream: bringing the air to rest raises T to T_t.
* π_r is the isentropic relation between total and static pressure.
* η_r is the shock loss from MIL-E-5008B. It is 1 for M₀ ≤ 1 (no shocks), so here π_d = π_d,max (wall friction only).

τ_r = 1 + 0.2(0.83)² = **1.13778**. π_r = 1.13778^3.5 = **1.57110**. η_r = 1, so π_d = **0.995**.

## Step 3. τ_λ

$$\tau_\lambda = \frac{c_{pt}\,T_{t4}}{c_{pc}\,T_0}$$

**Why:** this is burner-exit enthalpy over ambient enthalpy. It is the "how hot is the cycle" parameter that every later equation uses.

τ_λ = (1235 × 1560)/(1004 × 216.65) = **8.85726**

## Step 4. Compressor and fan

$$\tau_c = \pi_c^{(\gamma_c - 1)/(\gamma_c e_c)}, \qquad \tau_f = \pi_f^{(\gamma_c - 1)/(\gamma_c e_f)}$$

**Why:** polytropic efficiency e means every small stage has efficiency e. Integrating c_p dT/T = (R/e) dP/P gives the exponent (γ−1)/(γe). Dividing by e makes τ larger than the isentropic value: it costs more work to reach the same π.

τ_c = 36^(0.4/1.26) = **3.11935**. τ_f = 1.8^(0.4/1.246) = **1.20767**.
Optional isentropic efficiencies: η_c = (π_c^0.2857 − 1)/(τ_c − 1) = **0.8417** and η_f = **0.8805**.

*Check:* η_c < e_c and η_f < e_f, and the gap is bigger for the higher pressure ratio. ✓

## Step 5. Fuel/air ratio

$$f = \frac{\tau_\lambda - \tau_r\tau_c}{\dfrac{\eta_b h_{PR}}{c_{pc}T_0} - \tau_\lambda}$$

**Why:** this is the first law across the burner per unit air mass, c_pc T_t3 + f η_b h_PR = (1+f) c_pt T_t4, divided by c_pc T₀. Note that T_t3/T₀ = τ_r τ_c.

h_PR/(c_pc T₀) = 42 800/(1.004 × 216.65) = 196.77. With η_b: 0.999 × 196.77 = 196.57.
f = (8.85726 − 1.13778 × 3.11935)/(196.57 − 8.85726) = 5.30814/187.71 = **0.028278**

*Check:* 0.02–0.035 is normal for this T_t4. ✓

## Step 6. Turbine (power balance)

$$\tau_t = 1 - \frac{1}{\eta_m(1+f)}\,\frac{\tau_r}{\tau_\lambda}\Big[(\tau_c - 1) + \alpha(\tau_f - 1)\Big], \qquad \pi_t = \tau_t^{\gamma_t/[(\gamma_t - 1)e_t]}$$

**Why:**
* The turbine drives the compressor (core air) and the fan (core + bypass air): η_m(1+f) c_pt T_t4(1−τ_t) = c_pc T₀ τ_r[(τ_c−1) + α(τ_f−1)].
* For expansion the polytropic efficiency multiplies the exponent (dh = e dh_s), so π_t = τ_t^(γ/((γ−1)e)).

Bracket: (3.11935 − 1) + 10(1.20767 − 1) = 2.11935 + 2.07670 = 4.19605.
τ_t = 1 − (1/(0.99 × 1.028278)) × (1.13778/8.85726) × 4.19605 = 1 − 0.98233 × 0.128457 × 4.19605 = **0.47051**
π_t = 0.47051^(1.3/(0.3 × 0.89)) = 0.47051^4.8689 = **0.025455**

*Check:* this is where most errors happen. The exponent is about 4.87, so carry τ_t to 5 figures. π_t ≈ 0.025 is a very small pressure ratio; the turbine extracts most of the core's pressure to drive a 10:1 bypass fan. Optional: η_t = (1 − τ_t)/(1 − τ_t^(1/e_t)) = 0.9267 > e_t. Turbines gain from "reheat" the way compressors lose. ✓

## Step 7. Core nozzle

$$\frac{P_{t9}}{P_9} = \frac{P_0}{P_9}\,\pi_r\,\pi_d\,\pi_c\,\pi_b\,\pi_t\,\pi_n, \qquad M_9 = \sqrt{\frac{2}{\gamma_t - 1}\Big[\big(P_{t9}/P_9\big)^{(\gamma_t - 1)/\gamma_t} - 1\Big]}$$

$$\frac{T_9}{T_0} = \frac{\tau_\lambda\,\tau_t}{\big(P_{t9}/P_9\big)^{(\gamma_t - 1)/\gamma_t}}\,\frac{c_{pc}}{c_{pt}}, \qquad \frac{V_9}{a_0} = M_9\sqrt{\frac{\gamma_t R_t T_9}{\gamma_c R_c T_0}}$$

**Why:**
* Follow the total pressure from P₀ through every component to station 9.
* T_t9/T₀ = T_t4 τ_t/T₀ = τ_λ τ_t c_pc/c_pt, because nozzles are adiabatic so T_t9 = T_t5.
* Expanding isentropically from total to static gives T₉.

P_t9/P₉ = 1.5711 × 0.995 × 36 × 0.96 × 0.025455 × 0.995 = **1.36833**
M₉ = **0.70734** (subsonic, so a convergent nozzle is fine).
T₉/T₀ = 8.85726 × 0.47051/1.36833^0.23077 × (1004/1235) = **3.15142**
V₉/a₀ = 0.70734 × √(1.3 × 285 × 3.15142/(1.4 × 286.86)) = **1.20609** (V₉ ≈ 356 m/s)

## Step 8. Fan nozzle

$$\frac{P_{t19}}{P_{19}} = \frac{P_0}{P_{19}}\,\pi_r\,\pi_d\,\pi_f\,\pi_{fn}, \qquad \frac{T_{19}}{T_0} = \frac{\tau_r\,\tau_f}{\big(P_{t19}/P_{19}\big)^{(\gamma_c - 1)/\gamma_c}}, \qquad \frac{V_{19}}{a_0} = M_{19}\sqrt{T_{19}/T_0}$$

P_t19/P₁₉ = 1.5711 × 0.995 × 1.8 × 0.995 = **2.79977**. M₁₉ = **1.30765**.
T₁₉/T₀ = 1.13778 × 1.20767/2.79977^0.28571 = **1.02390**. V₁₉/a₀ = **1.32319** (V₁₉ ≈ 390 m/s).

*Note:* M₁₉ > 1 because P_t19/P₁₉ = 2.8 is above the critical 1.893. Full expansion (P₁₉ = P₀) would need a convergent-divergent fan nozzle. A real convergent fan nozzle would choke and leave P₁₉ > P₀, with a small pressure-thrust term instead. P₀/P₁₉ = 1 is an assumption of the problem, not something the hardware gives for free.

## Step 9. Specific thrust and TSFC

$$\frac{F}{\dot m_0} = \frac{a_0}{1+\alpha}\Big[(1+f)\frac{V_9}{a_0} - M_0 + (1+f)\frac{R_t}{R_c}\frac{T_9/T_0}{V_9/a_0}\frac{1 - P_0/P_9}{\gamma_c}\Big] + \frac{a_0\,\alpha}{1+\alpha}\Big[\frac{V_{19}}{a_0} - M_0 + \frac{T_{19}/T_0}{V_{19}/a_0}\frac{1 - P_0/P_{19}}{\gamma_c}\Big]$$

$$S = \frac{f}{(1+\alpha)\,F/\dot m_0}$$

**Why:**
* This is momentum out minus momentum in, plus pressure thrust, for each stream.
* It is divided by the **total** air flow (1+α)ṁ_c, which is why the 1/(1+α) and α/(1+α) weights appear.
* The pressure terms vanish because P₀/P₉ = P₀/P₁₉ = 1.

Core bracket: 1.028278 × 1.20609 − 0.83 = **0.41020**. Fan bracket: 1.32319 − 0.83 = **0.49319**.
F/ṁ₀ = (294.97/11)(0.41020) + (294.97 × 10/11)(0.49319) = 11.00 + 132.25 = **143.25 N/(kg/s)**
S = 0.028278/(11 × 143.25) = 1.7946 × 10⁻⁵ kg/(N·s) = **17.95 mg/(N·s)**

*Check:* the fan stream produces about 92 % of the thrust (thrust ratio FR = 0.083). That is what a 10:1 bypass ratio should do. ✓

## Step 10. Efficiencies (optional, but good to show)

$$\eta_T = \frac{\tfrac{a_0^2}{2}\big[(1+f)(V_9/a_0)^2 + \alpha(V_{19}/a_0)^2 - (1+\alpha)M_0^2\big]}{f\,h_{PR}}, \qquad \eta_P = \frac{(1+\alpha)\,(F/\dot m_0)\,V_0}{\text{same numerator}}$$

η_T = **0.4107**, η_P = **0.7761**, η_O = η_T η_P = **0.3187**.

## 7.D1 answer

| | Value |
|---|---|
| **F/ṁ₀** | **143.25 N/(kg/s)** |
| **S** | **17.95 mg/(N·s)** (= 0.0646 kg/(N·h)) |

`python scripts/hp1_design_problem.py` and `matlab/hp1_carpet.m` act as PARA for the comparison.

---

# Part 2. Problem 7.D2, HF-1 fighter (mixed flow)

This one is done in the book's US units. Two things are new compared with 7.D1:

1. **π_f is not given.** It is set by the mixer design rule P_t16 = P_t6 (Problem 7.27). You iterate for it.
2. **There is a mixer.** It gets its own set of steps (7a–7e).

## Step 0. Givens

| Quantity | Value | Source |
|---|---|---|
| M₀ | 1.6 | problem |
| T₀ | 390 °R (40 kft) | problem / atmosphere table |
| T_t4 | 3250 °R | problem |
| π_c / α | 25 / 0.5 | problem |
| γ_c, c_pc | 1.4, 0.240 Btu/(lbm·°R) | problem |
| γ_t, c_pt | 1.3, 0.296 Btu/(lbm·°R) | problem |
| h_PR | 18 400 Btu/lbm | problem |
| η_m | 0.99 | problem |
| M₆ | 0.4 | problem |
| π_M,max | 0.96 | problem |
| π_d,max | 0.96 | Table 6.2, diffuser type C |
| π_b, η_b | 0.96, 0.999 | Table 6.2 |
| e_c, e_f, e_t | 0.90, 0.89, 0.89 | Table 6.2, cooled turbine |
| π_n | 0.97 | Table 6.2, nozzle type F |
| P₀/P₉ | 1 | assume full expansion |
| g_c | 32.174 lbm·ft/(lbf·s²) | |
| J | 778.16 ft·lbf/Btu | |

## Step 1. Gas constants and speed of sound

R_c = (0.4/1.4)(0.240) = **0.068571 Btu/(lbm·°R)** (= 53.36 ft·lbf/(lbm·°R)).
R_t = (0.3/1.3)(0.296) = **0.068308 Btu/(lbm·°R)**.

$$a_0 = \sqrt{\gamma_c R_c\,J\,g_c\,T_0} = \sqrt{1.4 \times 0.068571 \times 778.16 \times 32.174 \times 390} = 968.18\ \text{ft/s}$$

*Units trap:* R in Btu needs J to get to ft·lbf, and g_c to get from lbm to slugs.

## Step 2. Free stream and inlet (supersonic now)

τ_r = 1 + 0.2(1.6)² = **1.512**. π_r = 1.512^3.5 = **4.25041**.

$$\eta_r = 1 - 0.075(M_0 - 1)^{1.35} = 1 - 0.075(0.6)^{1.35} = 0.96237$$

π_d = 0.96 × 0.96237 = **0.92387**

**Why:** above Mach 1 the inlet shock system loses total pressure. MIL-E-5008B is the standard curve for it. η_r = 0.962 at M 1.6 is a 3.8 % loss before the engine even starts.

## Step 3. τ_λ

τ_λ = (0.296 × 3250)/(0.240 × 390) = **10.27778**

## Step 4. Compressor

τ_c = 25^(0.4/1.26) = **2.77837** (η_c = 0.8482).

**Modelling note:** Mattingly applies e_c to the **whole** core compression (fan core side + HPC). A two-spool model that runs core air through the fan (e_f) and then the HPC (e_c) gives a slightly different τ_c. The Python package does the latter, and its test sets e_f = e_c to compare like with like.

## Step 5. Fuel/air ratio

h_PR/(c_pc T₀) = 18 400/(0.240 × 390) = 196.581.

f = (10.27778 − 1.512 × 2.77837)/(0.999 × 196.581 − 10.27778) = **0.032653**

## Step 6. Fan pressure ratio from P_t16 = P_t6 (Problem 7.27)

The design rule for a mixer is that both streams arrive with the same total pressure:

$$P_{t16} = P_{t6} \;\Rightarrow\; \pi_r\pi_d\pi_f = \pi_r\pi_d\pi_c\pi_b\pi_t \;\Rightarrow\; \pi_f = \pi_c\,\pi_b\,\pi_t$$

But π_t depends on π_f through the power balance (a harder-working fan needs more turbine work), so you iterate:

```
guess π_f
  τ_f = π_f^((γc−1)/(γc·e_f))
  τ_t = 1 − [1/(η_m(1+f))]·(τ_r/τ_λ)·[(τ_c − 1) + α(τ_f − 1)]
  π_t = τ_t^(γt/((γt−1)·e_t))
  π_f,new = π_c·π_b·π_t
repeat until π_f stops changing
```

Starting from π_f = 4.0:

| k | π_f | τ_f | τ_t | π_t | → π_f,new |
|---|---|---|---|---|---|
| 1 | 4.00000 | 1.56055 | 0.70376 | 0.18077 | 4.33842 |
| 2 | 4.33842 | 1.60177 | 0.70079 | 0.17709 | 4.25012 |
| 3 | 4.25012 | 1.59123 | 0.70155 | 0.17802 | 4.27256 |
| 4 | 4.27256 | 1.59393 | 0.70136 | 0.17778 | 4.26682 |
| 5 | 4.26682 | 1.59324 | 0.70141 | 0.17785 | 4.26828 |
| 6 | 4.26828 | 1.59341 | 0.70139 | 0.17783 | 4.26791 |

The iteration oscillates and converges. Averaging the last two values each time speeds it up on paper. Converged:

| | Value |
|---|---|
| **π_f** | **4.2680** |
| τ_f | 1.59338 |
| τ_t | 0.70140 |
| π_t | 0.17783 |
| η_f | 0.8659 |
| η_t | 0.9085 |

*Check:* π_c π_b π_t = 25 × 0.96 × 0.17783 = 4.268 = π_f. ✓

## Step 7. Mixer

**7a. Mixed-gas properties.** Bypass air per unit core *gas* is α' = α/(1+f) = 0.5/1.032653 = **0.48419**. The mixed stream's properties are mass-weighted averages:

$$c_{p6A} = \frac{c_{pt} + \alpha' c_{pc}}{1 + \alpha'} = 0.27773, \qquad R_{6A} = \frac{R_t + \alpha' R_c}{1 + \alpha'} = 0.068394, \qquad \gamma_{6A} = \frac{c_{p6A}}{c_{p6A} - R_{6A}} = 1.32672$$

**7b. Energy → τ_M = T_t6A/T_t6.**

$$\frac{T_{t16}}{T_{t6}} = \frac{\tau_r\tau_f}{\tau_\lambda\tau_t}\,\frac{c_{pt}}{c_{pc}} = 0.41218, \qquad \tau_M = \frac{c_{pt}}{c_{p6A}}\,\frac{1 + \alpha'\,(c_{pc}/c_{pt})\,(T_{t16}/T_{t6})}{1 + \alpha'} = 0.83429$$

**Why:** (1+α') c_p6A T_t6A = c_pt T_t6 + α' c_pc T_t16. Hot core gas and cold bypass air average out.

**7c. Equal static pressure → M₁₆, then area ratio.** The splitter can't sustain a pressure difference, so P₁₆ = P₆. With P_t16 = P_t6 this gives:

$$\frac{P_{t6}}{P_6} = \Big(1 + \frac{\gamma_t - 1}{2}M_6^2\Big)^{\gamma_t/(\gamma_t - 1)} = 1.10824 \;\Rightarrow\; M_{16} = \sqrt{\frac{2}{\gamma_c - 1}\Big[1.10824^{(\gamma_c - 1)/\gamma_c} - 1\Big]} = 0.38599$$

T₁₆/T₆ = (T_t16/T_t6)(1 + 0.15 M₆²)/(1 + 0.2 M₁₆²) = **0.40986**.

$$\frac{A_{16}}{A_6} = \alpha'\sqrt{\frac{R_c\,\gamma_t\,T_{16}}{R_t\,\gamma_c\,T_6}}\;\frac{M_6}{M_{16}} = 0.31014$$

(This is continuity, ṁ = P A M √(γ/RT), at equal P.)

**7d. Momentum → M₆A.** Mattingly's mixer function packs mass, momentum and energy into one number:

$$\Phi(M, \gamma) = \frac{M^2\big(1 + \frac{\gamma - 1}{2}M^2\big)}{(1 + \gamma M^2)^2}$$

Φ₆ = Φ(0.4, 1.3) = **0.11228**. Φ₁₆ = Φ(0.38599, 1.4) = **0.10504**.

$$\Phi_{6A} = \left[\frac{1 + \alpha'}{\dfrac{1}{\sqrt{\Phi_6}} + \alpha'\sqrt{\dfrac{R_c\gamma_t}{R_t\gamma_c}\dfrac{T_{t16}/T_{t6}}{\Phi_{16}}}}\right]^2 \frac{R_{6A}\,\gamma_t\,\tau_M}{R_t\,\gamma_{6A}} = 0.11791$$

$$M_{6A} = \sqrt{\frac{2\Phi_{6A}}{1 - 2\gamma_{6A}\Phi_{6A} + \sqrt{1 - 2(\gamma_{6A} + 1)\Phi_{6A}}}} = 0.41658$$

**Why:** √Φ is (mass flow × √(R T_t/γ)) / (impulse P A(1+γM²)). Impulse is conserved in a constant-area mixer with no wall force, so Φ of the mixed flow follows from the two inlet streams. Take the subsonic root.

**7e. Total-pressure ratio.** From continuity, ṁ = P_t A · MFP(M)/√T_t, with MFP = M√(γ/R)(1 + (γ−1)M²/2)^(−(γ+1)/(2(γ−1))):

$$\pi_{M,ideal} = (1 + \alpha')\sqrt{\tau_M}\;\frac{1}{1 + A_{16}/A_6}\;\frac{MFP(M_6)}{MFP(M_{6A})} = 1.48419 \times 0.91340 \times \frac{1}{1.31014} \times 0.95932 = 0.99265$$

π_M = π_M,max π_M,ideal = 0.96 × 0.99265 = **0.95294**

*Note:* the ideal mixing loss (0.7 %) is the irreversibility of mixing a ~1270 K core stream (T_t6 = τ_t T_t4) with ~520 K bypass air. The mixer still pays off because the mixed jet has a better propulsive efficiency than two separate jets would.

## Step 8. Nozzle

P_t9/P₉ = π_r π_d π_c π_b π_t π_M π_n = 4.25041 × 0.92387 × 25 × 0.96 × 0.17783 × 0.95294 × 0.97 = **15.4919**

M₉ = √(2/0.32672 × (15.4919^(0.32672/1.32672) − 1)) = **2.42884**. That is strongly supersonic, so the nozzle must be convergent-divergent (type F).

T_t9/T₀ = τ_M τ_t τ_λ c_pc/c_pt = **4.87639**. T₉/T₀ = 4.87639/15.4919^0.24626 = **2.48327**.

V₉/a₀ = M₉ √(γ₆A R₆A T₉/(γ_c R_c T₀)) = **3.72112** (V₉ ≈ 3603 ft/s)

## Step 9. Specific thrust and TSFC

The fuel per unit **total** air is f₀ = f/(1+α) = 0.032653/1.5 = **0.021768**.

$$\frac{F}{\dot m_0} = \frac{a_0}{g_c}\Big[(1+f_0)\frac{V_9}{a_0} - M_0 + (1+f_0)\frac{R_{6A}}{R_c}\frac{T_9/T_0}{V_9/a_0}\frac{1 - P_0/P_9}{\gamma_c}\Big]$$

F/ṁ₀ = (968.18/32.174)(1.021768 × 3.72112 − 1.6) = 30.092 × 2.20214 = **66.27 lbf/(lbm/s)**

S = f₀/(F/ṁ₀) = 0.021768/66.27 = 3.2849 × 10⁻⁴ lbm/(lbf·s) × 3600 = **1.1826 (lbm/h)/lbf**

## Step 10. Efficiencies

η_T = **0.5416**, η_P = **0.6081**, η_O = **0.3293**.

## 7.D2 answer

| | Value |
|---|---|
| **π_f** | **4.268** |
| **F/ṁ₀** | **66.27 lbf/(lbm/s)** (650 N/(kg/s)) |
| **S** | **1.183 (lbm/h)/lbf** (33.5 mg/(N·s)) |

---

# Part 3. The computer and aircraft parts of each problem

**Sweeps (the "PARA" part).** These repeat Steps 1–9 in a loop:
* `matlab/hp1_carpet.m`: π_c = 20–40 for each α, plus the TSFC-optimal α found with `fminbnd`.
* `matlab/hf1_carpet.m`: α = 0.10–1.00 for each π_c.

**Minimum uninstalled specific thrust:**
1. Drag from the polar: q = ½γP₀M₀², C_L = W/(qS), C_D = K₁C_L² + K₂C_L + C_D0, D = qSC_D.
2. Required installed thrust (Eq. 1.28): T = D + W·P_s/V₀.
3. Uninstalled thrust per engine: F = T/[(1 − φ_inlet − φ_nozzle) × 2 engines].
4. Largest airflow the inlet can capture: ṁ_max = ρ₀ A V₀.
5. F/ṁ₀,min = F/ṁ_max. A bigger inlet takes more air, so it needs less thrust per unit of air, which allows a higher-bypass, more efficient engine.

> The weights, drag polars and S_max in the scripts are **placeholders**. Replace them with your Chapter 1/2 HP-1 and HF-1 numbers before reading anything off the dashed lines.

**Reading the carpet plot.** Engines to the right of the F/ṁ_min line give enough thrust for that inlet. Engines below the S_max line burn little enough fuel. The promising cycles sit in the region that satisfies both.
* For HP-1 the optimum-α curve sits at nearly constant F/ṁ₀ ≈ 147 N/(kg/s). At the optimum, the core and fan jets leave at nearly the same velocity.
* For HF-1, a higher α and π_c lower S but also lower F/ṁ₀.

---

# Common mistakes

1. **Using γ_c where γ_t belongs.** Everything after the burner (turbine, core nozzle) uses γ_t and c_pt. For the mixed flow, everything after the mixer uses γ₆A.
2. **Rounding τ_t early.** π_t = τ_t^4.87: an error of 0.001 in τ_t moves π_t by about 1 %.
3. **f vs f₀.** In the separate-flow thrust equation the (1+f) terms use f per core air. In the mixed-flow equation f₀ = f/(1+α) is per total air.
4. **Forgetting η_r above Mach 1.** It is only 1 for M₀ ≤ 1.
5. **US units.** a₀ needs both J and g_c, and F/ṁ₀ needs 1/g_c. S comes out in lbm/(lbf·s); multiply by 3600 for (lbm/h)/lbf.
6. **Wrong root for M₆A.** Take the subsonic root; the supersonic one is a different (shocked) flow.
