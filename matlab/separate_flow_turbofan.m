function r = separate_flow_turbofan(in)
%SEPARATE_FLOW_TURBOFAN  Mattingly ch. 7 parametric cycle analysis of a
%turbofan with separate exhaust streams and losses ("PARA" equivalent).
%
%   r = separate_flow_turbofan(in)
%
%   in is a struct with fields (any consistent units; SI shown):
%     M0, T0 [K], Tt4 [K], pi_c, pi_f, alpha
%     gamma_c, cp_c [J/kg/K], gamma_t, cp_t [J/kg/K], h_PR [J/kg]
%     pi_d_max, pi_b, pi_n, pi_fn, e_c, e_f, e_t, eta_b, eta_m
%     P0_P9, P0_P19   (1 = fully expanded nozzles)
%     gc              (1 for SI; 32.174 for lbm-lbf-ft-s units)
%     J               (1 for SI; 778.16 ft lbf/Btu for US units; only used
%                      for the efficiencies)
%
%   r returns every intermediate ratio plus
%     r.F_mdot   uninstalled specific thrust   [N/(kg/s)] or [lbf/(lbm/s)]
%     r.S        thrust specific fuel consumption  [kg/(N s)] or [lbm/(lbf s)]
%
%   The step numbers match docs/hand_calcs.md.

% Step 1 - gas constants: R = cp (gamma-1)/gamma, from cp - cv = R
r.R_c = (in.gamma_c - 1)/in.gamma_c*in.cp_c;
r.R_t = (in.gamma_t - 1)/in.gamma_t*in.cp_t;
% speed of sound in the free stream.  US units: R in Btu/(lbm R) needs
% J = 778.16 ft lbf/Btu and gc = 32.174 lbm ft/(lbf s^2) to give ft/s
J = 1;  if isfield(in, 'J'), J = in.J; end
r.a0 = sqrt(in.gamma_c*r.R_c*in.gc*J*in.T0);

% Step 2 - free stream ("ram") ratios and inlet recovery
gc_ = in.gamma_c;  gt = in.gamma_t;
r.tau_r = 1 + (gc_ - 1)/2*in.M0^2;            % Tt0/T0 (energy equation)
r.pi_r  = r.tau_r^(gc_/(gc_ - 1));            % Pt0/P0 (isentropic)
r.eta_r = ram_recovery(in.M0);                % MIL-E-5008B shock loss
r.pi_d  = in.pi_d_max*r.eta_r;                % diffuser total-pressure ratio

% Step 3 - burner exit enthalpy over ambient enthalpy
r.tau_lambda = in.cp_t*in.Tt4/(in.cp_c*in.T0);

% Step 4 - compressor and fan: polytropic efficiency e -> tau = pi^((g-1)/(g e))
r.tau_c = in.pi_c^((gc_ - 1)/(gc_*in.e_c));
r.eta_c = (in.pi_c^((gc_ - 1)/gc_) - 1)/(r.tau_c - 1);
r.tau_f = in.pi_f^((gc_ - 1)/(gc_*in.e_f));
r.eta_f = (in.pi_f^((gc_ - 1)/gc_) - 1)/(r.tau_f - 1);

% Step 5 - burner energy balance: cp_c Tt3 + f eta_b h_PR = (1+f) cp_t Tt4
r.f = (r.tau_lambda - r.tau_r*r.tau_c)/(in.eta_b*in.h_PR/(in.cp_c*in.T0) - r.tau_lambda);

% Step 6 - turbine from the power balance (turbine drives HPC + fan)
%   eta_m (1+f) cp_t Tt4 (1 - tau_t) = cp_c T0 tau_r [(tau_c - 1) + alpha (tau_f - 1)]
r.tau_t = 1 - 1/(in.eta_m*(1 + r.f))*r.tau_r/r.tau_lambda*(r.tau_c - 1 + in.alpha*(r.tau_f - 1));
r.pi_t  = r.tau_t^(gt/((gt - 1)*in.e_t));
r.eta_t = (1 - r.tau_t)/(1 - r.tau_t^(1/in.e_t));

% Step 7 - core nozzle: chain total-pressure ratios from P0 to Pt9
r.Pt9_P9 = in.P0_P9*r.pi_r*r.pi_d*in.pi_c*in.pi_b*r.pi_t*in.pi_n;
r.M9     = exit_mach(r.Pt9_P9, gt);
r.T9_T0  = r.tau_lambda*r.tau_t/r.Pt9_P9^((gt - 1)/gt)*in.cp_c/in.cp_t;
r.V9_a0  = r.M9*sqrt(gt*r.R_t*r.T9_T0/(gc_*r.R_c));

% Step 8 - fan nozzle
r.Pt19_P19 = in.P0_P19*r.pi_r*r.pi_d*in.pi_f*in.pi_fn;
r.M19      = exit_mach(r.Pt19_P19, gc_);
r.T19_T0   = r.tau_r*r.tau_f/r.Pt19_P19^((gc_ - 1)/gc_);
r.V19_a0   = r.M19*sqrt(r.T19_T0);

% Step 9 - thrust per unit TOTAL air flow, and TSFC
core = (1 + r.f)*r.V9_a0 - in.M0 + (1 + r.f)*r.R_t/r.R_c*r.T9_T0/r.V9_a0*(1 - in.P0_P9)/gc_;
fan  = r.V19_a0 - in.M0 + r.T19_T0/r.V19_a0*(1 - in.P0_P19)/gc_;
r.F_mdot = r.a0/in.gc/(1 + in.alpha)*core + r.a0/in.gc*in.alpha/(1 + in.alpha)*fan;
r.S      = r.f/((1 + in.alpha)*r.F_mdot);
r.FR     = core/(in.alpha*fan);                     % core thrust / fan thrust

% Step 10 - efficiencies (per unit core air flow)
dKE     = r.a0^2/(2*in.gc)*((1 + r.f)*r.V9_a0^2 + in.alpha*r.V19_a0^2 - (1 + in.alpha)*in.M0^2);
r.eta_T = dKE/(r.f*in.h_PR*J);                     % jet kinetic-energy gain / fuel energy
r.eta_P = (1 + in.alpha)*r.F_mdot*in.M0*r.a0/dKE;  % thrust power / jet kinetic-energy gain
r.eta_O = r.eta_P*r.eta_T;
end

function M = exit_mach(Pt_P, g)
% isentropic Pt/P = (1 + (g-1)/2 M^2)^(g/(g-1)) solved for M
if Pt_P < 1
    error('turbofan:nojet', 'nozzle pressure ratio %.4f < 1: turbine cannot drive this much bypass flow', Pt_P);
end
M = sqrt(2/(g - 1)*(Pt_P^((g - 1)/g) - 1));
end
