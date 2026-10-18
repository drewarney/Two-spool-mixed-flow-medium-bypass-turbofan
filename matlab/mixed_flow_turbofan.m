function r = mixed_flow_turbofan(in)
%MIXED_FLOW_TURBOFAN  Mattingly ch. 7 parametric cycle analysis of a
%mixed-flow turbofan with losses, with the fan pressure ratio found from
%the mixer condition Pt16 = Pt6 (the Problem 7.27 system).
%
%   r = mixed_flow_turbofan(in)
%
%   in fields: M0, T0, Tt4, pi_c, alpha, gamma_c, cp_c, gamma_t, cp_t, h_PR,
%              pi_d_max, pi_b, pi_n, e_c, e_f, e_t, eta_b, eta_m,
%              M6, pi_M_max, P0_P9, gc (1 SI, 32.174 US), J (1 SI, 778.16 US)
%   optional:  pi_f  (skip the Pt16 = Pt6 solve and use this value)
%
%   The step numbers match docs/hand_calcs.md (problem 7.D2).

gc_ = in.gamma_c;  gt = in.gamma_t;
% Step 1 - gas constants and speed of sound
r.R_c = (gc_ - 1)/gc_*in.cp_c;
r.R_t = (gt - 1)/gt*in.cp_t;
J = 1;  if isfield(in, 'J'), J = in.J; end
r.a0  = sqrt(gc_*r.R_c*in.gc*J*in.T0);   % J: Btu -> ft lbf in US units
% Step 2 - free stream and inlet
r.tau_r = 1 + (gc_ - 1)/2*in.M0^2;
r.pi_r  = r.tau_r^(gc_/(gc_ - 1));
r.eta_r = ram_recovery(in.M0);
r.pi_d  = in.pi_d_max*r.eta_r;
% Step 3 - burner exit enthalpy ratio
r.tau_lambda = in.cp_t*in.Tt4/(in.cp_c*in.T0);
% Step 4 - compressor (whole core compression at e_c, Mattingly's lumping)
r.tau_c = in.pi_c^((gc_ - 1)/(gc_*in.e_c));
r.eta_c = (in.pi_c^((gc_ - 1)/gc_) - 1)/(r.tau_c - 1);
% Step 5 - burner
r.f = (r.tau_lambda - r.tau_r*r.tau_c)/(in.eta_b*in.h_PR/(in.cp_c*in.T0) - r.tau_lambda);

% Step 6 - fan and turbine together.  Turbine power balance gives tau_t(tau_f);
% the mixer design rule Pt16 = Pt6 says pi_f = pi_c pi_b pi_t.
tau_t_of = @(tau_f) 1 - 1/(in.eta_m*(1 + r.f))*r.tau_r/r.tau_lambda*(r.tau_c - 1 + in.alpha*(tau_f - 1));
pi_t_of  = @(tau_f) tau_t_of(tau_f)^(gt/((gt - 1)*in.e_t));
pi_f_of  = @(tau_f) tau_f^(gc_*in.e_f/(gc_ - 1));
if isfield(in, 'pi_f') && ~isempty(in.pi_f)
    r.tau_f = in.pi_f^((gc_ - 1)/(gc_*in.e_f));
else
    resid = @(tau_f) pi_f_of(tau_f) - in.pi_c*in.pi_b*pi_t_of(tau_f);
    % upper bracket: largest tau_f that leaves the turbine positive work room
    tf_hi = r.tau_c;
    while tau_t_of(tf_hi) <= 0
        tf_hi = 1 + 0.9*(tf_hi - 1);
    end
    r.tau_f = fzero(resid, [1 + 1e-9, tf_hi]);
end
r.pi_f  = pi_f_of(r.tau_f);
r.eta_f = (r.pi_f^((gc_ - 1)/gc_) - 1)/(r.tau_f - 1);
r.tau_t = tau_t_of(r.tau_f);
r.pi_t  = pi_t_of(r.tau_f);
r.eta_t = (1 - r.tau_t)/(1 - r.tau_t^(1/in.e_t));

% Step 7 - mixer
Pt16_Pt6  = r.pi_f/(in.pi_c*in.pi_b*r.pi_t);                  % = 1 at design
Pt6_P6    = (1 + (gt - 1)/2*in.M6^2)^(gt/(gt - 1));
r.M16     = sqrt(2/(gc_ - 1)*((Pt16_Pt6*Pt6_P6)^((gc_ - 1)/gc_) - 1));  % from P16 = P6
r.alpha_p = in.alpha/(1 + r.f);                               % bypass per unit core GAS
r.cp_6A   = (in.cp_t + r.alpha_p*in.cp_c)/(1 + r.alpha_p);    % mass-weighted
r.R_6A    = (r.R_t + r.alpha_p*r.R_c)/(1 + r.alpha_p);
r.gamma_6A = r.cp_6A/(r.cp_6A - r.R_6A);
Tt16_Tt6  = r.tau_r*r.tau_f/(r.tau_lambda*r.tau_t)*in.cp_t/in.cp_c;
r.tau_M   = in.cp_t/r.cp_6A*(1 + r.alpha_p*in.cp_c/in.cp_t*Tt16_Tt6)/(1 + r.alpha_p);   % energy
T16_T6    = Tt16_Tt6*(1 + (gt - 1)/2*in.M6^2)/(1 + (gc_ - 1)/2*r.M16^2);
r.A16_A6  = r.alpha_p*sqrt(r.R_c*gt*T16_T6/(r.R_t*gc_))*in.M6/r.M16;              % continuity
% momentum + mass + energy through the mixer -> Phi of the mixed flow
Phi = @(M, g) M^2*(1 + (g - 1)/2*M^2)/(1 + g*M^2)^2;
r.Phi_6  = Phi(in.M6, gt);
r.Phi_16 = Phi(r.M16, gc_);
r.Phi_6A = ((1 + r.alpha_p)/(1/sqrt(r.Phi_6) + r.alpha_p*sqrt(r.R_c*gt/(r.R_t*gc_)*Tt16_Tt6/r.Phi_16)))^2 ...
           *r.R_6A*gt*r.tau_M/(r.R_t*r.gamma_6A);
g6A = r.gamma_6A;
r.M6A = sqrt(2*r.Phi_6A/(1 - 2*g6A*r.Phi_6A + sqrt(1 - 2*(g6A + 1)*r.Phi_6A)));    % subsonic root
MFP = @(M, g, R) M*sqrt(g/R)*(1 + (g - 1)/2*M^2)^(-(g + 1)/(2*(g - 1)));
r.pi_M_ideal = (1 + r.alpha_p)*sqrt(r.tau_M)/(1 + r.A16_A6)*MFP(in.M6, gt, r.R_t)/MFP(r.M6A, g6A, r.R_6A);
r.pi_M = in.pi_M_max*r.pi_M_ideal;

% Step 8 - nozzle
r.Pt9_P9 = in.P0_P9*r.pi_r*r.pi_d*in.pi_c*in.pi_b*r.pi_t*r.pi_M*in.pi_n;
if r.Pt9_P9 < 1
    error('turbofan:nojet', 'nozzle pressure ratio < 1');
end
r.M9    = sqrt(2/(g6A - 1)*(r.Pt9_P9^((g6A - 1)/g6A) - 1));
r.Tt9_T0 = r.tau_M*r.tau_t*r.tau_lambda*in.cp_c/in.cp_t;
r.T9_T0 = r.Tt9_T0/r.Pt9_P9^((g6A - 1)/g6A);
r.V9_a0 = r.M9*sqrt(g6A*r.R_6A*r.T9_T0/(gc_*r.R_c));

% Step 9 - thrust and TSFC (f0 = fuel per unit TOTAL air)
r.f0     = r.f/(1 + in.alpha);
r.F_mdot = r.a0/in.gc*((1 + r.f0)*r.V9_a0 - in.M0 + (1 + r.f0)*r.R_6A/r.R_c*r.T9_T0/r.V9_a0*(1 - in.P0_P9)/gc_);
r.S      = r.f0/r.F_mdot;
end
