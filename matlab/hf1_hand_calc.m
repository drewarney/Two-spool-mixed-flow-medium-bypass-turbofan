%% HF1_HAND_CALC  Design problem 7.D2 (HF-1 fighter) - hand calculation
% Mixed-flow turbofan with losses, Mattingly ch. 7, in the book's US units.
% The fan pressure ratio is NOT given: it comes from the mixer design rule
% Pt16 = Pt6 (the system of equations in Problem 7.27).  Step 6 shows the
% fixed-point iteration you would do on paper, then checks it with fzero.
% Step numbers match docs/hand_calcs.md.
clear; clc;

%% Step 0 - givens (US units)
M0    = 1.6;        % supercruise Mach number
T0    = 390;        % degR, 40 kft
Tt4   = 3250;       % degR
pi_c  = 25;   alpha = 0.5;
gam_c = 1.4;  cp_c = 0.240;     % Btu/(lbm R)
gam_t = 1.3;  cp_t = 0.296;     % Btu/(lbm R)
h_PR  = 18400;      % Btu/lbm
eta_m = 0.99;
M6    = 0.4;        % mixer core-inlet Mach number
pi_M_max = 0.96;    % mixer wall-friction loss
gc = 32.174;        % lbm ft/(lbf s^2)
J  = 778.16;        % ft lbf/Btu
% Component technology level 4, Table 6.2 (CHECK against your book):
pi_d_max = 0.96;    % diffuser type C
pi_b  = 0.96;   eta_b = 0.999;
e_c   = 0.90;   e_f   = 0.89;   e_t = 0.89;   % cooled turbine
pi_n  = 0.97;       % nozzle type F
P0_P9 = 1;

fprintf('=== 7.D2 hand calculation ===\n');
%% Step 1 - gas constants and speed of sound
R_c = (gam_c - 1)/gam_c*cp_c;          % Btu/(lbm R)
R_t = (gam_t - 1)/gam_t*cp_t;
a0  = sqrt(gam_c*R_c*J*gc*T0);         % ft/s
fprintf('Step 1   R_c = %.6f  R_t = %.6f Btu/lbmR   a0 = %.2f ft/s\n', R_c, R_t, a0);

%% Step 2 - free stream and inlet (supersonic: MIL-E-5008B shock recovery)
tau_r = 1 + (gam_c - 1)/2*M0^2;
pi_r  = tau_r^(gam_c/(gam_c - 1));
eta_r = 1 - 0.075*(M0 - 1)^1.35;
pi_d  = pi_d_max*eta_r;
fprintf('Step 2   tau_r = %.5f   pi_r = %.5f   eta_r = %.5f   pi_d = %.5f\n', tau_r, pi_r, eta_r, pi_d);

%% Step 3 - tau_lambda
tau_lambda = cp_t*Tt4/(cp_c*T0);
fprintf('Step 3   tau_lambda = %.5f\n', tau_lambda);

%% Step 4 - compressor
tau_c = pi_c^((gam_c - 1)/(gam_c*e_c));
eta_c = (pi_c^((gam_c - 1)/gam_c) - 1)/(tau_c - 1);
fprintf('Step 4   tau_c = %.5f   eta_c = %.4f\n', tau_c, eta_c);

%% Step 5 - fuel/air ratio
f = (tau_lambda - tau_r*tau_c)/(eta_b*h_PR/(cp_c*T0) - tau_lambda);
fprintf('Step 5   h_PR/(cp_c T0) = %.4f   f = %.6f\n', h_PR/(cp_c*T0), f);

%% Step 6 - fan pressure ratio from Pt16 = Pt6 (Problem 7.27)
% Loop: guess pi_f -> tau_f -> tau_t (power balance) -> pi_t -> new pi_f = pi_c pi_b pi_t
fprintf('Step 6   fixed-point iteration on pi_f:\n');
pi_f = 4.0;                          % any sensible first guess (3 to 5)
for k = 1:8
    tau_f = pi_f^((gam_c - 1)/(gam_c*e_f));
    tau_t = 1 - 1/(eta_m*(1 + f))*tau_r/tau_lambda*(tau_c - 1 + alpha*(tau_f - 1));
    pi_t  = tau_t^(gam_t/((gam_t - 1)*e_t));
    pi_f_new = pi_c*pi_b*pi_t;
    fprintf('         k=%d  pi_f=%.5f  tau_f=%.5f  tau_t=%.5f  pi_t=%.5f  ->  pi_f=%.5f\n', ...
            k, pi_f, tau_f, tau_t, pi_t, pi_f_new);
    if abs(pi_f_new - pi_f) < 5e-5, pi_f = pi_f_new; break; end
    pi_f = pi_f_new;
end
% converged values (one more pass so tau_f, tau_t, pi_t all match the final pi_f)
tau_f = pi_f^((gam_c - 1)/(gam_c*e_f));
tau_t = 1 - 1/(eta_m*(1 + f))*tau_r/tau_lambda*(tau_c - 1 + alpha*(tau_f - 1));
pi_t  = tau_t^(gam_t/((gam_t - 1)*e_t));
eta_f = (pi_f^((gam_c - 1)/gam_c) - 1)/(tau_f - 1);
eta_t = (1 - tau_t)/(1 - tau_t^(1/e_t));
% exact root with fzero, for comparison
resid = @(tf) tf^(gam_c*e_f/(gam_c - 1)) - pi_c*pi_b*(1 - 1/(eta_m*(1 + f))*tau_r/tau_lambda*(tau_c - 1 + alpha*(tf - 1)))^(gam_t/((gam_t - 1)*e_t));
tau_f_exact = fzero(resid, [1.01, 2.0]);
fprintf('         converged: pi_f = %.5f  tau_f = %.5f (fzero: %.6f)  tau_t = %.5f  pi_t = %.5f\n', ...
        pi_f, tau_f, tau_f_exact, tau_t, pi_t);
fprintf('         eta_f = %.4f  eta_t = %.4f\n', eta_f, eta_t);

%% Step 7 - mixer
alpha_p  = alpha/(1 + f);                                  % bypass per lbm core gas
cp_6A    = (cp_t + alpha_p*cp_c)/(1 + alpha_p);
R_6A     = (R_t + alpha_p*R_c)/(1 + alpha_p);
gam_6A   = cp_6A/(cp_6A - R_6A);
Tt16_Tt6 = tau_r*tau_f/(tau_lambda*tau_t)*cp_t/cp_c;
tau_M    = cp_t/cp_6A*(1 + alpha_p*cp_c/cp_t*Tt16_Tt6)/(1 + alpha_p);
Pt6_P6   = (1 + (gam_t - 1)/2*M6^2)^(gam_t/(gam_t - 1));
Pt16_Pt6 = pi_f/(pi_c*pi_b*pi_t);                          % = 1 by construction
M16      = sqrt(2/(gam_c - 1)*((Pt16_Pt6*Pt6_P6)^((gam_c - 1)/gam_c) - 1));
T16_T6   = Tt16_Tt6*(1 + (gam_t - 1)/2*M6^2)/(1 + (gam_c - 1)/2*M16^2);
A16_A6   = alpha_p*sqrt(R_c*gam_t*T16_T6/(R_t*gam_c))*M6/M16;
Phi_6    = M6^2*(1 + (gam_t - 1)/2*M6^2)/(1 + gam_t*M6^2)^2;
Phi_16   = M16^2*(1 + (gam_c - 1)/2*M16^2)/(1 + gam_c*M16^2)^2;
Phi_6A   = ((1 + alpha_p)/(1/sqrt(Phi_6) + alpha_p*sqrt(R_c*gam_t/(R_t*gam_c)*Tt16_Tt6/Phi_16)))^2 ...
           *R_6A*gam_t*tau_M/(R_t*gam_6A);
M6A      = sqrt(2*Phi_6A/(1 - 2*gam_6A*Phi_6A + sqrt(1 - 2*(gam_6A + 1)*Phi_6A)));
MFP_6    = M6*sqrt(gam_t/R_t)*(1 + (gam_t - 1)/2*M6^2)^(-(gam_t + 1)/(2*(gam_t - 1)));
MFP_6A   = M6A*sqrt(gam_6A/R_6A)*(1 + (gam_6A - 1)/2*M6A^2)^(-(gam_6A + 1)/(2*(gam_6A - 1)));
pi_M_ideal = (1 + alpha_p)*sqrt(tau_M)/(1 + A16_A6)*MFP_6/MFP_6A;
pi_M     = pi_M_max*pi_M_ideal;
fprintf('Step 7a  alpha'' = %.5f  cp_6A = %.5f  R_6A = %.6f  gamma_6A = %.5f\n', alpha_p, cp_6A, R_6A, gam_6A);
fprintf('Step 7b  Tt16/Tt6 = %.5f  tau_M = %.5f\n', Tt16_Tt6, tau_M);
fprintf('Step 7c  Pt6/P6 = %.5f  M16 = %.5f  T16/T6 = %.5f  A16/A6 = %.5f\n', Pt6_P6, M16, T16_T6, A16_A6);
fprintf('Step 7d  Phi_6 = %.5f  Phi_16 = %.5f  Phi_6A = %.5f  M6A = %.5f\n', Phi_6, Phi_16, Phi_6A, M6A);
fprintf('Step 7e  MFP_6/MFP_6A = %.5f  pi_M,ideal = %.5f  pi_M = %.5f\n', MFP_6/MFP_6A, pi_M_ideal, pi_M);

%% Step 8 - nozzle
Pt9_P9 = P0_P9*pi_r*pi_d*pi_c*pi_b*pi_t*pi_M*pi_n;
M9     = sqrt(2/(gam_6A - 1)*(Pt9_P9^((gam_6A - 1)/gam_6A) - 1));
Tt9_T0 = tau_M*tau_t*tau_lambda*cp_c/cp_t;
T9_T0  = Tt9_T0/Pt9_P9^((gam_6A - 1)/gam_6A);
V9_a0  = M9*sqrt(gam_6A*R_6A*T9_T0/(gam_c*R_c));
fprintf('Step 8   Pt9/P9 = %.4f  M9 = %.5f  Tt9/T0 = %.5f  T9/T0 = %.5f  V9/a0 = %.5f\n', ...
        Pt9_P9, M9, Tt9_T0, T9_T0, V9_a0);

%% Step 9 - specific thrust and TSFC
f0     = f/(1 + alpha);
F_mdot = a0/gc*((1 + f0)*V9_a0 - M0 + (1 + f0)*R_6A/R_c*T9_T0/V9_a0*(1 - P0_P9)/gam_c);
S      = f0/F_mdot*3600;                                   % (lbm/h)/lbf
fprintf('Step 9   f0 = %.6f   F/mdot0 = %.3f lbf/(lbm/s)   S = %.4f (lbm/h)/lbf\n', f0, F_mdot, S);

%% Step 10 - efficiencies
dKE   = a0^2/(2*gc)*((1 + f0)*V9_a0^2 - M0^2);             % ft lbf per lbm of air
eta_T = dKE/(f0*h_PR*J);
eta_P = F_mdot*M0*a0/dKE;
fprintf('Step 10  eta_T = %.4f   eta_P = %.4f   eta_O = %.4f\n', eta_T, eta_P, eta_T*eta_P);

%% Self-check
in = struct('M0',M0,'T0',T0,'Tt4',Tt4,'pi_c',pi_c,'alpha',alpha,'gamma_c',gam_c,'cp_c',cp_c, ...
    'gamma_t',gam_t,'cp_t',cp_t,'h_PR',h_PR,'pi_d_max',pi_d_max,'pi_b',pi_b,'pi_n',pi_n, ...
    'e_c',e_c,'e_f',e_f,'e_t',e_t,'eta_b',eta_b,'eta_m',eta_m,'M6',M6,'pi_M_max',pi_M_max, ...
    'P0_P9',P0_P9,'gc',gc,'J',J);
r = mixed_flow_turbofan(in);
assert(abs(r.F_mdot/F_mdot - 1) < 1e-4 && abs(r.S*3600/S - 1) < 1e-4);   % 1e-4: hand iteration stopped at 5e-5
% Python package (turbofan.parametric): 66.27 lbf/(lbm/s), 1.1826 (lbm/h)/lbf, pi_f = 4.268
assert(abs(F_mdot - 66.266) < 0.01 && abs(S - 1.1826) < 2e-4 && abs(pi_f - 4.268) < 1e-3);
fprintf('Self-check passed (matches mixed_flow_turbofan.m and the Python package).\n');
