%% HP1_HAND_CALC  Design problem 7.D1 (HP-1 airliner) - hand calculation
% Separate-exhaust turbofan with losses, Mattingly ch. 7.
% Every step is written out the way you would do it on paper; the step
% numbers match docs/hand_calcs.md.  Run it next to your paper copy and
% check each printed line against yours.
clear; clc;

%% Step 0 - givens (SI units)
M0    = 0.83;       % flight Mach number
T0    = 216.65;     % K, ambient temperature at 11 km (standard atmosphere)
Tt4   = 1560;       % K, burner exit total temperature
pi_c  = 36;         % compressor pressure ratio (overall)
pi_f  = 1.8;        % fan pressure ratio
alpha = 10;         % bypass ratio
gam_c = 1.4;   cp_c = 1004;     % J/(kg K), cold section
gam_t = 1.3;   cp_t = 1235;     % J/(kg K), hot section
h_PR  = 42800e3;    % J/kg, fuel heating value
eta_m = 0.99;       % mechanical efficiency
% Component technology level 4, Table 6.2 (CHECK against your book):
pi_d_max = 0.995;   % diffuser type A
pi_b  = 0.96;   eta_b = 0.999;     % burner
e_c   = 0.90;   e_f   = 0.89;      % compressor, fan polytropic efficiency
e_t   = 0.89;                      % cooled turbine
pi_n  = 0.995;  pi_fn = 0.995;     % nozzle type D, core and fan
P0_P9 = 1;      P0_P19 = 1;        % assume fully expanded exhausts

fprintf('=== 7.D1 hand calculation ===\n');
%% Step 1 - gas constants and speed of sound
R_c = (gam_c - 1)/gam_c*cp_c;          % cp - cv = R and gamma = cp/cv
R_t = (gam_t - 1)/gam_t*cp_t;
a0  = sqrt(gam_c*R_c*T0);
fprintf('Step 1   R_c = %.2f J/kgK   R_t = %.2f J/kgK   a0 = %.2f m/s\n', R_c, R_t, a0);

%% Step 2 - free stream and inlet
tau_r = 1 + (gam_c - 1)/2*M0^2;        % Tt0/T0, energy equation
pi_r  = tau_r^(gam_c/(gam_c - 1));     % Pt0/P0, isentropic
eta_r = 1;                             % M0 < 1: no inlet shock loss
pi_d  = pi_d_max*eta_r;
fprintf('Step 2   tau_r = %.5f   pi_r = %.5f   eta_r = %.4f   pi_d = %.4f\n', tau_r, pi_r, eta_r, pi_d);

%% Step 3 - burner exit enthalpy over ambient enthalpy
tau_lambda = cp_t*Tt4/(cp_c*T0);
fprintf('Step 3   tau_lambda = %.5f\n', tau_lambda);

%% Step 4 - compressor and fan
tau_c = pi_c^((gam_c - 1)/(gam_c*e_c));
eta_c = (pi_c^((gam_c - 1)/gam_c) - 1)/(tau_c - 1);
tau_f = pi_f^((gam_c - 1)/(gam_c*e_f));
eta_f = (pi_f^((gam_c - 1)/gam_c) - 1)/(tau_f - 1);
fprintf('Step 4   tau_c = %.5f (eta_c = %.4f)   tau_f = %.5f (eta_f = %.4f)\n', tau_c, eta_c, tau_f, eta_f);

%% Step 5 - fuel/air ratio
f = (tau_lambda - tau_r*tau_c)/(eta_b*h_PR/(cp_c*T0) - tau_lambda);
fprintf('Step 5   h_PR/(cp_c T0) = %.3f   f = %.6f\n', h_PR/(cp_c*T0), f);

%% Step 6 - turbine (power balance)
tau_t = 1 - 1/(eta_m*(1 + f))*tau_r/tau_lambda*(tau_c - 1 + alpha*(tau_f - 1));
pi_t  = tau_t^(gam_t/((gam_t - 1)*e_t));
eta_t = (1 - tau_t)/(1 - tau_t^(1/e_t));
fprintf('Step 6   tau_t = %.5f   pi_t = %.6f   eta_t = %.4f\n', tau_t, pi_t, eta_t);

%% Step 7 - core nozzle
Pt9_P9 = P0_P9*pi_r*pi_d*pi_c*pi_b*pi_t*pi_n;
M9     = sqrt(2/(gam_t - 1)*(Pt9_P9^((gam_t - 1)/gam_t) - 1));
T9_T0  = tau_lambda*tau_t/Pt9_P9^((gam_t - 1)/gam_t)*cp_c/cp_t;
V9_a0  = M9*sqrt(gam_t*R_t*T9_T0/(gam_c*R_c));
fprintf('Step 7   Pt9/P9 = %.5f   M9 = %.5f   T9/T0 = %.5f   V9/a0 = %.5f\n', Pt9_P9, M9, T9_T0, V9_a0);

%% Step 8 - fan nozzle
Pt19_P19 = P0_P19*pi_r*pi_d*pi_f*pi_fn;
M19      = sqrt(2/(gam_c - 1)*(Pt19_P19^((gam_c - 1)/gam_c) - 1));
T19_T0   = tau_r*tau_f/Pt19_P19^((gam_c - 1)/gam_c);
V19_a0   = M19*sqrt(T19_T0);
fprintf('Step 8   Pt19/P19 = %.5f   M19 = %.5f   T19/T0 = %.5f   V19/a0 = %.5f\n', Pt19_P19, M19, T19_T0, V19_a0);

%% Step 9 - specific thrust and TSFC
core = (1 + f)*V9_a0 - M0 + (1 + f)*R_t/R_c*T9_T0/V9_a0*(1 - P0_P9)/gam_c;
fan  = V19_a0 - M0 + T19_T0/V19_a0*(1 - P0_P19)/gam_c;
F_mdot = a0/(1 + alpha)*core + a0*alpha/(1 + alpha)*fan;
S = f/((1 + alpha)*F_mdot);
fprintf('Step 9   core term = %.5f   fan term = %.5f\n', core, fan);
fprintf('         F/mdot0 = %.2f N/(kg/s)   S = %.4f mg/(N s)   FR = %.4f\n', F_mdot, S*1e6, core/(alpha*fan));

%% Step 10 - efficiencies
dKE   = a0^2/2*((1 + f)*V9_a0^2 + alpha*V19_a0^2 - (1 + alpha)*M0^2);
eta_T = dKE/(f*h_PR);
eta_P = (1 + alpha)*F_mdot*M0*a0/dKE;
fprintf('Step 10  eta_T = %.4f   eta_P = %.4f   eta_O = %.4f\n', eta_T, eta_P, eta_T*eta_P);

%% Self-check against the reusable function (the "PARA" equivalent)
in = struct('M0',M0,'T0',T0,'Tt4',Tt4,'pi_c',pi_c,'pi_f',pi_f,'alpha',alpha, ...
    'gamma_c',gam_c,'cp_c',cp_c,'gamma_t',gam_t,'cp_t',cp_t,'h_PR',h_PR, ...
    'pi_d_max',pi_d_max,'pi_b',pi_b,'pi_n',pi_n,'pi_fn',pi_fn,'e_c',e_c,'e_f',e_f, ...
    'e_t',e_t,'eta_b',eta_b,'eta_m',eta_m,'P0_P9',P0_P9,'P0_P19',P0_P19,'gc',1);
r = separate_flow_turbofan(in);
assert(abs(r.F_mdot/F_mdot - 1) < 1e-12 && abs(r.S/S - 1) < 1e-12);
% and against the Python package (turbofan.parametric): 143.25 N/(kg/s), 17.946 mg/(N s)
assert(abs(F_mdot - 143.249) < 0.01 && abs(S*1e6 - 17.9458) < 0.001);
fprintf('Self-check passed (matches separate_flow_turbofan.m and the Python package).\n');
