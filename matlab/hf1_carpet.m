%% HF1_CARPET  Design problem 7.D2 (HF-1) - computer part
%  1. sweep alpha = 0.10:0.05:1.00 for pi_c = 16, 18, 20, 22, 24, 28
%  2. minimum uninstalled specific thrust (5 ft^2 inlet)
%  3. carpet plot S vs F/mdot0
%
% AIRCRAFT DATA BELOW ARE PLACEHOLDERS - replace them with your HF-1 values.
clear; clc;

base = struct('M0',1.6,'T0',390,'Tt4',3250,'pi_c',25,'alpha',0.5,'gamma_c',1.4,'cp_c',0.240, ...
    'gamma_t',1.3,'cp_t',0.296,'h_PR',18400,'pi_d_max',0.96,'pi_b',0.96,'pi_n',0.97, ...
    'e_c',0.90,'e_f',0.89,'e_t',0.89,'eta_b',0.999,'eta_m',0.99,'M6',0.4,'pi_M_max',0.96, ...
    'P0_P9',1,'gc',32.174,'J',778.16);

%% Aircraft inputs (PLACEHOLDERS)
W_TO  = 40000;   % lbf takeoff weight
W     = 0.92*W_TO;
S_w   = 375;     % ft^2
CD0   = 0.028;  K1 = 0.20;  K2 = 0;
phi   = 0.05;    % phi_inlet + phi_nozzle
S_max = 1.094;   % (lbm/h)/lbf, from Chapter 1
mdot  = 146.3;   % lbm/s, given in the problem (5 ft^2 inlet)

%% 1. Sweep
alphas = 0.10:0.05:1.00;
pi_cs  = [16 18 20 22 24 28];
F = nan(numel(alphas), numel(pi_cs));  S = F;  PF = F;
for i = 1:numel(alphas)
    for j = 1:numel(pi_cs)
        in = base;  in.alpha = alphas(i);  in.pi_c = pi_cs(j);
        r = mixed_flow_turbofan(in);
        F(i, j) = r.F_mdot;  S(i, j) = r.S*3600;  PF(i, j) = r.pi_f;
    end
end
fprintf('F/mdot0 [lbf/(lbm/s)] / S [(lbm/h)/lbf] / pi_f\n');
fprintf('alpha   '); fprintf('   pi_c=%-12d', pi_cs); fprintf('\n');
for i = 1:numel(alphas)
    fprintf('%5.2f  ', alphas(i));
    fprintf('  %5.2f/%5.3f/%4.2f', [F(i, :); S(i, :); PF(i, :)]);
    fprintf('\n');
end

%% 2. Minimum specific thrust
P0 = 393.1;  gam = 1.4;  M0 = 1.6;             % lbf/ft^2 at 40 kft
q  = 0.5*gam*P0*M0^2;
CL = W/(q*S_w);
D  = q*S_w*(K1*CL^2 + K2*CL + CD0);
F_req = D/(1 - phi)/2;                         % uninstalled, per engine (P_s = 0 at cruise)
F_min = F_req/mdot;
fprintf('\nRequired uninstalled thrust/engine = %.0f lbf -> F/mdot_min = %.2f lbf/(lbm/s) (PLACEHOLDER aircraft data)\n', F_req, F_min);

%% 3. Carpet plot
fig = figure('visible', 'off');  hold on;  grid on;
for i = 1:numel(alphas)
    plot(F(i, :), S(i, :), 'b-', 'LineWidth', 1.2);
    if mod(i - 1, 2) == 0
        text(F(i, end) + 0.3, S(i, end), sprintf('\\alpha=%.2f', alphas(i)), 'FontSize', 7);
    end
end
for j = 1:numel(pi_cs)
    plot(F(:, j), S(:, j), 'r-', 'LineWidth', 0.8);
    text(F(1, j) + 0.3, S(1, j), sprintf('\\pi_c=%d', pi_cs(j)), 'FontSize', 7);
end
yl = ylim;  xl = xlim;
plot([F_min F_min], yl, 'k--');
plot([min(xl(1), F_min) - 2, xl(2)], [S_max S_max], 'k--');
xlabel('Specific thrust F/m_0  [lbf/(lbm/s)]');  ylabel('TSFC  [(lbm/h)/lbf]');
title('HF-1 (7.D2) supercruise M_0 = 1.6, 40 kft, mixed flow');
print(fig, 'carpet_hf1_matlab.png', '-dpng', '-r150');
fprintf('Saved carpet_hf1_matlab.png\n');
