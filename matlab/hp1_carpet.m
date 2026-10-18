%% HP1_CARPET  Design problem 7.D1 (HP-1) - computer part
%  1. sweep pi_c = 20:2:40 for alpha = 4, 6, 8, 10, 12 and the optimum alpha
%  2. minimum uninstalled specific thrust for inlet diameters 2.2 - 3.5 m
%  3. carpet plot S vs F/mdot0 with the S_max and F/mdot_min lines
%
% AIRCRAFT DATA BELOW ARE PLACEHOLDERS - replace them with your Chapter 1
% HP-1 values (weight, wing area, drag polar, S_max) before using the lines.
clear; clc;

%% Engine inputs (same as hp1_hand_calc.m)
base = struct('M0',0.83,'T0',216.65,'Tt4',1560,'pi_c',36,'pi_f',1.8,'alpha',10, ...
    'gamma_c',1.4,'cp_c',1004,'gamma_t',1.3,'cp_t',1235,'h_PR',42800e3, ...
    'pi_d_max',0.995,'pi_b',0.96,'pi_n',0.995,'pi_fn',0.995,'e_c',0.90,'e_f',0.89, ...
    'e_t',0.89,'eta_b',0.999,'eta_m',0.99,'P0_P9',1,'P0_P19',1,'gc',1);
pi_f_assigned = 1.8;            % <-- your assigned fan pressure ratio
base.pi_f = pi_f_assigned;

%% Aircraft inputs (PLACEHOLDERS)
W     = 0.95*1.6e6;   % N, weight at start of cruise
S_w   = 343;          % m^2, wing area
CD0   = 0.016;  K1 = 0.056;  K2 = -0.008;   % drag polar C_D = K1 C_L^2 + K2 C_L + CD0
Ps    = 1.5;          % m/s, required excess power
phi   = 0.02;         % phi_inlet + phi_nozzle
S_max = 17.0;         % mg/(N s), from Chapter 1
D_in  = [2.2 2.5 2.75 3.0 3.25 3.5];   % m, inlet diameters

%% 1. Sweep
pi_cs = 20:2:40;
alphas = [4 6 8 10 12];
F = nan(numel(alphas) + 1, numel(pi_cs));   % last row = optimum alpha
S = F;  a_opt = nan(size(pi_cs));
for j = 1:numel(pi_cs)
    in = base;  in.pi_c = pi_cs(j);
    for i = 1:numel(alphas)
        in.alpha = alphas(i);
        try
            r = separate_flow_turbofan(in);
            F(i, j) = r.F_mdot;  S(i, j) = r.S*1e6;
        catch
            % turbine cannot drive that much bypass air at this pi_c: leave NaN
        end
    end
    % optimum alpha minimises S.  Upper bound: largest alpha with a core jet.
    a_hi = 1;
    while true
        in.alpha = a_hi*1.05;
        try, separate_flow_turbofan(in); a_hi = a_hi*1.05; catch, break; end
    end
    sfun = @(a) tsfc_of(in, a);
    a_opt(j) = fminbnd(sfun, 0.5, a_hi);
    in.alpha = a_opt(j);  r = separate_flow_turbofan(in);
    F(end, j) = r.F_mdot;  S(end, j) = r.S*1e6;
end
fprintf('Optimum bypass ratio at pi_f = %.2f\n', pi_f_assigned);
fprintf('  pi_c   alpha*   F/mdot0 [N/(kg/s)]   S [mg/(N s)]\n');
fprintf('  %4d   %6.3f   %8.2f            %7.3f\n', [pi_cs; a_opt; F(end, :); S(end, :)]);

%% 2. Minimum specific thrust
T0 = 216.65;  P0 = 22632;  R = 287.05;  gam = 1.4;  M0 = 0.83;
rho0 = P0/(R*T0);  a0 = sqrt(gam*R*T0);  V0 = M0*a0;
q  = 0.5*gam*P0*M0^2;                         % dynamic pressure
CL = W/(q*S_w);
D  = q*S_w*(K1*CL^2 + K2*CL + CD0);
T_inst = D + W*Ps/V0;                         % Eq. 1.28: P_s = V (T - D)/W
F_req  = T_inst/(1 - phi)/2;                  % uninstalled, per engine
fprintf('\nRequired uninstalled thrust per engine = %.2f kN (PLACEHOLDER aircraft data)\n', F_req/1e3);
F_min = zeros(size(D_in));
for k = 1:numel(D_in)
    mdot = rho0*pi*D_in(k)^2/4*V0;            % largest free-stream tube the inlet captures
    F_min(k) = F_req/mdot;
    fprintf('  D = %4.2f m   mdot_max = %6.1f kg/s   F/mdot_min = %6.2f N/(kg/s)\n', D_in(k), mdot, F_min(k));
end

%% 3. Carpet plot
fig = figure('visible', 'off');  hold on;  grid on;
for i = 1:size(F, 1)
    if i <= numel(alphas)
        k = find(~isnan(F(i, :)), 1, 'last');
        if isempty(k)
            fprintf('alpha = %g: no feasible pi_c (the core jet has Pt9/P9 < 1) - not plotted\n', alphas(i));
            continue
        end
        plot(F(i, :), S(i, :), 'b-', 'LineWidth', 1.5);
        text(F(i, k) + 1, S(i, k), sprintf('\\alpha=%g', alphas(i)), 'FontSize', 8);
    else
        plot(F(i, :), S(i, :), 'b--', 'LineWidth', 1.5);
        text(F(i, 1) + 1, S(i, 1), '\alpha_{opt}', 'FontSize', 8);
    end
end
for j = 1:numel(pi_cs)
    ok = ~isnan(F(:, j));
    plot(F(ok, j), S(ok, j), 'r-', 'LineWidth', 0.8);
    text(F(1, j) + 1, S(1, j), sprintf('\\pi_c=%d', pi_cs(j)), 'FontSize', 7);
end
yl = ylim;  xl = xlim;
for k = 1:numel(F_min)
    plot([F_min(k) F_min(k)], yl, 'k--');
    text(F_min(k) + 1, yl(1) + 0.05*diff(yl), sprintf('D=%gm', D_in(k)), 'Rotation', 90, 'FontSize', 7);
end
plot([min(xl(1), min(F_min)) - 5, xl(2)], [S_max S_max], 'k--');
xlabel('Specific thrust F/m_0  [N/(kg/s)]');  ylabel('TSFC  [mg/(N s)]');
title(sprintf('HP-1 (7.D1) cruise M_0 = 0.83, 11 km, \\pi_f = %g', pi_f_assigned));
print(fig, 'carpet_hp1_matlab.png', '-dpng', '-r150');
fprintf('\nSaved carpet_hp1_matlab.png\n');
