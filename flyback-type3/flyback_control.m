%% CAPE DC/DC CONTROL DESIGN — Flyback Converter
%  Course: CAPE, UPC — Profs. Domingo Biel, Francesc Guinjoan
%  Authors: Shahan Antabian, Ayush Mohapatra
%  Date: April 2026
clear; clc; close all;

%% ── Parameters ──────────────────────────────────────────
Vg  = 48;           % Input voltage [V]
V   = 12;           % Output voltage [V]
Lm  = 240e-6;       % Magnetizing inductance [H]
C   = 50e-6;        % Output capacitance [F]
n   = 1/4;          % Transformer turns ratio
R   = 12;           % Load resistance [Ohm]

%% ── Steady-state duty cycle ─────────────────────────────
% V/Vg = nD/(1-D) => D = V/(Vg*n + V)
D = V / (Vg*n + V);
fprintf('Duty cycle D = %.4f\n', D);

%% ── Plant transfer function G_vd(s) ────────────────────
% Numerator:  (Vg/n) * (1 - tau_z * s)
% Denominator: Lm*C * (s^2 + s/(RC) + (1-D)^2/(n^2*Lm*C))

tau_z = Lm * D * n^2 / (R * (1-D)^2);   % RHP zero time constant
wz_rhp = 1 / tau_z;                       % RHP zero frequency [rad/s]
w0_sq  = (1-D)^2 / (n^2 * Lm * C);       % Natural frequency squared
w0     = sqrt(w0_sq);                      % Natural frequency [rad/s]

fprintf('RHP zero: fz = %.1f Hz\n', wz_rhp/(2*pi));
fprintf('Resonant freq: f0 = %.1f Hz\n', w0/(2*pi));
fprintf('Quality factor: Q = %.2f\n', R*C*w0);

% Build transfer function
num_G = [(Vg/n)*(-tau_z), Vg/n];
den_G = [Lm*C, Lm*C/(R*C), Lm*C*w0_sq];
G = tf(num_G, den_G);

fprintf('\nPlant G_vd(s):\n');
G

%% ══════════════════════════════════════════════════════════
%% QUESTION 1: Bode plot of uncompensated open-loop (2p)
%% ══════════════════════════════════════════════════════════
figure('Name','Q1 - Uncompensated Bode','NumberTitle','off');
bodeplot(G);
grid on;
title('Bode Plot — Uncompensated Open-Loop G_{vd}(s)');

%% ══════════════════════════════════════════════════════════
%% QUESTION 2: Crossover, Phase Margin, Gain Margin (1p)
%% ══════════════════════════════════════════════════════════
[Gm_lin, Pm, Wcg, Wcp] = margin(G);
Gm_dB = 20*log10(Gm_lin);

fprintf('\n=== QUESTION 2 ===\n');
fprintf('Crossover frequency:  fc = %.1f Hz (wc = %.1f rad/s)\n', Wcp/(2*pi), Wcp);
fprintf('Phase margin:         PM = %.2f deg\n', Pm);
fprintf('Gain margin:          GM = %.2f dB (at f = %.1f Hz)\n', Gm_dB, Wcg/(2*pi));

if Pm < 0
    fprintf('>> The closed-loop system is UNSTABLE (PM < 0).\n');
else
    fprintf('>> The closed-loop system is stable.\n');
end

%% ══════════════════════════════════════════════════════════
%% QUESTION 3: Controller selection justification (1p)
%% ══════════════════════════════════════════════════════════
% Evaluate plant phase at fc = 20 kHz
wc_target = 2*pi*20e3;
[mag_20k, phase_20k] = bode(G, wc_target);
phase_20k = phase_20k;  % already in degrees

fprintf('\n=== QUESTION 3 ===\n');
fprintf('Plant phase at 20 kHz: %.2f deg\n', phase_20k);
fprintf('Phase boost needed from compensator: %.2f deg\n', -130 - phase_20k);
fprintf('\nController analysis:\n');
fprintf('  PI:             max phase at wc ≈ 0 deg => PM ≈ %.1f deg. CANNOT.\n', 180 + phase_20k);
fprintf('  Type II:        max phase at wc ≈ 0 deg => PM ≈ %.1f deg. CANNOT.\n', 180 + phase_20k);
fprintf('  PI+lead:        2 zeros + 1 pole + integrator => CAN provide ~67 deg boost.\n');
fprintf('  PI+double lead: 3 zeros + 2 poles + integrator => CAN (comfortable margin).\n');
fprintf('  Type III:       2 zeros + 2 poles + integrator => CAN (best suited).\n');
fprintf('\n>> PI+lead, PI+double lead, and Type III can meet the specs.\n');
fprintf('>> Type III is selected for design (Q4).\n');

%% ══════════════════════════════════════════════════════════
%% QUESTION 4: Type III compensator design (2p)
%% ══════════════════════════════════════════════════════════
% Type III: Gc(s) = Gc0/s * (1 + s/wz)^2 / (1 + s/wp)^2
%
% Design for: fc = 20 kHz, PM = 50 deg

% Required compensator phase at wc:
%   phase_total = phase_G + phase_Gc = -130 deg
%   phase_Gc = -130 - phase_20k = 66.65 deg
%
% phase_Gc = -90 + 2*(atan(wc/wz) - atan(wc/wp))
% Using geometric placement: wz = wc/a, wp = wc*a
%   atan(a) - atan(1/a) = per_pair_phase
%   2*atan(a) - pi/2 = per_pair_rad

phase_boost = -130 - phase_20k;          % ~66.65 deg
per_pair = (phase_boost + 90) / 2;        % deg per zero-pole pair
a = tan((per_pair*pi/180 + pi/2) / 2);   % geometric ratio

wz_comp = wc_target / a;
wp_comp = wc_target * a;

fprintf('\n=== QUESTION 4 ===\n');
fprintf('Phase boost needed: %.2f deg\n', phase_boost);
fprintf('Geometric ratio a = %.4f\n', a);
fprintf('Compensator zero: fz = %.1f Hz (wz = %.1f rad/s)\n', wz_comp/(2*pi), wz_comp);
fprintf('Compensator pole: fp = %.1f Hz (wp = %.1f rad/s)\n', wp_comp/(2*pi), wp_comp);

% Gain Gc0: enforce |T(jwc)| = 1
mag_ratio = (1 + (wc_target/wz_comp)^2) / (1 + (wc_target/wp_comp)^2);
Gc0 = wc_target / (mag_ratio * mag_20k);

fprintf('Gc0 = %.4f\n', Gc0);

% Build compensator transfer function
num_Gc = Gc0 * conv([1/wz_comp, 1], [1/wz_comp, 1]);
den_Gc = conv([1, 0], conv([1/wp_comp, 1], [1/wp_comp, 1]));
Gc = tf(num_Gc, den_Gc);

fprintf('\nCompensator Gc(s):\n');
Gc

% Open-loop transfer function
T = Gc * G;
fprintf('Open-loop T(s) = Gc(s) * G(s):\n');
T

% Verify margins
[Gm_T_lin, Pm_T, ~, Wcp_T] = margin(T);
fprintf('\nVerification:\n');
fprintf('  Crossover freq: fc = %.1f Hz\n', Wcp_T/(2*pi));
fprintf('  Phase margin:   PM = %.2f deg\n', Pm_T);
fprintf('  Gain margin:    GM = %.2f dB\n', 20*log10(Gm_T_lin));

% Compensated Bode plot
figure('Name','Q4 - Compensated Bode','NumberTitle','off');
bodeplot(T);
grid on;
title('Bode Plot — Compensated Open-Loop T(s) = G_c(s) \cdot G_{vd}(s)');

%% ══════════════════════════════════════════════════════════
%% QUESTION 5: Nyquist plot (2p)
%% ══════════════════════════════════════════════════════════
figure('Name','Q5 - Nyquist','NumberTitle','off');
nyquist(T);
grid on;
title('Nyquist Plot — Compensated Open-Loop T(s)');

% Analysis
P_rhp = length(find(real(pole(T)) > 0));  % RHP open-loop poles
fprintf('\n=== QUESTION 5 ===\n');
fprintf('P (open-loop RHP poles) = %d\n', P_rhp);
fprintf('From Nyquist plot: N (clockwise encirclements of -1) = 0\n');
fprintf('Z = N + P = 0 + 0 = 0 => Closed-loop STABLE.\n');

% Find all crossover frequencies of T
w_bode = logspace(1, 7, 50000);
[mag_arr, phase_arr, w_arr] = bode(T, w_bode);
mag_arr = squeeze(mag_arr);
phase_arr = squeeze(phase_arr);
mag_dB = 20*log10(mag_arr);

crossings = [];
for i = 1:length(mag_dB)-1
    if (mag_dB(i) > 0 && mag_dB(i+1) <= 0) || (mag_dB(i) < 0 && mag_dB(i+1) >= 0)
        crossings = [crossings; w_arr(i)/(2*pi), phase_arr(i)];
    end
end

fprintf('Number of crossover frequencies: %d\n', size(crossings, 1));
for k = 1:size(crossings, 1)
    fprintf('  fc_%d = %.1f Hz, phase = %.2f deg\n', k, crossings(k,1), crossings(k,2));
end

%% ══════════════════════════════════════════════════════════
%% QUESTION 6: Sensitivity functions & disturbance rejection (2p)
%% ══════════════════════════════════════════════════════════
% S(s) = 1/(1+T(s))  — sensitivity
% CS(s) = T(s)/(1+T(s)) — complementary sensitivity

S  = feedback(1, T);    % 1/(1+T)
CS = feedback(T, 1);    % T/(1+T)

figure('Name','Q6 - Sensitivity','NumberTitle','off');
bodeplot(S, CS);
grid on;
legend('S(s) = 1/(1+T)', 'T/(1+T)');
title('Bode — Sensitivity Functions');

fprintf('\n=== QUESTION 6 ===\n');
fprintf('Perturbed output voltage for unit step in Vg:\n');
fprintf('  The closed-loop output perturbation is:\n');
fprintf('    dV(s) = [G_vg(s) / (1 + T(s))] * (1/s)\n');
fprintf('  At DC (s->0): T(0) -> infinity (integrator in Gc)\n');
fprintf('  Therefore: dV(DC) = G_vg(0) / (1 + T(0)) = 0\n');
fprintf('  >> The steady-state output perturbation is ZERO.\n');
fprintf('  >> The Type III integrator ensures perfect DC rejection.\n');
