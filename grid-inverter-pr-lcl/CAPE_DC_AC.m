%% CAPE DC/AC Control Design - Grid-connected inverter with LCL filter
%  Master in Electronic Engineering - UPC Barcelona
%  Authors: Shahan Antabian, Ayush Mohapatra
%  -------------------------------------------------------------------
%  This script answers questions 1 to 5 of the delivery.
%  Run section by section (Ctrl+Enter) or all at once (F5).
%  -------------------------------------------------------------------

clear; clc; close all;

%% Filter parameters
L1 = 2e-3;     % H, inverter-side inductor
C  = 10e-6;    % F, filter capacitor
L2 = 1e-3;     % H, grid-side inductor

s = tf('s');

%% =================================================================
%  Q1 - Uncompensated open-loop transfer function and Bode plot
%  =================================================================
%
%  For a non-losses LCL filter, the TF from the inverter voltage Vi
%  to the grid-side current I2 (the controlled variable) is:
%
%             I2(s)              1
%      G(s) = ----- = ---------------------------
%             Vi(s)    s^3*L1*L2*C + s*(L1+L2)
%
G = 1 / ( s^3*L1*L2*C + s*(L1+L2) );

% Resonant frequency (analytical)
wres = sqrt( (L1+L2) / (L1*L2*C) );    % rad/s
fres = wres/(2*pi);                    % Hz
fprintf('Q1) wres = %.2f rad/s  ->  fres = %.2f Hz\n', wres, fres);

% Bode plot
figure('Name','Q1: Uncompensated LCL');
bode(G, {1e1, 1e6}); grid on;
title('Bode Diagram - Uncompensated LCL filter G(s) = I_2(s)/V_i(s)');

% Peak (MATLAB returns a sampling-grid value; the ideal lossless system
% has an infinite peak at exactly w = wres)
[gpeak, ~] = getPeakGain(G);
fprintf('   Resonant peak (sampled): %.2f dB (theoretically +inf)\n\n', 20*log10(gpeak));

%% =================================================================
%  Q2 - Passive damping: find Rc such that the resonant peak is -20 dB
%  =================================================================
%
%  Adding Rc in series with C gives:
%
%             I2(s)        1 + s*Rc*C
%   G_d(s) = ----- = ----------------------------------------------
%             Vi(s)   s^3*L1*L2*C + s^2*Rc*C*(L1+L2) + s*(L1+L2)
%

% Sweep Rc and look for the resonant peak around wres
Rc_vec = linspace(0.1, 10, 1000);
peak_dB = zeros(size(Rc_vec));
w_eval  = logspace(2, 6, 4000);

for k = 1:length(Rc_vec)
    Rc = Rc_vec(k);
    Gd = (1 + s*Rc*C) / ( s^3*L1*L2*C + s^2*Rc*C*(L1+L2) + s*(L1+L2) );
    [mag,~] = bode(Gd, w_eval);
    mag = squeeze(mag);
    % Look only in the resonance band
    mask = (w_eval > 0.3*wres) & (w_eval < 3*wres);
    peak_dB(k) = 20*log10(max(mag(mask)));
end

% Find Rc that gives ~ -20 dB peak (smallest Rc that satisfies it)
[~, idx] = min(abs(peak_dB - (-20)));
Rc_opt = Rc_vec(idx);
fprintf('Q2) Rc that yields a -20 dB peak: Rc = %.2f ohm\n\n', Rc_opt);

% Plot peak vs Rc
figure('Name','Q2: Peak vs Rc');
plot(Rc_vec, peak_dB, 'b','LineWidth',1.4); hold on; grid on;
yline(-20,'r--','LineWidth',1.2,'Label','-20 dB target');
xline(Rc_opt,'g--','LineWidth',1.2,'Label',sprintf('R_C = %.2f \\Omega', Rc_opt));
xlabel('R_C (\Omega)'); ylabel('Resonant peak (dB)');
title('Resonant peak vs damping resistance R_C');

% Compare Bode: undamped vs damped
Gd_opt = (1 + s*Rc_opt*C) / ( s^3*L1*L2*C + s^2*Rc_opt*C*(L1+L2) + s*(L1+L2) );

figure('Name','Q2: Damped Bode');
bode(G, Gd_opt, {1e1, 1e6}); grid on;
legend('R_C = 0 (undamped)', ...
       sprintf('R_C = %.2f \\Omega (designed, -20 dB peak)', Rc_opt), ...
       'Location','SouthWest');
title('Bode Diagram - LCL filter with passive damping');

%% =================================================================
%  Q3 - Justification: why use a Proportional + Resonant controller
%  =================================================================
%
%  The reference current i_ref(t) is sinusoidal at the grid frequency
%  wg = 2*pi*50 rad/s. A standard PI controller (in the stationary
%  alpha-beta frame) cannot give zero steady-state error for AC
%  references because its infinite gain is located at s = 0 only.
%
%  The PR controller:
%       Gc(s) = kp + ki*s/(s^2 + wg^2)
%
%  places an infinite (unbounded) gain at s = +-j*wg, i.e. exactly
%  at the grid frequency. By the Internal Model Principle, this
%  guarantees zero steady-state error for sinusoidal references and
%  rejects sinusoidal disturbances of the same frequency.
%
%  Additional advantages versus PI in dq:
%   - No Park/Clarke transformation required (works in alpha-beta).
%   - No coupling between d/q axes.
%   - Lower computational cost and simpler implementation.
%
%  (See written justification in the PDF.)

%% =================================================================
%  Q4 - Design of kp
%  =================================================================
%
%  Above the grid frequency and well below the resonance, the PR
%  controller behaves as Gc(s) ~ kp, and the damped LCL plant behaves
%  asymptotically as G_d(s) ~ 1/(s*(L1+L2)).  The loop gain is then:
%
%       |L(jw)| = kp / ( w*(L1+L2) )
%
%  Forcing |L(jw_c)| = 1 at a desired crossover frequency w_c gives:
%
%       kp = w_c * (L1 + L2)
%
%  We choose w_c well below the resonance (a typical choice is
%  w_c = wres/10) so the closed loop is robust against the resonant
%  dynamics and has comfortable phase margin.

wc_target = wres/10;                 % desired crossover (rad/s)
fc_target = wc_target/(2*pi);        % Hz
kp = wc_target * (L1 + L2);
fprintf('Q4) Designed kp = %.3f  (target fc = %.1f Hz, target wc = %.1f rad/s)\n', ...
        kp, fc_target, wc_target);

% Loop gain with the designed kp (PR resonant part has no effect here
% because the grid frequency is far from the crossover region)
L_loop = kp * Gd_opt;

figure('Name','Q4: Loop gain');
margin(L_loop); grid on;
title(sprintf('Loop gain k_p*G_d(s) with k_p = %.2f, R_C = %.2f \\Omega', kp, Rc_opt));

%% =================================================================
%  Q5 - PR controller: which ki satisfies |Gc| >= 50 dB on 48-52 Hz?
%  =================================================================
%
%  Gc(s) = kp + ki*s/(s^2 + wg^2),   kp = 1
%
%  At s = j*w:
%   |Gc(jw)| = sqrt( kp^2 + ( ki*w / (wg^2 - w^2) )^2 )

wg = 2*pi*50;        % grid frequency in rad/s (nominal)
ki_list = [50, 250, 1000, 5000];

% Reproduce the figure (kp = 1)
figure('Name','Q5: PR controller');
hold on; grid on;
f_axis = linspace(45, 55, 4000);
w_axis = 2*pi*f_axis;
colors = lines(4);
for i = 1:length(ki_list)
    ki = ki_list(i);
    Gc = 1 + ki*s/(s^2 + wg^2);
    [magc,~] = bode(Gc, w_axis);
    magc = squeeze(magc);
    plot(f_axis, 20*log10(magc), 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('k_i = %d', ki));
end
yline(50,'k--','LineWidth',1.2,'DisplayName','50 dB threshold');
xline(48,':','Color',[0.5 0.5 0.5]); xline(52,':','Color',[0.5 0.5 0.5]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Bode Diagram - PR controller (k_p = 1)');
legend('Location','NorthEast');
ylim([0 120]);

% Check the worst-case gain (at f = 48 or 52 Hz)
fprintf('\nQ5) Worst-case gains at the edges of [48,52] Hz (kp = 1):\n');
fprintf('    ki   |  |Gc(48Hz)|  |  |Gc(52Hz)|  |  meets >=50 dB?\n');
fprintf('   ------+--------------+---------------+------------------\n');
for ki = ki_list
    g48 = sqrt(1 + (ki*2*pi*48/(wg^2 - (2*pi*48)^2))^2);
    g52 = sqrt(1 + (ki*2*pi*52/(wg^2 - (2*pi*52)^2))^2);
    g48_dB = 20*log10(g48);
    g52_dB = 20*log10(g52);
    meets = "NO ";
    if min(g48_dB, g52_dB) >= 50, meets = "YES"; end
    fprintf('  %5d  |   %6.2f dB  |    %6.2f dB |  %s\n', ki, g48_dB, g52_dB, meets);
end

%% =================================================================
%  Summary
%  =================================================================
fprintf('\n=================== FINAL DESIGN VALUES ===================\n');
fprintf(' Resonant frequency :  wres = %.0f rad/s  =  fres = %.2f Hz\n', wres, fres);
fprintf(' Damping resistor   :  Rc   = %.2f ohm\n', Rc_opt);
fprintf(' Proportional gain  :  kp   = %.2f   (fc ~ fres/10 = %.0f Hz)\n', kp, fc_target);
fprintf(' Best ki of {50,250,1000,5000} for the >=50 dB spec: 5000\n');
fprintf('===========================================================\n');
