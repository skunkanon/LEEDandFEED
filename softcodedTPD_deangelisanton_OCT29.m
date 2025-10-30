% === Simulated TPD from crystal temperature trace ===
% --- Cap the crystal temperature at 280 K ---

T_cap = 280;                               % target max temperature
idx_cap = find(T_a >= T_cap, 1);           % first index where T_a hits/exceeds 280
if ~isempty(idx_cap)
    T_a(idx_cap:end) = T_cap;              % keep temperature constant afterwards
end



Rgas = 8.314462618;          % J/mol/K

% ---- Kinetic parameters (edit these) ----
E_des_kJmol = 20 * 4.184;            % kJ/mol  (example barrier) 31.1 KCAL ORIGINAL 
nu_pref     = 7.6 * 10^14;          % s^-1    (1st-order prefactor)
order_n     = 1;             
theta0      = 1.0;           % initial coverage (dimensionless ML basis)

% Interpolate crystal temperature vs time
Tfun = @(tt) interp1(t, T_a, tt, 'pchip', 'extrap');

% Coverage ODE: dtheta/dt = -nu * theta^n * exp(-E/RT(t))
E_des = E_des_kJmol * 1e3;   % J/mol
rate_fun = @(tt,th) - nu_pref * (max(th,0)).^order_n .* exp(-E_des ./ (Rgas * Tfun(tt)));

opts_theta = odeset('RelTol',1e-9,'AbsTol',1e-12);
[t_theta, theta] = ode45(rate_fun, [t(1) t(end)], theta0, opts_theta);

% Desorption rate (positive, "MS signal")
T_theta = Tfun(t_theta);
rate = -rate_fun(t_theta, theta);  % s^-1 * (dimensionless coverage units)

% Heating rate beta(t) for reference
beta_t = gradient(T_a, t);   % K/s on the original grid

% ---- Plots ----
thetaACTUAL = trapz(t_theta, rate);



figure(15); clf;
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% (1) T(t) and beta(t)
nexttile(1);
yyaxis left; plot(t, T_a, 'LineWidth', 1.8); ylabel('T_{crystal} (K)');
yyaxis right; plot(t, beta_t, '--', 'LineWidth', 1.4); ylabel('\beta(t) (K s^{-1})');
xlabel('time (s)'); title('Crystal temperature and heating rate'); grid on; legend('T','\beta(t)','Location','best');

% (2–3) Combined: Desorption rate and coverage vs time
nexttile(2);
yyaxis left
h1 = plot(t_theta, rate, 'LineWidth', 1.8);
ylabel('Desorption rate (ML/sec)');
ylim([0, max(rate)*1.05]);

yyaxis right
h2 = plot(t_theta, theta, '--', 'LineWidth', 1.6);
ylabel('\theta (ML)');
ylim([0, max(theta)*1.05]);
ax = gca;
ax.YAxis(2).TickLabelFormat = '%.2f';

xlabel('time (s)');
title('TPD signal and surface coverage vs time');
grid on;

% --- Primary legend for plotted data ---
legend([h1 h2], {'Desorption rate','\theta'}, 'Location','northwest');

% --- Secondary legend for kinetic parameters ---
hold on
kin_text = sprintf(['E_{des} = %.1f kJ/mol\n' ...
                    '\\nu = %.1e s^{-1}\n' ...
                    'n = %d\n' ...
                    '\\theta_0 = %.1f ML'], ...
                    E_des_kJmol, nu_pref, order_n, theta0);
h_dummy = plot(nan, nan, 'w');  % invisible placeholder
legend(h_dummy, kin_text, 'Location','northeast', 'Box','on');

% (4) TPD signal vs temperature (often how TPD is shown)
% nexttile(4);
% plot(T_theta, rate, 'LineWidth', 1.8);
% xlabel('T (K)'); ylabel('desorption rate (arb)'); title('Simulated TPD (vs T)'); grid on;

% Optional: normalize to peak = 1 for MS-like display
%rate = rate ./ max(rate);
