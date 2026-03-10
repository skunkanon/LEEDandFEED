% === Simulated TPD from crystal temperature trace ===
% --- Cap the crystal temperature at 280 K ---

T_hold = 1200;                         % target hold temperature (K)

T_a_hold = T_a + 900;                        % copy original crystal trace
idx_hold = find(T_a_hold >= T_hold, 1, 'first');

if ~isempty(idx_hold)
    T_a_hold(idx_hold:end) = T_hold;   % hold at 1200 K after first hit
end

Rgas = 8.314462618;          % J/mol/K

% ---- Kinetic parameters (edit these) ----
E_des_kJmol = 160;            % kJ/mol  (example barrier) 31.1 KCAL ORIGINAL 
nu_pref     = 2 * 10^7;          % s^-1    (1st-order prefactor)
order_n     = 1;             
theta0      = 1.0;           % initial coverage (dimensionless ML basis)

% Interpolate crystal temperature vs time
Tfun = @(tt) interp1(t, T_a_hold, tt, 'pchip', 'extrap');

% Coverage ODE: dtheta/dt = -nu * theta^n * exp(-E/RT(t))
E_des = E_des_kJmol * 1e3;   % J/mol
rate_fun = @(tt,th) - nu_pref * (max(th,0)).^order_n .* exp(-E_des ./ (Rgas * Tfun(tt)));

opts_theta = odeset('RelTol',1e-9,'AbsTol',1e-12);
[t_theta, theta] = ode45(rate_fun, [t(1) t(end)], theta0, opts_theta);

% Desorption rate (positive, "MS signal")
T_theta = Tfun(t_theta);
rate = -rate_fun(t_theta, theta);  % s^-1 * (dimensionless coverage units)

% Heating rate beta(t) for reference
beta_t = gradient(T_a_hold, t);   % K/s on the original grid

% ---- Plots ----
thetaACTUAL = trapz(t_theta, rate);



figure(1); clf;
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% (1) T(t) and beta(t)
nexttile(1);
yyaxis left; plot(t, T_a_hold, 'LineWidth', 1.8); ylabel('T_{crystal} (K)');
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



%% === FIT DESCENDING BRANCH ONLY (after rate maximum) ===

% Use the whole simulated/measured trace here
t_fit_all    = t_theta(:);
rate_fit_all = rate(:);
T_fit_all    = Tfun(t_fit_all);

use_sigma = false;
if use_sigma
    sigma_fit_all = sigma_fit_all(:);
else
    sigma_fit_all = [];
end

% settings for peak detection / tail fit
smooth_pts = 9;          % moving-average width for peak detection
skip_pts   = 2;          % skip a couple points after peak before fitting
min_tail_pts = 8;        % require at least this many points in tail

[k_fit, theta_tail_fit, tailfit] = fit_pw_tail_after_peak( ...
    t_fit_all, rate_fit_all, order_n, sigma_fit_all, smooth_pts, skip_pts, min_tail_pts);

% effective temperature for the fitted tail
T_tail_eff = median(Tfun(tailfit.t_tail));

% convert apparent k -> apparent E using assumed nu_pref
E_app_kJmol = (Rgas * T_tail_eff / 1e3) * log(nu_pref / k_fit);

fprintf('\n=== DESCENDING-BRANCH FIT ===\n');
fprintf('peak time          = %.6f s\n', tailfit.t_peak);
fprintf('fit starts at      = %.6f s\n', tailfit.t_tail(1));
fprintf('T_tail_eff         = %.6f K\n', T_tail_eff);
fprintf('order_n            = %.4g\n', order_n);
fprintf('theta_tail_fit     = %.6e ML\n', theta_tail_fit);
%fprintf('k_fit              = %.6e s^-1\n', k_fit);
fprintf('E_app_kJmol        = %.6f kJ/mol  (using nu_pref = %.6e s^-1)\n', ...
        E_app_kJmol, nu_pref);
fprintf('squared residual sum = %.6e\n', tailfit.sse);

% fitted model on the tail
rate_tail_model  = pw_isothermal_rate_closedform(tailfit.t_tail, theta_tail_fit, order_n, k_fit);
theta_tail_model = pw_isothermal_theta_closedform(tailfit.t_tail, theta_tail_fit, order_n, k_fit);

% optional: use apparent E in your printed kinetic text
E_des_kJmol = E_app_kJmol;

%% --- diagnostic plot ---
figure(1); 
%tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile(3);
plot(t_fit_all, rate_fit_all, 'k-', 'LineWidth', 1.2); hold on
plot(t_fit_all, tailfit.r_smooth, 'b--', 'LineWidth', 1.0);
xline(tailfit.t_peak, ':', 'LineWidth', 1.3);
xline(tailfit.t_tail(1), '--', 'LineWidth', 1.3);
plot(tailfit.t_tail, rate_tail_model, 'r-', 'LineWidth', 1.8);
xlabel('time (s)');
ylabel('rate (ML/s)');
title('Peak detection and descending-branch fit');
legend('raw rate','smoothed rate','peak','fit start','PW tail fit','Location','best');


kin_text = sprintf('E_{app} = %.2f kJ/mol\n',  E_app_kJmol);

text(0.98, 0.95, kin_text, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k');

grid on;

function [k_best, theta_tail_best, out] = fit_pw_tail_after_peak(t, r_obs, n, sigma, smooth_pts, skip_pts, min_tail_pts)
% Fit only the descending branch after the rate maximum.
% Parameters fitted:
%   k               effective isothermal rate constant
%   theta_tail_best coverage remaining at start of fitted tail
%
% Inputs:
%   t, r_obs        full time/rate vectors
%   n               fixed desorption order
%   sigma           optional uncertainties, [] for unweighted
%   smooth_pts      smoothing width for peak detection
%   skip_pts        number of points to skip after peak
%   min_tail_pts    minimum number of points required in tail

    t = t(:);
    r_obs = r_obs(:);

    if nargin < 4 || isempty(sigma)
        sigma = ones(size(r_obs));
    else
        sigma = sigma(:);
        sigma = max(sigma, eps);
    end

    if nargin < 5 || isempty(smooth_pts),   smooth_pts = 9; end
    if nargin < 6 || isempty(skip_pts),     skip_pts = 2;   end
    if nargin < 7 || isempty(min_tail_pts), min_tail_pts = 8; end

    % keep only finite points
    good = isfinite(t) & isfinite(r_obs) & isfinite(sigma);
    t = t(good);
    r_obs = r_obs(good);
    sigma = sigma(good);

    if numel(t) < min_tail_pts
        error('Not enough valid data points for tail fitting.');
    end

    % smooth only for peak detection
    smooth_pts = max(3, 2*floor(smooth_pts/2)+1);   % force odd >= 3
    r_smooth = smoothdata(r_obs, 'movmean', smooth_pts);

    % rate maximum
    [~, i_peak] = max(r_smooth);

    % start fit after peak
    i_start = min(i_peak + skip_pts, numel(t) - min_tail_pts + 1);
    if i_start < 1
        i_start = 1;
    end

    t_tail = t(i_start:end);
    r_tail = r_obs(i_start:end);
    s_tail = sigma(i_start:end);

    if numel(t_tail) < min_tail_pts
        error('Tail segment too short after peak detection.');
    end

    % initial guesses
    theta_tail0_guess = max(trapz(t_tail, max(r_tail,0)), 1e-12);
    k0 = max(r_tail(1) / max(theta_tail0_guess.^n, eps), 1e-12);

    % for first-order tail, refine k0 from log-slope if possible
    if abs(n - 1) < 1e-12
        pos = r_tail > 0;
        if nnz(pos) >= 3
            tau0 = t_tail(pos) - t_tail(find(pos,1,'first'));
            p = polyfit(tau0, log(r_tail(pos)), 1);
            if isfinite(p(1)) && (-p(1) > 0)
                k0 = max(-p(1), 1e-12);
                A0 = exp(p(2));
                theta_tail0_guess = max(A0 / k0, 1e-12);
            end
        end
    end

    x0 = log([k0; theta_tail0_guess]);   % fit log(k), log(theta_tail0)

    obj = @(x) pw_tail_objective(x, t_tail, r_tail, n, s_tail);

    opts = optimset('Display','off', ...
                    'TolX',1e-12, ...
                    'TolFun',1e-12, ...
                    'MaxIter',5000, ...
                    'MaxFunEvals',10000);

    xbest = fminsearch(obj, x0, opts);

    k_best = exp(xbest(1));
    theta_tail_best = exp(xbest(2));

    r_fit = pw_isothermal_rate_closedform(t_tail, theta_tail_best, n, k_best);
    resid = (r_fit - r_tail) ./ s_tail;

    out.t_peak   = t(i_peak);
    out.i_peak   = i_peak;
    out.i_start  = i_start;
    out.t_tail   = t_tail;
    out.r_tail   = r_tail;
    out.r_smooth = r_smooth;
    out.r_fit    = r_fit;
    out.residuals = resid;
    out.sse      = sum(resid.^2);
end


function sse = pw_tail_objective(x, t_tail, r_tail, n, s_tail)
    k = exp(x(1));
    theta_tail0 = exp(x(2));

    r_fit = pw_isothermal_rate_closedform(t_tail, theta_tail0, n, k);
    resid = (r_fit - r_tail) ./ s_tail;
    sse = sum(resid.^2);
end


function theta = pw_isothermal_theta_closedform(t_abs, theta0, n, k)
% Closed-form solution of dtheta/dt = -k*theta^n

    t_abs = t_abs(:);
    tau = t_abs - t_abs(1);

    if theta0 <= 0
        theta = zeros(size(tau));
        return
    end

    if abs(n - 1) < 1e-12
        theta = theta0 .* exp(-k .* tau);
    else
        base = theta0.^(1-n) + (n-1) .* k .* tau;
        theta = zeros(size(base));
        pos = base > 0;
        theta(pos) = base(pos).^(1/(1-n));
    end

    theta = max(theta, 0);
end


function r = pw_isothermal_rate_closedform(t_abs, theta0, n, k)
    theta = pw_isothermal_theta_closedform(t_abs, theta0, n, k);
    r = k .* theta.^n;
end