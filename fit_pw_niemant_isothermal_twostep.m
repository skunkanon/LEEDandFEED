function [best1, best2, fit_error, hess1, hess2, cov1, cov2, se1, se2] = ...
    fit_pw_niemant_isothermal_twostep(exp_data, init1, init2, des_order)

% exp_data{i} = {t, r, N0_total, sigma, meta}
% meta.t1 = [t1_start t1_end], meta.T1 = hold temperature (K)
% meta.t2 = [t2_start t2_end], meta.T2 = hold temperature (K)

% --- scaling (same idea as your code; tune if needed)
sf = [1e5, 1e4, 1e7, 1e3];   % [Ea, w, nu, Tc] scaling (Tc not 1e10 here—keep sane)
x01 = init1 ./ sf;
x02 = init2 ./ sf;

% --- bounds (edit to your regime; these are just reasonable placeholders)
lb = [ 50e3, -200e3, 1e8,   0] ./ sf;
ub = [300e3,  200e3, 1e18, 5000] ./ sf;

opts = optimoptions('fmincon', ...
    'Display','iter', ...
    'MaxIterations',300, ...
    'MaxFunctionEvaluations',1200, ...
    'OptimalityTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'FunctionTolerance',1e-12, ...
    'Algorithm','interior-point');

% --- Fit step 1
obj1 = @(x) err_one_step(x, exp_data, 1, sf, des_order);
[xb1, f1, ~, ~, ~, ~, Hx1] = fmincon(obj1, x01, [], [], [], [], lb, ub, [], opts);
best1 = xb1 .* sf;

% --- Fit step 2
obj2 = @(x) err_one_step(x, exp_data, 2, sf, des_order);
[xb2, f2, ~, ~, ~, ~, Hx2] = fmincon(obj2, x02, [], [], [], [], lb, ub, [], opts);
best2 = xb2 .* sf;

fit_error = f1 + f2;

% --- uncertainty estimates (Hessian transform: Hp = S^{-T} Hx S^{-1})
S = diag(sf);

[hess1, cov1, se1] = hess_to_cov(Hx1, S, f1, exp_data, 1);
[hess2, cov2, se2] = hess_to_cov(Hx2, S, f2, exp_data, 2);

% --- optional plot
plot_two_step_fit(exp_data, best1, best2, des_order);

end


function err = err_one_step(x_scaled, exp_data, which_step, sf, n)
p = x_scaled .* sf;
Ea = p(1); w = p(2); nu = p(3); Tc = p(4);

total = 0;

for i = 1:numel(exp_data)
    t  = exp_data{i}{1}(:);
    r  = exp_data{i}{2}(:);
    N0 = exp_data{i}{3};
    s  = exp_data{i}{4}(:);
    meta = exp_data{i}{5};

    % pick plateau window + temperature
    if which_step == 1
        win = meta.t1;  Thold = meta.T1;
        t0 = win(1);
        % coverage at start of step 1 plateau = N0 minus desorbed BEFORE t0
        Nstart = max(N0 - trapz(t(t<=t0), r(t<=t0)), 0);
    else
        win = meta.t2;  Thold = meta.T2;
        t0 = win(1);
        % coverage at start of step 2 plateau = N0 minus desorbed BEFORE t0
        Nstart = max(N0 - trapz(t(t<=t0), r(t<=t0)), 0);
    end

    idx = (t >= win(1)) & (t <= win(2));
    if nnz(idx) < 5
        continue
    end

    t_seg = t(idx);
    r_seg = r(idx);
    s_seg = s(idx);

    % simulate isothermal decay for this plateau
    [t_m, r_m] = pw_niemant_isothermal(t_seg, Thold, Ea, w, nu, Tc, Nstart, n);

    % weighted residuals
    wts = 1 ./ max(s_seg.^2, eps);
    res = (r_m - r_seg) .* sqrt(wts);
    total = total + sum(res.^2);
end

err = total;
end


function [t_out, r_out] = pw_niemant_isothermal(t_query, T, Ea, w, nu, Tc, N0, n)
% Numerically integrate dN/dt = -nu * N^n * exp(-E(N,T)/(R*T))
% then return r(t) = -dN/dt evaluated at t_query

R = 8.314462618;  % J/mol/K
t0 = t_query(1);
tf = t_query(end);
tspan = [t0 tf];

% integrate coverage
ode = @(t,N) -nu * max(N,0).^n .* exp( -E_des_niemant(max(N,0), T, Ea, w, Tc) ./ (R*T) );

% stiff-safe is often helpful
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
[t_sol, N_sol] = ode15s(ode, tspan, N0, opts);

% evaluate N(t) on requested grid
Nq = interp1(t_sol, N_sol, t_query, 'pchip', 'extrap');

% compute rate on that grid
Eq = E_des_niemant(max(Nq,0), T, Ea, w, Tc);
r_out = nu .* max(Nq,0).^n .* exp( -Eq ./ (R*T) );

t_out = t_query;
end


function E = E_des_niemant(N, T, Ea, w, Tc)
% --- REPLACE THIS with your exact Niemantsverdriet coverage/ordering model ---
% This placeholder just implements a "turn-off" of interactions above Tc:
if Tc > 0
    f = max(0, 1 - T./Tc);   % 1 below Tc, 0 above Tc (crude)
else
    f = 1;
end
E = Ea + (w*f).*N;  % sign convention: w<0 => barrier decreases with coverage
end


function [Hp, Cov, SE] = hess_to_cov(Hx, S, fit_err, exp_data, which_step)
% Hp = S^{-T} Hx S^{-1}  (S diagonal)
if isempty(Hx) || any(size(Hx) ~= [4 4]) || any(isnan(Hx(:)))
    Hp = nan(4); Cov = nan(4); SE = nan(4,1);
    return
end

invS = diag(1./diag(S));
Hp = invS * Hx * invS;

% count points used in this step
n_data = 0;
for i = 1:numel(exp_data)
    meta = exp_data{i}{5};
    t = exp_data{i}{1}(:);
    if which_step == 1
        win = meta.t1;
    else
        win = meta.t2;
    end
    n_data = n_data + nnz((t>=win(1)) & (t<=win(2)));
end

n_params = 4;
resvar = fit_err / max(n_data - n_params, 1);

Cov = resvar * inv(Hp);
SE = sqrt(diag(Cov));
end


function plot_two_step_fit(exp_data, p1, p2, n)
figure; clf; hold on;

for i = 1:numel(exp_data)
    t  = exp_data{i}{1}(:);
    r  = exp_data{i}{2}(:);
    N0 = exp_data{i}{3};
    meta = exp_data{i}{5};

    plot(t, r, 'k-', 'HandleVisibility','off'); % raw trace in black

    % step 1
    win = meta.t1; T = meta.T1;
    Nstart = max(N0 - trapz(t(t<=win(1)), r(t<=win(1))), 0);
    idx = (t>=win(1)) & (t<=win(2));
    if nnz(idx)>5
        [~, r1] = pw_niemant_isothermal(t(idx), T, p1(1), p1(2), p1(3), p1(4), Nstart, n);
        plot(t(idx), r1, 'r-', 'LineWidth', 1.5, 'HandleVisibility','off');
    end

    % step 2
    win = meta.t2; T = meta.T2;
    Nstart = max(N0 - trapz(t(t<=win(1)), r(t<=win(1))), 0);
    idx = (t>=win(1)) & (t<=win(2));
    if nnz(idx)>5
        [~, r2] = pw_niemant_isothermal(t(idx), T, p2(1), p2(2), p2(3), p2(4), Nstart, n);
        plot(t(idx), r2, 'b-', 'LineWidth', 1.5, 'HandleVisibility','off');
    end
end

xlabel('Time (s)'); ylabel('Rate (ML/s)');
title('Two-step isothermal fits: step1=red, step2=blue');
end