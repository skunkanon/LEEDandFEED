%% wire_1D_spatial_with_dime_REVISED.m
% 1D spatial Ta-wire model (method of lines) + lumped copper "dime/crystal" node
%
% What this script is doing (in plain English):
% - Discretizes ONE representative Ta wire into N nodes along x (0..L).
% - Solves transient heat balance at each node:
%     (thermal mass) * dT/dt =
%         + axial conduction to neighbors
%         + Joule heating in that segment (I^2 * Rseg(T))
%         - radiative losses to environment
%         - optional convective losses to environment
%         - end heat sinking to a bath at both ends (clamp/leads)
%         - plus a special coupling term at one "attachment node" where the wire dumps heat into the crystal/dime.
% - The crystal/dime is a single lumped temperature Ta(t) with its own heat capacity and radiative/convective losses.
% - "n_wires" multiplies the wire→crystal heat flow (because many identical wires feed the same crystal).
%
% Key knobs:
% - x_attach: where along the modeled wire the crystal is attached.
%       If you want a symmetric wire temperature profile with a peak at the midpoint,
%       you generally want x_attach = L/2 AND identical end sinks at both ends.
%       (Your previous asymmetry was largely from attaching at the end and/or disabling one end sink.)
% - endSinkScale: multiplies end heat sinking strength (how hard the ends are clamped to Tbath).
% - spotweldscale: scales wire<->crystal coupling conductance.
%       spotweldscale = 1  means "use the naive conductance scale k*A/Lwire"
%       spotweldscale < 1  weaker weld/contact (less coupling)
%       spotweldscale > 1  stronger effective coupling
%
% Experimental data:
% - span_IRL is YOUR measured crystal/dime temperature vs time for this exact configuration.
% - This file keeps your wire diameter, number of wires, and current total fixed.

clear; clc;

%% 12/28; ACTUAL DATA (crystal/dime temperature vs time)
span_IRL = [
303, 0
313, 7
323, 13
333, 19
343, 24
353, 31
363, 33
373, 42
383, 51
393, 58
403, 69
413, 74
423, 83
433, 93
443, 103
453, 116
463, 129
473, 141
483, 159
493, 175
503, 188
];
timespan_IRL = span_IRL(:,2)';   % s
tempspan_IRL = span_IRL(:,1)';   % K

%% =========================
% Universal constants / environment
sigma = 5.670374419e-8;     % W/m^2/K^4
T_env = 300;                % K (radiation/ambient reference)
T_b   = 300;                % K (bath / clamp sink temperature)

%% =========================
% Tantalum wire properties (baseline constants)
molarheatcap_Ta = 25.36;    % J/mol/K
atomicmass_Ta   = 180.947;  % g/mol
density_Ta      = 16.678;   % g/cm^3
eps_w           = 0.90;     % emissivity of Ta wire (effective)

% Default electrical resistivity baseline (ohm*m) if no Ta_scan provided
rho_elec_Ta_300K = 1.38e-7;      % ohm*m around 300 K (order-of-magnitude)
alpha_rho_lin    = 3.5e-3;       % 1/K crude linear slope (fallback)

% Default thermal conductivity baseline (W/m/K) if no Ta_scan_CONDUCT provided
conduct_Ta_300K = 57.5;          % W/m/K around 300 K

%% =========================
% Wire geometry + multi-wire setup (KEEP THESE AS YOUR CONFIG)
length_wire    = 8e-3;           % m
diameter_wire  = 0.25e-3;        % m
area_wire      = pi*(diameter_wire/2)^2;      % m^2
n_wires        = 6;              % identical wires

% lateral surface area per unit length = pi*d
% later we use SA_seg = pi*d*dx

%% =========================
% Electrical drive (KEEP AS YOUR CONFIG)
I_total = 20;                    % A total
I_amp   = I_total / (n_wires/2);     % A per wire
% Current gate: on for t < t_on, off after
%t_on  = max(timespan_IRL);       % s (default: match experiment duration)
t_on = 1000;
t_end = t_on + 60;               % s (some cooldown window)
I_gate = @(t) double(t < t_on);  % 0/1

%% =========================
% Copper "dime/crystal" properties
dime_radius   = 17.91e-3/2;      % m (US dime diameter 17.91 mm)
dime_thick    = 1.35e-3;         % m
cu_density    = 8.96e3;          % kg/m^3
C_cu          = 24.47;           % J/mol/K
atomicmass_Cu = 63.546;          % g/mol
eps_a         = 0.25;            % emissivity of crystal/dime (effective)

% Derived crystal mass + heat capacity
mass_crystal = pi*dime_radius^2 * dime_thick * cu_density;       % kg
heatcap_crystal = C_cu * (mass_crystal*1000) / atomicmass_Cu;    % J/K

% If your advisor wants a hard-coded effective crystal heat capacity:
% (keep this line if you want to force the fit knob)
heatcap_crystal = 0.235;  % J/K DOC HARDCODE 

% Crystal surface area for radiation/conv
A_crys = 2*pi*dime_radius^2 + 2*pi*dime_radius*dime_thick;       % m^2

%% =========================
% Spatial discretization + "symmetry" knobs
N_nodes     = 41;                 % number of wire nodes along x
x           = linspace(0, length_wire, N_nodes).';
dx          = x(2) - x(1);

% Where is the crystal attached along the modeled wire?
% If you want midpoint hottest and symmetric gradient, default to the midpoint:
x_attach = length_wire/2;         % <-- CHANGE THIS if your physical attachment is at an end
[~, idx_attach] = min(abs(x - x_attach));
idx_mid = round((N_nodes+1)/2);

% Optional: choose time at which to plot profile/gradient
t_grad_plot = t_on;

%% =========================
% Heat sinking / gas cooling knobs
endSinkScale  = 5;     % multiplies end sinking strength (both ends)
h_wire        = 0.0;   % W/m^2/K (wire convection, set 0 for vacuum)
h_crys        = 0.0;   % W/m^2/K (crystal convection, set 0 for vacuum)

% Wire<->crystal coupling knob
spotweldscale = 1.0;   % scales the wire↔crystal conductance

%% =========================
% Derived wire segment properties
% Convert Cp(molar) to Cp(mass):
cp_Ta      = molarheatcap_Ta / (atomicmass_Ta/1000);  % J/kg/K
rho_Ta_kgm = density_Ta * 1000;                        % kg/m^3

% Segment mass, heat cap, surface area
m_seg   = rho_Ta_kgm * area_wire * dx;                 % kg
C_seg   = m_seg * cp_Ta;                               % J/K per node/segment
SA_seg  = pi * diameter_wire * dx;                     % m^2 lateral area per segment

%% =========================
% rho(T) and k(T) functions
% If you have Ta_scan arrays in your workspace, we use them.
% Otherwise we fall back to simple constant/linear models.

% --- Resistivity rho(T) [ohm*m] ---
if exist('Ta_scan','var') && size(Ta_scan,2) >= 2
    Ta_temp = Ta_scan(:,1).';
    Ta_rho  = Ta_scan(:,2).';        % whatever units your Ta_scan stores
    p_rho   = polyfit(Ta_temp, Ta_rho, 3);

    % Your prior code multiplied by 1e-8; keep that convention if Ta_scan is in (nOhm*cm) etc.
    Ta_rho_FUNCT = @(T) max(polyval(p_rho, T) * 1e-8, 1e-12);
else
    % fallback: crude linear model around 300 K
    Ta_rho_FUNCT = @(T) max(rho_elec_Ta_300K .* (1 + alpha_rho_lin .* (T - 300)), 1e-12);
end

% --- Thermal conductivity k(T) [W/m/K] ---
if exist('Ta_scan_CONDUCT','var') && size(Ta_scan_CONDUCT,2) >= 2
    Ta_tempK = Ta_scan_CONDUCT(:,1).';
    Ta_k     = Ta_scan_CONDUCT(:,2).';
    p_k      = polyfit(Ta_tempK, Ta_k, 2);
    Ta_k_FUNCT = @(T) max(polyval(p_k, T), 1e-6);
else
    Ta_k_FUNCT = @(T) conduct_Ta_300K + 0*T;  % constant fallback
end

%% =========================
% Pack parameters for ODE
p = struct();
p.N             = N_nodes;
p.idx_attach    = idx_attach;
p.dx            = dx;
p.A             = area_wire;
p.SAseg         = SA_seg;
p.Lwire         = length_wire;

p.Cseg          = C_seg;
p.Ccrys         = heatcap_crystal;

p.rhoFUN        = Ta_rho_FUNCT;
p.kFUN          = Ta_k_FUNCT;

p.I_gate        = I_gate;
p.I_amp         = I_amp;

p.epsw          = eps_w;
p.epsa          = eps_a;
p.sigma         = sigma;

p.Tenv          = T_env;
p.Tbath         = T_b;

p.Acrys         = A_crys;
p.hwire         = h_wire;
p.hcrys         = h_crys;

p.nwires        = n_wires;
p.spotweldscale = spotweldscale;
p.endSinkScale  = endSinkScale;

% numeric safety clamps
p.k_min   = 1e-6;
p.rho_min = 1e-12;

%% =========================
% Solve ODE
y0 = [T_env*ones(N_nodes,1); T_env];
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);

odefun = @(t,y) wireCrystalODE(t,y,p);
[t,y] = ode15s(odefun, [0 t_end], y0, opts);

Twire = y(:, 1:N_nodes);
T_a   = y(:, N_nodes+1);
T_w_mid = Twire(:, idx_mid);

% Crystal heating rate (model): dT/dt [K/s]
dTcrys_dt = gradient(T_a, t);

% Crystal heating rate (data): dT/dt [K/s]  (optional overlay)
t_IRL = timespan_IRL(:);
T_IRL = tempspan_IRL(:);
dTcrys_dt_IRL = gradient(T_IRL, t_IRL);

% Optional light smoothing (uncomment if you want less jaggedness)
% dTcrys_dt = movmean(dTcrys_dt, 5);
% dTcrys_dt_IRL = movmean(dTcrys_dt_IRL, 3);


%% =========================
% Diagnostics: Joule power + losses
It  = arrayfun(I_gate, t) * I_amp;   % per-wire current
It2 = It.^2;

R_seg_mat = Ta_rho_FUNCT(Twire) .* (dx ./ area_wire);      % [Nt x N_nodes]
Pjoule_perWire = It2 .* sum(R_seg_mat, 2);                  % W, one wire
Pjoule_total   = n_wires * Pjoule_perWire;

Prad_wire_perWire = eps_w * sigma * sum(SA_seg * (Twire.^4 - T_env^4), 2);
Prad_wire_total   = n_wires * Prad_wire_perWire;
Prad_crys          = eps_a * sigma * A_crys .* (T_a.^4 - T_env^4);

%% =========================
% Gradient plot at chosen time
t_target = min(t_grad_plot, t(end));
[~, itg] = min(abs(t - t_target));
Tprof = Twire(itg,:).';
dTdx  = gradient(Tprof, x);  % K/m

%% =========================
% PLOTS (4x2 layout)
figure(14); clf;
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');

% (1) current command
nexttile(1);
stairs(t, arrayfun(I_gate, t), 'LineWidth', 1.8); ylim([-0.05 1.05]);
xlabel('time (s)'); ylabel('I(t) gate'); title('Current command (0/1)');
grid on;

% (2) wire temp profile at t_target
nexttile(2);
plot(x*1e3, Tprof, 'LineWidth', 1.8);
xlabel('x along wire (mm)'); ylabel('T_{wire}(x) (K)');
title(sprintf('Wire profile at t = %.1f s', t(itg)));
grid on;

% (3) midpoint wire temperature vs time
nexttile(3);
plot(t, T_w_mid, 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('T_{wire,mid} (K)');
title('Wire temperature (midpoint node)');
grid on;

% (4) Joule power (total)
nexttile(4);
plot(t, Pjoule_total, 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('Power (W)');
title('Joule power (total, all wires)');
grid on;

% (5) crystal temperature + your data
nexttile(5);
hold on;
plot(t, T_a, 'LineWidth', 1.8, 'DisplayName','Model');
scatter(timespan_IRL, tempspan_IRL, 30, 'LineWidth', 1.5, 'DisplayName','span\_IRL');
xlabel('time (s)'); ylabel('T_{crystal} (K)');
title('Crystal/dime temperature');
grid on; legend('Location','best'); hold off;

% (6) spatial gradient dT/dx
nexttile(6);
plot(x*1e3, dTdx, 'LineWidth', 1.8);
xlabel('x along wire (mm)'); ylabel('dT/dx (K/m)');
title(sprintf('Wire gradient at t = %.1f s', t(itg)));
grid on;

% (7) radiative losses
nexttile(7);
plot(t, Prad_wire_total, 'LineWidth', 1.8, 'DisplayName','Wire (total)');
hold on;
plot(t, Prad_crys, 'LineWidth', 1.8, 'DisplayName','Crystal');
plot(t, Prad_wire_total + Prad_crys, '--', 'LineWidth', 2, 'DisplayName','Total');
xlabel('time (s)'); ylabel('Radiative loss (W)');
title('Radiative heat loss');
legend('Location','best'); grid on; hold off;

% (8) crystal heating rate
nexttile(8);
hold on;
plot(t, dTcrys_dt, 'LineWidth', 1.8, 'DisplayName','Model dT/dt');
plot(t_IRL, dTcrys_dt_IRL, '--', 'LineWidth', 1.8, 'DisplayName','Data dT/dt');
xlabel('time (s)'); ylabel('dT_{crystal}/dt (K/s)');
title('Crystal heating rate');
legend('Location','best'); grid on; hold off;

%% =========================
% Local ODE function
function dydt = wireCrystalODE(t, y, p)
    % State
    N  = p.N;
    Tw = y(1:N);
    Ta = y(N+1);

    % Current per wire
    I = p.I_gate(t) * p.I_amp;

    % Electrical heating
    rho_e = p.rhoFUN(Tw);
    rho_e = max(rho_e, p.rho_min);
    Rseg  = rho_e .* (p.dx ./ p.A);     % ohm per segment
    Pjseg = (I^2) .* Rseg;              % W per segment (per wire)

    % Losses along wire
    Prad  = p.epsw * p.sigma * p.SAseg .* (Tw.^4 - p.Tenv^4);
    Pconv = p.hwire * p.SAseg .* (Tw - p.Tenv);

    % Axial conduction
    k_i = p.kFUN(Tw);
    k_i = max(k_i, p.k_min);

    k_half = 0.5 * (k_i(1:end-1) + k_i(2:end));
    G_half = (k_half .* p.A) ./ p.dx;              % W/K between nodes

    % End sinks (BOTH ends, symmetric unless you choose otherwise)
    G_end_1 = p.endSinkScale * (k_i(1) * p.A / p.dx);
    G_end_N = p.endSinkScale * (k_i(N) * p.A / p.dx);

    % Wire <-> crystal coupling at attachment node (per wire)
    ia  = p.idx_attach;
    G_wc = p.spotweldscale * (k_i(ia) * p.A / p.Lwire);  % W/K per wire
    Qwc  = G_wc * (Tw(ia) - Ta);                         % W leaving wire into crystal (per wire)

    % Node balances
    dTdt = zeros(N,1);

    % i=1
    Q_from_2 = G_half(1) * (Tw(2) - Tw(1));
    Q_sink_1 = G_end_1   * (Tw(1) - p.Tbath);
    dTdt(1)  = (Q_from_2 + Pjseg(1) - Prad(1) - Pconv(1) - Q_sink_1) / p.Cseg;

    % interior
    for i = 2:N-1
        Q_left  = G_half(i-1) * (Tw(i-1) - Tw(i));
        Q_right = G_half(i)   * (Tw(i+1) - Tw(i));
        dTdt(i) = (Q_left + Q_right + Pjseg(i) - Prad(i) - Pconv(i)) / p.Cseg;
    end

    % i=N
    Q_from_Nm1 = G_half(end) * (Tw(N-1) - Tw(N));
    Q_sink_N   = G_end_N     * (Tw(N) - p.Tbath);
    dTdt(N)    = (Q_from_Nm1 + Pjseg(N) - Prad(N) - Pconv(N) - Q_sink_N) / p.Cseg;

    % Apply coupling at attachment node
    dTdt(ia) = dTdt(ia) - Qwc / p.Cseg;

    % Crystal balance (lumped)
    Prad_a  = p.epsa * p.sigma * p.Acrys * (Ta^4 - p.Tenv^4);
    Pconv_a = p.hcrys * p.Acrys * (Ta - p.Tenv);
    dTa_dt  = (p.nwires * Qwc - Prad_a - Pconv_a) / p.Ccrys;

    dydt = [dTdt; dTa_dt];
end
