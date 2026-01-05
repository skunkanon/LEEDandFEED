%% wire_1D_spatial_with_dime.m
% 1D spatial wire model (method of lines) + lumped crystal/dime node
% Recycles your parameters; adds a few spatial knobs for gradient plotting.


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


% ---- Universal constants ----
sigma = 5.670374419e-8;        % W/m^2/K^4
T_env = 300;                   % K, ambient
T_b   = 300;                   % K, bath (clamps/manipulator sink)

% ---- Tantalum wire properties ----
molarheatcap_Ta = 25.36;       % J/mol·K
atomicmass_Ta   = 180.947;     % g/mol
density_Ta      = 16.678;      % g/cm^3
conduct_Ta      = 57.5;        % W/m·K  (thermal conductivity)
eps_w           = 0.90;        % emissivity (Ta wire)

% ---- Wire geometry ----
length_wire    = 8e-3;         % m
diameter_wire  = 0.25e-3;      % m
area_wire      = pi * (diameter_wire/2)^2;          % m^2
SA_wire_single = pi * diameter_wire * length_wire;  % m^2 (lateral area)

% ---- Multi-wire setup ----
n_wires = 6;                   % identical wires

% ---- Electrical parameters ----
rho_elec_Ta = 1.38e-7;         % ohm*m target at ~300 K (used to scale Ta_scan fit)
I_max    = 20;                 % A total available current
I_max_single = I_max / (n_wires); 

% ---- Copper dime properties (as your "crystal") ----
dime_radius   = 17.91e-3/2;    % m
dime_thick    = 1.35e-3;       % m
cu_density    = 8.96e3;        % kg/m^3
C_cu          = 24.47;         % J/mol·K
atomicmass_Cu = 63.546;        % g/mol
eps_a         = 0.25;          % emissivity of dime (effective)

spotweldscale = 1;             % NEW VARIABLE (scales wire<->dime coupling)

% ---- Derived dime properties ----
mass_crystal = pi * dime_radius^2 * dime_thick * cu_density;           % kg
heatcap_crystal = C_cu * (mass_crystal * 1000) / atomicmass_Cu;        % J/K
heatcap_crystal = 0.235; % DOC HARDCODE @ 0.235, only ~1.1 works (lowering spotweldscale doesn't) 
A_crys = 2*pi*dime_radius^2 + 2*pi*dime_radius*dime_thick;             % m^2

% ---- Timing ----
t_on  = 1090.0;   % s, heating pulse duration
t_end = 1260.0;   % s, total simulation horizon

%% ==================================

N_nodes       = 41;              % number of wire nodes along x 
x_attach      = length_wire;   % where the dime is welded (m). midpoint default
t_grad_plot   = t_on;            % time at which to plot dT/dx profile (s)

endSinkScale  = 5;               % clamp/bath heat-sinking strength multiplier at each end
h_wire        = 0.0;             % W/m^2/K (optional gas cooling along wire; set 0 to ignore)
h_crys        = 0.0;             % W/m^2/K (optional gas cooling of crystal; set 0 to ignore)

%% ========================
% Resistivity fit -> rho(T) from Ta_scan
% =======================
Ta_temp = Ta_scan(:,1).';
Ta_rho  = Ta_scan(:,2).';

p = polyfit(Ta_temp, Ta_rho, 3);



Ta_rho_FUNCT = @(T) polyval(p, T) * 1e-8;    % ohm*m


%% =======================================
% Conductivity fit 
% ========

Ta_temp_CONDUCT = Ta_scan_CONDUCT(:,1)';
Ta_CONDUCT = Ta_scan_CONDUCT(:,2)';

p2 = polyfit(Ta_temp_CONDUCT, Ta_CONDUCT, 2);



Ta_CONDUCT_FUNCT = @(T) polyval(p2, T);             % watts/ meters K 



%% ==================================
% Build 1-D wire discretization
% =============================
x  = linspace(0, length_wire, N_nodes).';
dx = x(2) - x(1);

% Mass-based heat capacity of Ta (J/kg/K) from molar Cp
cp_Ta      = molarheatcap_Ta / (atomicmass_Ta/1000);     % J/kg/K
rho_Ta_kgm = density_Ta * 1000;                          % kg/m^3   (g/cm^3 -> kg/m^3)

% Segment properties (per node/segment)
m_seg   = rho_Ta_kgm * area_wire * dx;                   % kg
C_seg   = m_seg * cp_Ta;                                 % J/K
SA_seg  = pi * diameter_wire * dx;                       % m^2 lateral area per segment

% Axial thermal conductance between neighboring nodes (W/K)
G_ax = conduct_Ta * area_wire / dx;

% End sinking to bath (each end)
G_end = endSinkScale * G_ax;

%{

% Wire<->dime weld conductance (W/K)
G_wc_base = conduct_Ta * area_wire / length_wire;        
G_wc = spotweldscale * G_wc_base;
%} 

% Choose attachment node index
[~, idx_attach] = min(abs(x - x_attach));
idx_mid = round((N_nodes+1)/2);

%% ==================================
% Current gate
% ==================================
I_gate = @(t) double(t < t_on);        % 0/1
I_amp  = I_max_single;                % per-wire current amplitude (your convention)

%% ==================================
% ODE system: y = [Twire_nodes ; Tcrystal]
% ==================================
y0 = [T_env*ones(N_nodes,1); T_env];

opts = odeset('RelTol',1e-7,'AbsTol',1e-7);  % start reasonable; tighten later if you want

odefun = @(t,y) wireSpatialODE(t,y, ...
    N_nodes, idx_attach, ...
    C_seg, ...
    area_wire, dx, SA_seg, length_wire, ...
    Ta_rho_FUNCT, Ta_CONDUCT_FUNCT, ...
    spotweldscale, endSinkScale, ...
    eps_w, eps_a, sigma, ...
    T_env, T_b, ...
    A_crys, heatcap_crystal, ...
    h_wire, h_crys, ...
    n_wires, I_gate, I_amp);


% PDE discretizations often become stiff as N grows -> ode15s is a good default here.
[t,y] = ode15s(odefun, [0 t_end], y0, opts);

Twire = y(:, 1:N_nodes);          % [Nt x N_nodes]
T_a   = y(:, N_nodes+1);          % [Nt x 1]
T_w_mid = Twire(:, idx_mid);

%% ==================================
% Diagnostics: total Joule + losses (all wires)
% ==================================
It = arrayfun(I_gate, t) * I_amp;                    % per-wire current vs time
It2 = It.^2;

% per-segment resistance at each time: R_seg = rho(T)*dx/A
R_seg_mat = Ta_rho_FUNCT(Twire) .* (dx ./ area_wire);                % [Nt x N_nodes]
Pjoule_perWire = It2 .* sum(R_seg_mat, 2);                            % W, one wire
Pjoule_total   = n_wires * Pjoule_perWire;

Prad_wire_perWire = eps_w * sigma * sum(SA_seg * (Twire.^4 - T_env^4), 2);
Prad_wire_total   = n_wires * Prad_wire_perWire;

Prad_crys = eps_a * sigma * A_crys .* (T_a.^4 - T_env^4);

%% ==================================
% Gradient plot at chosen time
% ==================================
t_target = min(t_grad_plot, t(end));
[~, itg] = min(abs(t - t_target));
Tprof = Twire(itg, :).';
dTdx  = gradient(Tprof, x);  % K/m

%% ==================================
% PLOTS (reuse your 4x2 layout, swap one tile for dT/dx)
% ==================================
figure(14); clf;
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');

% (1) current command
nexttile(1);
stairs(t, arrayfun(I_gate, t), 'LineWidth', 1.8); ylim([-0.05 1.05]);
xlabel('time (s)'); ylabel('I(t) (scaled)'); title('Current command'); grid on;

% (2) wire temp profile at t_target
nexttile(2);
plot(x*1e3, Tprof, 'LineWidth', 1.8);
xlabel('x along wire (mm)'); ylabel('T_{wire}(x) (K)');
title(sprintf('Wire profile at t = %.1f s', t(itg)));
grid on;

% (3) mid-wire temperature vs time
nexttile(3);
plot(t, T_w_mid, 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('T_{wire,mid} (K)');
title('Wire temperature (midpoint)'); grid on;

% (4) Joule power (total)
nexttile(4);
plot(t, Pjoule_total, 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('Power (W)'); title('Joule power (total)'); grid on;

% (5) crystal temperature + data
nexttile(5);
hold on;
plot(t, T_a, 'LineWidth', 1.8);
scatter(timespan_IRL, tempspan_IRL, 30, 'LineWidth', 1.5);
xlabel('time (s)'); ylabel('T_{crystal} (K)');
title('Crystal/dime temperature'); grid on; hold off;

% (6) spatial gradient dT/dx (this is your new subplot)
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

% (8) resistance sanity: total wire R(t) per wire (sum of segment R)
nexttile(8);
plot(t, sum(R_seg_mat,2), 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('R_{wire} (\Omega)');
title('Per-wire R(t) = sum R_{seg}(t)'); grid on;

%% ==================================
% Local function
% ==================================
function dydt = wireSpatialODE(t,y, ...
    N, idx_attach, ...
    Cseg, ...
    A, dx, SAseg, Lwire, ...
    rhoFUN, kFUN, ...
    spotweldscale, endSinkScale, ...
    epsw, epsa, sigma, ...
    Tenv, Tbath, ...
    Acrys, Ccrys, ...
    hwire, hcrys, ...
    nwires, I_gate, I_amp)

    % State
    Tw = y(1:N);
    Ta = y(N+1);

    % Current (per wire)
    I = I_gate(t) * I_amp;

    % Electrical: segment resistance + Joule heat (per wire)
    rho_e = rhoFUN(Tw);                % ohm*m
    Rseg  = rho_e .* (dx ./ A);        % ohm
    Pjseg = (I^2) .* Rseg;             % W per segment (per wire)

    % Losses (per wire)
    Prad  = epsw * sigma * SAseg .* (Tw.^4 - Tenv^4);   % W
    Pconv = hwire * SAseg .* (Tw - Tenv);               % W

    % Thermal conductivity k(T) and interface conductances
    k_i = kFUN(Tw);                    % W/m/K at nodes
    k_i = max(k_i, 1e-6);              % safety clamp against negative/zero polyfit weirdness

    k_half = 0.5 * (k_i(1:end-1) + k_i(2:end));   % W/m/K at interfaces
    G_half = (k_half .* A) ./ dx;                 % W/K at interfaces (N-1)

    % End sinks to bath (use local k at ends)
    G_end_1 = endSinkScale * (k_i(1) * A / dx);
    G_end_N = endSinkScale * (k_i(N) * A / dx);

    % Wire <-> crystal weld conductance (use local k at attach node)
    G_wc = spotweldscale * (k_i(idx_attach) * A / Lwire);   % W/K (per wire)

    % Axial conduction + node balances
    dTdt = zeros(N,1);

    % Left boundary node (i=1)
    Q_from_2 = G_half(1) * (Tw(2) - Tw(1));
    Q_sink_1 = G_end_1   * (Tw(1) - Tbath);
    dTdt(1)  = (Q_from_2 + Pjseg(1) - Prad(1) - Pconv(1) - Q_sink_1) / Cseg;

    % Interior nodes
    for i = 2:N-1
        Q_left  = G_half(i-1) * (Tw(i-1) - Tw(i));
        Q_right = G_half(i)   * (Tw(i+1) - Tw(i));
        dTdt(i) = (Q_left + Q_right + Pjseg(i) - Prad(i) - Pconv(i)) / Cseg;
    end

    % Right boundary node (i=N)
    Q_from_Nm1 = G_half(end) * (Tw(N-1) - Tw(N));
    Q_sink_N   = G_end_N     * (Tw(N) - Tbath);
    dTdt(N)    = (Q_from_Nm1 + Pjseg(N) - Prad(N) - Pconv(N) - (Q_sink_N * 0) ) / Cseg;

    % Wire <-> crystal coupling at attachment node (per wire)
    Qwc = G_wc * (Tw(idx_attach) - Ta);            % W leaving wire into crystal (per wire)
    dTdt(idx_attach) = dTdt(idx_attach) - Qwc / Cseg;

    % Crystal (lumped)
    Prad_a  = epsa * sigma * Acrys * (Ta^4 - Tenv^4);
    Pconv_a = hcrys * Acrys * (Ta - Tenv);
    dTa_dt  = (nwires * Qwc - Prad_a - Pconv_a) / Ccrys;

    % OUTPUT (this line is what your error says is missing)
    dydt = [dTdt; dTa_dt];
end
