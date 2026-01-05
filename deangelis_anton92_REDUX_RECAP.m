%% 12/28; ACTUAL DATA

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
timespan_IRL = span_IRL(:,2)';
tempspan_IRL = span_IRL(:,1)';


%%

% DEFINE VARIABLES


% ---- Universal constants ----
sigma = 5.670374419e-8;        % W/m^2/K^4, Stefan–Boltzmann constant
T_env = 300;                   % K, ambient environment
T_b   = 300;                    % K, bath temperature

% ---- Tantalum wire properties ----
molarheatcap_Ta = 25.36;       % J/mol·K
atomicmass_Ta   = 180.947;     % g/mol
density_Ta      = 16.678;      % g/cm^3
conduct_Ta      = 57.5;        % W/m·K
eps_w           = 0.30;        % emissivity (Ta wire)

% ---- Wire geometry ----
length_wire   = 8e-3;         % m
diameter_wire = 0.25e-3;        % m
area_wire     = pi * (diameter_wire/2)^2;          % m^2 (≈2.83e-7)
SA_wire_single = pi * diameter_wire * length_wire;  % m^2 (lateral area)

% ---- Derived wire properties ----
mass_wire = density_Ta * length_wire * area_wire * 1e6;  % g → convert cm^3→m^3
heatcap_wire = (mass_wire / atomicmass_Ta) * molarheatcap_Ta; % J/K (per wire)
conduct_wire = conduct_Ta * (area_wire / length_wire);         % W/K


% ---- Multi-wire setup ----
n_wires = 6;                   % number of identical wires
SA_wire_total = n_wires * SA_wire_single;  % total radiative area (m^2)

% ---- Electrical parameters ----
rho_elec_Ta = 2.2e-7;          % Ω·m (approx @300 K)
R_single = rho_elec_Ta * length_wire / area_wire;  % Ω, CONSTANT 
I_max    = 20;                 % A, total available current
I_max_single = I_max / (n_wires/2);  % current per active span (your convention)

% ---- Nickel crystal properties ----
ni_radius  = 0.005;            % m
ni_thick   = 0.0002;           % m
ni_density = 8.9e3;            % kg/m^3
C_ni       = 26.07;            % J/mol·K
atomicmass_Ni = 58.69;         % g/mol
eps_a      = 0.05;             % emissivity (Ni crystal)

% ---- Copper dime properties (geometry of a US dime; material assumed pure Cu) ----
dime_radius  = 17.91e-3/2;      % m  (US dime diameter = 17.91 mm)
dime_thick   = 1.35e-3;         % m  (US dime thickness = 1.35 mm)
cu_density   = 8.96e3;          % kg/m^3 (copper)
C_cu         = 24.47;           % J/mol·K (molar Cp of Cu at ~25 °C)
atomicmass_Cu = 63.546;         % g/mol
eps_a        = 0.45;            % emissivity (polished copper 0.1 to 0.2; oxidized can be much higher), doesnt affect ramp rate though 

spotweldscale = 1; %NEW VARIABLE, 12/28

ni_radius = dime_radius;
ni_thick = dime_thick;
ni_density = cu_density;
C_ni = C_cu;
atomicmass_Ni = atomicmass_Cu;



%C_ni = 24.98; %RHODIUM
%atomicmass_Ni = 102.91;
%eps_a = ;



% ---- Derived Ni properties ----
mass_crystal = pi * ni_radius^2 * ni_thick * ni_density; % kg
heatcap_crystal = C_ni * (mass_crystal * 1000) / atomicmass_Ni; % J/K
A_crys = 2*pi*ni_radius^2 + 2*pi*ni_radius*ni_thick;     % m^2, total area


% ---- Timing ----
t_on  = 1590.0;   % s, heating pulse
t_end = 1760.0;   % s, total simulation horizon


 %=== Resistivity fit → R(T) for the Ta span ===
% Expect Ta_scan(:,1) = T [K], Ta_scan(:,2) = rho data (see units note below)

Ta_temp = Ta_scan(:,1).';
Ta_rho  = Ta_scan(:,2).';

% Fit a cubic (you already tested this elsewhere)
p = polyfit(Ta_temp, Ta_rho, 3);

% Units note: Desai published as 10^-8
scale_BASE = 1e-8;
scale = scale_BASE * 2.2/1.63; %makes rho @ 300 K the same as RECAP's 

% Resistivity vs T [ohm·m]
Ta_rho_FUNCT = @(T) ( ( (p(1)*T + p(2)).*T + p(3) ).*T + p(4) ) * scale;

% Geometry for THIS page's single Ta span
R_wire_FUNCT = @(T) Ta_rho_FUNCT(T) .* (length_wire ./ area_wire);  % [ohm]

%% === ODE with temperature-dependent R and proper radiating area ===
I = @(t) (t < t_on);                 % 0/1 gate
I_amp_single = I_max_single;         % amps 
%{ 
%CONSTANT RESISTIVITY

odefun = @(t,y) [ ...
    % --- wire node (y1) ---
    ( (I(t) * I_amp_single).^2 .* R_single ...   % Joule heat, W, R_wire_FUNCT(y(1))
      -  (conduct_wire*spotweldscale) * (y(1) - y(2)) ...                  % to crystal
      -  conduct_wire * (y(1) - T_b) ...                   % to bath
      -  eps_w * sigma * SA_wire_total * ( y(1).^4 - T_env^4 ) ... % radiation (all wires)
    ) / heatcap_wire ; ...
    % --- crystal node (y2) ---
    ( (conduct_wire * spotweldscale) * (y(1) - y(2)) * n_wires ...           % from all wires
      - eps_a * sigma * A_crys * ( y(2).^4 - T_env^4 ) ... % crystal radiative loss
    ) / heatcap_crystal ...
];
%}


odefun = @(t,y) [ ...
    % --- wire node (y1) ---
    ( (I(t) * I_amp_single).^2 .* R_wire_FUNCT(y(1)) ...   % Joule heat, W, R_wire_FUNCT(y(1))
      -  (conduct_wire*spotweldscale) * (y(1) - y(2)) ...                  % to crystal
      -  (conduct_wire*spotweldscale) * (y(1) - T_b) ...                   % to bath
      -  eps_w * sigma * SA_wire_total * ( y(1).^4 - T_env^4 ) ... % radiation (all wires)
    ) / heatcap_wire ; ...
    % --- crystal node (y2) ---
    ( (conduct_wire * spotweldscale) * (y(1) - y(2)) * n_wires ...           % from all wires
      - eps_a * sigma * A_crys * ( y(2).^4 - T_env^4 ) ... % crystal radiative loss
    ) / heatcap_crystal ...
];

% ---- Integrate and plot ----
y0 = [300; 300];                 % K
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(odefun,[0 t_end],y0,opts);
T_w = y(:,1); T_a = y(:,2);

tI  = linspace(0, t_end, 800);
ItI = arrayfun(I, tI);

% Crystal heating rate (model): dT/dt [K/s]
dTcrys_dt = gradient(T_a, t);

% Crystal heating rate (data): dT/dt [K/s]  (optional overlay)
t_IRL = timespan_IRL(:);
T_IRL = tempspan_IRL(:);
dTcrys_dt_IRL = gradient(T_IRL, t_IRL);

% Optional light smoothing (uncomment if you want less jaggedness)
% dTcrys_dt = movmean(dTcrys_dt, 5);
% dTcrys_dt_IRL = movmean(dTcrys_dt_IRL, 3);

%% DIAGNOSTICS

% --- Radiative losses (needed for diagnostics & plots)
radiative_loss_wire    = eps_w * sigma * SA_wire_total * (T_w.^4 - T_env^4);  % W
radiative_loss_crystal = eps_a * sigma * A_crys        * (T_a.^4 - T_env^4);  % W
total_radiative_loss   = radiative_loss_wire + radiative_loss_crystal;        % W



% --- CURRENT used in ODE ---
% If single wire: 
%I_inst = I(t) * I_max_single;                   % amps (single-wire convention)

% If multiple wires
 I_inst = I(t) * (I_max / n_wires);            % amps in each span

% Temperature-dependent resistance of the Ta span
R_t = R_wire_FUNCT(T_w);                        % ohm

%If single wire: 
%P_joule = (I_inst.^2) .* R_t;                   % W

%If multiple wires: 
 P_joule = n_wires * (I_inst.^2) .* R_t;      % W

% Conduction terms (note: your ODE used conduct_wire twice from wire→crystal and wire→bath)
Q_cond_wire_to_crys =  conduct_wire * (T_w - T_a);   % W
Q_cond_wire_to_bath =  conduct_wire * (T_w - T_b);   % W

% Radiation already computed above for the wire node: radiative_loss_wire (W)

% Energy balance residual at the wire node (should match heatcap_wire * dT_w/dt)
P_balance = P_joule - Q_cond_wire_to_crys - Q_cond_wire_to_bath - radiative_loss_wire;

% Optional: numerical dT_w/dt to check consistency
dT_w_dt = gradient(T_w, t);
lhs_wire = heatcap_wire * dT_w_dt;   % predicted by ODE right-hand side


%% EVERYTHING PLOT
figure(14); clf;
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');

% ============ LEFT COLUMN (thermal evolution) ============

% Row 1: current command
nexttile(1);
stairs(tI, ItI, 'LineWidth', 1.8); ylim([-0.05 1.05]);
xlabel('time (s)'); ylabel('I(t) (scaled)');
title('Current command');
grid on;

% Row 2: wire temperature
nexttile(3);
plot(t, T_w, 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('T_{wire} (K)');
title('Wire temperature');
grid on;

% Row 3: crystal temperature
nexttile(5);
hold on;
plot(t, T_a, 'LineWidth', 1.8);
scatter(timespan_IRL, tempspan_IRL, 'LineWidth', 2);
xlabel('time (s)'); ylabel('T_{Ni} (K)');
title('Crystal temperature');
hold off;
grid on;

% Row 4: radiative losses
nexttile(7);
plot(t, radiative_loss_wire, 'LineWidth', 1.8, 'DisplayName','Wire');
hold on;
plot(t, radiative_loss_crystal, 'LineWidth', 1.8, 'DisplayName','Crystal');
plot(t, total_radiative_loss, '--', 'LineWidth', 2, 'DisplayName','Total');
xlabel('time (s)'); ylabel('Radiative loss (W)');
title('Radiative heat loss');
legend('Location','best'); grid on; hold off;


% ============ RIGHT COLUMN (diagnostics) ============

% Row 1: resistance vs time
nexttile(2);
plot(t, R_t, 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('R_{wire} (\Omega)');
title('Temperature-dependent R(T) of wire');
grid on;

% Row 2: Joule power
nexttile(4);
plot(t, P_joule, 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('Power (W)');
title('Joule power');
grid on;
% 
% % Row 3: conduction + radiation
nexttile(6);
plot(t, Q_cond_wire_to_crys, 'LineWidth',1.8,'DisplayName','Cond → crystal'); hold on;
plot(t, Q_cond_wire_to_bath, 'LineWidth',1.8,'DisplayName','Cond → bath');
plot(t, radiative_loss_wire, '--','LineWidth',1.8,'DisplayName','Radiation');
xlabel('time (s)'); ylabel('W');
title('Wire heat-flow channels');
legend('Location','best'); grid on; hold off;


nexttile(8); hold on;

% Data + fit

% (8) crystal heating rate
nexttile(8);
hold on;
plot(t, dTcrys_dt, 'LineWidth', 1.8, 'DisplayName','Model dT/dt');
plot(t_IRL, dTcrys_dt_IRL, '--', 'LineWidth', 1.8, 'DisplayName','Data dT/dt');
xlabel('time (s)'); ylabel('dT_{crystal}/dt (K/s)');
title('Crystal heating rate');
legend('Location','best'); grid on; hold off;


%% GRADIENT PLOT 1/2 
 figure(15); clf;

 % === Define spatial grid along wire (in mm for clarity) ===
 L_wire = length_wire;            % meters (already defined in script)
 L_mm = L_wire * 1000;            % convert to mm for plotting
 x_m  = linspace(0, L_wire, 51);  % spatial grid in meters
 x_mm = x_m * 1000;               % same grid in mm for the plots

 % === Select several time points to sample profiles ===
 % (more contours in the top subplot for a richer time sweep)
 numProfiles = 10;
 idx_sample = round(linspace(1, length(t), numProfiles));
 t_sample = t(idx_sample);
 T_sample = T_a(idx_sample);                 % crystal temps at those times (right end)

 % === Helper: parabolic profile from uniform Joule heating in a rod ===
 % T(x) = T_b + (T_R - T_b)*(x/L) + (q'''/(2*k))*(L*x - x^2)
 % where q''' is volumetric heating [W/m^3] and k = conduct_Ta [W/m/K].
 volume_wire = L_wire * area_wire;              % per wire
 P_wire_sample = P_joule(idx_sample) / n_wires; % W in a single wire at samples
 q_vol_sample = P_wire_sample ./ volume_wire;   % W/m^3

 % Precompute the linear and quadratic coefficients for the analytic gradient
 % dT/dx = (T_R - T_b)/L + (q'''/(2*k))*(L - 2x)

 % Distinct style helpers so overlapping contours stay visible
 colors = lines(numProfiles);
 lineStyles = {'-','--',':','-.'};
 markers = {'o','s','d','^','v','>','<','p','h','x'};
 markerStride = max(2, floor(numel(x_mm) / 10));

 % === Subplot 1: T(x) profiles at selected times ===
 subplot(2,1,1); hold on;
 for k = 1:numProfiles
     T_left = T_b;                  % bath end (x = 0)
     T_right = T_sample(k);         % crystal end (x = L)
     T_x = T_left ...
         + (T_right - T_left) * (x_m / L_wire) ...
         + (q_vol_sample(k) / (2 * conduct_Ta)) .* (L_wire * x_m - x_m.^2);
     plot(x_mm, T_x, 'LineWidth', 1.6, ...
          'Color', colors(k,:), 'LineStyle', lineStyles{mod(k-1, numel(lineStyles))+1}, ...
          'Marker', markers{mod(k-1, numel(markers))+1}, 'MarkerIndices', 1:markerStride:numel(x_mm), ...
          'DisplayName', sprintf('t = %.1f s', t_sample(k)));
 end
 xlabel('Position along wire (mm)');
 ylabel('Temperature (K)');
 title('Wire Temperature Profile with Uniform Joule Heating');
 legend('Location','best');
 grid on; hold off;

 % === Subplot 2: Spatial gradient dT/dx for the same snapshots ===
 subplot(2,1,2); hold on;
 for k = 1:numProfiles
     T_left = T_b;
     T_right = T_sample(k);
     dTdx_analytic = (T_right - T_left) / L_wire ...
         + (q_vol_sample(k) / (2 * conduct_Ta)) .* (L_wire - 2 * x_m);
     plot(x_mm, dTdx_analytic, 'LineWidth', 1.6, ...
          'Color', colors(k,:), 'LineStyle', lineStyles{mod(k-1, numel(lineStyles))+1}, ...
          'Marker', markers{mod(k-1, numel(markers))+1}, 'MarkerIndices', 1:markerStride:numel(x_mm), ...
          'DisplayName', sprintf('t = %.1f s', t_sample(k)));
 end
 xlabel('Position along wire (mm)');
 ylabel('dT/dx (K/m)');
 title('Spatial Temperature Gradient Along Wire (heat flux \propto -k dT/dx)');
 legend('Location','best');
 grid on; hold off;

