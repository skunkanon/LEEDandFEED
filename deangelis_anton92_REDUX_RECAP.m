%OCT 22, new DEANGELIS
%run TANT.m first

% === Resistivity fit → R(T) for the Ta span ===
% Expect Ta_scan(:,1) = T [K], Ta_scan(:,2) = rho data (see units note below)

Ta_temp = Ta_scan(:,1).';
Ta_rho  = Ta_scan(:,2).';

% Fit a cubic (you already tested this elsewhere)
p = polyfit(Ta_temp, Ta_rho, 3);

% Units note: Desai published as 10^-8
scale = 1e-8;  

% Resistivity vs T [ohm·m]
Ta_rho_FUNCT = @(T) ( ( (p(1)*T + p(2)).*T + p(3) ).*T + p(4) ) * scale;

% Geometry for THIS page's single Ta span
R_wire_FUNCT = @(T) Ta_rho_FUNCT(T) .* (length_wire ./ area_wire);  % [ohm]
%% DEFINE VARIABLES


% ---- Universal constants ----
sigma = 5.670374419e-8;        % W/m^2/K^4, Stefan–Boltzmann constant
T_env = 300;                   % K, ambient environment
T_b   = 80;                    % K, bath temperature

% ---- Tantalum wire properties ----
molarheatcap_Ta = 25.36;       % J/mol·K
atomicmass_Ta   = 180.947;     % g/mol
density_Ta      = 16.678;      % g/cm^3
conduct_Ta      = 57.5;        % W/m·K
eps_w           = 0.30;        % emissivity (Ta wire)

% ---- Wire geometry ----
length_wire   = 0.005;         % m
diameter_wire = 0.0006;        % m
area_wire     = pi * (diameter_wire/2)^2;          % m^2 (≈2.83e-7)
SA_wire_single = pi * diameter_wire * length_wire;  % m^2 (lateral area)

% ---- Derived wire properties ----
mass_wire = density_Ta * length_wire * area_wire * 1e6;  % g → convert cm^3→m^3
heatcap_wire = (mass_wire / atomicmass_Ta) * molarheatcap_Ta; % J/K (per wire)
conduct_wire = conduct_Ta * (area_wire / length_wire);         % W/K


% ---- Multi-wire setup ----
n_wires = 4;                   % number of identical wires
SA_wire_total = n_wires * SA_wire_single;  % total radiative area (m^2)

% ---- Electrical parameters ----
rho_elec_Ta = 2.2e-7;          % Ω·m (approx @300 K)
R_single = rho_elec_Ta * length_wire / area_wire;  % Ω
R_total  = R_single / 2;       % effective for two-wire pair (legacy)
I_max    = 40;                 % A, total available current
I_max_single = I_max / (n_wires/2);  % current per active span (your convention)

% ---- Nickel crystal properties ----
ni_radius  = 0.005;            % m
ni_thick   = 0.0002;           % m
ni_density = 8.9e3;            % kg/m^3
C_ni       = 26.07;            % J/mol·K
atomicmass_Ni = 58.69;         % g/mol
eps_a      = 0.15;             % emissivity (Ni crystal)

% ---- Derived Ni properties ----
mass_crystal = pi * ni_radius^2 * ni_thick * ni_density; % kg
heatcap_crystal = C_ni * (mass_crystal * 1000) / atomicmass_Ni; % J/K
A_crys = 2*pi*ni_radius^2 + 2*pi*ni_radius*ni_thick;     % m^2, total area


% ---- Timing ----
t_on  = 2.0;   % s, heating pulse
t_end = 8.0;   % s, total simulation horizon




%% === ODE with temperature-dependent R and proper radiating area ===
I = @(t) (t < t_on);                 % 0/1 gate
I_amp_single = I_max_single;         % amps 

odefun = @(t,y) [ ...
    % --- wire node (y1) ---
    ( (I(t) * I_amp_single).^2 .* R_wire_FUNCT(y(1)) ...   % Joule heat, W
      -  conduct_wire * (y(1) - y(2)) ...                  % to crystal
      -  conduct_wire * (y(1) - T_b) ...                   % to bath
      -  eps_w * sigma * SA_wire_total * ( y(1).^4 - T_env^4 ) ... % radiation (all wires)
    ) / heatcap_wire ; ...
    % --- crystal node (y2) ---
    ( conduct_wire * (y(1) - y(2)) * n_wires ...           % from all wires
      - eps_a * sigma * A_crys * ( y(2).^4 - T_env^4 ) ... % crystal radiative loss
    ) / heatcap_crystal ...
];

% ---- Integrate and plot ----
y0 = [200; 200];                 % K
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(odefun,[0 t_end],y0,opts);
T_w = y(:,1); T_a = y(:,2);

tI  = linspace(0, t_end, 800);
ItI = arrayfun(I, tI);


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
figure(1); clf;
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
plot(t, T_a, 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('T_{Ni} (K)');
title('Crystal temperature');
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
title('Temperature-dependent R(T)');
grid on;

% Row 2: Joule power
nexttile(4);
plot(t, P_joule, 'LineWidth', 1.8);
xlabel('time (s)'); ylabel('Power (W)');
title('Joule power');
grid on;

% Row 3: conduction + radiation
nexttile(6);
plot(t, Q_cond_wire_to_crys, 'LineWidth',1.8,'DisplayName','Cond → crystal'); hold on;
plot(t, Q_cond_wire_to_bath, 'LineWidth',1.8,'DisplayName','Cond → bath');
plot(t, radiative_loss_wire, '--','LineWidth',1.8,'DisplayName','Radiation');
xlabel('time (s)'); ylabel('W');
title('Wire heat-flow channels');
legend('Location','best'); grid on; hold off;

% Row 4: energy balance
nexttile(8);
plot(Ta_temp, Ta_rho * scale, 'ko', 'MarkerFaceColor',[0.2 0.2 0.2], ...
    'DisplayName','Data');
T_fit = linspace(min(Ta_temp), max(Ta_temp), 400);
rho_fit = Ta_rho_FUNCT(T_fit);
plot(T_fit, rho_fit, 'r-', 'LineWidth', 1.8, 'DisplayName','Cubic fit');

T_min_sim = min(T_w);
T_max_sim = max(T_w);
y_limits = ylim;
patch([T_min_sim T_max_sim T_max_sim T_min_sim], ...
      [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
     'b', 'EdgeColor','none', 'FaceAlpha',0.3, ...
      'DisplayName','Simulated T-range');

xlabel('Temperature (K)');
ylabel('\rho_{Ta} (ohm·m)');
title('Tantalum resistivity vs temperature');
legend('Location','northwest');
% overall title
sgtitle('Temperature evolution (left) vs Electrical/thermal diagnostics (right)');
