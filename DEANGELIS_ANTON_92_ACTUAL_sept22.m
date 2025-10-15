% 9/22 - ORGANIZING DEANGELIS_ANTON_AUG27. STANDALONE SECTION
molarheatcap_Ta = 25.36; %J/mol K 
atomicmass_Ta = 180.947; %g/mol
length_wire = 0.005; %m, original 0.005 
diameter_wire = 0.0006; %m, original 0.0006
area_wire = pi() * (diameter_wire/2)^2; %2.8 * 10^-7 m^2 

conduct_Ta = 57.5; %W/ m*K
conduct_wire = conduct_Ta * (area_wire/length_wire); %W/K, ~3.3 x 10^-3

density_Ta = 16.678; %g/cm^3
mass_wire = density_Ta * length_wire * area_wire * 10^6;
heatcap_wire = (mass_wire / atomicmass_Ta) * molarheatcap_Ta; 

ni_radius = 0.005; %m
ni_thick = 0.0002; %m 
ni_density = 8.9 * 10^3; %kg / m^3
mass_crystal = pi() * ni_radius^2 * ni_thick * ni_density; %1.398 * 10^-4 kg
C_ni = 26.07; %J/mol K
atomicmass_Ni = 58.69; %g/mol
heatcap_crystal = C_ni * (mass_crystal*1000)/atomicmass_Ni; %0.0621 J/K

T_b = 80; %K


% 'ORIGINAL' CODE: 9/19 - MULTIPLE WIRES + RADIATIVE LOSS TEST, 2D
% WORKS - Ta's actual resistivity + considering multiple wires WRT 6a, 6b
% check out with Fig 3 
fprintf('NEW INSTANCE \n');

n_wires = 4;  % set number of identical wires in parallel

% Electrical: R_total = R_single / n, Heat capacity: Cw_total = n * Cw_single
rho_elec_Ta  = 1.38e-7;                              % ohm*m (Ta @ ~300 K)
R_single     = rho_elec_Ta * length_wire / area_wire;
R_total      = R_single /2 ; %'of two wires'

I_max_single  = 40 /  (n_wires/2);                                  % total available current (A)
I_max = 40;
% Heating parameter with fixed TOTAL current budget:
%rho = (I_max_single^2 * R_total) / Cw_total;         % = I_max^2 * R_single / (n^2 * Cw_single)

% Radiative area scales with n (lateral area only, same geometry as before)
SA_wire_single = pi() * diameter_wire *length_wire;   % m^2
SA_wire_total  = n_wires * SA_wire_single;            % m^2

%fprintf('n_wires=%d  |  kappa_w/rho = %4e \n', n_wires, kap_w/rho);

% ---- ODE with radiation (wire uses total area & total heat capacity) ----
sigma = 5.670374419e-8;  % W/m^2/K^4
T_env = 300;             % K, roomtemp
eps_w = 0.30;            % Ta wire emissivity
eps_a = 0.15;            % Ni crystal emissivity

A_crys = 2*pi*ni_radius^2 + 2*pi*ni_radius*ni_thick; % m^2
t_on   = 2.0;            % s   heating pulse length
t_end  = 8.0;            % s   simulation horizon

I = @(t) (t < t_on);
%{
odefun = @(t,y) [ ...
    rho*I(t) ...
  - kap_w*(y(1) - y(2)) - kap_w*(y(1) - T_b) ...
  - eps_w*sigma*A_wire_total * ( y(1).^4 - T_env^4 ) / Cw_total ;           % wire radiation (n wires)
    kap_a*(y(1) - y(2)) * n_wires ...
  - eps_a*sigma*A_crys * ( y(2).^4 - T_env^4 ) / heatcap_crystal            % crystal radiation
];
%}

odefun = @(t,y) [ ...
    ( (R_total*I(t) * ((I_max )^2)...
  -  (conduct_wire) * (y(1) - y(2)) -  (conduct_wire) * (y(1) - T_b)) ) / (heatcap_wire)  ...
  - eps_w*sigma*SA_wire_total * ( y(1).^4 - T_env^4 ) / heatcap_wire  ;          % wire radiation (n wires)
    (conduct_wire *(y(1) - y(2)) * n_wires)/ heatcap_crystal ...
  - eps_a*sigma*A_crys * ( y(2).^4 - T_env^4 ) / heatcap_crystal       % crystal radiation
];



% ---- Integrate and plot (unchanged below) ----
figure(4); hold on;
y0 = [200; 200];                 % K
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(odefun,[0 t_end],y0,opts);
T_w = y(:,1); T_a = y(:,2);

tI  = linspace(0, t_end, 800);
ItI = arrayfun(I, tI);
tiledlayout(3,1);
nexttile; stairs(tI, ItI, 'LineWidth', 1.8); ylim([-0.05 1.05]);
xlabel('time (s)'); ylabel('I(t) (scaled)'); title('Current command'); grid on;

nexttile; plot(t,T_w,'LineWidth',1.8);
xlabel('time (s)'); ylabel('T_{wire} (K)'); title('T_{wire}(t)'); grid on;

nexttile; plot(t,T_a,'LineWidth',1.8);
xlabel('time (s)'); ylabel('T_{Ni} (K)'); title('T_{Ni}(t)'); grid on; hold off;


% Calculate the heating rate (dT_a/dt) by taking the numerical derivative
dT_a_dt = gradient(T_a, t);

% Find the maximum heating rate
max_heating_rate = max(dT_a_dt);
fprintf('Maximum heating rate of crystal (T_a): %.4f K/s\n', max(dT_a_dt));

% Calculate radiative losses for both wire and crystal
radiative_loss_wire = eps_w * sigma * SA_wire_total * (T_w.^4 - T_env^4);  % W
radiative_loss_crystal = eps_a * sigma * A_crys * (T_a.^4 - T_env^4);     % W
total_radiative_loss = radiative_loss_wire + radiative_loss_crystal;      % W

% Modify the existing plotting section to include 4 subplots
figure(4); 
tiledlayout(4,1);  % Changed from 3,1 to 4,1

nexttile; stairs(tI, ItI, 'LineWidth', 1.8); ylim([-0.05 1.05]);
xlabel('time (s)'); ylabel('I(t) (scaled)'); title('Current command'); grid on;

nexttile; plot(t,T_w,'LineWidth',1.8);
xlabel('time (s)'); ylabel('T_{wire} (K)'); title('T_{wire}(t)'); grid on;

nexttile; plot(t,T_a,'LineWidth',1.8);
xlabel('time (s)'); ylabel('T_{Ni} (K)'); title('T_{Ni}(t)'); grid on;

% New fourth subplot for radiative losses
nexttile; 
plot(t, radiative_loss_wire, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Wire Radiation');
hold on;
plot(t, radiative_loss_crystal, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Crystal Radiation');
plot(t, total_radiative_loss, 'k--', 'LineWidth', 2, 'DisplayName', 'Total Radiation');
xlabel('time (s)'); 
ylabel('Radiative Loss (W)'); 
title('Radiative Heat Loss vs Time');
legend('Location', 'best');
grid on;
hold off;

