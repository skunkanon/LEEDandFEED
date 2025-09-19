%% 9/9 - FIG 3 DERIVE (eq 8 visualization as fct of radius and length of wire; add radiative loss later) 
%9/11 - really just replication of doc's calculations, adjustable wire
%dimensions

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


kapp_w = conduct_wire/heatcap_wire;

%% 9/18 - 3D PLOT OF KAPP_W 

length_wireSPAN = linspace(3/1000, 6/1000, 15);    % wire length [m]
diameter_wireSPAN = linspace(0.1/1000,1.2/1000, 15); % wire diameter [m]
%diameter_wireSPAN = 0.0006;
[length_wireSPAN, diameter_wireSPAN] = meshgrid(length_wireSPAN, diameter_wireSPAN);
area_wireSPAN = pi*(diameter_wireSPAN/2).^2;                 % m^2 cross-sectional area

conduct_wireSPAN = conduct_Ta .* (area_wireSPAN ./ length_wireSPAN);                   % W/K
mass_wireSPAN  = density_Ta .* (length_wireSPAN .* area_wireSPAN) .* 1e6;            % g  (m^3 -> cm^3)
heatcap_wireSPAN = (mass_wireSPAN ./ atomicmass_Ta) .* molarheatcap_Ta; % J/K

% Inverse time constant
kapp_wSPAN = conduct_wireSPAN ./ heatcap_wireSPAN;                  % 1/s

% --- Plot ---
figure('Color','w');
surf(length_wireSPAN, diameter_wireSPAN, kapp_wSPAN);
shading interp; box on; grid on; hold on;
contour3(length_wireSPAN, diameter_wireSPAN, kapp_wSPAN, 12, 'k:', 'LineWidth', 0.8); % iso-kappa lines

xlabel('Length L (mm)', 'FontSize', 12);
ylabel('Diameter D (mm)', 'FontSize', 12);
zlabel('\kappa_w (s^{-1})', 'FontSize', 12);
title('\kappa_w vs. Wire Length and Diameter (Ta)', 'FontSize', 13);
colorbar; view(-35, 28);

% Optional: show Z on log scale if the L-range is wide
% set(gca,'ZScale','log');

% --- (Optional) sanity check: diameter-independence ---
% For any fixed L, kappa_w should be (numerically) constant vs D.
% rel_range = max(kappa_w) - min(kappa_w) along D, normalized by mean
% rr = max(range(kappa_w,1) ./ mean(kappa_w,1));
% fprintf('Max relative variation along D at fixed L: %.2e\n', rr);


%% 9/18 - KAPP_A

ni_radius = 0.005; %m
ni_thick = 0.0002; %m 
ni_density = 8.9 * 10^3; %kg / m^3
mass_crystal = pi() * ni_radius^2 * ni_thick * ni_density; %1.398 * 10^-4 kg
C_ni = 26.07; %J/mol K
atomicmass_Ni = 58.69; %g/mol
heatcap_crystal = C_ni * (mass_crystal*1000)/atomicmass_Ni; %0.0621 J/K

kapp_a = conduct_wire/heatcap_crystal; %0.052 s^-1

%% 9/18 3D PLOT OF KAPP_A
% Inverse time constant
kapp_ASPAN = conduct_wireSPAN ./ heatcap_crystal;                  % 1/s

% --- Plot ---
figure('Color','w');
surf(length_wireSPAN, diameter_wireSPAN, kapp_ASPAN);
shading interp; box on; grid on; hold on;
contour3(length_wireSPAN, diameter_wireSPAN, kapp_ASPAN, 12, 'k:', 'LineWidth', 0.8); % iso-kappa lines

xlabel('Length L (mm)', 'FontSize', 12);
ylabel('Diameter D (mm)', 'FontSize', 12);
zlabel('\kappa_A (s^{-1})', 'FontSize', 12);
title('\kappa_A vs. Wire Length and Diameter (Ta)', 'FontSize', 13);
colorbar; view(-35, 28);



%% SEPT04 RE-DERIVATION OF SECT 3, GIVES FIG 3 

R_gas = 8.314; % J/mol K
T_b = 80; %cooling block temp, K 
%kap_a = 0.049; % k / m_a * c_a, slow inverse time constant, s^-1 
kap_a = kapp_a; % SLOW INVERSE TIME CONSTANT CALCULATED FROM WIRE DIMENSIONS
kap_w = kapp_w; %kappa_w derived from dimensions as opposed to fit
%kap_w = 0.81; % k / m_w * c_w, fast inverse time constant, s^-1 GIVEN 




% GEOMETRY MOD, DISPLACING RHO = 1900 INDEPENDENT OF WIRE DIAMETER
rho_elec_Ta  = 2.2e-7;           % ohm * m (for Ta @ 300K) 1.38 
R     = rho_elec_Ta * length_wire / area_wire;  % ohm 
I_max = 40;                 
rho   = (I_max^2 * R) / heatcap_wire;    % K/s when input=1



%rho = 1900; % GIVEN 'effective wire resistivity for rapid heating, i_max^2 * R / m_w * c_w
fprintf('kappa_w / rho = %4e \n', kap_w/rho);

% time and input
t_on   = 2.0;            % s   heating pulse length
t_end  = 8.0;            % s   simulation horizon
I = @(t) (t < t_on);     % **scaled** current: 1 while on, 0 after


% ODE:  y(1)=T_w (wire), y(2)=T_a (adsorbent/crystal)
odefun = @(t,y) [ ...
    rho*I(t) - kap_w*(y(1) - y(2)) - kap_w*(y(1) - T_b) ;  % dT_w/dt
    kap_a*(y(1) - y(2))                                    % dT_a/dt
];

figure(2);
hold on;

y0 = [200; 200];                 % K  initial wire & crystal temps
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(odefun,[0 t_end],y0,opts);
T_w = y(:,1);
T_a = y(:,2);
tI  = linspace(0, t_end, 800);
ItI = arrayfun(I, tI); 
tiledlayout(3,1);  % 1) current, 2) T_wire, 3) T_crystal

% 1) current vs time
nexttile;
stairs(tI, ItI, 'LineWidth', 1.8);
ylim([-0.05 1.05]);
xlabel('time (s)'); ylabel('I(t) (scaled)');
title('Current command'); grid on;

nexttile;
plot(t,T_w,'LineWidth',1.8);
xlabel('time (s)'); ylabel('T_{wire} (K)');
title('T_{wire}(t)'); grid on;

nexttile;
plot(t,T_a,'LineWidth',1.8);
xlabel('time (s)'); ylabel('T_{Ni} (K)');
title('T_{Ni}(t)'); grid on;
hold off;

% (Optional) overlay the input to sanity-check the timing:
% figure; plot(t, arrayfun(I,t)); ylim([-0.1 1.1]); title('Scaled current I(t)');

t_NOLOSS = t;

% ^ CURRENT IS SCALED AS (CURRENT_ACTUAL/CURRENT_MAX)^2

%% 9/16 - OVERLAYING CRYSTAL TEMP VS TIME AND DIAMETER AND ALL 

figure(4);
hold on;
plot(t, T_a, 'LineWidth', 2, ...
     'DisplayName', sprintf('Length = %.2e, Diameter = %.2e', length_wire, diameter_wire));
xlabel('time (s)');
ylabel('T_{Ni} (K)');
title('T_{Ni}(t)');
grid on;
legend show;



%% 9/18 3D PLOT OF MAX RATES 


RSPAN   = rho_elec_Ta .* (length_wireSPAN ./ area_wireSPAN);   % ohm
rhoSPAN = (I_max.^2 .* RSPAN) ./ heatcap_wireSPAN;             % K/s (u=1)

% --- Loop over grid using rhoSPAN values ---
u = @(t) double(t < t_on);

max_dTa = zeros(size(kapp_ASPAN));
for j = 1:numel(kapp_ASPAN)
    kw  = kapp_wSPAN(j);
    ka  = kapp_ASPAN(j);
    rho = rhoSPAN(j);   % geometry-dependent rho

    odefun_j = @(t,y) [ ...
        rho*u(t) - kw*(y(1)-y(2)) - kw*(y(1)-T_b);  % dT_w/dt
        ka*(y(1)-y(2))                               % dT_a/dt
    ];

    [t, y] = ode45(odefun_j, [0 t_end], y0, opts);
    dTa = ka .* (y(:,1) - y(:,2));
    max_dTa(j) = max(dTa);
end

% --- Plot ---
figure('Color','w');
surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, max_dTa);
shading interp; box on; grid on; hold on;
contour3(length_wireSPAN*1e3, diameter_wireSPAN*1e3, max_dTa, 12, 'k:', 'LineWidth', 0.7);
xlabel('Wire Length (mm)'); ylabel('Wire Diameter (mm)'); zlabel('max dT_{Ni}/dt (K/s)');
title('Max crystal heating rate vs. wire length & diameter (geom. \rho)');
colorbar; view(-35, 28);

%% 9/19 - RADIATIVE LOSS TEST, 2D 
length_wire = 0.01;
diameter_wire = 0.0005;
area_wire = (diameter_wire/2)^2 * pi();

% SEPT04 RE-DERIVATION OF SECT 3, GIVES FIG 3 

R_gas = 8.314; % J/mol K
T_b = 80; %cooling block temp, K 
%kap_a = 0.049; % k / m_a * c_a, slow inverse time constant, s^-1 
kap_a = kapp_a; % SLOW INVERSE TIME CONSTANT CALCULATED FROM WIRE DIMENSIONS
kap_w = kapp_w; %kappa_w derived from dimensions as opposed to fit
%kap_w = 0.81; % k / m_w * c_w, fast inverse time constant, s^-1 GIVEN 

% GEOMETRY MOD, DISPLACING RHO = 1900 INDEPENDENT OF WIRE DIAMETER
rho_elec_Ta  = 2.2e-7;                 % ohm*m (Ta @ ~300 K) DOC ORIGINAL 2.2 * 10^-7
R     = rho_elec_Ta * length_wire / area_wire;   % ohm 
I_max = 40;                 
rho   = (I_max^2 * R) / heatcap_wire;            % K/s when input=1

% ^ CURRENT IS SCALED AS (CURRENT_ACTUAL/CURRENT_MAX)^2
fprintf('kappa_w / rho = %4e \n', kap_w/rho);

% time and input
t_on   = 2.0;            % s   heating pulse length
t_end  = 8.0;            % s   simulation horizon
I = @(t) (t < t_on);     % **scaled** current: 1 while on, 0 after

% ---------- Radiative loss (wire + crystal) ----------
sigma = 5.670374419e-8;  % W/m^2/K^4
T_env = 300;             % K, chamber/walls
eps_w = 0.30;            % wire emissivity (Ta, hot/oxidized ~0.2–0.4)
eps_a = 0.15;            % crystal emissivity (Ni, finish-dependent)

% Crystal area (disk: two faces + rim)
A_crys = 2*pi*ni_radius^2 + 2*pi*ni_radius*ni_thick; % m^2
% Wire lateral area using only length_wire & area_wire:
% D = 2*sqrt(area/pi) => A_wire = circumference*L = pi*D*L = 2*sqrt(pi*area)*L
A_wire = 2*sqrt(pi*area_wire)*length_wire;           % m^2
% -----------------------------------------------------

% ODE:  y(1)=T_w (wire), y(2)=T_a (adsorbent/crystal)
odefun = @(t,y) [ ...
    rho*I(t) ...
  - kap_w*(y(1) - y(2)) - kap_w*(y(1) - T_b) ...
  - eps_w*sigma*A_wire * ( y(1).^4 - T_env^4 ) / heatcap_wire ;          % wire radiation
    kap_a*(y(1) - y(2)) ...
  - eps_a*sigma*A_crys * ( y(2).^4 - T_env^4 ) / heatcap_crystal         % crystal radiation
];

figure(2);
hold on;

y0 = [200; 200];                 % K  initial wire & crystal temps
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t_LOSS,y] = ode45(odefun,[0 t_end],y0,opts);
T_w_LOSS = y(:,1);
T_a_LOSS = y(:,2);
tI  = linspace(0, t_end, 800);
ItI = arrayfun(I, tI); 
tiledlayout(3,1);  % 1) current, 2) T_wire, 3) T_crystal

% 1) current vs time
nexttile;
stairs(tI, ItI, 'LineWidth', 1.8);
ylim([-0.05 1.05]);
xlabel('time (s)'); ylabel('I(t) (scaled)');
title('Current command'); grid on;

nexttile;
plot(t_LOSS,T_w_LOSS,'LineWidth',1.8);
xlabel('time (s)'); ylabel('T_{wire} (K)');
title('T_{wire}(t)'); grid on;

nexttile;
plot(t_LOSS,T_a_LOSS,'LineWidth',1.8);
xlabel('time (s)'); ylabel('T_{Ni} (K)');
title('T_{Ni}(t)'); grid on;
hold off;

figure(3); hold on;
plot(t_LOSS, T_a_LOSS);



%% 9/19 - MULTIPLE WIRES + RADIATIVE LOSS TEST, 2D
% WORKS - Ta's actual resistivity + considering multiple wires WRT 6a, 6b
% check out with Fig 3 

n_wires = 4;  % set number of identical wires in parallel

% Electrical: R_total = R_single / n, Heat capacity: Cw_total = n * Cw_single
rho_elec_Ta  = 1.38e-7;                              % ohm*m (Ta @ ~300 K)
R_single     = rho_elec_Ta * length_wire / area_wire;
R_total      = R_single;

I_max_total  = 40 /  (n_wires/2);                                  % total available current (A)
Cw_total     =  heatcap_wire;

% Heating parameter with fixed TOTAL current budget:
rho = (I_max_total^2 * R_total) / Cw_total;         % = I_max^2 * R_single / (n^2 * Cw_single)

% Radiative area scales with n (lateral area only, same geometry as before)
A_wire_single = 2*sqrt(pi*area_wire)*length_wire;   % m^2
A_wire_total  = n_wires * A_wire_single;            % m^2

fprintf('n_wires=%d  |  kappa_w/rho = %4e \n', n_wires, kap_w/rho);

% ---- ODE with radiation (wire uses total area & total heat capacity) ----
sigma = 5.670374419e-8;  % W/m^2/K^4
T_env = 300;             % K
eps_w = 0.30;            % Ta wire emissivity
eps_a = 0.15;            % Ni crystal emissivity

A_crys = 2*pi*ni_radius^2 + 2*pi*ni_radius*ni_thick; % m^2

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
    (R*I(t) * (20 ^2/2)...
  - conduct_wire*(y(1) - y(2)) - conduct_wire*(y(1) - T_b)) / heatcap_wire ...
  - eps_w*sigma*A_wire_total * ( y(1).^4 - T_env^4 ) / Cw_total  ;          % wire radiation (n wires)
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

%%

T_DIFF_WIRE_REAL = T_w.^4 - 300^4;
plot(t, T_DIFF_WIRE_REAL * eps_w * sigma*A_wire_total);

%% 9/19 RADIATIVE LOSS TEST 3D, ONE WIRE 



RSPAN   = rho_elec_Ta .* (length_wireSPAN ./ area_wireSPAN);   % ohm
rhoSPAN = (I_max.^2 .* RSPAN) ./ heatcap_wireSPAN;             % K/s (u=1)

% --- Loop over grid using rhoSPAN values ---
u = @(t) double(t < t_on);

max_dTa = zeros(size(kapp_ASPAN));
for j = 1:numel(kapp_ASPAN)
    kw  = kapp_wSPAN(j);
    ka  = kapp_ASPAN(j);
    rho = rhoSPAN(j);   % geometry-dependent rho

    odefun_j = @(t,y) [ ...
        rho*u(t) - kw*(y(1)-y(2)) - kw*(y(1)-T_b);  % dT_w/dt
        ka*(y(1)-y(2))                               % dT_a/dt
    ];

    [t, y] = ode45(odefun_j, [0 t_end], y0, opts);
    dTa = ka .* (y(:,1) - y(:,2));
    max_dTa(j) = max(dTa);
end

% --- Plot ---
figure('Color','w');
surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, max_dTa);
shading interp; box on; grid on; hold on;
contour3(length_wireSPAN*1e3, diameter_wireSPAN*1e3, max_dTa, 12, 'k:', 'LineWidth', 0.7);
xlabel('Wire Length (mm)'); ylabel('Wire Diameter (mm)'); zlabel('max dT_{Ni}/dt (K/s)');
title('Max crystal heating rate vs. wire length & diameter (geom. \rho)');
colorbar; view(-35, 28);

%% POWER TEST

T_w_DIFF4 = (T_w_LOSS .^4 - 200^4) .* diameter_wire * pi() * length_wire * eps_w * sigma;

%T_w_DIFF4 = (T_w_LOSS .^4 - 300^4) .* area_wire * eps_w * sigma;

T_a_DIFF4 = (T_a_LOSS .^ 4 - 200^4) .* A_crys * eps_a * sigma;
