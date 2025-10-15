%% PAGE 1 
% SEPT04 RE-DERIVATION OF SECT 3, GIVES FIG 3 
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

R_gas = 8.314; % J/mol K
T_b = 80; %cooling block temp, K 
kap_a = 0.049; % k / m_a * c_a, slow inverse time constant, s^-1 
%kap_a = kapp_a; % SLOW INVERSE TIME CONSTANT CALCULATED FROM WIRE DIMENSIONS
%kap_w = kapp_w; %kappa_w derived from dimensions as opposed to fit
kap_w = 0.81; % k / m_w * c_w, fast inverse time constant, s^-1 GIVEN 


%length_wire = 19/1000; %m
%area_wire = ((.18 / 1000) / 2)^2 * pi();

% GEOMETRY MOD, DISPLACING RHO = 1900 INDEPENDENT OF WIRE DIAMETER
rho_elec_Ta  = 2.2e-7;           % ohm * m (for Ta @ 300K) 1.38 
R     = rho_elec_Ta * length_wire / area_wire;  % ohm 
I_max = 40;                 
rho   = (I_max^2 * R) / heatcap_wire;    % K/s when input=1
%rho = 1900;


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

%% PAGE 3

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

T_b = 300; %K, original 80 K 


% 'ORIGINAL' CODE: 9/19 - MULTIPLE WIRES + RADIATIVE LOSS TEST, 2D
% WORKS - Ta's actual resistivity + considering multiple wires WRT 6a, 6b
% check out with Fig 3 
fprintf('NEW INSTANCE \n');

n_wires = 4;  % set number of identical wires in parallel

% Electrical: R_total = R_single / n, Heat capacity: Cw_total = n * Cw_single
rho_elec_Ta  = 2.2e-7;                              % ohm*m (Ta @ ~300 K)
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
t_on   = 15.0;            % s   heating pulse length
t_end  = 30.0;            % s   simulation horizon

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
    ( (R_single*I(t) * ((I_max_single )^2)...
  -  (conduct_wire) * (y(1) - y(2)) -  (conduct_wire) * (y(1) - T_b)) ) / (heatcap_wire)  ...
  - eps_w*sigma*SA_wire_single * ( y(1).^4 - T_env^4 ) / heatcap_wire  ;          % wire radiation (n wires)
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

%%  9/25 - EQUATION 8 FOR MAX RATE CONTOUR 


% Sweeps
length_wireSPAN   = linspace(5e-3, 20e-3, 150);       % m
diameter_wireSPAN = linspace(0.5e-3, 1.2e-3, 150);    % m
[length_wireSPAN, diameter_wireSPAN] = meshgrid(length_wireSPAN, diameter_wireSPAN);

% Geometry-dependent terms
area_wireSPAN      = pi*(diameter_wireSPAN/2).^2;                          % m^2
conduct_wireSPAN   = conduct_Ta .* (area_wireSPAN ./ length_wireSPAN);     % W/K
mass_wireSPAN_g    = density_Ta .* (length_wireSPAN .* area_wireSPAN) .* 1e6; % g
heatcap_wireSPAN   = (mass_wireSPAN_g ./ atomicmass_Ta) .* molarheatcap_Ta;   % J/K

% Inverse time constants
kapp_wSPAN = conduct_wireSPAN ./ heatcap_wireSPAN;         % s^-1
kapp_aSPAN = conduct_wireSPAN ./ heatcap_crystal;          % s^-1

% Drive tied to geometry
rhoSPAN      = (I_max.^2 .* RSPAN) ./ heatcap_wireSPAN;                    % K/s
Ta0 = 200; % K 
Tb = 80; 
% Closed-form max heating rate
dTa_dt_MAX_SPAN = (kapp_aSPAN .* rhoSPAN) ./ (2 .* kapp_wSPAN) ...
                  - (kapp_aSPAN./2) .* (Ta0 - Tb);                         % K/s

% Plot
figure('Color','w');
hSurf = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, dTa_dt_MAX_SPAN);
shading interp; box on; grid on; hold on;

% Make the surface semi-transparent
set(hSurf, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.3);

% Overlay contour lines
contour3(length_wireSPAN*1e3, diameter_wireSPAN*1e3, dTa_dt_MAX_SPAN, 25, 'LineWidth', 1.5);

xlabel('Wire length (mm)');
ylabel('Wire diameter (mm)');
zlabel('dT_Ni/dt_{MAX} (K/s)');
title('Max crystal heating rate vs wire length & diameter');
view(135,30);

%%  9/26 - VOLTAGE CONTOUR 


% Sweeps
length_wireSPAN   = linspace(1e-3, 20e-3, 150);       % m
diameter_wireSPAN = linspace(0.1e-3, 1.2e-3, 150);    % m
[length_wireSPAN, diameter_wireSPAN] = meshgrid(length_wireSPAN, diameter_wireSPAN);

% Geometry-dependent terms
area_wireSPAN      = pi*(diameter_wireSPAN/2).^2;                          % m^2
conduct_wireSPAN   = conduct_Ta .* (area_wireSPAN ./ length_wireSPAN);     % W/K
mass_wireSPAN_g    = density_Ta .* (length_wireSPAN .* area_wireSPAN) .* 1e6; % g
heatcap_wireSPAN   = (mass_wireSPAN_g ./ atomicmass_Ta) .* molarheatcap_Ta;   % J/K

% Inverse time constants
kapp_wSPAN = conduct_wireSPAN ./ heatcap_wireSPAN;         % s^-1
kapp_aSPAN = conduct_wireSPAN ./ heatcap_crystal;          % s^-1

% Drive tied to geometry
RSPAN        = rho_elec_Ta .* (length_wireSPAN ./ area_wireSPAN);          % ohm

% Closed-form max heating rate
volt_SPAN = I_max_single .* RSPAN; 

% Plot
figure('Color','w');
hSurf = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, volt_SPAN);
shading interp; box on; grid on; hold on;

% Make the surface semi-transparent
set(hSurf, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.3);

% Overlay contour lines
contour3(length_wireSPAN*1e3, diameter_wireSPAN*1e3, volt_SPAN, 25, 'LineWidth', 1.5);

xlabel('Wire length (mm)');
ylabel('Wire diameter (mm)');
zlabel('Voltage');

title('Voltage vs wire length & diameter');
view(135,30);
%% RESISTANCE VS DIMENSIONS, fit test 9/30 
Ta_temp = Ta_scan(:,1)';
Ta_rho = Ta_scan(:,2)';
figure;
hold on;
plot(Ta_temp, Ta_rho);

p = polyfit(Ta_temp, Ta_rho, 3);
Ta_rho_fit = polyval(p, Ta_temp);
plot(Ta_temp, Ta_rho_fit);
%hold off;
fprintf('\n rho(T) = %.4e*T^3 + %.4e*T^2 + %.4e*T + %.4e\n', p);

Ta_rho_FUNCT = @(T) (p(1) * T.^3 + p(2) * T.^2 + p(3) * T + p(4)) * 10^-8;% ohm * m

RSPAN_FUNCT = @(T) Ta_rho_FUNCT(T) .* (length_wireSPAN ./ area_wireSPAN);

RSPAN_FUNCTTEST = RSPAN_FUNCT(300);

figure('Color','w');
hSurf = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, RSPAN_FUNCTTEST);
shading interp; box on; grid on; hold on;

% Make the surface semi-transparent
set(hSurf, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.3);

% Overlay contour lines
contour3(length_wireSPAN*1e3, diameter_wireSPAN*1e3, RSPAN_FUNCTTEST, 25, 'LineWidth', 1.5);

xlabel('Wire length (mm)');
ylabel('Wire diameter (mm)');
zlabel('Resistance (ohm)');

title('Resistance vs wire length & diameter, multiple temperatures');
view(135,30);


%% RESISTANCE VS TEMP AND DIMENSIONS, 10/1 



% Plot surfaces for multiple temperatures (0:100:2000 K)
T_list = 0:20:500;                      % K
figure('Color','w'); hold on; box on; grid on;

% Choose a nice colormap and map T to colors
cmap = parula(numel(T_list));
clim = [min(T_list) max(T_list)];

for k = 1:numel(T_list)
    Tk = T_list(k);
    Rk = RSPAN_FUNCT(Tk);                % same size as your grids
    h = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, Rk);
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.35);    % semi-transparent
    % Color this surface by its temperature
    set(h, 'CData', Tk*ones(size(Rk)), 'FaceColor','interp');
end

% Axes & labels
xlabel('Wire length (mm)');
ylabel('Wire diameter (mm)');
zlabel('Resistance (Ω)');
title('Resistance vs wire length & diameter across temperature (0–2000 K)');

% Colorbar encodes temperature
colormap(parula);
caxis(clim);
zlim([0,3]);
cb = colorbar; cb.Label.String = 'Temperature (K)';

view(135,30);

%% RESISTANCE, BLOCKING OUT MINIMUM 


% Plot surfaces for multiple temperatures (0:200:2000 K)
T_list = 0:300:1500;                      % K
figure('Color','w'); hold on; box on; grid on;

cmap = parula(numel(T_list));
clim = [min(T_list) max(T_list)];

for k = 1:numel(T_list)
    Tk = T_list(k);
    Rk = RSPAN_FUNCT(Tk);

    % mask regions
    mask_low  = (Rk < 0.01);
    mask_high = ~mask_low;

    % normal part (>= 0.01 Ω)
    h1 = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, Rk);
    set(h1, 'EdgeColor','none', 'FaceAlpha',0.6);
    set(h1, 'CData', Tk*ones(size(Rk)), 'FaceColor','interp');
    % hide the low values
    Rk_masked = Rk; Rk_masked(mask_low) = NaN;
    set(h1, 'ZData', Rk_masked);

    % gray overlay for < 10 mΩ
    Rk_gray = Rk; Rk_gray(mask_high) = NaN;
    h2 = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, Rk_gray);
    set(h2, 'EdgeColor','none', 'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.3);
end

xlabel('Wire length (mm)');
ylabel('Wire diameter (mm)');
zlabel('Resistance (Ω)');
title('Resistance vs wire length & diameter across temperature (0–2000 K)');

colormap(parula);
caxis(clim);
zlim([0,3]);
cb = colorbar; cb.Label.String = 'Temperature (K)';

view(135,30);
%%

% === Heating rate 3D plot with masks: V>3 V or R<10 mΩ shown in gray ===
set(gcf,'Renderer','opengl');  % smoother transparency

% Build a single mask from your voltage & resistance grids
mask_gray = (volt_SPAN > 3) | (RSPAN < 0.01);   % gray-out condition
mask_main = ~mask_gray;

Z_main = dTa_dt_MAX_SPAN;    Z_main(mask_gray) = NaN;   % keep valid region colored
Z_gray = dTa_dt_MAX_SPAN;    Z_gray(mask_main) = NaN;   % show masked region in gray

figure('Color','w'); hold on; box on; grid on;

% Main colored surface (valid region)
hMain = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, Z_main, ...
    'EdgeColor','none', 'FaceAlpha',0.6, 'FaceColor','interp');

% Gray overlay (invalid region: V>3 or R<10 mΩ)
hGray = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, Z_gray, ...
    'EdgeColor','none', 'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.35);

% Contours for the main (optional—skip if you want maximum speed)
contour3(length_wireSPAN*1e3, diameter_wireSPAN*1e3, Z_main, 25, 'LineWidth', 1.2);

% Labels & colorbar (color encodes heating-rate magnitude)
xlabel('Wire length (mm)');
ylabel('Wire diameter (mm)');
zlabel('dT_a/dt_{MAX} (K/s)');
title('Max crystal heating rate vs wire length & diameter');

colormap(parula);
cb = colorbar; cb.Label.String = 'dT_a/dt_{MAX} (K/s)';

view(135,30);

% (Optional) Add a tiny legend entry for the gray region:
hFake = plot3(NaN,NaN,NaN,'s','MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerEdgeColor','none','MarkerSize',10);
legend([hFake],{'V > 3 V or R < 10 mΩ'},'Location','northeast');
%%

% === Max heating rate vs geometry across temperature, with gray masks ===
T_list = 0:500:1500;                       % K
figure('Color','w'); hold on; box on; grid on;
set(gcf,'Renderer','opengl'); opengl hardware;

clim = [min(T_list) max(T_list)];
colormap(parula);

for k = 1:numel(T_list)
    Tk   = T_list(k);

    % Temperature-dependent resistance on the same (L,D) grid
    Rk   = RSPAN_FUNCT(Tk);                % Ω (same size as length/diameter grids)
    Vk   = I_max .* Rk;                    % V (driving current I_max)

    % Geometry-dependent terms (precomputed above in your script):
    % heatcap_wireSPAN  (J/K), kapp_wSPAN (s^-1), kapp_aSPAN (s^-1)

    % Joule input per K (K/s) term that ties drive to geometry
    rhoSPAN_k = (I_max.^2 .* Rk) ./ heatcap_wireSPAN;   % K/s

    % Closed-form max heating rate at this temperature
    Zk = (kapp_aSPAN .* rhoSPAN_k) ./ (2 .* kapp_wSPAN) ...
         - (kapp_aSPAN./2) .* (Ta0 - Tb);               % K/s

    % Build masks: gray-out where V>3 V OR R<10 mΩ
    mask_gray = (Vk > 3) | (Rk < 0.01);
    mask_main = ~mask_gray;

    Z_main = Zk;  Z_main(mask_gray) = NaN;
    Z_gray = Zk;  Z_gray(mask_main) = NaN;

    % ---- Main colored surface (valid region) ----
    hMain = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, Z_main, ...
        'EdgeColor','none', 'FaceAlpha',0.60, 'FaceColor','interp');
    % color by temperature (uniform per surface)
    set(hMain,'CData', Tk*ones(size(Z_main)));

    % ---- Gray overlay (masked-out region) ----
    hGray = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, Z_gray, ...
        'EdgeColor','none', 'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.35);

    % Optional contours for the valid region (comment out for speed)
    % contour3(length_wireSPAN*1e3, diameter_wireSPAN*1e3, Z_main, 20, 'LineWidth',1.0);
end

% Axes, labels, and colorbar
xlabel('Wire length (mm)');
ylabel('Wire diameter (mm)');
zlabel('dT_a/dt_{MAX} (K/s)');
title('Max crystal heating rate vs wire length & diameter (T-colored)');
caxis(clim);
cb = colorbar; cb.Label.String = 'Temperature (K)';

view(135,30);
% If ranges are wild, you can clamp Z (optional):
% zlim([max(0, min(Zk(:))), max(Zk(:))]);

% Legend entry for gray condition
hFake = plot3(NaN,NaN,NaN,'s','MarkerFaceColor',[0.5 0.5 0.5], ...
    'MarkerEdgeColor','none','MarkerSize',10);
legend([hFake], {'V > 3 V or R < 10 mΩ'}, 'Location','northeast');

%%

%wire time constant contour, 10/1
% 9/26 - VOLTAGE CONTOUR 


% Sweeps
length_wireSPAN   = linspace(1e-3, 20e-3, 150);       % m
diameter_wireSPAN = linspace(0.1e-3, 1.2e-3, 150);    % m
[length_wireSPAN, diameter_wireSPAN] = meshgrid(length_wireSPAN, diameter_wireSPAN);

% Geometry-dependent terms
area_wireSPAN      = pi*(diameter_wireSPAN/2).^2;                          % m^2
conduct_wireSPAN   = conduct_Ta .* (area_wireSPAN ./ length_wireSPAN);     % W/K
mass_wireSPAN_g    = density_Ta .* (length_wireSPAN .* area_wireSPAN) .* 1e6; % g
heatcap_wireSPAN   = (mass_wireSPAN_g ./ atomicmass_Ta) .* molarheatcap_Ta;   % J/K

% Inverse time constants
kapp_wSPAN = conduct_wireSPAN ./ heatcap_wireSPAN;         % s^-1
kapp_aSPAN = conduct_wireSPAN ./ heatcap_crystal;          % s^-1

tau_wSPAN = heatcap_wireSPAN ./ conduct_wireSPAN;
% Drive tied to geometry
RSPAN        = rho_elec_Ta .* (length_wireSPAN ./ area_wireSPAN);          % ohm

% Closed-form max heating rate
volt_SPAN = I_max_single .* RSPAN; 

% Plot
figure('Color','w');
hSurf = surf(length_wireSPAN*1e3, diameter_wireSPAN*1e3, tau_wSPAN);
shading interp; box on; grid on; hold on;

% Make the surface semi-transparent
set(hSurf, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.3);

% Overlay contour lines
contour3(length_wireSPAN*1e3, diameter_wireSPAN*1e3, tau_wSPAN, 25, 'LineWidth', 1.5);

xlabel('Wire length (mm)');
ylabel('Wire diameter (mm)');
zlabel('Time constant (s)');

title('Wire time constant vs wire length & diameter');
view(135,30);

%%

% --- 2D plot: tau_w vs wire length at a chosen diameter ---
target_d_mm = 0.6;                       % choose diameter in mm (your original)
[~, idx] = min(abs(diameter_wireSPAN(:,1)*1e3 - target_d_mm));  % nearest row

Lmm = length_wireSPAN(1,:)*1e3;          % mm (x-axis)
tau_line = tau_wSPAN(idx,:);             % s  (y-axis)

figure('Color','w'); 
plot(Lmm, tau_line, 'LineWidth', 2);
grid on; box on;
xlabel('Wire length (mm)');
ylabel('\tau_w (s)');
title(sprintf('\\tau_w vs wire length', diameter_wireSPAN(idx,1)*1e3));

% Optional: if tau spans decades, use log y-scale:
% set(gca, 'YScale','log');
