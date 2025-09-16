%% 9/9 - FIG 3 DERIVE (eq 8 visualization as fct of radius and length of wire; add radiative loss later) 
%9/11 - really just replication of doc's calculations, adjustable wire
%dimensions

molarheatcap_Ta = 25.36; %J/mol K 
atomicmass_Ta = 180.947; %g/mol
length_wire = 0.005; %m, original 0.005 
diameter_wire = 0.0006 * 10; %m, original 0.0006
area_wire = pi() * (diameter_wire/2)^2; %2.8 * 10^-7 m^2 

conduct_Ta = 57.5; %W/ m*K
conduct_wire = conduct_Ta * (area_wire/length_wire); %1/ohm, ~3.3 x 10^-3

density_Ta = 16.678; %g/cm^3
mass_wire = density_Ta * length_wire * area_wire * 10^6;
heatcap_wire = (mass_wire / atomicmass_Ta) * molarheatcap_Ta; 


kapp_w = conduct_wire/heatcap_wire;


%% SEPT04 RE-DERIVATION OF SECT 3, GIVES FIG 3 

R_gas = 8.314; % J/mol K
T_b = 80; %cooling block temp, K 
kap_a = 0.049; % k / m_a * c_a, slow inverse time constant, s^-1 
kap_w = kapp_w; %kappa_w derived from dimensions as opposed to fit
%kap_w = 0.81; % k / m_w * c_w, fast inverse time constant, s^-1 GIVEN 
rho = 1900; % 'effective wire resistivity for rapid heating, i_max^2 * R / m_w * c_w
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

figure;
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


% ^ CURRENT IS SCALED AS (CURRENT_ACTUAL/CURRENT_MAX)^2


%% 9/15 - LENGTH OR AREA VS MAX HEATING RATE 

% Length span (m)
L_vec = linspace(1e-3, 50e-3, 40);      % 1 mm -> 200 mm

% Diameters
d_vec = linspace(0.0006, 0.001, 10); ;    

% Storage: rows = diameters, cols = lengths
max_dTadt_mat = zeros(numel(d_vec), numel(L_vec));

for id = 1:numel(d_vec)
    d = d_vec(id);
    A = pi*(d/2)^2;                     % m^2

    for k = 1:numel(L_vec)
        L = L_vec(k);

        % Geometry-dependent kappa_w for this L, d
        conduct_wire = conduct_Ta * (A / L);                         % W/K
        mass_wire_g  = density_Ta * L * A * 1e6;                    % g
        heatcap_wire = (mass_wire_g / atomicmass_Ta) * molarheatcap_Ta; % J/K
        kap_w_local  = conduct_wire / heatcap_wire;                  % s^-1

        % ODE with this geometry's kap_w; same structure as your second section
        odefun_geo = @(t,y) [ ...
            rho*I(t) - kap_w_local*(y(1) - y(2)) - kap_w_local*(y(1) - T_b); ... % dT_w/dt
            kap_a*(y(1) - y(2))                                                   % dT_a/dt (crystal)
        ];

        % Solve
        [tL, yL] = ode45(odefun_geo, [0 t_end], y0, opts);

        % Compute instantaneous crystal heating rate along the trajectory
        dTadt = zeros(size(tL));
        for j = 1:numel(tL)
            dy = odefun_geo(tL(j), yL(j,:).');   % dy = [dT_w/dt; dT_a/dt]
            dTadt(j) = dy(2);
        end

        % Max crystal heating rate while current is ON
        on_idx = tL <= (t_on + 1e-9);
        max_dTadt_mat(id, k) = max(dTadt(on_idx));
    end
end

% Plot
figure; hold on; grid on;
for id = 1:numel(d_vec)
    plot(L_vec*1e3, max_dTadt_mat(id,:), 'LineWidth', 1.8);
end
xlabel('Wire length L (mm)');
ylabel('Max crystal heating rate (K/s)');
legtxt = arrayfun(@(d) sprintf('d = %.2f mm', d*1e3), d_vec, 'UniformOutput', false);
legend(legtxt, 'Location', 'northeast');
title('Max crystal heating rate vs wire dimensions');
