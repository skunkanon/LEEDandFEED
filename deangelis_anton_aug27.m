%% 9/9 - FIG 3 DERIVE (eq 8 visualization as fct of radius and length of wire; add radiative loss later) 
%9/11 - really just replication of doc's calculations, adjustable wire
%dimensions

molarheatcap_Ta = 25.36; %J/mol K 
atomicmass_Ta = 180.947; %g/mol
length_wire = 0.01 ; %m, original 0.005 
diameter_wire = 0.0003; %m, original 0.0006
area_wire = pi() * (diameter_wire/2)^2; %2.8 * 10^-7 m^2 

conduct_Ta = 57.5; %W/ m*K
conduct_wire = conduct_Ta * (area_wire/length_wire); %W/K, ~3.3 x 10^-3

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




% GEOMETRY MOD, DISPLACING RHO = 1900 INDEPENDENT OF WIRE DIAMETER
rho_elec_Ta  = 1.348e-7;           % ohm * m (for Ta @ 300K)
R     = rho_elec_Ta * length_wire / area_wire;  % ohm 
I_max = 40;                 
rho   = (I_max^2 * R) / heatcap_wire;    % K/s when input=1



%rho = 1900; % 'effective wire resistivity for rapid heating, i_max^2 * R / m_w * c_w
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

figure(1);
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

%% 9/16 - OVERLAYING CRYSTAL TEMP VS TIME AND DIAMETER AND ALL 

figure(2);
hold on;
plot(t, T_a, 'LineWidth', 2, ...
     'DisplayName', sprintf('Length = %.2e, Diameter = %.2e', length_wire, diameter_wire));
xlabel('time (s)');
ylabel('T_{Ni} (K)');
title('T_{Ni}(t)');
grid on;
legend show;



%%
