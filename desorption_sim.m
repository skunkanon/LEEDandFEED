function [tmp_span,rate_span] = desorption_sim(N_0)

%{
2/21 notes:
ultimately should be able to simulate spectra as per falconer-maddix
method with the g(c) and everything; take real data, find activation E and
pre-exponential factor. 

2/23: 

pretty straightforward to get index of max dNdt for peak
temperature. dont know why it starts exploding at higher beta. paper said
35 K/s is max but can't handle even 15. 

%}

%some constants
v = 10^13; %pre-exponential
R = 8.314; %gas constant in J/(K*mol)
%N_0 = 0.6; %initial surface 
kcal_to_j = 4184; %conversion factor
Ea_0 = 60 * 10^3; %initial activation E in J/mol

y_E = 4 * kcal_to_j; %4kcal/mol into j/mol, eq) 14 

%defining experimental conditions
beta = 9; % Kelvin/s
init_tmp =70; %initial temperature
ds_time = 30; %define desorption duration in seconds
t_span = linspace(0,ds_time,1000); %array of times from start through whole duration

%temperature as function of time 
tmp = @(t) beta*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

%coverage-dependent activation energy
Ea = @(N) Ea_0 - N*y_E;

%solve coverage differential equation 
dNdt = @(t, N) -v * N * min(exp(-Ea(N) / (R * tmp(t))), 2);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); % Adjust tolerances
[t, N] = ode15s(dNdt, t_span, N_0, options);

%make coords for plot
rate_span = arrayfun(@(t, N) -dNdt(t, N), t, N);
%Ea_span = arrayfun(@(N) Ea(N),N);


[M,I] = max(rate_span);



end